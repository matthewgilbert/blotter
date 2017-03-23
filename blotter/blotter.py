import pandas as pd
import numpy as np
import json
from collections import namedtuple
from array import array
from . import marketdata


class _Event():
    # this class manages the actions which are performed on the Holdings class
    # and encapsulates all the data required for this action. The class
    # proviedes a consistent interface regardless of whether the event is
    # reconstituted from a string representation in a trading log using
    # Blotter.read_log or created from within the context of Blotter. The class
    # also manages the implementation details of turning events into a string
    # representation for appending to an event log.

    def __init__(self, event_type, data):
        # data is dictionary of parameters for corresponding Holdings function
        self._type = event_type
        self._data = data

    @classmethod
    def fromstring(cls, event_str):
        parts = event_str.split("|")
        event_type = parts[0]
        data = json.loads(parts[1])
        data['timestamp'] = pd.Timestamp(data['timestamp'])
        if 'prices' in data:
            data['prices'] = pd.Series(data['prices'])
        return cls(event_type, data)

    @property
    def type(self):
        """
        Returns the event type
        """
        return self._type

    @property
    def data(self):
        """
        Returns the data associated to the event
        """
        return self._data

    def __repr__(self):
        # sorted dictionary representation except with timestamp as first entry
        str_repr = (self.type +
                    '|{"timestamp": ' + json.dumps(str(self.data['timestamp']))
                    )
        keys = list(self.data.keys())
        keys.sort()
        for key in keys:
            if key == 'timestamp':
                continue

            jkey = json.dumps(key)
            if key == 'prices':
                str_repr = (str_repr + ", " + jkey + ": " +
                            self.data[key].to_json())
            else:
                str_repr += ', ' + jkey + ': ' + json.dumps(self.data[key])

        str_repr += '}'

        return str_repr


_metadata = namedtuple('metadata', ['ccy', 'margin', 'multiplier',
                                    'commission', 'isFX'])


class Blotter():
    """
    This is a financial blotter which is used for maintaing positions and PnL
    in conjunction with backtesting a strategy historically.
    The main purpose is for calculating historical PnL for both open and
    closed trades as well as maintaining marked to market holdings.

    This class maintains market pricing data for marking holdings to market as
    well as interest rate data for charging interest and margin on open
    positions. The class manages interest and margin charges based of a user
    defined time of day and also provides functionality for repatriating closed
    PnL to the user defined base currency on a daily user defined time.
    """

    def __init__(self,
                 prices=None,
                 interest_rates=None,
                 accrual_time=pd.Timedelta(16, unit="h"),
                 eod_time=pd.Timedelta(16, unit="h"),
                 sweep_time=pd.Timedelta(16, unit="h"),
                 base_ccy="USD",
                 margin_charge=0.015):
        """
        Parameters
        ----------
        prices: str
            path to folder of data for all traded instruments. Refer to
            blotter.MarketData for more information on file format. Names for
            FX instruments should be in the form 'XXXYYY' where XXX is the
            first currency and YYY is the second currency, e.g. AUDUSD, USDCAD
        interest_rates: str
            Path to csv of data for all traded interest bearing instruments.
            These rates should be daily annualized rates. Refer to
            blotter.MarketData for more information on file format.
        accrual_time: pandas.Timedelta
            Time of day which interest is charged/paid for interest bearing
            instruments (FX) as well as margin costs, if no automatic charges
            desired set to None
        eod_time: pandas.Timedelta
            End of day time used for automatic PnL calculation, if no automatic
            PnL calculation desired set to None
        sweep_time: pandas.Timedelta
            Automatic end of day time used for sweeping PnL calculation, if no
            automatic sweeping is desired set to None
        base_ccy: str
            Base currency of blotter, used when sweeping pnl to base currency
        margin_charge: float
            Interest rate spread above daily base currency interest rate which
            is paid on margin, e.g. if daily interest rate is 0.5% annualized,
            margin_charge=0.015 implies daily balance paid on margin is
            (0.005 + 0.015)/365
        """

        actions = []
        if accrual_time is not None:
            actions.append((accrual_time, "INTEREST"))
            actions.append((accrual_time, "MARGIN"))
        if eod_time is not None:
            actions.append((eod_time, "PNL"))
        if sweep_time is not None:
            actions.append((sweep_time, "PNL_SWEEP"))
        self._actions = actions

        self._base_ccy = base_ccy
        self._margin_charge = margin_charge
        self._event_log = []
        # dictionaries of instrument level data
        self._gnrc_meta = dict()
        self._instr_map = dict()

        self._prices = prices
        self._rates = interest_rates

        self._holdings = Holdings()
        self.get_holdings_history = self._holdings.get_holdings_history
        self.get_instrument_pnl_history = self._holdings.get_instrument_pnl_history  # NOQA
        self.get_pnl_history = self._holdings.get_pnl_history  # NOQA

    def connect_market_data(self):
        """
        Initialize MarketData class, should be called before calling trade()
        """
        self._mdata = marketdata.MarketData(prices=self._prices,
                                            rates=self._rates)

    @property
    def event_log(self):
        """
        Returns the event log of events which have acted on the Blotter
        """
        return self._event_log

    def define_generic(self, generic, ccy=None, margin=0, multiplier=1,
                       commission=0, isFX=False):
        """
        Define meta data for a tradeable instruments associated with a generic.

        Parameters
        ----------
        generic: str
            Name for the instrument type used for looking up meta data, e.g. we
            would define 'CL' and the associated meta data for these type of
            contracts
        ccy: str
            Currency that contract is traded in, default is base currency of
            blotter
        margin: float
            Amount of margin required for contract
        multiplier: int
            The multiplier to multiply the price by to get the notional
            amount of the instrument
        commission: float
            Commission charged for trading the instrument
        isFX: boolean
            Indicate if this instrument is an FX instrument. Affects
            whether cash balances are updated for calculating payable
            interest.

        """
        if ccy is None:
            ccy = self._base_ccy
        self._gnrc_meta[generic] = _metadata(ccy, margin, multiplier,
                                             commission, isFX)

    def map_instrument(self, generic, instrument):
        """
        Define a mapping between tradeable instruments and generics, used for
        looking up meta data on instruments. Note in the case of a single
        instrument such as a currency pair the generic and the instrument
        can be the same value.

        Parameters
        ----------
        generic: str
            Name for the instrument type used for looking up meta data, e.g. we
            would define 'CL' and the associated meta data for these type of
            contracts
        instrument: str
            Tradeable instrument name
        """
        self._instr_map[instrument] = generic

    def trade(self, timestamp, instrument, quantity, price):
        """
        Record an instrument trade in the Blotter. This will also make a
        call to automatic_events to trigger all automatic events up to the time
        of this trade.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of trade
        instrument: str
            Tradeable instrument name
        quantity: int
            Number of instrument traded
        price: float
            Price of trade
        """

        # side effects since trade() also manages time state for when to sweep
        # pnl and charge/pay interest/margin
        self.automatic_events(timestamp)
        self._trade(timestamp, instrument, quantity, price)

    def _trade(self, timestamp, instrument, quantity, price):
        # create and dispatch trade events
        events = self._create_trade(timestamp, instrument, quantity, price)
        self.dispatch_events(events)

    def _create_trade(self, timestamp, instrument, quantity, price):
        # implements trade event logic and updates the cash balances for FX
        # instruments where applicable, returns events for these actions

        # returns empty trade list when 0 quantity is traded, this is done for
        # convenience so user can call trade() method without first validating
        # input is non 0 and alternatively calling automatic_events()
        if quantity == 0:
            return []
        generic = self._instr_map[instrument]
        metadata = self._gnrc_meta[generic]
        com = metadata.commission
        ccy = metadata.ccy
        quantity = quantity * metadata.multiplier

        events = [_Event("TRADE", {"timestamp": timestamp,
                                   "instrument": instrument,
                                   "price": price, "quantity": quantity,
                                   "commission": com, "ccy": ccy})]
        ev_cashs = []
        if metadata.isFX:
            counter_quantity = -quantity * price
            cash_counter = _Event("CASH", {"timestamp": timestamp, "ccy": ccy,
                                           "quantity": counter_quantity})
            base_ccy = instrument[:3]
            cash_base = _Event("CASH", {"timestamp": timestamp,
                                        "ccy": base_ccy, "quantity": quantity})
            ev_cashs = [cash_counter, cash_base]

        events.extend(ev_cashs)
        return events

    def automatic_events(self, timestamp):
        """
        Update the current time of the Blotter, triggering all scheduled events
        between previous clock time and new clock time such as interest
        charges, margin charges, PnL calculations and PnL sweeps. See
        create_events() for more information on the type of events.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time to update clock to and tigger internal events up until
        """

        current_time = self._holdings.timestamp
        # first event so there is nothing automatic that needs to be done
        if current_time is pd.NaT:
            return
        actions = self._get_actions(current_time, timestamp, self._actions)
        for ts, action in actions.iteritems():
            events = self.create_events(ts, action)
            self.dispatch_events(events)

    @staticmethod
    def _get_actions(old_ts, new_ts, action_times):
        # calculates the actions between two datetimes and returns them as
        # ordered pandas.Series, filters out weekends since assumption is
        # nothing happens here. This could be extended to allow more advanced
        # user defined filtering based on things such as holiday calendars
        # action_times is a list of tuples with Timedelta and string for action
        # type

        if not action_times:
            return pd.Series([])

        timestamps = pd.date_range(old_ts, new_ts, normalize=True)
        wknds = (timestamps.dayofweek == 5) + (timestamps.dayofweek == 6)
        timestamps = timestamps[~wknds]
        actions = []
        for ts, ac_type in action_times:
            ac_ts = timestamps + ts
            ac_ts = ac_ts[ac_ts > old_ts]
            ac_ts = ac_ts[ac_ts <= new_ts]
            # this will give an empty DataFrame is ac_ts is an empty
            # DateTimeIndex resulting in no actions as desired
            actions.append(pd.Series(ac_type, index=ac_ts))

        actions = pd.concat(actions, axis=0)
        actions.sort_index(inplace=True)

        return actions

    def create_events(self, timestamp, action):
        """
        Create internal event for updating Holdings class contained within
        Blotter instance. Manages creation of INTEREST, MARGIN, PNL and
        PNL_SWEEP events based on internal Blotter data.

        This method is exposed to allow users greater flexibility in calling
        internal events however by default this is automatically called through
        automatic_events() and best not called unless user understands what
        they are doing.

        MARGIN event charges interest in the base currency based on the margin
        required for the current open positions at a rate equal to the base
        currency interest rate + the margin_charge.

        INTEREST events charges interest on the outstandings cash balances in
        different currencies based on the current interest rates.

        PNL event calculates and saves the PNL based on current market prices
        for all open positions.

        PNL_SWEEP event repatriates closed PnL for non base currencies to the
        base currency based on the current FX rates.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time to create event for
        action: str
            Type of event to create, supports INTEREST, MARGIN, PNL and
            PNL_SWEEP

        Returns
        -------
        A list of events for being dispatched using dispatch_events()

        """
        events = []
        if action is "INTEREST":
            cashs = self._holdings.get_cash_balances()
            if not cashs.empty:
                rates = self._mdata.rates.loc[timestamp, cashs.index]
                rates = self._adjusted_rates(timestamp, rates)
                interests = cashs * rates
                for ccy, qty in interests.iteritems():
                    ev = _Event("INTEREST", {"timestamp": timestamp,
                                             "ccy": ccy,
                                             "quantity": qty})
                    events.append(ev)
        elif action is "MARGIN":
            # calculate total margin charge
            base_hlds_value = np.abs(self.get_holdings_value(timestamp))
            int_rate = self._mdata.rates.loc[timestamp, self._base_ccy]
            mrate = int_rate + self._margin_charge
            mrate = self._adjusted_rates(timestamp, mrate)
            charge = 0
            for instr, value in base_hlds_value.iteritems():
                metadata = self._gnrc_meta[self._instr_map[instr]]
                charge += mrate * metadata.margin * value
            if charge:
                ev = _Event("INTEREST", {"timestamp": timestamp,
                                         "ccy": self._base_ccy,
                                         "quantity": charge})
                events.append(ev)
        elif action is "PNL":
            assets = self._holdings.get_assets()
            if assets:
                prices = self._get_prices(timestamp, assets)
            else:
                prices = pd.Series([])
            ev = _Event("PNL", {"timestamp": timestamp, "prices": prices})
            events.append(ev)
        elif action is "PNL_SWEEP":
            assets = self._holdings.get_assets()
            if assets:
                prices = self._get_prices(timestamp, assets)
            else:
                prices = None
            pnls = self._holdings.get_pnl(timestamp, prices, cache=False)
            pnl_sweep = pnls.loc[:, 'closed pnl']
            for ccy, pnl in pnl_sweep.iteritems():
                if ccy is self._base_ccy:
                    continue
                if pnl != 0:
                    conv_rate = self._get_fx_conversion(timestamp, ccy)
                    base_pnl = pnl * conv_rate
                    ev = _Event("PNL_SWEEP", {"timestamp": timestamp,
                                              "ccy1": ccy, "quantity1": -pnl,
                                              "ccy2": self._base_ccy,
                                              "quantity2": base_pnl})
                    events.append(ev)
        else:
            raise NotImplementedError("Unknown event type")

        return events

    @staticmethod
    def _adjusted_rates(timestamp, interest_rates):
        # adjust rates to actual daily payable amount
        interest_rates = interest_rates / 365
        if timestamp.dayofweek == 4:
            # pay interest for Friday, Saturday, Sunday
            interest_rates = interest_rates * 3
        return interest_rates

    def dispatch_events(self, events):
        """
        Update Blotter._holdings based on event. See create_events() for the
        type of events supported. This method is best not called directly
        unless user understands what is going on.

        Parameters
        ----------
        events: list
            list of _Event to dispatch
        """

        for event in events:
            if event.type == "TRADE":
                self._holdings.record_trade(**event.data)
            elif event.type == "CASH":
                self._holdings.update_cash(**event.data)
            elif event.type == "INTEREST":
                self._holdings.charge_interest(**event.data)
            elif event.type == "PNL":
                self._holdings.get_instrument_pnl(**event.data)
            elif event.type == "PNL_SWEEP":
                self._holdings.sweep_pnl(**event.data)
            else:
                raise NotImplementedError("Unknown event type")

            self._event_log.append(str(event))

    def _get_prices(self, timestamp, instruments):
        prices = []
        for instr in instruments:
            prices.append(self._mdata.prices[instr].loc[timestamp])
        prices = pd.concat(prices, axis=0)
        return prices

    def get_holdings_value(self, timestamp):
        """
        Return pandas.Series of values of holdings converted to Blotter base
        currency sorted by index name. Note that for each currency for which
        instruments are traded in, FX rates must be available for the given
        timestamp in order to convert. E.g. if Blotter base ccy is USD, and an
        instrument traded is in AUD, then AUDUSD must be available in the
        prices data folder.

        Parameters
        ----------
        timestamp: pandas.Timestamp which corresponds to the time for
        marking to market blotter holdings
        """
        if self._holdings.timestamp > timestamp:
            raise ValueError('Must mark to market holdings after'
                             'Holdings.timestamp')

        hlds = self._holdings.get_holdings()
        if not hlds:
            return pd.Series()
        base_hlds_value = []
        for ccy in hlds:
            prices_ccy = self._get_prices(timestamp, hlds[ccy].index)
            conv_rate = self._get_fx_conversion(timestamp, ccy)
            value = hlds[ccy] * prices_ccy * conv_rate
            base_hlds_value.append(value)
        base_hlds_value = pd.concat(base_hlds_value, axis=0)
        base_hlds_value.sort_index(inplace=True)
        return base_hlds_value

    def get_instruments(self):
        """
        Return pandas.Series of the number of instruments held.
        """

        hlds = self._holdings.get_holdings()
        if not hlds:
            return pd.Series()

        instr_nums = []
        for ccy in hlds:
            instr_num = hlds[ccy]
            for ast in instr_num.index:
                gnrc = self._instr_map[ast]
                multiplier = self._gnrc_meta[gnrc].multiplier
                instr_num.loc[ast] = instr_num.loc[ast] / multiplier
            instr_nums.append(instr_num)
        instr_nums = pd.concat(instr_nums, axis=0)
        instr_nums.sort_index(inplace=True)
        return instr_nums

    def _get_fx_conversion(self, timestamp, ccy):
        # return rate to multiply through be to convert given ccy
        # to base cccy
        ccy_pair1 = ccy + self._base_ccy
        ccy_pair2 = self._base_ccy + ccy
        if ccy is self._base_ccy:
            conv_rate = 1
        elif ccy_pair1 in self._mdata.prices:
            conv_rate = self._mdata.prices[ccy_pair1].loc[timestamp]
            conv_rate = conv_rate.values
        elif ccy_pair2 in self._mdata.prices:
            conv_rate = 1 / self._mdata.prices[ccy_pair2].loc[timestamp]
            conv_rate = conv_rate.values
        else:
            raise(KeyError(ccy_pair1, ccy_pair2))

        return float(conv_rate)

    def write_log(self, fp):
        """
        Write log of blotter events to file. This can be used for
        reconstituting blotter.

        Parameters
        ----------
        fp: str
            path to write log to
        """
        with open(fp, 'w') as thefile:
            for line in self._event_log:
                thefile.write("%s\n" % line)

    def read_log(self, fp):
        """
        Reconstitute a Blotter object from an event log. Note that this will
        only replay all the events, meta data and market data sources will
        need to be reloaded as well.

        Parameters
        ----------
        fp: str
            path to read log from
        """
        events = self._create_log_events(fp)
        self.dispatch_events(events)

    @staticmethod
    def _create_log_events(fp):
        events = []
        with open(fp, 'r') as thefile:
            for line in thefile:
                events.append(_Event.fromstring(line))
        return events

    def write_meta(self, fp):
        """
        Write meta data of associated instruments in a Blotter to a file. This
        can be used later to reconstitute a Blotter.

        Parameters
        ----------
        fp: str
            path to write meta data
        """
        # https://stackoverflow.com/questions/483666/python-reverse-invert-a-mapping#485368  # NOQA
        inv_map = {}
        for k, v in self._instr_map.items():
            inv_map[v] = inv_map.get(v, [])
            inv_map[v].append(k)

        for key in inv_map:
            inv_map[key].sort()

        keys = list(self._gnrc_meta.keys())
        keys.sort()
        with open(fp, 'w') as myfile:
            for key in keys:
                meta_data_str = json.dumps(self._gnrc_meta[key]._asdict())
                map_str = '{"' + str(key) + '": ' + json.dumps(inv_map[key]) + '}'  # NOQA
                line = meta_data_str + "|" + map_str + "\n"
                myfile.write(line)

    def read_meta(self, fp):
        """
        Reconstitute the meta data of a Blotter from a file. Reads as input
        files output by write_meta(). File formats should be of the following
        form

        Parameters
        ----------
        fp: str
            Path to file. File should have the following format

           {"ccy": "CAD", "margin": 0.1, "multiplier": 100, "commission": 2.5,"isFX": false}|{"CL": ["CLU16", "CLZ16"]} 
           {"ccy": "CAD", "margin": 0, "multiplier": 1, "commission": 2.5, "isFX": true}|{"USDCAD": ["USDCAD"]}

            ...
        """  # NOQA

        with open(fp, 'r') as thefile:
            for line in thefile:
                meta, mapping = line.split("|")
                meta_dict = json.loads(meta)
                mapping_dict = json.loads(mapping)
                generic = list(mapping_dict.keys())[0]
                meta_dict['generic'] = generic
                self.define_generic(**meta_dict)
                instrs = mapping_dict[generic]
                for instr in instrs:
                    self.map_instrument(generic, instr)


class Holdings():
    """
    The Holdings class is designed to manage holdings data and PnL data. The
    class stores instrument level holdings data on a per currency basis and
    calculates PnL on a per currency basis given instrument prices. The class
    is primarily designed to manage these aspects from within the context
    of the Blotter class however can also provide this functionality stand
    alone.

    The main features of the Holdings class include:

    - Store per currency per instrument holindgs
    - Calculate per currency per instrument PnL
    - Maintain interest payable cash balances per currency
    - Maintain charged/payed interest per currency
    - Provide functionality to sweep PnL from one currency to another
    - Return historical holdings
    - Return historical PnL

    Calculating PnL is done on a as of current holdings basis, there is no
    functionality for looking up historical holdings for calculating historic
    PnL.

    Note: For interest bearing instruments, when users are using the Holdings
    call standalone, users are responsible for calling charge_interest() at
    appropriate intervals and with appropriate interest rates to ensure that
    the PnL calculations are correct. This is handled by the Blotter class.

    All actions on the Holdings class must follow in time sequential order.

    """
    def __init__(self):
        self._position_data_per_ccy = {}
        self._cash = {}
        self._interest = {}
        self._pnl_sweep = {}
        self._pnl_data = {}
        self._timestamp = pd.NaT

    @staticmethod
    def _make_empty_holding():
        holding = namedtuple('holding', ['timestamp', 'trade', 'position',
                                         'avg_pos_price', 'fees',
                                         'avg_sell_price', 'total_sell',
                                         'avg_buy_price', 'total_buy'])
        return holding(array('d'), array('d'), array('d'), array('d'),
                       array('d'), array('d'), array('d'), array('d'),
                       array('d'))

    @staticmethod
    def _make_empty_qty():
        cash = namedtuple('cash', ['timestamp', 'amount'])
        return cash(array('d'), array('d'))

    @staticmethod
    def _make_empty_hist_pnl():
        pnl_hist = namedtuple('hist_pnl', ['time', 'pnl'])
        return pnl_hist([], [])

    @property
    def timestamp(self):
        """
        Returns the current timestamp of the Holdings
        """
        return self._timestamp

    def get_holdings(self):
        """
        Get the current instrument holdings. This includes any multiplier
        associated with the instrument.

        Returns
        -------
        dictionary
            Dictionary with currencies as keys and pandas.Series as values
            where that Series contain the most recent holdings for each of the
            holdings in a given currency
        """

        pos_data = self._position_data_per_ccy
        positions = dict()
        for ccy in pos_data:
            ccy_pos_data = pos_data[ccy]
            idx = list(ccy_pos_data)
            idx.sort()
            h = pd.Series(index=idx)
            for asst in ccy_pos_data:
                h.loc[asst] = ccy_pos_data[asst].position[-1]
            # filter closed positions
            h = h.loc[h != 0]
            if not h.empty:
                positions[ccy] = h

        return positions

    def get_holdings_history(self):
        """
        Get the full history of holdings for each instrument traded.
        This includes any multiplier associated with the instrument.

        Returns
        -------
        dictionary
            Dictionary with currencies as keys and dictionary of pandas.Series
            as values where the keys of the nested dictionary are instrument
            names and the pandas.Series is a timeseries of holdings
        """
        pos_data = self._position_data_per_ccy
        positions = dict()
        for ccy in pos_data:
            ccy_pos_data = pos_data[ccy]
            ccy_positions = dict()
            for asst in ccy_pos_data:
                pos_array = ccy_pos_data[asst]
                ts = self._to_timestamp(pos_array.timestamp)
                pos = pd.Series(pos_array.position, index=ts, copy=True)
                ccy_positions[asst] = pos
            positions[ccy] = ccy_positions

        return positions

    @staticmethod
    def _to_timestamp(array):
        # convert array of floats representing POSIX timestamps to a
        # list of Timestamps
        return [pd.Timestamp.fromtimestamp(i) for i in array]

    def get_assets(self):
        """
        Get the names of instruments held.

        Returns
        -------
        list
            Sorted list of strings of current assets which have holdings
        """
        pos_data = self._position_data_per_ccy
        asts = []
        for ccy in pos_data:
            ccy_pos_data = pos_data[ccy]
            for asst in ccy_pos_data:
                if ccy_pos_data[asst].position[-1] != 0:
                    asts.append(asst)

        asts.sort()
        return asts

    def record_trade(self, timestamp, instrument, price, quantity, commission,
                     ccy):
        """
        Record an instrument trade in Holdings. Trades must be time ordered.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of trade
        instrument: str
            Tradeable instrument name
        price: float
            Price of trade
        quantity: int
            Number of instruments traded. This should include any multiplier
            associated with the contract. E.g. for trading an ES contract,
            this has a multipler of 50, therefore 1 ES contract should
            correspond to a quantity of 50 * 1
        commission: float
            total commission for the trade
        ccy: str
            currency of instrument denomination
        """

        if quantity == 0:
            raise ValueError("Cannot trade 0 quantity of an instrument")

        if np.isnan(quantity):
            raise ValueError("Cannot trade nan quantity of an instrument")

        if quantity > 0:
            price_attr = "avg_buy_price"
            total_attr = "total_buy"
        elif quantity < 0:
            price_attr = "avg_sell_price"
            total_attr = "total_sell"

        if ccy in self._position_data_per_ccy:
            ccy_holdings = self._position_data_per_ccy[ccy]
        else:
            ccy_holdings = {}
            self._position_data_per_ccy[ccy] = ccy_holdings

        if instrument in ccy_holdings:
            holdings = ccy_holdings[instrument]
        else:
            holdings = self._make_empty_holding()
            ccy_holdings[instrument] = holdings

        # deals with first access being non existent
        prev_hldings = self._get_last(holdings, 'position')
        avg_price = self._get_last(holdings, price_attr)
        total = self._get_last(holdings, total_attr)

        if self._timestamp > timestamp:
            raise ValueError('Operations on Holdings must follow in time'
                             ' sequential order')
        holdings.timestamp.append(timestamp.timestamp())
        holdings.position.append(prev_hldings + quantity)
        holdings.trade.append(quantity)
        self._timestamp = timestamp

        fees = self._get_last(holdings, "fees", default=0)
        holdings.fees.append(commission + fees)

        aqntity = np.abs(quantity)
        new_price = (total * avg_price + aqntity * price) / (total + aqntity)
        getattr(holdings, price_attr).append(new_price)
        getattr(holdings, total_attr).append(total + aqntity)

        # when adding to position or flipping position sign update
        # average price
        ADDING = np.sign(quantity) == np.sign(prev_hldings)
        NEW_POS = np.sign(quantity + prev_hldings) not in {np.sign(prev_hldings), 0}  # NOQA
        if ADDING:
            a_price = holdings.avg_pos_price[-1]
            new_pos_price = (a_price * prev_hldings + price * quantity) / (prev_hldings + quantity)  # NOQA
            holdings.avg_pos_price.append(new_pos_price)
        elif NEW_POS:
            holdings.avg_pos_price.append(price)
        else:
            holdings.avg_pos_price.append(holdings.avg_pos_price[-1])

    def _get_last(self, obj, attr, default=0):
        try:
            value = getattr(obj, attr)[-1]
        except IndexError:
            value = default
        return value

    def update_cash(self, timestamp, ccy, quantity):
        """
        Update the amount of cash in a certain type of currency, used for
        charging interest on that balance.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of trade
        ccy: str
            currency of cash balance
        quantity: float
            Amount of cash
        """
        self._update_property(timestamp, ccy, quantity, '_cash')

    def charge_interest(self, timestamp, ccy, quantity):
        """
        Update the amount of interest charged in the account of a currency.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of trade
        ccy: str
            currency of interest charge/payment
        quantity: float
            Amount of interest
        """
        self._update_property(timestamp, ccy, quantity, '_interest')

    def sweep_pnl(self, timestamp, ccy1, quantity1, ccy2, quantity2):
        """
        Convert PnL from one currency to another. The user is
        responsible for ensuring that the implicit FX rates used are sensible.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of trade
        ccy1: str
            currency of first leg of sweep
        quantity1: float
            Amount of currency from first leg of sweep
        ccy2: str
            currency of second leg of sweep
        quantity2: float
            Amount of currency from second leg of sweep

        Examples
        --------
        >>> ts = pd.Timestamp('2016-12-01T10:00:00')
        aud = 5000
        usd = 5000 * 0.80
        holder.sweep_pnl(ts, 'AUD', -aud, 'USD', usd)
        """

        self._update_property(timestamp, ccy1, quantity1, '_pnl_sweep')
        self._update_property(timestamp, ccy2, quantity2, '_pnl_sweep')

    def _update_property(self, timestamp, ccy, quantity, attr):
        if self._timestamp > timestamp:
            raise ValueError('Operations on Holdings must follow in time'
                             ' sequential order')

        attr_dict = getattr(self, attr)
        if ccy in attr_dict:
            field = attr_dict[ccy]
        else:
            field = self._make_empty_qty()
            attr_dict[ccy] = field

        prev_amnt = self._get_last(field, "amount", default=0)
        self._timestamp = timestamp
        field.amount.append(prev_amnt + quantity)
        field.timestamp.append(timestamp.timestamp())

    def get_cash_balances(self):
        """
        Return a pandas.Series of the cash balances for each currency
        """

        currencies = list(self._cash)
        currencies.sort()
        cashs = pd.Series(index=currencies)
        for ccy in self._cash:
            cashs.loc[ccy] = self._cash[ccy].amount[-1]
        cashs = cashs.loc[cashs != 0]
        return cashs

    def get_instrument_pnl(self, timestamp, prices=None, cache=True):
        """
        Calculate and return pnl, closed pnl and open pnl for traded
        instruments in each currency.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of PnL calculation, used for caching the result
        prices: pandas.Series
            series of instrument prices for current holdings
        cache: boolean
            Cache this result for later retrieval and advance internal Holdings
            event clock

        Returns
        -------
        dictionary
            Dictionary with currencies as keys and pandas.DataFrame as values
            where the DataFrame contains columns
            ['pnl', 'closed pnl', 'open pnl'] and the index is the set of
            holdings of current instruments
        """

        # allows PnL calculation without having to pass dummy series of prices
        # when all positions are closed
        if prices is None:
            prices = pd.Series()

        if self._timestamp > timestamp:
            raise ValueError('Operations on Holdings must follow in time'
                             ' sequential order')

        pos_data = self._position_data_per_ccy
        pnls = dict()
        for ccy in pos_data:
            ccy_pos_data = pos_data[ccy]
            asts = list(ccy_pos_data)
            asts.sort()
            pos = pd.Series(index=asts)
            fees = pd.Series(index=asts)
            avg_buy_price = pd.Series(index=asts)
            tot_buy = pd.Series(index=asts)
            avg_sell_price = pd.Series(index=asts)
            tot_sell = pd.Series(index=asts)
            avg_pos_price = pd.Series(index=asts)
            for asst in ccy_pos_data:
                ast_dat = ccy_pos_data[asst]
                pos.loc[asst] = self._get_last(ast_dat, 'position')
                fees.loc[asst] = self._get_last(ast_dat, 'fees')
                avg_buy_price[asst] = self._get_last(ast_dat, 'avg_buy_price')
                tot_buy[asst] = self._get_last(ast_dat, 'total_buy')
                avg_sell_price[asst] = self._get_last(ast_dat,
                                                      'avg_sell_price')
                tot_sell[asst] = self._get_last(ast_dat, 'total_sell')
                avg_pos_price[asst] = self._get_last(ast_dat, 'avg_pos_price')

            # this is required to avoid needing to pass in prices for
            # instruments with 0 current holdings but holdings historically
            asts_not0 = pos.loc[pos != 0].index
            prices_ccy = prices.loc[asts_not0]

            if len(asts_not0) == 0:
                pos_value = 0.0
                ccy_open_pnl = pd.Series(0.0, index=asts)
            else:
                pos_value = pos.loc[asts_not0].mul(prices_ccy)
                ccy_open_pnl = pos.loc[asts_not0].mul(prices_ccy - avg_pos_price.loc[asts_not0])  # NOQA    

            ccy_pnl = tot_sell * avg_sell_price + pos_value - avg_buy_price * tot_buy - fees  # NOQA
            ccy_closed_pnl = ccy_pnl - ccy_open_pnl
            df_pnl = pd.concat([ccy_pnl, ccy_closed_pnl, ccy_open_pnl], axis=1)
            df_pnl.columns = ['pnl', 'closed pnl', 'open pnl']
            pnls[ccy] = df_pnl

        if cache:
            for ccy in pnls:
                instr_pnls = pnls[ccy]
                for instr in instr_pnls.index:
                    instr_pnl = instr_pnls.loc[instr, :].tolist()
                    if ccy in self._pnl_data:
                        ccy_pnl_datas = self._pnl_data[ccy]
                    else:
                        ccy_pnl_datas = {}
                        self._pnl_data[ccy] = ccy_pnl_datas
                    if instr in ccy_pnl_datas:
                        instr_pnl_data = ccy_pnl_datas[instr]
                    else:
                        instr_pnl_data = self._make_empty_hist_pnl()
                        ccy_pnl_datas[instr] = instr_pnl_data
                    instr_pnl_data.time.append(timestamp)
                    instr_pnl_data.pnl.append(instr_pnl)
            self._timestamp = timestamp

        return pnls

    def get_pnl(self, timestamp, prices=None, cache=True):
        """
        Calculate open, closed and total pnl in each currency where instruments
        are traded based on given prices.

        Parameters
        ----------
        timestamp: pandas.Timestamp
            Time of PnL calculation
        prices: pandas.Series
            series of instrument prices
        cache: boolean
            Cache this result for later retrieval and advance internal Holdings
            event clock

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns ['pnl', 'closed pnl', 'open pnl'] and an
            index of currencies of instrument denominations. Note that this
            will return a row for each currency that an instrument has ever
            been traded in, even if the current PnL in the currency is all
            0's due to sweeps.
        """

        # allows PnL calculation without having to pass dummy series of prices
        # when all positions are closed
        if prices is None:
            prices = pd.Series()

        pnls = self.get_instrument_pnl(timestamp, prices, cache)

        ccys = list(set().union(pnls, self._interest, self._pnl_sweep))
        ccys.sort()
        ccy_pnls = pd.DataFrame(index=ccys,
                                columns=['pnl', 'closed pnl', 'open pnl'],
                                dtype='float64')
        for ccy in ccys:
            try:
                pnl_sums = pnls[ccy].sum()
            except KeyError:
                pnl_sums = pd.Series(0, index=['pnl', 'closed pnl',
                                               'open pnl'])
            if ccy in self._interest:
                interest = self._get_last(self._interest[ccy], 'amount')
            else:
                interest = 0
            if ccy in self._pnl_sweep:
                swept_pnl = self._get_last(self._pnl_sweep[ccy], 'amount')
            else:
                swept_pnl = 0
            pnl_sums.loc['pnl'] = pnl_sums.loc['pnl'] + interest + swept_pnl
            pnl_sums.loc['closed pnl'] = (pnl_sums.loc['closed pnl'] +
                                          interest + swept_pnl)
            ccy_pnls.loc[ccy] = pnl_sums

        return ccy_pnls

    def get_pnl_history(self):
        """
        Return open, closed and total PnL in each currency where instruments
        are traded based on cached values from previous calls to
        get_instrument_pnl

        Returns
        -------
        dictionary
            Dictionary of pandas.DataFrames where keys are currencies and the
            DataFrames have columns ['pnl', 'closed pnl', 'open pnl'] and
            index of timestamps
        """

        ccy_pnls = self.get_instrument_pnl_history()
        ccys = list(set().union(ccy_pnls, self._interest, self._pnl_sweep))
        ccys.sort()
        hist_pnls = dict()
        PNL_COLS = ['pnl', 'closed pnl', 'open pnl']

        def reindex(df, index):
            df = df.reindex(index, method='ffill')
            df = df.fillna(value=0)
            return df

        for ccy in ccys:
            try:
                instr_pnls = ccy_pnls[ccy]
                instr_idx = pd.DatetimeIndex([])
                instrs = list(instr_pnls.keys())
                instrs.sort()
                for instr in instrs:
                    instr_idx = instr_idx.union(instr_pnls[instr].index)

                instr_pnl_sum = reindex(instr_pnls[instrs[0]], instr_idx)
                for instr in instrs[1:]:
                    pnl = reindex(instr_pnls[instr], instr_idx)
                    instr_pnl_sum = instr_pnl_sum + pnl
            except KeyError:
                instr_pnl_sum = pd.DataFrame([], columns=PNL_COLS)

            try:
                interest_data = self._interest[ccy]
                dts = self._to_timestamp(interest_data.timestamp)
                interest = pd.DataFrame(0, index=dts, columns=PNL_COLS)
                interest.loc[:, 'closed pnl'] = interest_data.amount
                interest.loc[:, 'pnl'] = interest_data.amount
                interest = interest.groupby(interest.index).last()
            except KeyError:
                interest = pd.DataFrame([], columns=PNL_COLS)

            try:
                sweep_data = self._pnl_sweep[ccy]
                dts = self._to_timestamp(sweep_data.timestamp)
                sweep = pd.DataFrame(0, index=dts, columns=PNL_COLS)
                sweep.loc[:, 'closed pnl'] = sweep_data.amount
                sweep.loc[:, 'pnl'] = sweep_data.amount
                # multiple sweeps can happen at same time which all build on
                # each other so only last one is relevant
                sweep = sweep.groupby(sweep.index).last()
            except KeyError:
                sweep = pd.DataFrame([], columns=PNL_COLS)

            idx = instr_pnl_sum.index.union(interest.index).union(sweep.index)
            pnl_ccy = (reindex(instr_pnl_sum, idx) + reindex(sweep, idx) +
                       reindex(interest, idx))
            hist_pnls[ccy] = pnl_ccy

        return hist_pnls

    def get_instrument_pnl_history(self):
        """
        Return open, closed and total PnL in each currency for each traded
        instrumentbased on cached values from previous calls to
        get_instrument_pnl

        Returns
        -------
        dictionary
            Dictionary of dictionaries where to top level dictionary contains
            keys for each currency where there has been PnL historically and
            the nested dictionaries contain keys for each instrument and values
            which are pandas.DataFrame with columns
            ['pnl', 'closed pnl', 'open pnl'] and index of timestamps
        """

        pnl_data = self._pnl_data
        hist_pnl = dict()
        for ccy in pnl_data:
            pnl_data_ccy = pnl_data[ccy]
            hist_pnl_ccy = dict()
            for instr in pnl_data_ccy:
                ts = pnl_data_ccy[instr].time
                instr_pnl = pnl_data_ccy[instr].pnl
                instr_pnl = pd.DataFrame(instr_pnl, index=ts,
                                         columns=['pnl', 'closed pnl',
                                                  'open pnl'])
                hist_pnl_ccy[instr] = instr_pnl
            hist_pnl[ccy] = hist_pnl_ccy

        return hist_pnl
