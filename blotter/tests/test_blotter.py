import unittest
import os
import tempfile
from collections import namedtuple
from blotter import blotter
from pandas.util.testing import assert_frame_equal, assert_series_equal, \
    assert_dict_equal
import pandas as pd
import numpy as np


class TestBlotter(unittest.TestCase):

    def setUp(self):
        cdir = os.path.dirname(__file__)
        self.prices = os.path.join(cdir, 'data/prices')
        self.rates = os.path.join(cdir, 'data/rates/daily_interest_rates.csv')
        self.log = os.path.join(cdir, 'data/events.log')
        self.meta_log = os.path.join(cdir, 'data/meta_data.log')

    def tearDown(self):
        pass

    def assertEventsEqual(self, evs1, evs2):
        if len(evs1) != len(evs2):
            raise(ValueError("Event lists length mismatch"))
        for ev1, ev2 in zip(evs1, evs2):
            self.assertEqual(ev1.type, ev2.type)
            assert_dict_equal(ev1.data, ev2.data)

    def assertEventTypes(self, evs1, evs2):
        msg = "Event lists length mismatch\n\nLeft:\n%s \nRight:\n%s"
        left_msg = ""
        for ev in evs1:
            left_msg += str(ev) + "\n"
        right_msg = ""
        for ev in evs2:
            right_msg += ev.type + "\n"
        msg = msg % (left_msg, right_msg)

        if len(evs1) != len(evs2):
            raise(ValueError(msg))
        for ev1, ev2 in zip(evs1, evs2):
            if ev1.type is not ev2.type:
                raise(ValueError(msg))

    def assertDictDataFrameEqual(self, dict1, dict2):
        self.assertEqual(dict1.keys(), dict2.keys())
        for key in dict1.keys():
            try:
                assert_frame_equal(dict1[key], dict2[key])
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def make_blotter(self):
        blt = blotter.Blotter(self.prices, self.rates)
        return blt

    def test_get_actions(self):
        actions = [(pd.Timedelta("16h"), "PNL"),
                   (pd.Timedelta("16h"), "INTEREST")]
        old_ts = pd.Timestamp("2017-01-04T10:30")
        new_ts = pd.Timestamp("2017-01-06T10:30")
        ac_ts = blotter.Blotter._get_actions(old_ts, new_ts, actions)

        idx = pd.DatetimeIndex([pd.Timestamp("2017-01-04T16:00"),
                                pd.Timestamp("2017-01-04T16:00"),
                                pd.Timestamp("2017-01-05T16:00"),
                                pd.Timestamp("2017-01-05T16:00")])
        ac_ts_ex = pd.Series(["PNL", "INTEREST", "PNL", "INTEREST"], index=idx)
        assert_series_equal(ac_ts, ac_ts_ex)

    def test_get_actions_weekend_filter(self):
        actions = [(pd.Timedelta("16h"), "PNL"),
                   (pd.Timedelta("16h"), "INTEREST")]
        old_ts = pd.Timestamp("2017-01-06T10:30")
        new_ts = pd.Timestamp("2017-01-09T16:30")
        ac_ts = blotter.Blotter._get_actions(old_ts, new_ts, actions)

        idx = pd.DatetimeIndex([pd.Timestamp("2017-01-06T16:00"),
                                pd.Timestamp("2017-01-06T16:00"),
                                pd.Timestamp("2017-01-09T16:00"),
                                pd.Timestamp("2017-01-09T16:00")])
        ac_ts_ex = pd.Series(["PNL", "INTEREST", "PNL", "INTEREST"], index=idx)
        assert_series_equal(ac_ts, ac_ts_ex)

    def test_trade_undefined_instrument(self):
        blt = self.make_blotter()
        ts = pd.Timestamp('2016-12-10T08:30:00')
        instr = 'CLZ6'
        qty = 1
        price = 48.56

        def make_trade():
            blt._trade(ts, instr, qty, price)

        self.assertRaises(KeyError, make_trade)

    def test_get_meta_data(self):
        blt = self.make_blotter()
        blt.define_generic("CL", "USD", 0.1, 100, 2.5, False)

        meta = namedtuple('metadata', ['ccy', 'margin', 'multiplier',
                                       'commission', 'isFX'])
        metadata_exp = meta("USD", 0.1, 100, 2.5, False)
        metadata = blt._gnrc_meta["CL"]
        self.assertEqual(metadata, metadata_exp)

    def test_get_holdings_empty(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        hlds = blt.get_holdings_value(ts)
        assert_series_equal(hlds, pd.Series())

    def test_get_holdings_timestamp_before(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-05T00:00:00')
        instr = 'ESZ15'
        qty = 1
        price = 2081
        blt.define_generic("ES", "USD", 0.1, 100, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, instr, qty, price)
        ts = pd.Timestamp('2015-08-04T00:00:00')

        def get_holdings():
            blt.get_holdings_value(ts)

        self.assertRaises(ValueError, get_holdings)

    def test_get_holdings_base_ccy(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        instr = 'ESZ15'
        qty = 1
        price = 2081
        blt.define_generic("ES", "USD", 0.1, 100, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, instr, qty, price)
        ts = pd.Timestamp('2015-08-05T00:00:00')
        hlds = blt.get_holdings_value(ts)

        hlds_exp = pd.Series([2082.73 * 100], index=['ESZ15'])
        assert_series_equal(hlds, hlds_exp)

    def test_get_holds_AUD_instr_AUDUSD_fxrate(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        instr = 'APZ15'
        qty = 1
        price = 5200
        blt.define_generic("AP", "AUD", 0.1, 1, 2.5)
        blt.map_instrument("AP", "APZ15")
        blt._trade(ts, instr, qty, price)
        ts = pd.Timestamp('2015-08-05T00:00:00')
        hlds = blt.get_holdings_value(ts)

        hlds_exp = pd.Series([5283 * 0.73457], index=['APZ15'])
        assert_series_equal(hlds, hlds_exp)

    def test_get_holds_CAD_instr_USDCAD_fxrate(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        instr = 'SXMZ15'
        qty = 1
        price = 802.52
        blt.define_generic("SXM", "CAD", 0.1, 1, 2.5)
        blt.map_instrument("SXM", "SXMZ15")
        blt._trade(ts, instr, qty, price)
        ts = pd.Timestamp('2015-08-05T00:00:00')
        hlds = blt.get_holdings_value(ts)

        hlds_exp = pd.Series([795.95 / 1.3183], index=['SXMZ15'])
        assert_series_equal(hlds, hlds_exp)

    def test_get_instruments_empty(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        instrs = blt.get_instruments()
        assert_series_equal(instrs, pd.Series())

    def test_get_instruments_multiplier(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        instr = 'ESZ15'
        qty = 1.0
        price = 2081
        blt.define_generic("ES", "USD", 0.1, 100, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, instr, qty, price)
        instrs = blt.get_instruments()

        instrs_exp = pd.Series([qty], index=['ESZ15'])
        assert_series_equal(instrs, instrs_exp)

    def test_get_instruments_two_ccy(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        instr1 = 'ESZ15'
        instr2 = 'CLZ15'
        qty = 1.0
        price = 2081
        blt.define_generic("ES", "USD", 0.1, 100, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt.define_generic("CL", "CAD", 0.1, 1, 2.5)
        blt.map_instrument("CL", "CLZ15")
        blt._trade(ts, instr1, qty, price)
        blt._trade(ts, instr2, qty, price)
        instrs = blt.get_instruments()

        instrs_exp = pd.Series([qty, qty], index=['CLZ15', 'ESZ15'])
        assert_series_equal(instrs, instrs_exp)

    def test_create_interest_event(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-03T00:00:00')
        blt._holdings.update_cash(ts, "AUD", 1000000)
        blt._holdings.update_cash(ts, "JPY", 1000000)
        ts = pd.Timestamp('2015-08-04T00:00:00')
        evs = blt.create_events(ts, "INTEREST")
        irates = pd.read_csv(self.rates, index_col=0, parse_dates=True)
        aud_int = irates.loc[ts, "AUD"] / 365 * 1000000
        jpy_int = irates.loc[ts, "JPY"] / 365 * 1000000

        evs_exp = [blotter._Event("INTEREST", {"timestamp": ts, "ccy": "AUD",
                                               "quantity": aud_int}),
                   blotter._Event("INTEREST", {"timestamp": ts, "ccy": "JPY",
                                               "quantity": jpy_int})]
        self.assertEventsEqual(evs, evs_exp)

    def test_create_interest_weekend_event(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-06T00:00:00')
        blt._holdings.update_cash(ts, "AUD", 1000000)
        blt._holdings.update_cash(ts, "JPY", 1000000)
        ts = pd.Timestamp('2015-08-07T00:00:00')
        evs = blt.create_events(ts, "INTEREST")
        irates = pd.read_csv(self.rates, index_col=0, parse_dates=True)
        aud_int = irates.loc[ts, "AUD"] / 365 * 3 * 1000000
        jpy_int = irates.loc[ts, "JPY"] / 365 * 3 * 1000000

        evs_exp = [blotter._Event("INTEREST", {"timestamp": ts, "ccy": "AUD",
                                               "quantity": aud_int}),
                   blotter._Event("INTEREST", {"timestamp": ts, "ccy": "JPY",
                                               "quantity": jpy_int})]
        self.assertEventsEqual(evs, evs_exp)

    def test_create_margin_event(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD",
                              margin_charge=0.015)
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        qty = 1
        price = 0
        blt.define_generic("SXM", "CAD", 0.1, 1, 2.5)
        blt.map_instrument("SXM", "SXMZ15")
        blt.define_generic("ES", "USD", 0.05, 1, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, 'SXMZ15', qty, price)
        blt._trade(ts, "ESZ15", qty, price)

        ts = pd.Timestamp('2015-08-05T00:00:00')
        ev = blt.create_events(ts, "MARGIN")

        rates = pd.read_csv(self.rates, index_col=0, parse_dates=True)
        es_fp = os.path.join(self.prices, 'ESZ15.csv')
        es = pd.read_csv(es_fp, index_col=0, parse_dates=True)
        sxm_fp = os.path.join(self.prices, 'SXMZ15.csv')
        sxm = pd.read_csv(sxm_fp, index_col=0, parse_dates=True)
        usdcad_fp = os.path.join(self.prices, 'USDCAD.csv')
        usdcad = pd.read_csv(usdcad_fp, index_col=0, parse_dates=True)

        es_notional = es.loc[ts].values * qty * 0.05
        sxm_notional = sxm.loc[ts].values * qty * 0.1 / usdcad.loc[ts].values
        notnl = float(es_notional + sxm_notional)
        quantity = notnl * (rates.loc[ts, "USD"] + 0.015) / 365
        ev_exp = [blotter._Event("INTEREST", {"timestamp": ts, "ccy": "USD",
                                              "quantity": quantity})]
        self.assertEventsEqual(ev, ev_exp)

    def test_create_short_margin_event(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD",
                              margin_charge=0.015)
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        qty = -1
        price = 0
        blt.define_generic("ES", "USD", 0.05, 1, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, "ESZ15", qty, price)

        ts = pd.Timestamp('2015-08-05T00:00:00')
        ev = blt.create_events(ts, "MARGIN")

        rates = pd.read_csv(self.rates, index_col=0, parse_dates=True)
        es_fp = os.path.join(self.prices, 'ESZ15.csv')
        es = pd.read_csv(es_fp, index_col=0, parse_dates=True)

        es_notional = float(es.loc[ts].values * np.abs(qty) * 0.05)
        quantity = es_notional * (rates.loc[ts, "USD"] + 0.015) / 365
        ev_exp = [blotter._Event("INTEREST", {"timestamp": ts, "ccy": "USD",
                                              "quantity": quantity})]
        self.assertEventsEqual(ev, ev_exp)

    def test_create_pnl_event(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        qty = 1
        price = 0
        blt.define_generic("SXM", "CAD", 0.1, 1, 2.5)
        blt.map_instrument("SXM", "SXMZ15")
        blt.define_generic("ES", "USD", 0.05, 1, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, 'SXMZ15', qty, price)
        blt._trade(ts, "ESZ15", qty, price)

        ts = pd.Timestamp('2015-08-05T00:00:00')
        ev = blt.create_events(ts, "PNL")

        es_fp = os.path.join(self.prices, 'ESZ15.csv')
        es = pd.read_csv(es_fp, index_col=0, parse_dates=True)
        sxm_fp = os.path.join(self.prices, 'SXMZ15.csv')
        sxm = pd.read_csv(sxm_fp, index_col=0, parse_dates=True)
        prices = pd.concat([es.loc[ts], sxm.loc[ts]], axis=0)
        ev_exp = [blotter._Event("PNL", {"timestamp": ts, "prices": prices})]
        self.assertEventsEqual(ev, ev_exp)

    def test_closed_position_pnl_event(self):
        blt = self.make_blotter()
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-04T00:00:00')
        qty = 1
        price = 0
        blt.define_generic("ES", "USD", 0.05, 1, 2.5)
        blt.map_instrument("ES", "ESZ15")
        blt._trade(ts, "ESZ15", qty, price)
        ts = pd.Timestamp('2015-08-05T00:00:00')
        blt._trade(ts, "ESZ15", -qty, price)

        ts = pd.Timestamp('2015-08-06T00:00:00')
        ev = blt.create_events(ts, "PNL")
        ev_exp = [blotter._Event("PNL", {"timestamp": ts,
                                         "prices": pd.Series([])})]
        self.assertEventsEqual(ev, ev_exp)

    def test_create_pnl_sweep_event_closed_pnl(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-03T12:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 50.50, 1, 0, "CAD")
        ts = pd.Timestamp('2015-08-03T14:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 51.50, -1, 0, "CAD")
        ts = pd.Timestamp('2015-08-04T00:00:00')
        evs = blt.create_events(ts, "PNL_SWEEP")
        evs_exp = [blotter._Event("PNL_SWEEP", {"timestamp": ts, "ccy1": "CAD",
                                                "quantity1": -1.00,
                                                "ccy2": "USD",
                                                "quantity2": 1/1.3125})]

        self.assertEventsEqual(evs, evs_exp)

    def test_create_pnl_sweep_no_event_open_pnl_only(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD",
                              sweep_time=None,
                              accrual_time=pd.Timedelta("0h"),
                              eod_time=pd.Timedelta("0h"))
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-03T12:00:00')
        pos = 1
        blt.define_generic("SXM", "CAD", 0, 1, 0)
        blt.map_instrument("SXM", "SXMZ15")
        blt.trade(ts, 'SXMZ15', pos, 800)
        ts = pd.Timestamp('2015-08-04T00:00:00')
        blt.automatic_events(ts)
        evs = blt.create_events(ts, "PNL_SWEEP")
        evs_exp = []
        self.assertEventsEqual(evs, evs_exp)

    def test_create_pnl_sweep_no_event_base(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-03T12:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 50.50, 1, 0, "USD")
        ts = pd.Timestamp('2015-08-03T14:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 51.50, -1, 0, "USD")
        ts = pd.Timestamp('2015-08-04T00:00:00')
        evs = blt.create_events(ts, "PNL_SWEEP")
        evs_exp = []
        self.assertEqual(evs, evs_exp)

    def test_create_pnl_sweep_no_event_pnl_already_swept(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        ts = pd.Timestamp('2015-08-03T12:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 50.50, 1, 0, "CAD")
        ts = pd.Timestamp('2015-08-03T14:00:00')
        blt._holdings.record_trade(ts, 'CLZ15', 51.50, -1, 0, "CAD")
        ts = pd.Timestamp('2015-08-04T00:00:00')
        evs = blt.create_events(ts, "PNL_SWEEP")
        blt.dispatch_events(evs)
        ts = pd.Timestamp('2015-08-05T00:00:00')
        evs = blt.create_events(ts, "PNL_SWEEP")
        self.assertEqual(evs, [])

    def test_create_trade_fx_AUDUSD(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        blt.define_generic("AUDUSD", "USD", 0, 1, 0, True)
        blt.map_instrument("AUDUSD", "AUDUSD")
        ts = pd.Timestamp('2015-08-03T12:00:00')
        evs = blt._create_trade(ts, "AUDUSD", 1000, 0.80)

        ev_exp = [blotter._Event("TRADE", {"timestamp": ts,
                                           "instrument": "AUDUSD",
                                           "ccy": "USD", "price": 0.80,
                                           "quantity": 1000, "commission": 0}),
                  blotter._Event("CASH", {"timestamp": ts, "ccy": "USD",
                                          "quantity": -1000 * 0.80}),
                  blotter._Event("CASH", {"timestamp": ts, "ccy": "AUD",
                                          "quantity": 1000})]
        self.assertEventsEqual(evs, ev_exp)

    def test_create_trade_fx_USDCAD(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        blt.define_generic("USDCAD", "CAD", 0, 1, 0, True)
        blt.map_instrument("USDCAD", "USDCAD")
        ts = pd.Timestamp('2015-08-03T12:00:00')
        evs = blt._create_trade(ts, "USDCAD", 1000, 1.31)

        ev_exp = [blotter._Event("TRADE", {"timestamp": ts,
                                           "instrument": "USDCAD",
                                           "ccy": "CAD", "price": 1.31,
                                           "quantity": 1000, "commission": 0}),
                  blotter._Event("CASH", {"timestamp": ts, "ccy": "CAD",
                                          "quantity": -1000 * 1.31}),
                  blotter._Event("CASH", {"timestamp": ts, "ccy": "USD",
                                          "quantity": 1000})]
        self.assertEventsEqual(evs, ev_exp)

    def test_create_trade_future(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        blt.define_generic("ES", "USD", 0, 1, 0, False)
        blt.map_instrument("ES", "ESZ15")
        ts = pd.Timestamp('2015-08-03T12:00:00')
        evs = blt._create_trade(ts, "ESZ15", 1, 1800)

        ev_exp = [blotter._Event("TRADE", {"timestamp": ts,
                                           "instrument": "ESZ15",
                                           "ccy": "USD", "price": 1800,
                                           "quantity": 1, "commission": 0})]

        self.assertEventsEqual(evs, ev_exp)

    def test_create_trade_0_quantity(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.connect_market_data()
        blt.define_generic("ES", "USD", 0, 1, 0, False)
        blt.map_instrument("ES", "ESZ15")
        ts = pd.Timestamp('2015-08-03T12:00:00')
        evs = blt._create_trade(ts, "ESZ15", 0, 1800)
        self.assertEqual(evs, [])

    def test_create_read_log(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        evs = blt._create_log_events(self.log)

        ts1 = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-02T10:00:00')
        exp_evs = [blotter._Event("TRADE", {"timestamp": ts1,
                                            "instrument": "CLZ16",
                                            "ccy": "USD", "price": 53.46,
                                            "quantity": 100,
                                            "commission": 2.50}),
                   blotter._Event("TRADE", {"timestamp": ts2,
                                            "instrument": "CLZ16",
                                            "ccy": "USD", "price": 55.32,
                                            "quantity": 100,
                                            "commission": 2.50})]
        self.assertEventsEqual(evs, exp_evs)

    def test_write_log(self):
        blt = blotter.Blotter(self.prices, self.rates, base_ccy="USD")
        blt.define_generic("CL", "USD", 0.1, 1, 2.50)
        blt.map_instrument("CL", "CLZ16")
        ts1 = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-02T10:00:00')
        blt._trade(ts1, "CLZ16", 100, 53.46)
        blt._trade(ts2, "CLZ16", 100, 55.32)
        tmp_file = tempfile.mktemp()
        blt.write_log(tmp_file)

        with open(tmp_file, 'r') as fp:
            output_str = fp.read()

        with open(self.log, 'r') as fp:
            expected_output_str = fp.read()

        self.assertMultiLineEqual(expected_output_str, output_str)

    def test_automatic_events_future_type_creation(self):
        blt = blotter.Blotter(self.prices, self.rates,
                              accrual_time=pd.Timedelta(0, unit='h'),
                              eod_time=pd.Timedelta(0, unit='h'),
                              sweep_time=pd.Timedelta(0, unit='h'))
        blt.connect_market_data()
        blt.define_generic("ES", "USD", 0.1, 100, 2.50)
        blt.map_instrument("ES", "ESZ15")
        ts = pd.Timestamp('2015-08-04T10:00:00')
        number_instr = 1
        blt._trade(ts, "ESZ15", number_instr, 2000)
        blt.automatic_events(pd.Timestamp('2015-08-05T10:00:00'))

        ev_types = []
        for ev in blt.event_log:
            ev_types.append(ev.split("|")[0])

        ev_types_exp = ["TRADE", "INTEREST", "PNL"]
        self.assertEqual(ev_types, ev_types_exp)

    def test_automatic_events_fx_type_creation(self):
        blt = blotter.Blotter(self.prices, self.rates,
                              accrual_time=pd.Timedelta(0, unit='h'),
                              eod_time=pd.Timedelta(0, unit='h'),
                              sweep_time=pd.Timedelta(0, unit='h'))
        blt.connect_market_data()

        blt.define_generic("AUDUSD", "USD", 0, 1, 2.50, isFX=True)
        blt.map_instrument("AUDUSD", "AUDUSD")
        ts = pd.Timestamp('2015-08-04T10:00:00')
        number_instr = 1000000
        blt._trade(ts, "AUDUSD", number_instr, 0.80)
        blt.automatic_events(pd.Timestamp('2015-08-05T10:00:00'))

        ev_types = []
        for ev in blt.event_log:
            ev_types.append(ev.split("|")[0])

        ev_types_exp = ["TRADE", "CASH", "CASH", "INTEREST", "INTEREST",
                        "PNL", "PNL_SWEEP"]
        self.assertEqual(ev_types, ev_types_exp)

    def test_automatic_events_closed_pnl_mark(self):
        blt = blotter.Blotter(self.prices, self.rates,
                              accrual_time=pd.Timedelta(0, unit='h'),
                              eod_time=pd.Timedelta(0, unit='h'),
                              sweep_time=pd.Timedelta(0, unit='h'))
        blt.connect_market_data()
        blt.define_generic("ES", ccy="USD", margin=0, multiplier=1,
                           commission=0, isFX=False)
        blt.map_instrument("ES", "ESZ15")

        ts = pd.Timestamp("2015-08-04T11:00:00")
        blt.trade(ts, "ESZ15", 1, 2000)

        ts = pd.Timestamp("2015-08-04T12:00:00")
        hlds = blt.get_instruments()
        for instr, qty in hlds.iteritems():
            blt.trade(ts, instr, -qty, 2001)

        ts = pd.Timestamp("2015-08-05T00:00:00")
        blt.automatic_events(ts)
        pnl_history = blt.get_pnl_history()
        usd = pd.DataFrame([[1.0, 1.0, 0.0]], index=[ts],
                           columns=["pnl", "closed pnl", "open pnl"])
        pnl_history_exp = {"USD": usd}
        self.assertDictDataFrameEqual(pnl_history, pnl_history_exp)

    def test_empty_automatic_events(self):
        blt = blotter.Blotter(self.prices, self.rates,
                              accrual_time=None, eod_time=None,
                              sweep_time=None)
        blt.connect_market_data()
        blt.define_generic("ES", ccy="USD", margin=0, multiplier=1,
                           commission=0, isFX=False)
        blt.map_instrument("ES", "ESZ15")

        ts = pd.Timestamp("2015-08-04T11:00:00")
        blt.trade(ts, "ESZ15", 1, 2000)

        ts = pd.Timestamp("2015-08-05T00:00:00")
        blt.automatic_events(ts)

    def test_read_meta(self):
        blt = blotter.Blotter()
        blt.read_meta(self.meta_log)
        blt.define_generic(**{'generic': 'CL', 'ccy': 'CAD', 'margin': 0.1,
                              'multiplier': 100, 'commission': 2.5,
                              'isFX': False})
        blt.map_instrument('CL', 'CLU16')
        blt.map_instrument('CL', 'CLZ16')

        blt.define_generic(**{'generic': 'USDCAD', 'ccy': 'CAD',
                              'margin': 0,  'multiplier': 1,
                              'commission': 2.5, 'isFX': True})
        blt.map_instrument('USDCAD', 'USDCAD')

        gen_meta_exp = {}
        gen_meta_exp["CL"] = blotter._metadata(**{'ccy': 'CAD', 'margin': 0.1,
                                                  'multiplier': 100,
                                                  'commission': 2.5,
                                                  'isFX': False})
        gen_meta_exp["USDCAD"] = blotter._metadata(**{'ccy': 'CAD',
                                                      'margin': 0,
                                                      'multiplier': 1,
                                                      'commission': 2.5,
                                                      'isFX': True})

        map_instr_exp = {"USDCAD": "USDCAD", "CLU16": "CL", "CLZ16": "CL"}
        self.assertEqual(blt._gnrc_meta, gen_meta_exp)
        self.assertEqual(blt._instr_map, map_instr_exp)

    def test_write_meta(self):
        blt = blotter.Blotter()

        blt.define_generic(**{'generic': 'CL', 'ccy': 'CAD', 'margin': 0.1,
                              'multiplier': 100, 'commission': 2.5,
                              'isFX': False})
        blt.map_instrument('CL', 'CLU16')
        blt.map_instrument('CL', 'CLZ16')

        blt.define_generic(**{'generic': 'USDCAD', 'ccy': 'CAD',
                              'margin': 0,  'multiplier': 1,
                              'commission': 2.5, 'isFX': True})
        blt.map_instrument('USDCAD', 'USDCAD')

        tmp_file = tempfile.mktemp()
        blt.write_meta(tmp_file)

        with open(tmp_file, 'r') as fp:
            output_str = fp.read()

        with open(self.meta_log, 'r') as fp:
            expected_output_str = fp.read()

        self.assertMultiLineEqual(expected_output_str, output_str)
