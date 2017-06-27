.. toctree::
   :maxdepth: 2

Tutorial on using blotter
=========================

This tutorial provides some simple use cases for ``blotter``. To start with,
import the library and create a ``Blotter()`` object

.. ipython:: python

    from blotter import blotter
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    pd.options.display.max_rows = 10
    blt = blotter.Blotter(
               prices="data/prices",
               interest_rates="data/daily_interest_rates.csv",
               accrual_time=pd.Timedelta(16, unit="h"),
               eod_time=pd.Timedelta(16, unit="h"),
               sweep_time=pd.Timedelta(16, unit="h"),
               base_ccy="USD",
               margin_charge=0.015
              )
    blt.connect_market_data()

In the above we provide a path to a folder of csv files for the prices of the
instruments we intend to trade as well as a file containing interest rates
for each of the currencies transacted in. The ``accrual_time`` corresponds to
the time of day which interest and margin are accrued. The interest rates
provided must include an entry for this time of day. ``eod_time`` is the time
at which daily PnL is calculated. The csv files of prices provided must
included a price for this time of day for each day when an instrument has an
open position excluding weekends. Note that this includes market holidays where
there may be no trading, thus the user is responsible for appropriately filling
the data. ``sweep_time`` defines the time of day when closed PnL in
currencies other than the base currency are repatriated to the base currency.
Daily FX rates corresponding to the ``sweep_time`` must be provided in the csv
files of prices for each currency in the non base currency, used for conversion
to the base.

Now we define the meta data for the instruments we plan to trade.

.. ipython:: python

    blt.define_generic("AUDUSD", ccy="USD", margin=0, multiplier=1,
                           commission=2.5, isFX=True)
    blt.define_generic("USDJPY", ccy="JPY", margin=0, multiplier=1,
                           commission=2.5, isFX=True)
    blt.map_instrument(generic="AUDUSD", instrument="AUDUSD")
    blt.map_instrument("USDJPY", "USDJPY")
    blt.define_generic("CL", ccy="USD", margin=0.1, multiplier=1000,
                          commission=2.5, isFX=False)
    blt.map_instrument("CL", "CLZ2008")

We use a vary simple trading strategy which buys and holds Crude and which
uses 5 day momentum for deciding whether to buy or sell AUDUSD and USDJPY.
For simplicity, we will also trade at the 4:00 P.M. and use the same prices as
used internally in the Blotter.

.. ipython:: python

    crude = pd.read_csv("data/prices/CLZ2008.csv", parse_dates=True, index_col=0)
    aud = pd.read_csv("data/prices/AUDUSD.csv", parse_dates=True, index_col=0)
    jpy = pd.read_csv("data/prices/USDJPY.csv", parse_dates=True, index_col=0)

.. ipython:: python

    timestamps = pd.date_range("2008-01-02", "2008-05-27", freq="b")
    timestamps = timestamps + pd.Timedelta("16h")
    timestamps = (timestamps.intersection(crude.index)
                  .intersection(aud.index).intersection(jpy.index))
    ts = timestamps[0]
    blt.trade(ts, "CLZ2008", 10, float(crude.loc[ts]))
    aud_alpha = np.sign(np.log(aud / aud.shift(5))).fillna(value=0)
    jpy_alpha = np.sign(np.log(jpy / jpy.shift(5))).fillna(value=0)


    signal = pd.concat([aud_alpha, jpy_alpha], axis=1)
    signal.columns = ["AUDUSD", "USDJPY"]
    @savefig alpha.png
    signal.plot(title="Alphas")

.. ipython:: python

    def calc_trade(position, alpha):
        if position.empty:
            return alpha
        else:
            return alpha*1000000 - position

.. ipython:: python

    for ts in timestamps:
        pos = blt.get_instruments()
        pos = pos.drop("CLZ2008")
        trds = calc_trade(pos, signal.loc[ts,:])
        aud_qty = float(trds.loc["AUDUSD"])
        aud_price = float(aud.loc[ts])
        jpy_qty = float(trds.loc["USDJPY"])
        jpy_price = float(jpy.loc[ts])
        blt.trade(ts, "AUDUSD", aud_qty, aud_price)
        blt.trade(ts, "USDJPY", jpy_qty, jpy_price)

If at any given time we want to get the current holdings for each instrument in
units of the instrument we can obtain these using

.. ipython:: python

    blt.get_instruments()

If we want to get the value of the current holindgs denominated in the base
currency, we can use

.. ipython:: python

    ts = pd.Timestamp("2008-05-28T16:00:00")
    blt.get_holdings_value(ts)

Note that in the call above we need to provide a timestamp since we are doing
an FX conversion under the hood. Finally we will close out all our positions
and calculate our closed PnL.

.. ipython:: python

    prices = {"AUDUSD": aud, "USDJPY": jpy, "CLZ2008": crude}
    instrs = blt.get_instruments()
    for instr, qty in instrs.iteritems():
        price = float(prices[instr].loc[ts])
        blt.trade(ts, instr, -qty, price)

If we want to look at a full history of positions we can use

.. ipython:: python

    blt.get_holdings_history()

Lastly we calculate all automatic events (interest charges, PnL sweeps and
margin charges).

.. ipython:: python

    ts = pd.Timestamp("2008-05-29T16:00:00")
    blt.automatic_events(ts)


Get holdings intraday, show errors when getting holdings without prices
Show event log of trades

.. ipython:: python

    pnls = blt.get_pnl_history()
    pnls
    @savefig pnl_plot_USD.png
    pnls["USD"].plot(title="USD PnL")
    @savefig pnl_plot_JPY.png
    pnls["JPY"].plot(title="JPY PnL")

Another helpful thing can be looking at the blotter event log, which contains
a record of all actions performed by the Blotter object on its instance of
Holdings.

.. ipython:: python

    blt.event_log[-10:]

You can also get a DataFrame of trades using

.. ipython:: python

    blt.get_trades()
