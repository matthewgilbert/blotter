# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import pandas as pd
import numpy as np
import os
from blotter import blotter


class TimeSuite:
    """
    Time basic functionality of Blotter
    """
    def setup(self):
        # TODO handle these paths in a better way

        ts = pd.date_range('2008-01-02', '2008-05-30', freq='b')
        ts = ts + pd.Timedelta("16h")
        self.ts_daily = ts
        cdir = os.path.dirname(__file__)
        updir = os.path.split(cdir)[0]
        rate_fp = os.path.join(updir, 'doc/data/daily_interest_rates.csv')
        prices_fp = os.path.join(updir, 'doc/data/prices')
        aud_fp = os.path.join(updir, 'doc/data/prices/AUDUSD.csv')
        aud = pd.read_csv(aud_fp, index_col=0, parse_dates=True)
        jpy_fp = os.path.join(updir, 'doc/data/prices/USDJPY.csv')
        jpy = pd.read_csv(jpy_fp, index_col=0, parse_dates=True)
        self.prices = pd.concat([aud, jpy], axis=1)

        self.blt = blotter.Blotter(prices=prices_fp, interest_rates=rate_fp)
        self.blt.connect_market_data()
        self.blt.define_generic("AUDUSD", ccy="USD", margin=0, multiplier=1,
                                   commission=2.5, isFX=True)
        self.blt.define_generic("USDJPY", ccy="JPY", margin=0, multiplier=1,
                                   commission=2.5, isFX=True)
        self.blt.map_instrument("AUDUSD", "AUDUSD")
        self.blt.map_instrument("USDJPY", "USDJPY")

        idx = pd.date_range("2008-01-02T00:00:00", "2008-01-02T23:59:59",
                            freq="5T")
        self.intraday_prices = pd.DataFrame(0.8, index=idx, columns=["AUDUSD"])

    def time_one_trade(self):        
        ts = pd.Timestamp("2008-01-02T16:00:00")
        aud_price = float(self.prices.loc[ts, "AUDUSD"])
        self.blt.trade(ts, "AUDUSD", 100000, aud_price)

    def time_6_months_2_asset_daily(self):
        for ts in self.ts_daily:
            aud_price = float(self.prices.loc[ts, "AUDUSD"])
            jpy_price = float(self.prices.loc[ts, "USDJPY"])
            self.blt.trade(ts, "AUDUSD", 100000, aud_price)
            self.blt.trade(ts, "USDJPY", 100000, jpy_price)

    def time_trade_5_min_intraday(self):
        for ts, p in self.intraday_prices.iterrows():
            self.blt.trade(ts, "AUDUSD", 100000, float(p))


#class MemSuite:
#    def mem_list(self):
#        return [0] * 256
#