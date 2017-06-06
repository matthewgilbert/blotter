import unittest
import os
from blotter import marketdata
from pandas.util.testing import assert_frame_equal
import pandas as pd


class TestMarketData(unittest.TestCase):

    def assertDictOfFrames(self, dict1):
        for key in dict1:
            try:
                assert isinstance(dict1[key], pd.DataFrame)
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def setUp(self):
        cdir = os.path.dirname(__file__)
        self.prices = os.path.join(cdir, 'data/prices')
        self.rates = os.path.join(cdir, 'data/rates/daily_interest_rates.csv')
        self.price_file = os.path.join(cdir, 'data/prices/EURUSD.csv')

    def tearDown(self):
        pass

    def test_read_price_data(self):
        mdata = marketdata.MarketData(prices=self.prices, rates=self.rates)

        exp_keys = {"ESZ15": "", "AUDUSD": "", "EURUSD": "", "USDCAD": "",
                    "SXMZ15": "", "APZ15": "", "USDZAR": ""}.keys()
        self.assertEqual(mdata.prices.keys(), exp_keys)
        self.assertDictOfFrames(mdata.prices)

        assert isinstance(mdata.rates, pd.DataFrame)

    def test_read_file(self):
        df = marketdata.MarketData._read_file(self.price_file)
        df_exp = pd.DataFrame(
            index=pd.date_range("2015-08-04", "2015-08-05"),
            columns=["EURUSD"],
            data=[1.09624, 1.0873]
        )
        df_exp.index.name = 'date'
        assert_frame_equal(df, df_exp)
