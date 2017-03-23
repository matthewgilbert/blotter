import unittest
import array
from blotter import blotter
from pandas.util.testing import assert_series_equal, assert_frame_equal
import pandas as pd
import numpy as np

PNL_COLS = ['pnl', 'closed pnl', 'open pnl']


class TestHoldings(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def assertDictSeriesEqual(self, dict1, dict2):
        self.assertEqual(dict1.keys(), dict2.keys())
        for key in dict1:
            try:
                assert_series_equal(dict1[key], dict2[key])
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def assertDictFrameEqual(self, dict1, dict2):
        self.assertEqual(dict1.keys(), dict2.keys())
        for key in dict1:
            try:
                assert_frame_equal(dict1[key], dict2[key])
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def assertNestedDictSeriesEqual(self, dict1, dict2):
        self.assertEqual(dict1.keys(), dict2.keys())
        for key in dict1:
            try:
                self.assertDictSeriesEqual(dict1[key], dict2[key])
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def assertNestedDictFrameEqual(self, dict1, dict2):
        self.assertEqual(dict1.keys(), dict2.keys())
        for key in dict1:
            try:
                self.assertDictFrameEqual(dict1[key], dict2[key])
            except AssertionError as e:
                e.args = (("\nfor key %s\n" % key) + e.args[0],)
                raise e

    def test_empty_asts(self):
        holder = blotter.Holdings()
        asts = holder.get_assets()
        self.assertEqual(asts, [])

    def test_empty_holdings(self):
        holder = blotter.Holdings()
        holdings = holder.get_holdings()
        holdings_hst = holder.get_holdings_history()
        self.assertEqual(holdings, {})
        self.assertEqual(holdings_hst, {})

    def test_timestamp_conversion(self):
        ts = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-01T11:00:00')
        timestamps = array.array('d', [ts.timestamp(), ts2.timestamp()])
        expected_ts = [ts, ts2]
        ret_ts = blotter.Holdings._to_timestamp(timestamps)
        self.assertEqual(expected_ts, ret_ts)

    def test_one_trade(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)

        asts = holder.get_assets()
        pos = holder.get_holdings()
        pos_hist = holder.get_holdings_history()

        self.assertEqual(asts, ['CLZ6'])
        self.assertDictSeriesEqual(pos,
                                   {"USD": pd.Series(1.0, index=['CLZ6'])})
        exp_pos_hst = {"USD": {"CLZ6": pd.Series(1.0, index=[ts])}}
        self.assertNestedDictSeriesEqual(pos_hist, exp_pos_hst)

    def test_trade_out_of_order(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T09:00:00')

        def book_out_of_order():
            holder.record_trade(ts, instr, price, quantity, commission, ccy)

        self.assertRaises(ValueError, book_out_of_order)

    def test_trade_0_quantity(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        quantity = 0
        commission = 2.50
        ccy = 'USD'

        def trade_0():
            holder.record_trade(ts, instr, price, quantity, commission, ccy)

        self.assertRaises(ValueError, trade_0)

    def test_trade_nan_quantity(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        quantity = np.NaN
        commission = 2.50
        ccy = 'USD'

        def trade_nan():
            holder.record_trade(ts, instr, price, quantity, commission, ccy)

        self.assertRaises(ValueError, trade_nan)

    def test_two_trades_same_ast(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-01T11:00:00')
        instr = 'CLZ6'
        price = 53.36
        price2 = 53.84
        quantity = 1
        quantity2 = 7
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts2, instr, price2, quantity2, commission, ccy)

        asts = holder.get_assets()
        pos = holder.get_holdings()
        pos_hst = holder.get_holdings_history()

        self.assertEqual(asts, ['CLZ6'])

        self.assertDictSeriesEqual(pos,
                                   {"USD": pd.Series(8.0, index=['CLZ6'])})

        exp_pos_hst = {"USD": {"CLZ6": pd.Series([1.0, 8.0], index=[ts, ts2])}}
        self.assertNestedDictSeriesEqual(pos_hst, exp_pos_hst)

    def test_closed_trade(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-01T11:00:00')
        instr = 'CLZ6'
        price = 53.36
        price2 = 53.84
        quantity = 1
        quantity2 = -1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts2, instr, price2, quantity2, commission, ccy)

        asts = holder.get_assets()
        pos = holder.get_holdings()
        pos_hst = holder.get_holdings_history()

        self.assertEqual(asts, [])
        self.assertDictSeriesEqual(pos, {})
        exp_pos_hst = {"USD": {"CLZ6": pd.Series([1.0, 0.0], index=[ts, ts2])}}
        self.assertNestedDictSeriesEqual(pos_hst, exp_pos_hst)

    def test_two_trades_two_asts_one_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        instr2 = 'COZ6'
        price2 = 53.84
        quantity = 1
        quantity2 = 7
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts2, instr2, price2, quantity2, commission, ccy)

        asts = holder.get_assets()
        pos = holder.get_holdings()
        pos_hst = holder.get_holdings_history()

        self.assertEqual(asts, ['CLZ6', 'COZ6'])
        s = pd.Series([1.0, 7.0], index=['CLZ6', 'COZ6'])
        self.assertDictSeriesEqual(pos, {"USD": s})
        exp_pos_hst = {'USD': {'CLZ6': pd.Series([1.0], index=[ts]),
                               'COZ6': pd.Series([7.0], index=[ts2])}}
        self.assertNestedDictSeriesEqual(pos_hst, exp_pos_hst)

    def test_two_trades_two_asts_two_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        ts2 = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        instr2 = 'COZ6'
        price2 = 53.84
        quantity = 1
        quantity2 = 7
        commission = 2.50
        ccy = 'USD'
        ccy2 = 'CAD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts2, instr2, price2, quantity2, commission, ccy2)

        asts = holder.get_assets()
        pos = holder.get_holdings()
        pos_hst = holder.get_holdings_history()

        self.assertEqual(asts, ['CLZ6', 'COZ6'])
        s1 = pd.Series([1.0], index=['CLZ6'])
        s2 = pd.Series([7.0], index=['COZ6'])
        self.assertDictSeriesEqual(pos, {"USD": s1, "CAD": s2})
        exp_pos_hst = {'USD': {'CLZ6': pd.Series([1.0], index=[ts])},
                       'CAD': {'COZ6': pd.Series([7.0], index=[ts2])}}
        self.assertNestedDictSeriesEqual(pos_hst, exp_pos_hst)

    def test_get_holdings_make_trade(self):
        # this test is necessary to test for a BufferError which arises when
        # assigning the output to get_holdings() and then making a trade.
        # Copying the data when creating a the pd.Series resolves
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 53.36
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)

        asts = holder.get_assets()  # NOQA
        pos = holder.get_holdings()  # NOQA
        pos_hst = holder.get_holdings_history()  # NOQA

        holder.record_trade(ts, instr, price, quantity, commission, ccy)

    def test_get_no_cash(self):
        holder = blotter.Holdings()
        cashs = holder.get_cash_balances()
        assert_series_equal(cashs, pd.Series())

    def test_get_USD_cash(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        holder.update_cash(ts, 'USD', 1000)
        cashs = holder.get_cash_balances()
        assert_series_equal(cashs, pd.Series(1000.0, index=['USD']))

    def test_get_USD_and_CAD_cash(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        holder.update_cash(ts, 'USD', 1000)
        holder.update_cash(ts, 'CAD', 1000)
        cashs = holder.get_cash_balances()
        s = pd.Series([1000.0, 1000.0], index=['CAD', 'USD'])
        assert_series_equal(cashs, s)

    def test_get_USD_closed_cash(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        holder.update_cash(ts, 'USD', 1000)
        ts = pd.Timestamp('2016-12-02T10:00:00')
        holder.update_cash(ts, 'USD', -1000)
        cashs = holder.get_cash_balances()
        s = pd.Series()
        assert_series_equal(cashs, s)

    def test_get_USD_interest(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        holder.charge_interest(ts, 'USD', 1000)
        interest = holder._interest['USD'].amount[-1]
        self.assertEqual(interest, 1000)

    # instrument level PnL tests

    def test_instrument_pnl_no_trades(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr_pnls = holder.get_instrument_pnl(ts)
        self.assertDictEqual(instr_pnls, {})

    def test_instrument_pnl_one_instrument_one_trade(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        pnls = holder.get_instrument_pnl(ts, pd.Series([56], index=[instr]))

        df_pnl = pd.concat([pd.Series([-1.5], index=[instr]),
                            pd.Series([-2.5], index=[instr]),
                            pd.Series([1.0], index=[instr])], axis=1)
        df_pnl.columns = PNL_COLS

        pnls_expected = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_expected)

    def test_instrument_pnl_one_instrument_one_trade_extra_prices(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        # add prices for irrelevant instruments, these should be ignored
        prices = pd.Series([56, 80, 0.80], index=[instr, 'COZ6', 'AUDUSD'])
        pnls = holder.get_instrument_pnl(ts, prices)

        df_pnl = pd.concat([pd.Series([-1.5], index=[instr]),
                            pd.Series([-2.5], index=[instr]),
                            pd.Series([1.0], index=[instr])], axis=1)
        df_pnl.columns = PNL_COLS

        pnls_expected = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_expected)

    def test_instrument_pnl_two_instrument_one_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        instr2 = 'COZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts, instr2, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57, 56], index=[instr, instr2])
        pnls = holder.get_instrument_pnl(ts, prices)

        df_pnl = pd.concat([pd.Series([-0.5, -1.5], index=[instr, instr2]),
                            pd.Series([-2.5, -2.5], index=[instr, instr2]),
                            pd.Series([2.0, 1.0], index=[instr, instr2])],
                           axis=1)
        df_pnl.columns = PNL_COLS
        pnls_exp = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_exp)

    def test_instrument_pnl_two_instrument_two_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        instr2 = 'COZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        ccy2 = 'CAD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        holder.record_trade(ts, instr2, price, quantity, commission, ccy2)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57, 56], index=[instr, instr2])
        pnls = holder.get_instrument_pnl(ts, prices)
        df_pnl1 = pd.concat([pd.Series([-1.5], index=[instr2]),
                            pd.Series([-2.5], index=[instr2]),
                            pd.Series([1.0], index=[instr2])], axis=1)
        df_pnl2 = pd.concat([pd.Series([-0.5], index=[instr]),
                            pd.Series([-2.5], index=[instr]),
                            pd.Series([2.0], index=[instr])], axis=1)

        df_pnl1.columns = PNL_COLS
        df_pnl2.columns = PNL_COLS
        pnls_exp = {'CAD': df_pnl1,
                    'USD': df_pnl2}

        self.assertDictFrameEqual(pnls, pnls_exp)

    def test_instrument_pnl_one_instrument_closed_no_price_needed(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        price = 65
        quantity = -1
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T12:00:00')
        # when there is no position in an instrument the price for that is
        # ignored so doesn't need to be passed
        pnls = holder.get_instrument_pnl(ts, pd.Series([]))

        df_pnl = pd.concat([pd.Series([5.0], index=[instr]),
                            pd.Series([5.0], index=[instr]),
                            pd.Series([0.0], index=[instr])], axis=1)
        df_pnl.columns = PNL_COLS
        pnls_exp = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_exp)

    def test_instrument_pnl_missing_price(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series()

        def no_price():
            holder.get_instrument_pnl(ts, prices)

        self.assertRaises(KeyError, no_price)

    def test_one_instrument_pnl_one_nan_price(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([pd.np.NaN], index=[instr])

        pnls = holder.get_instrument_pnl(ts, prices)

        df_pnl = pd.concat([pd.Series([pd.np.NaN], index=[instr]),
                            pd.Series([pd.np.NaN], index=[instr]),
                            pd.Series([pd.np.NaN], index=[instr])], axis=1)
        df_pnl.columns = PNL_COLS
        pnls_exp = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_exp)

    def test_two_instrument_pnl_one_nan_price(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr1 = 'CLZ6'
        instr2 = 'COZ6'
        price = 55
        quantity = 1
        commission = 0
        ccy = 'USD'
        holder.record_trade(ts, instr1, price, quantity, commission, ccy)
        holder.record_trade(ts, instr2, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([pd.np.NaN, price], index=[instr1, instr2])

        pnls = holder.get_instrument_pnl(ts, prices)

        df_pnl1 = pd.concat([pd.Series([pd.np.NaN], index=[instr1]),
                             pd.Series([pd.np.NaN], index=[instr1]),
                             pd.Series([pd.np.NaN], index=[instr1])], axis=1)
        df_pnl2 = pd.concat([pd.Series([0.0], index=[instr2]),
                             pd.Series([0.0], index=[instr2]),
                             pd.Series([0.0], index=[instr2])], axis=1)
        df_pnl = pd.concat([df_pnl1, df_pnl2], axis=0)
        df_pnl.columns = PNL_COLS
        pnls_exp = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_exp)

    # historical instrument pnl

    def test_instrument_pnl_hist_no_trades(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        prices = pd.Series([])
        holder.get_instrument_pnl(ts, prices)

        instr_pnls = holder.get_instrument_pnl_history()
        self.assertNestedDictFrameEqual(instr_pnls, {})

    def test_instrument_pnl_hist_one_instrument(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57], index=[instr])
        holder.get_instrument_pnl(ts, prices)

        pnls = holder.get_instrument_pnl_history()
        df_pnl = pd.DataFrame([[-0.5, -2.5, 2.0]], index=[ts],
                              columns=PNL_COLS)
        pnls_expected = {'USD': {'CLZ6': df_pnl}}
        self.assertNestedDictFrameEqual(pnls, pnls_expected)

    def test_instrument_pnl_hist_one_instrument_no_cache(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57], index=[instr])
        holder.get_instrument_pnl(ts, prices, cache=False)

        pnls = holder.get_instrument_pnl_history()
        pnls_expected = {}
        self.assertNestedDictFrameEqual(pnls, pnls_expected)

    def test_instrument_pnl_hist_two_instrument_one_ccy_diff_lengths(self):
        holder = blotter.Holdings()
        ts1 = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts1, instr, price, quantity, commission, ccy)

        ts2 = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([56], index=[instr])
        holder.get_instrument_pnl(ts2, prices)

        ts3 = pd.Timestamp('2016-12-01T12:00:00')
        instr2 = 'COZ6'
        price = 50
        holder.record_trade(ts3, instr2, price, quantity, commission, ccy)

        ts4 = pd.Timestamp('2016-12-01T13:00:00')
        prices = pd.Series([57, 49], index=[instr, instr2])
        holder.get_instrument_pnl(ts4, prices)

        pnls = holder.get_instrument_pnl_history()

        pnl_cz = [[-1.5, -2.5, 1.0], [-0.5, -2.5, 2.0]]
        pnl_co = [[-3.5, -2.5, -1.0]]

        df_pnl_cz = pd.DataFrame(pnl_cz, index=[ts2, ts4], columns=PNL_COLS)
        df_pnl_co = pd.DataFrame(pnl_co, index=[ts4], columns=PNL_COLS)

        df_pnl_cz.columns = PNL_COLS
        df_pnl_co.columns = PNL_COLS
        pnls_exp = {'USD': {'CLZ6': df_pnl_cz, 'COZ6': df_pnl_co}}

        self.assertNestedDictFrameEqual(pnls, pnls_exp)

    def test_pnl_calc_pre_trade(self):
        holder = blotter.Holdings()
        ts1 = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts1, instr, price, quantity, commission, ccy)

        ts2 = pd.Timestamp('2016-12-01T09:00:00')
        prices = pd.Series([56], index=[instr])

        def pnl_out_of_order():
            holder.get_instrument_pnl(ts2, prices)

        self.assertRaises(ValueError, pnl_out_of_order)

    def test_instrument_pnl_hist_two_instrument_two_ccy_two_marks(self):
        holder = blotter.Holdings()
        ts1 = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        instr2 = 'COZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        ccy2 = 'CAD'
        holder.record_trade(ts1, instr, price, quantity, commission, ccy)
        holder.record_trade(ts1, instr2, price, quantity, commission, ccy2)

        ts2 = pd.Timestamp('2016-12-01T10:30:00')
        prices = pd.Series([56, 57], index=[instr, instr2])
        holder.get_instrument_pnl(ts2, prices)

        ts3 = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57, 58], index=[instr, instr2])
        holder.get_instrument_pnl(ts3, prices)

        pnls = holder.get_instrument_pnl_history()

        pnl_cz = [[-1.5, -2.5, 1.0], [-0.5, -2.5, 2.0]]
        pnl_co = [[-0.5, -2.5, 2.0], [0.5, -2.5, 3.0]]

        df_pnl_cz = pd.DataFrame(pnl_cz, index=[ts2, ts3], columns=PNL_COLS)
        df_pnl_co = pd.DataFrame(pnl_co, index=[ts2, ts3], columns=PNL_COLS)

        df_pnl_cz.columns = PNL_COLS
        df_pnl_co.columns = PNL_COLS
        pnls_exp = {'USD': {'CLZ6': df_pnl_cz},
                    'CAD': {'COZ6': df_pnl_co}}

        self.assertNestedDictFrameEqual(pnls, pnls_exp)

    # historical aggregate level PnL tests

    def test_pnl_hist_no_trades(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        prices = pd.Series([])
        holder.get_instrument_pnl(ts, prices)

        pnls = holder.get_pnl_history()
        self.assertNestedDictFrameEqual(pnls, {})

    def test_pnl_hist_one_instrument(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57], index=[instr])
        holder.get_instrument_pnl(ts, prices)

        pnls = holder.get_pnl_history()
        df_pnl = pd.DataFrame([[-0.5, -2.5, 2.0]], index=[ts],
                              columns=PNL_COLS)
        pnls_expected = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_expected)

    def test_pnl_hist_two_instrument(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr1 = 'CLZ6'
        price1 = 55
        quantity1 = 1
        instr2 = 'COZ6'
        price2 = 54
        quantity2 = 2
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr1, price1, quantity1, commission, ccy)
        holder.record_trade(ts, instr2, price2, quantity2, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        prices = pd.Series([57, 53], index=[instr1, instr2])
        holder.get_instrument_pnl(ts, prices)

        pnls = holder.get_pnl_history()
        df_pnl = pd.DataFrame([[-5.0, -5.0, 0.0]], index=[ts],
                              columns=PNL_COLS)
        pnls_expected = {'USD': df_pnl}
        self.assertDictFrameEqual(pnls, pnls_expected)

    def test_pnl_hist_two_instr_two_ccy_sweep_interest(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')

        holder.record_trade(ts, 'CLZ6', 55, 1, 2.50, 'USD')
        holder.record_trade(ts, 'USDJPY', 110.00, 1000, 0, 'JPY')

        ts_eod1 = pd.Timestamp('2016-12-01T16:00:00')
        prices = pd.Series([56, 111.00], index=['CLZ6', 'USDJPY'])
        holder.charge_interest(ts_eod1, 'USD', 0.1)
        holder.charge_interest(ts_eod1, 'JPY', -1)
        holder.get_instrument_pnl(ts_eod1, prices)

        ts_eod2 = pd.Timestamp('2016-12-02T16:00:00')
        holder.record_trade(ts_eod2, 'USDJPY', 112.0, -1000, 0, 'JPY')
        usd = ((112 - 110) * 1000 - 1) / 112.0
        jpy = -((112 - 110) * 1000 - 1)
        holder.sweep_pnl(ts_eod2, 'USD', usd, 'JPY', jpy)
        prices = pd.Series([57], index=['CLZ6'])
        holder.get_instrument_pnl(ts_eod2, prices)

        pnls = holder.get_pnl_history()

        closed = 0.1 - 2.5 + usd
        pnl_usd = pd.DataFrame([[-1.4, -2.4, 1.0], [2 + closed, closed, 2.0]],
                               index=[ts_eod1, ts_eod2], columns=PNL_COLS)

        pnl_jpy = pd.DataFrame([[999.0, -1.0, 1000.0], [0.0, 0.0, 0.0]],
                               index=[ts_eod1, ts_eod2], columns=PNL_COLS)

        pnls_expected = {'USD': pnl_usd, 'JPY': pnl_jpy}
        self.assertDictFrameEqual(pnls, pnls_expected)

    def test_pnl_hist_multi_sweep_same_time(self):
        holder = blotter.Holdings()
        ts_eod1 = pd.Timestamp('2016-12-01T16:00:00')
        holder.sweep_pnl(ts_eod1, 'USD', 1.0, 'JPY', -1.0)
        holder.sweep_pnl(ts_eod1, 'USD', 1.0, 'JPY', -1.0)
        ts_eod2 = pd.Timestamp('2016-12-02T16:00:00')
        holder.charge_interest(ts_eod2, 'USD', 0.1)

        pnl_history = holder.get_pnl_history()
        jpy = pd.DataFrame([[-2.0, -2.0, 0]], index=[ts_eod1],
                           columns=PNL_COLS)
        usd = pd.DataFrame([[2.0, 2.0, 0.0], [2.1, 2.1, 0.0]],
                           index=pd.DatetimeIndex([ts_eod1, ts_eod2]),
                           columns=PNL_COLS)
        pnl_history_exp = {"JPY": jpy, "USD": usd}
        self.assertDictFrameEqual(pnl_history, pnl_history_exp)

    def test_pnl_hist_multi_interest_charge_same_time(self):
        holder = blotter.Holdings()
        ts_eod1 = pd.Timestamp('2016-12-01T16:00:00')
        holder.charge_interest(ts_eod1, 'USD', 1.0)
        holder.charge_interest(ts_eod1, 'USD', 5.0)

        ts_eod2 = pd.Timestamp('2016-12-02T16:00:00')
        holder.sweep_pnl(ts_eod2, 'USD', 1.0, 'JPY', -1)

        pnl_history = holder.get_pnl_history()
        jpy = pd.DataFrame([[-1, -1, 0]], index=[ts_eod2],
                           columns=PNL_COLS)
        usd = pd.DataFrame([[6.0, 6.0, 0.0], [7.0, 7.0, 0.0]],
                           index=pd.DatetimeIndex([ts_eod1, ts_eod2]),
                           columns=PNL_COLS)
        pnl_history_exp = {"JPY": jpy,
                           "USD": usd}

        self.assertDictFrameEqual(pnl_history, pnl_history_exp)

    # aggregate level PnL tests
    def test_pnl_no_trades(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        prices = pd.Series([])
        pnl = holder.get_pnl(ts, prices)
        exp_pnl = pd.DataFrame([], columns=PNL_COLS,
                               dtype='float64')
        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_one_instrument_one_trade_one_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        pnl = holder.get_pnl(ts, pd.Series([56], index=[instr]))

        pnl_tot = pd.Series([-1.50], index=['USD'])
        pnl_closed = pd.Series([-2.50], index=['USD'])
        pnl_open = pd.Series([1.0], index=['USD'])
        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_one_instrument_all_closed_one_ccy_dummy_prices(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        price = 65
        quantity = -1
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T12:00:00')
        pnl = holder.get_pnl(ts, pd.Series())

        pnl_tot = pd.Series([5.0], index=['USD'])
        pnl_closed = pd.Series([5.0], index=['USD'])
        pnl_open = pd.Series([0.0], index=['USD'])
        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_one_instrument_all_closed_one_ccy_no_prices(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 1
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        price = 65
        quantity = -1
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T12:00:00')
        pnl = holder.get_pnl(ts)
        pnl_tot = pd.Series([5.0], index=['USD'])
        pnl_closed = pd.Series([5.0], index=['USD'])
        pnl_open = pd.Series([0.0], index=['USD'])
        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_one_instrument_two_trades_one_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'CLZ6'
        price = 55
        quantity = 5
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T11:00:00')
        price = 65
        quantity = -1
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T12:00:00')
        new_price = pd.Series([70], index=[instr])
        pnl = holder.get_pnl(ts, new_price)

        pnl_tot = pd.Series([65.0], index=['USD'])
        pnl_closed = pd.Series([5.0], index=['USD'])
        pnl_open = pd.Series([60.0], index=['USD'])
        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_one_instrument_one_trade_with_interest_one_ccy(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        instr = 'AUDUSD'
        price = 0.80
        quantity = 1000000
        commission = 2.50
        ccy = 'USD'
        holder.record_trade(ts, instr, price, quantity, commission, ccy)
        ts = pd.Timestamp('2016-12-01T16:00:00')
        # accrue interest on long AUD position, pay on short USD position
        holder.charge_interest(ts, 'AUD', 55)
        holder.charge_interest(ts, 'USD', -15)
        new_price = pd.Series([0.81], index=[instr])
        pnl = holder.get_pnl(ts, new_price)

        pnl_tot = pd.Series([55, 10000.0 - 2.50 - 15], index=['AUD', 'USD'])
        pnl_closed = pd.Series([55, -2.5 - 15], index=['AUD', 'USD'])
        pnl_open = pd.Series([0, 10000.0], index=['AUD', 'USD'])

        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_sweep(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        aud = 5000
        usd = 5000 * 0.80
        holder.sweep_pnl(ts, 'AUD', -aud, 'USD', usd)
        # pass dummy prices since never used since no trades
        pnl = holder.get_pnl(ts, pd.Series([]))

        pnl_tot = pd.Series([-5000.0, 4000.0], index=['AUD', 'USD'])
        pnl_closed = pd.Series([-5000.0, 4000.0], index=['AUD', 'USD'])
        pnl_open = pd.Series([0.0, 0.0], index=['AUD', 'USD'])
        exp_pnl = pd.concat([pnl_tot, pnl_closed, pnl_open], axis=1)
        exp_pnl.columns = PNL_COLS

        assert_frame_equal(pnl, exp_pnl)

    def test_pnl_sweep_out_of_order(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')
        aud = 5000
        # just use filler exchange rate
        holder.sweep_pnl(ts, 'AUD', aud, 'USD', -aud)
        ts = pd.Timestamp('2016-12-01T09:00:00')

        def out_of_order_sweep():
            holder.sweep_pnl(ts, 'AUD', aud, 'USD', -aud)

        self.assertRaises(ValueError, out_of_order_sweep)

    def test_pnl_sweep_closed_position(self):
        holder = blotter.Holdings()
        ts = pd.Timestamp('2016-12-01T10:00:00')

        holder.record_trade(ts, 'CLZ6', 55, 1, 2.50, 'USD')
        holder.record_trade(ts, 'USDJPY', 110.00, 1000, 0, 'JPY')

        ts_eod1 = pd.Timestamp('2016-12-01T16:00:00')
        prices = pd.Series([56, 111.00], index=['CLZ6', 'USDJPY'])
        holder.charge_interest(ts_eod1, 'USD', 0.1)
        holder.charge_interest(ts_eod1, 'JPY', -1)
        holder.get_instrument_pnl(ts_eod1, prices)

        ts_eod2 = pd.Timestamp('2016-12-02T16:00:00')
        holder.record_trade(ts_eod2, 'USDJPY', 112.0, -1000, 0, 'JPY')
        usd = ((112 - 110) * 1000 - 1) / 112.0
        jpy = -((112 - 110) * 1000 - 1)
        holder.sweep_pnl(ts_eod2, 'USD', usd, 'JPY', jpy)
        prices = pd.Series([57], index=['CLZ6'])
        pnls = holder.get_pnl(ts_eod2, prices)

        closed = 0.1 - 2.5 + usd
        pnls_exp = pd.DataFrame([[0.0, 0.0, 0.0], [2 + closed, closed, 2.0]],
                                index=['JPY', 'USD'], columns=PNL_COLS)

        assert_frame_equal(pnls, pnls_exp)
