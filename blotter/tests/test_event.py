import unittest
import ast
import pandas as pd
from blotter import blotter
from pandas.util.testing import assert_dict_equal


class TestBlotter(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def assertEventEqual(self, ev1, ev2):
        self.assertEqual(ev1.type, ev2.type)
        assert_dict_equal(ev1.data, ev2.data)

    def test_no_timestamp_key(self):
        data = {'somevalue': pd.Timestamp('2016-12-01T10:00:00'),
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'commission': 2.50, 'ccy': 'USD'}

        def make_ev():
            blotter._Event('TRADE', data)

        self.assertRaises(ValueError, make_ev)

    def test_no_timestamp_value(self):
        data = {'timestamp': '2016-12-01T10:00:00',
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'commission': 2.50, 'ccy': 'USD'}

        def make_ev():
            blotter._Event('TRADE', data)

        self.assertRaises(ValueError, make_ev)

    def test_trade(self):
        blt = blotter.Blotter()
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'multiplier': 1, 'commission': 2.50, 'ccy': 'USD'}
        ev = [blotter._Event('TRADE', data)]
        blt.dispatch_events(ev)

    def test_trade_fromstring(self):
        trd_str = ('TRADE|{"timestamp":"2016-12-01 10:00:00",'
                   '"instrument":"CLZ6","price":53.46,"quantity":100,'
                   '"commission":2.50,"ccy":"USD"}')
        ev = blotter._Event.fromstring(trd_str)

        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'commission': 2.50, 'ccy': 'USD'}
        ev_exp = blotter._Event('TRADE', data)
        self.assertEventEqual(ev, ev_exp)

    def test_trade_tostr(self):
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'commission': 2.50, 'ccy': 'USD'}
        trd_str = str(blotter._Event('TRADE', data))

        exp_str = ('TRADE|{"timestamp": "2016-12-01 10:00:00", '
                   '"ccy": "USD", "commission": 2.5, "instrument": "CLZ6", '
                   '"price": 53.46, "quantity": 100}')

        self.assertEqual(trd_str, exp_str)

    def test_tostr_valid_dict(self):
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'instrument': 'CLZ6', 'price': 53.46, 'quantity': 100,
                'commission': 2.50, 'ccy': 'USD'}
        trd_str = str(blotter._Event('TRADE', data))
        ast.literal_eval(trd_str.split('|')[1])

    def test_pnl(self):
        blt = blotter.Blotter()
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'prices': pd.Series([53.46, 52], index=['CLZ6', 'COZ6'])}
        ev = [blotter._Event('PNL', data)]
        blt.dispatch_events(ev)

    def test_pnl_fromstring(self):
        trd_str = ('PNL|{"timestamp":"2016-12-01 10:00:00",'
                   '"prices":{"CLZ6": 53.46, "COZ6": 52}}')
        ev = blotter._Event.fromstring(trd_str)

        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'prices': pd.Series([53.46, 52], index=['CLZ6', 'COZ6'])}
        ev_exp = blotter._Event('PNL', data)
        self.assertEventEqual(ev, ev_exp)

    def test_pnl_tostring(self):
        trd_str = ('PNL|{"timestamp": "2016-12-01 10:00:00", '
                   '"prices": {"CLZ6":53.46,"COZ6":52.0}}')
        ev = blotter._Event.fromstring(trd_str)
        ev_str = str(ev)

        self.assertEqual(ev_str, trd_str)

    def test_interest(self):
        blt = blotter.Blotter()
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'ccy': 'USD', 'quantity': 1000}
        ev = [blotter._Event('INTEREST', data)]
        blt.dispatch_events(ev)

    def test_sweep(self):
        blt = blotter.Blotter()
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'ccy1': 'USD', 'quantity1': 1000,
                'ccy2': 'CAD', 'quantity2': -1300}
        ev = [blotter._Event('PNL_SWEEP', data)]
        blt.dispatch_events(ev)

    def test_cash(self):
        blt = blotter.Blotter()
        data = {'timestamp': pd.Timestamp('2016-12-01T10:00:00'),
                'ccy': 'USD', 'quantity': 1000}
        ev = [blotter._Event('CASH', data)]
        blt.dispatch_events(ev)
