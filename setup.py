from setuptools import setup
import re

# https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package#7071358
VERSIONFILE = "blotter/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." %
                       (VERSIONFILE,))

LONG_DESCRIPTION = """
**blotter** is a python object for representing a financial blotter. The main
goal of **blotter** is to provide a decoupled way to manage historical
holdings and track performance when backtesting a trading strategy which trades
FX and Futures. **blotter** is not intended for calculating alpha's, simulating
execution or performing portfolio optimization but rather serves as an endpoint
to book the price of executed trades and provide marked to market holdings and
PnL.

## Main Features

- Handle open and closed PnL calculations across multiple currencies
- Ability to specify base currency for marked to market holdings calculations,
automated PnL sweeps and margin charges
- Accounts for commissions
- Allows for trading instruments denominated in different currencies
- Automatically manage daily interest payments, margin charges, PnL
calculations and PnL sweeps or manage manually
- Get intraday mark to market holdings denominated in base currency
"""

setup(name='blotter',
      version=verstr,
      description='Financial blotter for Futures and FX',
      long_description=LONG_DESCRIPTION,
      url='https://github.com/MatthewGilbert/blotter',
      author='Matthew Gilbert',
      author_email='matthew.gilbert12@gmail.com',
      license='MIT',
      platforms='any',
      install_requires=['pandas>=0.18.0'],
      packages=['blotter', 'blotter.tests'],
      package_data = {'blotter.tests': ['data/*.log', 'data/prices/*.csv',
                                        'data/rates/*.csv']},
      test_suite='blotter.tests',
      zip_safe=False)
