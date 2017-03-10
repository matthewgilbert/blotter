# blotter: a financial blotter for Python

## Description

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

## Install

You can pip install this package from github, i.e.

```
pip install git+git://github.com/matthewgilbert/pdblp.git@master
```

the package is not yet available on Pypi.

## Requires

- `pandas` >= 0.18.1

For building the documentation `sphinx` is also required and performance tests
use `asv`

## Documentation

Online documentation can be found at https://matthewgilbert.github.io/blotter/
The documentation is built using Sphinx and its source is available in `doc/`

## Tests

Unit tests are located in `blotter/tests/` and can be run using `make test`
