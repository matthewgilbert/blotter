.. blotter documentation master file, created by
   sphinx-quickstart on Sat Dec 24 15:05:39 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

About blotter
=============

The ``blotter`` package aims to provide a mechanism for managing historical
instrument level positions for FX and futures and calculating PnL for
backtests in a manner which is decoupled from upstream trading logic and
downstream rebalance logic.

``blotter`` provides users with an interface to book trades, after which
``blotter`` manages margin and interest charges and allows the user to
obtain marked to market holdings as well as open and closed PnL at future
times. There is also functionality for specifying simple schedules for
sweeping closed PnL in foreign currencies into the base currency.
The two main classes for doing this are the Blotter class and the
Holdings class.

Blotter
*********

The Blotter class manages all the market data for the traded instruments as
well as instrument level meta data which includes commissions, margin
requirements, instrument multiplier and instrument currency. In addition,
the Blotter class allows users to define a time of day for charging interest
based on open positions, charging margin based on open positions,
calculating open and closed PnL and sweeping closed PnL back to the base
currency defined by the user. The Blotter class contains a set of Holdings
which maintains currency level instrument holdings and PnL data. A key thing to
note about the Blotter class is that it does not support trading FX forwards.
The functionality provided for trading FX equivalent to trading spot FX and
rolling positions daily, thus interest rates should be chosen appropriatly
to be representative of this.

After instantiating a Blotter and defining instruments, the main methods the
user uses to simulate a trading strategy are ``trade`` and
``automatic_events``. These methods allow the user to perform actions on the
Blotter in a time seequential manner. ``trade`` allows the user to trade an
instrument and also triggers a call to ``automatic_events``.
``automatic_events`` triggers daily PnL calculations, interest charges,
margin charges and PnL repatriation to the base currency of the Blotter on a
daily time schedule between the previous time state of the Blotter and the new
time state.

Holdings
**********

The Holdings class manages positions and PnL calculations. The class provides
both instrument level PnL calculations and currency level PnL calculations.
Instrument level PnL is given by


.. math::

   \text{Open PnL} &= \text{Position} \cdot (\text{Market Price} - \text{Avg. Transacted Price}) \\
   \text{PnL} &= \text{Total Sold} \cdot \text{Avg. Sell Price} -
   \text{Total Bought} \cdot \text{Avg. Purchase Price} \\
   & + \text{Position} \cdot \text{Market Price} - \text{Total Commission} \\
   \text{Closed PnL} &= \text{PnL} - \text{Open PnL}

Currency level PnL for each currency is given by

.. math::

   \text{Open PnL} &= \sum \text{Open Instrument PnL} \\
   \text{PnL} &= \sum \text{Instrument PnL} + \text{Interest} + \text{PnL Sweeps}\\
   \text{Closed PnL} &= \sum \text{Closed PnL} + \text{Interest} + \text{PnL Sweeps}\\

One thing to note in the above calculations is that interest is attributed at
the currency level, thus at present there is no way to get instrument level PnL
calculations which account for interest.

There is currently no support for margin netting as well as settling
contracts.

There is also a section with a simple tutorial

.. toctree::
   :hidden:

   self

.. toctree::
   :titlesonly:

   tutorial

as well as a section detailing the API

.. toctree::
   :titlesonly:

   api

Installation
============
Please refer to the
`README <https://github.com/matthewgilbert/blotter/blob/master/README.md>`_
for installation instructions


Contributions
=============
Contributions are welcome, the source repository is available
`here <https://github.com/matthewgilbert/blotter>`_
