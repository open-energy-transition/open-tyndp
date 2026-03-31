.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

###########################
Cost-Benefit Analysis (CBA)
###########################

CBA Workflow
============

At the high-level, the Scenario Building (SB) builds and solves the base network(s) for each scenario/planning horizon, optimizing investments and dispatch.
The CBA does not re‑optimize capacities -- it reuses the SB’s solved network as a starting point, makes modifications (such as fixing capacities), 
and then runs dispatch-only optimizations on the prepared networks.

A diagram of the workflow between SB and CBA is shown below:

.. image:: img/tyndp/SB-CBA-workflow-subsequent-h.png

The CBA workflow is described in detail in the ``doc/cba-rolling-horizon.rst`` document.


Calculate CBA Indicators
========================

After solving both the reference and project networks, several key indicators are calculated to assess benefits of each project and to determine whether the project provides a positive net benefit to the energy system.
The indicators calculated in the CBA are described in additional detail in the ``doc/cba-indicators.rst`` document.
The calculations of the CBA indicators are implemented in the ``scripts/cba/make_indicators.py`` script. 
The calculated CBA indicators for each project are saved in the CSV file: ``results/cba/{cba_method}/project_{cba_project}_{planning_horizons}.csv``.
