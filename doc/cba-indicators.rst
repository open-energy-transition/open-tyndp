.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

######################################
Cost-Benefit Analysis (CBA) Indicators
######################################

Within the Cost-Benefit Analysis (CBA), several key indicators are calculated to assess the economic and 
operational impact of projects on the energy system. The calculation of these indicators is based on the results of the solved 
CBA (reference and project) networks.

The following indicators are computed in the CBA workflow:
- B1: Social Economic Welfare (SEW)
- B2: Social costs of CO2 emissions
- B3: Renewable Energy Sources (RES) integration costs

B1: Social Economic Welfare (SEW)
=================================

This indicator quantifies the change in total system costs resulting from the implementation of a project.

**Computation**

1. Compute the total system costs for the reference network and the project network. This is done by summing together the capital expenditures (CAPEX) and 
operational expenditures (OPEX) of all components in each network (done separately for the reference network and the project network).

.. code-block:: python

   capex = n.statistics.capex().sum()
   opex = n.statistics.opex(aggregate_time="sum").sum()
   total = capex + opex

2. Calculate the difference in total system costs between the reference and project networks (depending on the method used, TOOT or PINT):

   - For TOOT: ``B1 = Total System Costs (Reference Network) - Total System Costs (Project Network)``
   - For PINT: ``B1 = Total System Costs (Project Network) - Total System Costs (Reference Network)``

**Outputs**

The following columns are saved in the CBA indicators CSV file (``results/cba/{cba_method}/indicators/cba_indicators_{planning_horizons}.csv``):

- ``B1_total_system_cost_change``
- ``cost_reference``
- ``capex_reference``
- ``opex_reference``
- ``cost_project``
- ``capex_project``
- ``opex_project``
- ``capex_change``
- ``opex_change``
- ``is_beneficial`` (TRUE/FALSE)
- ``interpretation``: a sentence stating whether the project is beneficial or not based on the B1 indicator.

B2: Social costs of CO2 emissions
=================================

This indicator quantifies the social cost impact of CO2 emissions resulting from the implementation of a project, 
assuming different societal cost values for CO2.

**Computation**

1. Extract the net CO2 from the final spashot of the CO2 store (done separately for the reference network and the project network):

..code-block:: python

    net_co2 = n.stores_t.e.T.groupby(n.stores.carrier).sum().T["co2"].iloc[-1]

2. Calculate the difference in net CO2 emissions between the reference and project networks (depending on the method used, TOOT or PINT):

   - For TOOT: ``Delta CO2 = Net CO2 (Reference Network) - Net CO2 (Project Network)``
   - For PINT: ``Delta CO2 = Net CO2 (Project Network) - Net CO2 (Reference Network)``

3. Calculate the social cost of CO2 emissions using the societal cost values ((``cba.co2_societal_cost``)) and CO2 ETS price (``costs.emission_prices.co2``):

   ``B2 = Delta CO2 * (Societal Cost of CO2 - CO2 ETS Price)``

   Note that the societal cost of CO2 is planning horizon specific (e.g., different values for 2030 and 2040). 
   Additionally, there are three different values for the societal cost of CO2: low, central, and high.
   The values for the societal cost of CO2 are taken from page 68 of the 2024 CBA Implementation Guidelines.

**Outputs**

The following columns are saved in the CBA indicators CSV file (``results/cba/{cba_method}/indicators/cba_indicators_{planning_horizons}.csv``):

- ``co2_diff``
- ``co2_ets_price``
- ``co2_societal_cost_low``
- ``co2_societal_cost_central``
- ``co2_societal_cost_high``
- ``B2_societal_cost_variation_low``
- ``B2_societal_cost_variation_central``
- ``B2_societal_cost_variation_high``

B3: Renewable Energy Sources (RES) integration costs
====================================================

This indicator is to quantify how the project improves the integration of renewable energy sources (RES) into the energy system.
RES technologies are defined in ``electricity.tyndp_renewable_carriers`` in the configuration file.

**Computation**

1. Compute the total capacities, generation, and curtailment of RES technologies in both the reference and project networks.

    - Capacities: ``n.statistics.optimal_capacity(...).groupby("carrier").sum()``
    - Generation: ``n.statistics.energy_balance(aggregate_time="sum", ...).groupby("carrier").sum()``
    - Curtailment: ``n.statistics.curtailment(...).groupby("carrier").sum()``

2. Calculate the change in RES capacities, generation, and curtailment between the reference and project networks 
(depending on the method used, TOOT or PINT):

   - For TOOT: ``Delta RES = RES (Reference Network) - RES (Project Network)``
   - For PINT: ``Delta RES = RES (Project Network) - RES (Reference Network)``

**Outputs**

The following columns are saved in the CBA indicators CSV file (``results/cba/{cba_method}/indicators/cba_indicators_{planning_horizons}.csv``):
- ``B3_res_capacity_change_mw``
- ``B3_res_generation_change_mwh``
- ``B3_annual_avoided_curtailment_mwh``