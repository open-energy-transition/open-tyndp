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

**Methodology**
1. Compute the total system costs for the reference network and the project network. This can be done by summing together the capital expenditures (CAPEX) and operational expenditures (OPEX) of all components in the network.
```python
capex = n.statistics.capex().sum()
opex = n.statistics.opex(aggregate_time="sum").sum()
total = capex + opex
```

2. Calculate the difference in total system costs between the reference and project networks (depending on the method used, TOOT or PINT):
   - For TOOT: B1 = Total System Costs (Reference Network) - Total System Costs (Project Network)
   - For PINT: B1 = Total System Costs (Project Network) - Total System Costs (Reference Network)

**Outputs**
The following columns are saved in the CBA indicators CSV file (`results/cba/{cba_method}/indicators/cba_indicators_{planning_horizons}.csv`):
- `B1_total_system_cost_change`
- `cost_reference`
- `capex_reference`
- `opex_reference`
- `cost_project`
- `capex_project`
- `opex_project`
- `capex_change`
- `opex_change`
- `is_beneficial` (TRUE/FALSE)
-  `interpretation``: a sentence stating whether the project is beneficial or not based on the B1 indicator.

B2: Social costs of CO2 emissions
=================================


B3: Renewable Energy Sources (RES) integration costs
====================================================
