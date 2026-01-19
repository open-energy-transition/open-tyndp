.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

####################################
Cost-Benefit Analysis (CBA) Workflow
####################################

At the high-level, the Scenario Building (SB) builds and solves the base network(s) for each scenario/planning horizon, optimizing investments and dispatch.
The CBA does not re‑optimize capacities -- it reuses the SB’s solved network as a starting point, makes modifications (such as fixing capacities), 
and then runs dispatch-only optimizations on the prepared networks.

TOOT and PINT
=============

Within the CBA, there are two methods to evaluate the cost-benefit impact of a project: TOOT (Take Out One at a Time) and PINT (Put IN at a Time).
The TOOT method evaluates the impact of removing a project from the reference network. 
Conversely, the PINT method evaluates the impact of adding a project to the reference network.
The method(s) are configured in the `cba.methods` section of the configuration file.

Build CBA Network
=================

Currently, the CBA network is built by pulling from a solved SB (Scenario Building) network (saved in `results/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc`). 
The CBA network is constructed by copying the SB network and then applying certain modifications to simplify and alter the SB network (in `scripts/cba/simplify_sb_network.py`):
- **Fix optimal capacities**: The CBA dispatch uses SB‑optimized capacities, so the SB’s capacity expansion decisions are embed into the CBA.
- **Disable volume limits for generators and links**: The volume limits (`e_sum_min`) of generators and links are set to negative infinity, as volume limits constrain minimum energy production over the optimization period.
- **Apply hurdle costs**: The hurdle costs from the `cba.hurdle_costs` configuration section are set as marginal costs for links.
- **Disable cycling for long-term storage**: The `e_cyclic` and `e_cyclic_per_period` for stores and the `cyclic_state_of_charge` for storage units are set to False. 

After the simiplified SB network is created, this simplified SB network is further adapted to build the CBA reference network. 
The CBA reference grid is then used to build the project network -- the project network refers to the network used to evaluate the cost-benefit impact of a (tranmission, storage, etc) project.
The method used to build the project network is dependent on the method used (TOOT vs PINT):
- **TOOT**: The project network is built by removing the project from the CBA reference network.
- **PINT**: The project network is built by adding the project to the CBA reference network.

Solve CBA Network
=================

The CBA solve optimizes the dispatch of generators, storage, and transmission in both the reference and project networks.

The CBA network is solved using a rolling horizon, which is configured in the
`cba.solving.options.horizon` and `cba.solving.options.overlap` sections of the configuration file. 
The rolling horizon approach splits the entire time horizon into smaller time windows (or "horizons"), 
which may overlap (if configured). Each horizon is solved sequentially.

Within the CBA solve, the storage state of charge (SOC) is carried over between rolling horizons: 
`stores.e_initial` and `storage_units.state_of_charge_initial` at each rolling horizon are updated using 
the previous rolling horizon. 