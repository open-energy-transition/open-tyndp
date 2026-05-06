.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

###########################
Cost-Benefit Analysis (CBA)
###########################

The Cost-Benefit Analysis (CBA) evaluates transmission and storage projects by comparing two dispatch-only simulations: a **reference network** (the baseline grid) and a **project network** (the grid with a specific project added or removed)[cite: 4, 75]. Unlike the Scenario Building (SB) phase, the CBA does not re-optimize capacities; it reuses the SB's solved network, fixes capacities, and runs dispatch-only optimizations[cite: 74, 75].

.. image:: img/tyndp/SB-CBA-workflow-subsequent-h.png
    :align: center
    :alt: Workflow between Scenario Building and CBA

CBA Workflow Methodology
========================

The workflow evaluates projects using a **rolling horizon** approach where the full year is divided into sequential weekly windows (168 hourly snapshots each, with an overlap of 1 snapshot)[cite: 5]. 

To resolve **myopia**—where the optimizer cannot see beyond the current week and makes suboptimal decisions for seasonal storage (H2, gas, large hydro)—the workflow uses Marginal Storage Values (MSV) derived from a full-year optimization[cite: 2, 6].

.. image:: img/tyndp/cba-rolling-horizon-pipeline.jpeg
    :width: 60%
    :align: center
    :alt: CBA rolling horizon pipeline diagram

Pipeline Stages
---------------

Stage 0 — Network Simplification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The SB network is transformed into a dispatch-ready CBA network[cite: 10]. 
* **Fixed Capacities:** Capacities are fixed via ``n.optimize.fix_optimal_capacities()``[cite: 9].
* **Hurdle Costs:** A cost of 0.01 €/MWh is applied to all DC links[cite: 11].
* **Fuel Capacities:** Primary fuel generator capacities (coal, gas, oil, nuclear) are set to infinity to prevent dispatch from being artificially restricted by fuel-supply limits during peak hours[cite: 12, 15].

Stage 1 — Reference Network
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The simplified network is extended to form the CBA reference baseline by adding all TOOT project capacities[cite: 17]. This ensures the reference and MSV extraction operate on the same topology[cite: 19].

Stage 2 — MSV Extraction
^^^^^^^^^^^^^^^^^^^^^^^^
The reference network is solved with **perfect foresight** (entire year, single LP)[cite: 25]. This exposes the shadow prices (``mu_energy_balance``) of energy balance constraints, representing the **Marginal Storage Value (MSV)**—the opportunity cost of stored energy at that moment[cite: 26, 27].

Stage 3 — Rolling Horizon Preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The reference network and MSV results are combined through five transformations:
* **(a) Initial storage state:** Seasonal components have their initial state set to the perfect foresight solution's last-snapshot value[cite: 33].
* **(b) Disable cyclicity:** Short-term storage (battery) keeps cyclicity, while seasonal units (H2, gas, hydro) have it disabled, guided instead by MSVs[cite: 35, 37].
* **(c) Remove global constraints:** Annual CO2 and biomass limits are removed as they cannot be enforced consistently in weekly windows[cite: 39].
* **(d) Disable annual volume limits:** Annual budgets for biomass/biogas are distributed proportionally to the perfect foresight dispatch[cite: 41, 42].
* **(e) Apply MSV:** The ``mu_energy_balance`` time series is written into the ``marginal_cost`` of storage components[cite: 43].
* **(f) Hydro Pinning:** Large hydro reservoirs are pinned to their perfect-foresight state-of-charge values at window boundaries to guide dispatch where duals are near-zero[cite: 46, 48, 68].

Stage 4 & 5 — Solve
^^^^^^^^^^^^^^^^^^^
The prepared network is solved for the **Reference** baseline and then for each **Project**[cite: 52, 60]. Projects are evaluated using either **TOOT** (Take Out One at a Time) or **PINT** (Put IN at a Time) methods[cite: 56, 58].

CBA Indicators
==============

Indicators are computed as the difference in system costs and emissions between the reference and project dispatch solutions[cite: 64, 77].

* **B1: Social Economic Welfare (SEW):** Quantifies the change in total system costs (CAPEX + OPEX)[cite: 126].
* **B2: Social costs of CO2 emissions:** Calculates the impact using societal cost assumptions (low, central, high)[cite: 131, 134].
* **B3: RES integration costs:** Tracks changes in renewable capacity, generation, and avoided curtailment[cite: 136, 138].
* **B4: Non-direct greenhouse emissions:** Quantifies pollutants (NOx, SO2, PM, etc.) using fuel consumption multipliers[cite: 140, 142].

Configuration Reference
=======================

CBA settings are defined in the ``cba`` section of the configuration file[cite: 66, 80].

Project Selection
-----------------
* ``planning_horizons``: Selects horizons (e.g., 2030, 2040)[cite: 82].
* ``projects``: Defines project identifiers (e.g., ``t1-t35``)[cite: 84].
* ``cba_scenario_input``: If ``use_presolved`` is true, the workflow retrieves pre-solved SB networks from an archive[cite: 87, 88].

Rolling Horizon Settings
------------------------
* ``storage.cyclic_carriers``: Carriers that remain cyclic within each weekly window[cite: 99].
* ``storage.soc_boundary_carriers``: Carriers pinned at window boundaries[cite: 101].
* ``msv_extraction.resolution``: Controls temporal resolution for the MSV solve (e.g., ``24H``)[cite: 102, 104].

Running Single vs Multiple Climate Years
========================================

Climate-year collections allow project benefits to be assessed across multiple weather years, consistent with the 2024 TYNDP implementation[cite: 112].

The CBA entry point ``snakemake -call cba`` expects a **collection scenario** that defines a list of child (climate years) scenarios under ``cba.scenarios``[cite: 113, 114].

Example Collection (``config/scenarios.tyndp.yaml``):

.. code-block:: yaml

    NT-cyears:
      cba:
        scenarios: [NT-cy2009, NT-cy2008, NT-cy1995]

Individual child scenarios (e.g., ``NT-cy2009``) must define their specific ``snapshots``, ``atlite.default_cutout``, and the ``cba.sb_scenario`` used as input[cite: 115, 116].

Running Multiple Years
----------------------
To run a collection like ``NT-cyears``, modify ``run.name`` in ``config/config.tyndp.yaml`` or override it via command line[cite: 118]:

.. code-block:: console

    $ snakemake -call cba --configfile config/config.tyndp.yaml --config run='{"name":"NT-cyears"}'

Running a Single Climate Year
-----------------------------
Because the CBA rule requires a list of scenarios, a single climate year must be wrapped in a **one-entry collection scenario**[cite: 121, 122]:

.. code-block:: yaml

    NT-cy2009-only:
      cba:
        scenarios: [NT-cy2009]