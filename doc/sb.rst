.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

######################
Scenario Building (SB)
######################

Scenario Building (SB) is the first phase of the Open-TYNDP workflow. It constructs a fully sector-coupled European energy system model for a given scenario and planning horizon and solves a least-cost capacity expansion and dispatch optimisation. The solved network produced by SB serves as the direct input to the `Cost-Benefit Analysis (CBA) <cba.html>`_ phase.

Open-TYNDP implements the TYNDP 2024 Scenario Building methodology as a soft-fork of `PyPSA-Eur <https://pypsa-eur.readthedocs.io/en/latest/>`_, inheriting its full modelling framework—optimisation structure, network representation, and sector-coupling capabilities—while replacing specific inputs and assumptions to match TYNDP 2024 reference data.

.. image:: img/tyndp/SB-CBA-workflow-subsequent-h.png
    :align: center
    :alt: Workflow between Scenario Building and CBA

SB Workflow
===========

The SB workflow transforms raw ENTSO-E input datasets into a solved, sector-coupled PyPSA network. The key stages are: integrating public input data, constructing the multi-sector network, applying TYNDP-specific constraints, and solving the capacity expansion optimisation. The solved network is archived and retrieved by the CBA workflow.

Public input data from ENTSO-E is used wherever available. Where publicly available data does not match the fixed values observed in the TYNDP 2024 Market Model output files, the output files are used as the reference for fixed input assumptions. This applies strictly to exogenous variables that are not part of the optimisation—primarily H₂ demand profiles, reference grid topologies, and generator maintenance profiles.

.. image:: img/tyndp/sb-workflow-overview.png
    :width: 70%
    :align: center
    :alt: Scenario Building workflow stages

Input Data Integration
^^^^^^^^^^^^^^^^^^^^^^

The following ENTSO-E datasets are ingested and transformed into PyPSA-compatible format by dedicated ``build_tyndp_*`` Snakemake rules:

* **PEMMDB 2.5:** Installed generation and storage capacities, must-run constraints, and per-unit cost assumptions by country and technology.
* **PECD 3.1:** Hourly capacity factor time series for wind (onshore/offshore) and solar PV, derived from ERA5 reanalysis.
* **Hydro Inflows:** Hourly inflow profiles for reservoir and run-of-river hydro plants.
* **Demand Profiles:** Hourly electricity and hydrogen demand profiles by country, interpolated to the target planning horizon.
* **Line Data:** Electricity and hydrogen transmission network topology; the base grid is taken from the TYNDP reference, with candidate lines applicable from 2030 onward.
* **Hydrogen datasets:** Hydrogen storage parameters, steam methane reforming (SMR) capacities, and import pipeline assumptions.
* **Investment Candidates:** Optional extendable transmission and storage assets for 2035 and 2040 planning horizons.

Network Construction
^^^^^^^^^^^^^^^^^^^^

The PyPSA network is built at **country-level resolution** for the electricity and hydrogen sectors. Buses represent national-level aggregations; AC lines and DC links represent cross-border interconnectors with capacities and impedances taken from the TYNDP Line Data.

Generator and storage components are attached to country buses using PEMMDB 2.5 capacity data. Which carriers are **extendable**—meaning the optimiser may invest in additional capacity beyond the fixed assumptions—is controlled by the ``electricity.extendable_carriers`` configuration key and varies by scenario and planning horizon.

Sector Coupling
^^^^^^^^^^^^^^^

Open-TYNDP models the electricity and hydrogen sectors as fully coupled. For the Distributed Energy (DE) scenario, heating sector links are included in addition. Cross-sector components include:

* **Electrolysers:** Convert electricity to hydrogen; capacity is either fixed per PEMMDB or left extendable depending on the scenario.
* **Fuel cells and back-pressure plants:** Reconvert hydrogen or gas to electricity.
* **Hydrogen network:** Dedicated H₂ pipelines between country buses are included for planning horizons from 2030 onward, using TYNDP Line Data. Zones with split H₂ grids (e.g., the Iberian Peninsula) are represented with separate H₂ buses.
* **Demand-side electrification:** Where the scenario specifies it, electricity demand incorporates direct electrification of heat and transport end-uses.

Capacity Optimisation
^^^^^^^^^^^^^^^^^^^^^

The SB optimisation minimises **total annualised system cost** (investment plus variable operating cost) subject to:

* Hourly supply–demand balance for each carrier at every bus.
* Transmission capacity constraints, with optional extendability for candidate lines.
* CO₂ emission budgets derived from the TYNDP 2024 scenario pathway.
* Minimum and maximum generation constraints from PEMMDB, including must-run levels and scheduled maintenance outages.
* Country-level annual hydrogen supply and demand balances.

The problem is formulated as a **linear programme (LP)** and solved with the configured solver (Gurobi by default; HiGHS is supported as an open-source fallback). Each planning horizon is solved independently using the capacity assumptions fixed for that horizon.

CBA Handoff
^^^^^^^^^^^

On completion, the solved ``network.nc`` file is written to the run archive directory. The CBA workflow retrieves this file as its starting point, fixing all optimised capacities before running dispatch-only simulations for project evaluation. See the `CBA documentation <cba.html>`_ for details.

Configuration
=============

SB settings are split across ``config/config.tyndp.yaml`` (run-level settings) and ``config/scenarios.tyndp.yaml`` (scenario-specific overrides).

Scenarios and Planning Horizons
--------------------------------

* ``scenario``: Selects the TYNDP 2024 scenario. Supported values are ``NT`` (National Trends) and ``DE`` (Distributed Energy).
* ``planning_horizons``: List of target years to solve (e.g., ``[2030, 2040]``). Each horizon is solved as an independent optimisation.
* ``run.name``: Identifies the run and determines the output directory; typically set to the scenario/climate-year identifier (e.g., ``NT-cy2009``).

Climate Years
-------------

Each scenario run is tied to a specific historical climate year, which determines the renewable generation and hydro inflow profiles:

* ``snapshots``: Defines the modelling time window, e.g.:

  .. code-block:: yaml

      snapshots:
        start: "2009-01-01"
        end:   "2009-12-31"
        inclusive: "left"

* ``atlite.default_cutout``: ERA5 reanalysis cutout used to compute PECD-compatible capacity factor profiles (e.g., ``europe-2009-era5``).

Solver Settings
---------------

* ``solving.solver.name``: Solver to use (``gurobi`` or ``highs``).
* ``solving.solver_options``: Solver-specific parameters such as optimality gap and memory limits.
* ``solving.options.linearized_unit_commitment``: Set to ``true`` to enable linearised unit commitment constraints for thermal plant cycling.

Running Scenario Building
=========================

Running a Single Scenario
^^^^^^^^^^^^^^^^^^^^^^^^^^

Set ``run.name`` in ``config/config.tyndp.yaml`` to the desired scenario identifier, then execute:

.. code-block:: console

    $ snakemake -call solve_sector_networks --configfile config/config.tyndp.yaml

This solves all ``planning_horizons`` defined in the active scenario sequentially and writes the solved networks to the run archive.

Running Multiple Climate Years
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Define each climate year as a named scenario in ``config/scenarios.tyndp.yaml``:

.. code-block:: yaml

    NT-cy2009:
      snapshots:
        start: "2009-01-01"
        end:   "2009-12-31"
        inclusive: "left"
      atlite:
        default_cutout: europe-2009-era5

    NT-cy2008:
      snapshots:
        start: "2008-01-01"
        end:   "2008-12-31"
        inclusive: "left"
      atlite:
        default_cutout: europe-2008-era5

Then override ``run.name`` on the command line for each climate year:

.. code-block:: console

    $ snakemake -call solve_sector_networks --configfile config/config.tyndp.yaml \
        --config run='{"name":"NT-cy2008"}'

Running Both Scenarios
^^^^^^^^^^^^^^^^^^^^^^^^

To run National Trends and Distributed Energy side by side, define both in ``config/scenarios.tyndp.yaml`` with their respective capacity and demand overrides, then invoke the workflow once per scenario name.
