.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _scenarios:

####################
Open TYNDP Scenarios
####################

The modelling for the reference, National Trends (NT), 
and the deviation, Distributed Energy (DE) and Global Ambition (GA), scenarios
differ in terms of assumptions, storyline and method (capacity expansion vs dispatch modelling).

Here we explain these modelling differences and these differences have been implemented into open-tyndp.
We discuss the relevant configuration settings and the implications of these implementation decisions.

Background information can be found in the report from ENTSO-E and ENTSO-G 
`TYNDP 2024 Scenarios Methdoology Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Methodology_Report_Final_Version_250128.pdf>`__.

A description of the storylines underlying the NT, DE and GA scenarios can be found in the
`TYNDP 2024 Scenarios Storyline Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2023/12/ENTSOs_TYNDP_2024_Scenarios_Storyline_Report_2023-12-04.pdf>`__

National Trends+
~~~~~~~~~~~~~~~~

The NT+ scenario was developed in alignment with the National Energy and Climate Plans that were available during the development
of TYNDP 2024 (e.g. during 2023). 

- Uses predefined final energy demands and generation capacity collected from transmission system operators.
- All energy carriers included in final energy demand, not just electricity and gas (as with TYDNP < 2024)
- Run for 2030 and 2040

Deviation Scenarios
~~~~~~~~~~~~~~~~~~~

- Begins with National Trends+ 2030 scenario results for grid topology, generation capacity and final energy demands
- Uses predefined final energy demands for 2040 and 2050 collected from transmission system operators and public consultation
- Capacity expansion model is used to identify where and when investment in generation capacity is required
- Exogenous constraints on generation investment are imposed to force:
    - no new nuclear power capacity
    - all existing fossil gas plants are decommissioned
- upper and lower bounds are also imposed to force trajectories for
    - the expansion of solar, wind, prosumer batteries and large scale batteries
    - H2 import potentials

Distributed Energy
------------------

.. epigraph::

   *This scenario pictures a pathway achieving EU27 carbon neutrality target by 2050 with
   higher European Economy. The scenario is driven by a willingness of the society to achieve
   high levels of independence in terms of energy supply and goods of strategic importance
   (e. g., industrial and agricultural produce). It translates into both a behavioural shift and
   strong decentralised drive towards decarbonisation through local initiatives by citizens,
   communities and businesses, supported by authorities*

   -- `TYNDP 2024 Scenarios Storyline Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2023/12/ENTSOs_TYNDP_2024_Scenarios_Storyline_Report_2023-12-04.pdf>`__

Global Ambition
---------------

.. epigraph::

   *This scenario pictures a pathway to achieving carbon neutrality by 2050, driven by a fast
   and global move towards the Paris Agreement targets. It translates into development of a
   very wide range of technologies (many being centralised) and the use of global energy
   trade as a tool to accelerate decarbonisation.*

   -- `TYNDP 2024 Scenarios Storyline Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2023/12/ENTSOs_TYNDP_2024_Scenarios_Storyline_Report_2023-12-04.pdf>`__

Open TYNDP Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

The TYNDP scenarios are defined in ``config/config.tyndp.yaml`` and ``config/scenarios.tyndp.yaml``.

The NT+ scenario is defined in full in ``config/config.tyndp.yaml`` and the modifications to NT+
to form the two deviation scenarios are defined in ``config/scenarios.tyndp.yaml``.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Config Key
     - Description
   * - run
     - Run naming/prefix and scenario file reference; scenarios enabled.
   * - foresight
     - Sets planning foresight to myopic.
   * - tyndp_scenario
     - Selects TYNDP scenario code (NT).
   * - scenario
     - Defines cluster set and planning horizons (2030, 2040).
   * - countries
     - Lists modeled countries/regions.
   * - snapshots
     - Time span for simulation (2009 calendar year).
   * - co2_budget
     - Placeholder (empty).
   * - electricity
     - Core electricity settings—network base, extendable/conventional/renewable carriers, TYNDP mappings, storage types, renewable capacity estimation toggle, PECD profiles (years/techs), PEMMDB hydro profiles and capacities (years/techs), transmission limit version.
   * - atlite
     - Cutout configuration for weather data (extent, resolution, time).
   * - links
     - Default link power limits (``p_max_pu``/``p_min_pu``).
   * - transmission_projects
     - Enables projects; all source sets disabled.
   * - load
     - Demand source and year availability; gap fill and adjustments disabled.
   * - pypsa_eur
     - Carrier-to-component mappings for imported PyPSA-Eur data.
   * - biomass
     - Sustainable/unsustainable biomass shares over time.
   * - sector
     - Toggles for sector coupling and detailed demand shares; transport/shipping/aviation settings; CO2 sequestration options; fuels and networks; biomass and e-fuels options; imports; offshore hubs limits.
   * - costs
     - Overwrites for lifetimes/efficiencies and CO2 price trajectory.
   * - clustering
     - Spatial/temporal clustering and network simplification settings.
   * - adjustments
     - Optional scaling factors for sector components; currently off.
   * - solving
     - Solver selection and option set (HiGHS).
   * - plotting
     - Thresholds, map projection, balance map settings and factors.
   * - benchmarking
     - Enables benchmarking.
   * - cba
     - Cost-benefit analysis settings (hurdle costs, horizons, methods, projects, solver options).

The base config is ``config.tyndp.yaml``. The scenario file ``scenarios.tyndp.yaml`` defines per-scenario override blocks (e.g., NT, DE, GA). 
When a scenario is selected, its keys are merged onto the base config: matching keys override the base values, 
and nested keys override only their sub-keys.

Examples:

- ``tyndp_scenario`` is overwritten by the scenario's value.
- In DE/GA, ``electricity.extendable_carriers.Generator`` replaces the base list for that path.
- In DE/GA, ``sector.land_transport_ice_share``, ``sector.h2_zones_tyndp``, ``sector.force_biomass_potential``, ``sector.force_biogas_potential``, 
  and ``sector.co2_sequestration_potential`` override the corresponding base entries.
- If a key is not present in the scenario block, the base config value remains unchanged
