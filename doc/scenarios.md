<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Open-TYNDP Scenarios {#scenarios}

The modelling for the reference scenario, National Trends+ (NT+),
and its economic variants, High Economic Variant (HEV) and Low Economic Variant (LEV),
differ in terms of price and demand assumptions.

!!! warning
    Not all HEV/LEV assumptions are included yet. HEV/LEV currently differ from NT+ only in fuel and CO2 prices. Still missing are:

    - demand adjustments (electricity market, prosumer demand, EV, H2, thermal energy demand)
    - hydrogen import prices
    - benchmarking

Here we explain these modelling differences and how these differences have been implemented into open-tyndp.
We discuss the relevant configuration settings and the implications of these implementation decisions.

Background information can be found in the report from ENTSO-E and ENTSO-G
[TYNDP 2026 Scenarios Methodology Report](https://2026.entsos-tyndp-scenarios.eu/scenarios-methodology-report/introduction/).

## National Trends+

The NT+ scenario is the central scenario of TYNDP 2026. It is derived from the National Trends (NT) dataset, based primarily on the
National Energy and Climate Plans and other official national strategies available at the data collection cut-off date,
to which a gap-filling methodology is applied to ensure consistency with selected EU-level targets.

- Uses predefined final energy demands and generation capacity collected from transmission system operators
- Modelled as a pure dispatch problem, with fixed generation capacities for 2030, 2035, 2040 and 2050
- All energy carriers included in final energy demand, not just electricity and gas (as with TYNDP < 2024)


## Open-TYNDP Implementation

The TYNDP scenarios are defined in `config/config.tyndp.yaml` and `config/scenarios.tyndp.yaml`.

The NT+ scenario is defined in full in `config/config.tyndp.yaml` and the modifications to NT+
to form the two economic variants are defined in `config/scenarios.tyndp.yaml`.

| Config Key | Description |
|---|---|
| `run` | Run naming/prefix and scenario file reference; scenarios enabled. |
| `foresight` | Sets planning foresight to myopic. |
| `tyndp_scenario` | Selects TYNDP scenario code (NT, HEV or LEV). |
| `scenario` | Defines cluster set and planning horizons (2030, 2035, 2040, 2050). |
| `countries` | Lists modeled countries/regions. |
| `snapshots` | Time span for simulation (2009 calendar year). |
| `co2_budget` | Placeholder (empty). |
| `electricity` | Core electricity settings—network base, extendable/conventional/renewable carriers, TYNDP mappings, storage types, renewable capacity estimation toggle, PECD profiles (years/techs), PECD hydro profiles, PEMMDB capacities (years/techs), transmission limit version. |
| `atlite` | Cutout configuration for weather data (extent, resolution, time). |
| `links` | Default link power limits (`p_max_pu`/`p_min_pu`). |
| `transmission_projects` | Enables projects; all source sets disabled. |
| `load` | Demand source and year availability; default gap fill and adjustments disabled. |
| `pypsa_eur` | Carrier-to-component mappings for imported PyPSA-Eur data. |
| `biomass` | Sustainable/unsustainable biomass shares over time. |
| `sector` | Toggles for sector coupling and detailed demand shares; transport/shipping/aviation settings; CO2 sequestration options; fuels and networks; biomass and e-fuels options; imports; offshore hubs limits. |
| `costs` | Overwrites for lifetimes/efficiencies and CO2 price trajectory. |
| `clustering` | Spatial/temporal clustering and network simplification settings. |
| `adjustments` | Optional scaling factors for components; used to scale the biogas and solid biomass potentials to TYNDP values. |
| `solving` | Solver selection and option set (uses HiGHS by default). |
| `plotting` | Thresholds, map projection, balance map settings and factors. |
| `benchmarking` | Enables benchmarking. |
| `cba` | Cost-benefit analysis settings (hurdle costs, horizons, projects, solver options). |

The base config is `config.tyndp.yaml`. The scenario file `scenarios.tyndp.yaml` defines per-scenario override blocks (e.g., NT, HEV, LEV).
When a scenario is selected, its keys are merged onto the base config: matching keys override the base values,
and nested keys override only their sub-keys.

Examples:

- `tyndp_scenario` is overwritten by the scenario's value.
- In HEV/LEV, `costs.emission_prices.co2` replaces the base CO2 price trajectory.
- If a key is not present in the scenario block, the base config value remains unchanged.
