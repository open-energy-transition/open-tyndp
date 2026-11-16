<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# TYNDP 2026 data
If you move from TYNDP 2024 data, see the [TYNDP 2024 to 2026 migration guide](#migration-tyndp-2026).


The TYNDP 2026 data is downloaded and extracted using the `retrieve_tyndp_2026` and `retrieve_tyndp_pecd` rules. Originally published by ENTSO-E and ENTSOG under
Creative Commons Attribution 4.0 International License (CC-BY 4.0).
The original data files can be found under <https://2026.entsos-tyndp-scenarios.eu/download/>.
The `{version}`, `{pecd_version}` and `{pecd_prebuilt_version}` placeholders refer to the dataset versions listed in `data/tyndp_versions.csv`.

## `Demand`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Demand.zip`
- **Description:** Contains hourly electricity, hydrogen, synthetic fuel, behind-the-meter, EV fleet and thermal (hydrogen and methane) demand profiles by node and weather scenario, for 2030, 2035, 2040 and 2050.

## `Electric-Vehicle`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Electric-Vehicle.zip`
- **Description:** Contains the passenger EV fleet inputs: EV stock and demand, market zone, flexible fleet and V2G shares, charging station and battery parameters, and availability, driving and minimum state of charge profiles.

## `Heat`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Heat.zip`
- **Description:** Contains coefficient of performance (COP) profiles per planning horizon and the hybrid heat pump parameters (heat pumps and methane or hydrogen boilers).

## `Hydro-Inflow`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Hydro-Inflow.zip`
- **Description:** Contains PECD v4.2 hydro inflow time series and hydro constraints per weather scenario, for 2030, 2035, 2040 and 2050.

## `Hydrogen`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Hydrogen.zip`
- **Description:** Contains hydrogen storage parameters, steam methane reforming (SMR and SMR+CCS) assets, hydrogen import generators and profiles, and synthetic fuel lines and generator properties.

## `Line-data`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Line-data.zip`
- **Description:** Contains the electricity and hydrogen reference grids and the wheeling charges.

## `Nodes`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Nodes.zip`
- **Description:** Contains a list of electricity and hydrogen nodes.

## `PEMMDB_2.X`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/PEMMDB_2.X.zip`
- **Description:** Contains the PEMMDB 2.X installed generation and storage capacities, must-run constraints and technical parameters per node, for 2030, 2035, 2040 and 2050.

## `PECD/PECD_{pecd_version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_version}.zip`
- **Description:** Contains hourly PECD v4.2 capacity factor time series for Solar PV rooftop,
  Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `PECD/PECD_{pecd_prebuilt_version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_prebuilt_version}.zip`
- **Description:** Contains pre-built hourly PECD v4.2 capacity factor time series for weather scenarios WS003, WS021, WS029, WS032, WS037, WS059, WS065, WS071, WS077, WS091, WS092 and WS106, covering Solar PV rooftop, Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `Supply-Tool`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/outputs/Supply-Tool.zip`
- **Description:** Contains the TYNDP 2026 Supply Tool, with one workbook per scenario (NT+, HEV and LEV).

## `NT+_TimeSeriesDashboard`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/outputs/NT+_TimeSeriesDashboard.zip`
- **Description:** Contains the hourly market model outputs of the NT+ scenario, one file per country, for each planning horizon and weather scenario. Used for benchmarking and to build the nuclear availability profiles: the hourly nuclear generation of the `E-Market` sheet is divided by the installed capacity of the `Installed Capacity` sheet, and the weekly maximum is retained to account for maintenance schedules and forced outages.

## `TYNDP_2026_Scenario_Report_Data_Figures.xlsx`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/additional/TYNDP_2026_Scenario_Report_Data_Figures.xlsx`
- **Description:** Contains the data used to produce the figures in the TYNDP 2026 Scenarios Report. Used for benchmarking.


# TYNDP Migration guide (2024 -> 2026) {#migration-tyndp-2026}

This guide is for users who ran Open-TYNDP on TYNDP 2024 data, up to release `v0.9`. It explains what changed for the TYNDP 2026 cycle and lists all renamed, added and removed names. Read the upstream [Migration Guide](migration.md) first, as it covers the workflow restructuring of the same release.

!!! danger "Removed and renamed keys are ignored, not rejected"
    Removed or renamed config keys are silently ignored, so their settings have no effect. Old values of choice lists, such as `tyndp_scenario: DE`, fail validation at start-up. List unknown keys outside `plotting`, `conventional` and `benchmarking` with:

    ```console
    $ pixi run -e open-tyndp validate-config config/config.yaml
    ```

## Part A: What changed and why

**Data cycle.** TYNDP 2026 is a new dataset under `data: tyndp_2026:`, with version `2026-11-11` ([#893](https://github.com/open-energy-transition/open-tyndp/pull/893), [#1148](https://github.com/open-energy-transition/open-tyndp/pull/1148)). PECD moves from `3.1` to `4.2.1`. The 2024 dataset is still used for the methane demand and the CBA project list. With the default data config, no change is needed. If you pin versions, keep both `data: tyndp:` and `data: tyndp_2026:`.

**Scenarios.** HEV and LEV replace the DE and GA scenarios ([#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066)). Following the TYNDP 2026 methodology, they add +10 % and -10 % to the NT+ fuel and CO2 prices. In Open-TYNDP, they differ from NT+ only in these prices. Replace `DE` or `GA` with `NT`, `HEV` or `LEV` in `tyndp_scenario` and your scenario files.

**Weather scenarios and full year.** Numbered weather scenarios per horizon replace the 2024 climate years, and `cyear` became `wscenario` ([#878](https://github.com/open-energy-transition/open-tyndp/pull/878)). TYNDP 2026 models the full year, so snapshots cover 8760 h instead of 52 weeks ([#1065](https://github.com/open-energy-transition/open-tyndp/pull/1065)). Set `wscenarios_tyndp` and rename `pecd_renewable_profiles.pre_built.cyears` to `wscenarios`. Set `snapshots.end` to the first day of the next year.

**Horizons and countries.** TYNDP 2026 provides data for 2030, 2035, 2040 and 2050, so the TYNDP config solves all four. The 2026 node list adds Moldova, Turkey and Ukraine ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). Set the top-level `planning_horizons`, as described in the [upstream guide](migration.md#scenario-block). Add `MD`, `TR` and `UA` to `countries`.

**Offshore.** 2026 offshore nodes are real substations with a real reference-grid topology ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). The 2024 offshore hubs are therefore removed, and offshore buses carry 2026 codes such as `BEO1_OFF`. Offshore wind becomes `offwind-r` and `offwind-oh` plus SRES and DRES variants ([#999](https://github.com/open-energy-transition/open-tyndp/pull/999)). Remove `sector.offshore_hubs_tyndp` and the old offshore technologies. Update analyses that use `...OH001` buses.

**Generation.** Capacities follow PEMMDB in all horizons without expansion, so trajectories and investment candidates are removed ([#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051)). SRES and DRES generation gets own buses to feed electrolysers ([#999](https://github.com/open-energy-transition/open-tyndp/pull/999)). Nuclear profiles come from the 2026 dashboards, and hydro inflows use PECD codes. Remove `tyndp_investment_candidates` and update `tyndp_renewable_carriers` to the 2026 names. Rename `electricity.pemmdb_hydro_profiles` to `pecd_hydro_profiles` and use the PECD codes.

**Electricity demand and distribution.** Demand splits into market, prosumer and BTM, with prosumer demand on a `RETE` bus ([#944](https://github.com/open-energy-transition/open-tyndp/pull/944), [#935](https://github.com/open-energy-transition/open-tyndp/pull/935)). The distribution grid uses wheeling charges, and batteries are no longer extendable. Keep `sector.electricity_distribution_grid` enabled, since only it adds the prosumer buses. Add `home battery` to `tyndp_stores` and do not list batteries as extendable stores. Remove `electricity distribution grid` from `sector.transmission_efficiency.enable`, as TYNDP demand includes grid losses.

**Hydrogen.** The 2026 node list always splits H2 buses into Zone 1 and Zone 2 ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). Imports become hourly generators, storage follows flexibility classes and SMR can follow ramp limits ([#931](https://github.com/open-energy-transition/open-tyndp/pull/931)). P2G shares split electrolysers between zones. Remove `h2_zones_tyndp`, `h2_reference_grid_source` and `h2_demand_patch_with_mm`. Rename the stores in `tyndp_stores` and update analyses that use ` H2 Z1` or ` H2 Z2` buses.

**Transport and synthetic fuels.** EV fleets and a flat EU27 oil demand replace the PyPSA-Eur land transport ([#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039)). The oil demand is flat for EU27, because TYNDP reports EU27 only. Keep `sector.transport` enabled, since fixed capacities would otherwise oversupply the system. Remove the `land_transport_*_share` and `bev_charge_efficiency` overrides. Synthetic fuels need `sector.synfuels_tyndp.enable` and `h2_topology_tyndp`, without `industry` or `methanation` ([#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098)).

**Grid.** Reference grids keep one topology in all horizons, and only the NTCs change ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). The network is therefore built once from `electricity.tyndp_reference_year`. Hurdle costs now also apply in Scenario Building, on `DC` and `H2 pipeline` links ([#935](https://github.com/open-energy-transition/open-tyndp/pull/935)). Move `cba.hurdle_costs` to the top-level `hurdle_costs`.

**Benchmarking.** Benchmarks now compare against the 2026 Time Series Dashboards and Scenario Report figures ([#960](https://github.com/open-energy-transition/open-tyndp/pull/960)). The `generation_profiles` table served only DE and GA and is removed ([#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066)). New tables cover prosumer demand, EV charging and synthetic fuels, and HEV and LEV are not benchmarked. Rename `market_out` to `dashboard_out` and remove `vp` and `generation_profiles`. Read outputs from `benchmarks/tyndp-2026/`.

**CBA.** The CBA still runs on the 2024 project list, and its reference grid corrections are now placeholders. Run names changed their prefix from `cy` to `ws` ([#878](https://github.com/open-energy-transition/open-tyndp/pull/878)). `t335` becomes `t225` and `t1106` is removed for now, as both use 2024 offshore hub buses ([#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079)). Rename runs such as `NT-cy2009` to `NT-ws2009`.

**Workflow.** TYNDP rules and files now use the `tyndp` prefix and 2026 names. One rule builds all TYNDP demand with the new `{demand_type}` wildcard ([#808](https://github.com/open-energy-transition/open-tyndp/pull/808)). The weather scenario test config and its pixi task were renamed. Update custom Snakemake targets, scripts that read resource files and calls of the test task.

## Part B: Full list of changes

"Before" names are those of release `v0.9`, "Now" names those of the current release.

### Config keys

Removed and renamed keys:

| Before | Now | PR |
|---|---|---|
| `tyndp_investment_candidates` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `sector.offshore_hubs_tyndp` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `sector.h2_zones_tyndp` | removed (zoning always on) | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `sector.h2_reference_grid_source` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `sector.h2_demand_patch_with_mm` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `load.patch_demand_with_mm` | removed | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944), [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `load.available_years_tyndp` | removed | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| `sector.land_transport_fuel_cell_share`, `sector.land_transport_electric_share`, `sector.land_transport_ice_share` (TYNDP config) | removed | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| `sector.bev_charge_efficiency` (TYNDP config) | removed | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| `data.tyndp_nuclear_profiles` | removed | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `data.tyndp_vis_plfm` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `data.countries_centroids` | removed | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| `sector.electricity_distribution_grid_cost_factor` | removed (upstream) | [PyPSA-Eur #2299](https://github.com/PyPSA/pypsa-eur/pull/2299) |
| `electricity.pemmdb_hydro_profiles` | `electricity.pecd_hydro_profiles` | [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `electricity.pecd_renewable_profiles.pre_built.cyears` | `electricity.pecd_renewable_profiles.pre_built.wscenarios` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `cba.hurdle_costs` | `hurdle_costs` | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |

Added keys:

| Now | PR |
|---|---|
| `data.tyndp_2026` | [#893](https://github.com/open-energy-transition/open-tyndp/pull/893) |
| `wscenarios_tyndp` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `electricity.tyndp_reference_year` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `sector.smr_ramp_limits_tyndp` | [#931](https://github.com/open-energy-transition/open-tyndp/pull/931) |
| `sector.land_transport_oil_demand_tyndp` | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| `sector.synfuels_tyndp.enable`, `sector.synfuels_tyndp.carriers` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| `sector.imports.carriers` gains `e-liquids` and `sng` (TYNDP config) | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| `cba.storage.cyclic_carriers` gains `home battery` (TYNDP config) | [#892](https://github.com/open-energy-transition/open-tyndp/pull/892) |

New validation rules. A TYNDP run fails at start-up if one of them is violated:

| Rule |
|---|
| `sector.electricity_distribution_grid` must be enabled for a TYNDP scenario. |
| `sector.transport` must be enabled for a TYNDP scenario. |
| `electricity distribution grid` must not be in `sector.transmission_efficiency.enable` for a TYNDP scenario. |
| `battery` and `home battery` must not be in both `electricity.tyndp_stores` and `electricity.extendable_carriers.Store`. |
| `sector.synfuels_tyndp.enable` requires `sector.h2_topology_tyndp`. |
| `sector.synfuels_tyndp.enable` cannot be combined with `sector.industry` or `sector.methanation`. |

Changed values: `sector.biogas_upgrading_cc` is now `true` in `config/config.tyndp.yaml` (default `false`, unchanged), see [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098).

Many values in `config/config.tyndp.yaml` were updated to TYNDP 2026 NT+, among them CO2 prices, CO2 sequestration potential, biomass and biogas values and the shipping oil share. See the [release notes](release_notes.md) for the details.

### Scenarios, horizons and weather scenarios

| Before | Now | PR |
|---|---|---|
| `tyndp_scenario`: `NT`, `DE`, `GA` | `tyndp_scenario`: `NT`, `HEV`, `LEV` | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| climate years `1995`, `2008`, `2009` | weather scenarios in `wscenarios_tyndp`, e.g. `2030: [3, 21, 29]` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| planning horizons 2030, 2040 | 2030, 2035, 2040, 2050 | |
| snapshots end `2009-12-31` (8736 h) | snapshots end `2010-01-01` (8760 h) | [#1065](https://github.com/open-energy-transition/open-tyndp/pull/1065) |
| no `MD`, `TR`, `UA` | `countries` adds `MD`, `TR`, `UA` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |

### Technology names

| Before | Now | PR |
|---|---|---|
| `offwind-ac-fb-r`, `offwind-ac-fl-r`, `offwind-dc-fb-r`, `offwind-dc-fl-r` | `offwind-r` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| `offwind-dc-fb-oh`, `offwind-dc-fl-oh` | `offwind-oh` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| `offwind-h2-fb-oh`, `offwind-h2-fl-oh` | removed | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| `offwind` | removed | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| none | `onwind-sres`, `onwind-dres`, `solar-pv-utility-sres`, `solar-pv-utility-dres`, `offwind-sres-r`, `offwind-sres-oh`, `offwind-dres-r`, `offwind-dres-oh` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| none | per H2 zone: `onwind-sres-h2z1`, `solar-pv-utility-sres-h2z2`, `offwind-sres-h2z1-r`, `offwind-sres-h2z2-oh` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| none | `H2 Electrolysis Z1`, `H2 Electrolysis Z2`, `H2 Electrolysis SRES Z1`, `H2 Electrolysis DRES Z2` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999), [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
| none | `h2-ocgt`; index names `other-non-res-h2-ocgt`, `other-non-res-h2-ccgt` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| PECD `LFSolarPVRooftop`, `LFSolarPVUtility` | `Solar PV Rooftop`, `Solar PV Utility` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887) |
| none | PECD `CSP_noStorage_0h_dispatched` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887) |
| hydro `Run_of_River`, `Pondage`, `Reservoir`, `PS_Open`, `PS_Closed` | `HRR`, `HPI`, `HRI`, `HOL`, `HCL` | [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `tyndp_stores`: `h2_cavern`, `h2_tank` | `H2 cavern-storage`, `H2 tank-storage` | [#931](https://github.com/open-energy-transition/open-tyndp/pull/931) |
| none | `tyndp_stores`: `home battery` | [#892](https://github.com/open-energy-transition/open-tyndp/pull/892) |
| PEMMDB `Other Non-RES Oil shale` | `Other Non-RES Shale oil` | [#865](https://github.com/open-energy-transition/open-tyndp/pull/865) |
| none | PEMMDB `Other Non-RES Hydrogen` | [#865](https://github.com/open-energy-transition/open-tyndp/pull/865) |

Columns of `data/tyndp_technology_map.csv`:

| Before | Now | PR |
|---|---|---|
| `pemmdb_hydro_inflows` | `pecd_hydro_inflows` | [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `tyndp_output_carrier`, `tyndp_vp_carrier` | `tyndp_dashboard_carrier` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `benchmarking_generation_profiles` | removed | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| `investment_dataset_carrier` | removed | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| none | `benchmarking_synfuel_supply` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |

### Carriers

| Before | Now | PR |
|---|---|---|
| `AC_OH`, `DC_OH`, `H2_OH`, `H2 pipeline OH` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| none | bus carrier `AC_SRES`, link `SRES grid connection` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| none | loads `electricity prosumer`, `electricity prosumer btm` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| one extendable two-way `electricity distribution grid` link | two fixed one-way `electricity distribution grid` links with wheeling charges | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| `home battery charger`, `home battery discharger` (PyPSA-Eur) | same carriers, now for residential batteries at the `RETE` bus | [#892](https://github.com/open-energy-transition/open-tyndp/pull/892) |
| H2 import links | generators `H2 import Pipeline`, `H2 import LH2` with hourly profiles | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| none | EV fleets `EV Street Fixed`, `EV Street Commuter`, `EV Street Non-Commuter`, `EV Prosumer Fixed`, `EV Prosumer Commuter`, `EV Prosumer Non-Commuter`; links `BEV charger`, `V2G` | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| none | `e-liquids`, `sng`, `e-liquids import`, `sng import`, `H2 to e-liquids`, `H2 to sng`, `e-liquids to oil`, `sng to gas`, `biogas to gas CC` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |

### Bus names

| Before | Now | PR |
|---|---|---|
| offshore hubs, e.g. `BEOH001`, `DEOH002` | 2026 offshore nodes, e.g. `BEO1_OFF`, `DEKF_OFF` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| H2 buses `{country} H2 Z1`, `{country} H2 Z2`, or `{country} H2` without zoning | 2026 node codes, e.g. `DEh2Z1` (Z1), `DEh2` (Z2), `FRh2N`, plus a `category` tag | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| prosumer bus `{node} low voltage` | `{node}RETE` | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| none | SRES buses with carrier `AC_SRES` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |

### Data bundles and versions

| Before | Now | PR |
|---|---|---|
| none | dataset `tyndp_2026`, version `2026-11-11` | [#1148](https://github.com/open-energy-transition/open-tyndp/pull/1148) |
| `tyndp_pecd` `3.1+pre-built.0.1` | `4.2.1+pre-built.0.1` | [#1148](https://github.com/open-energy-transition/open-tyndp/pull/1148) |
| `Hydro-Inflows.zip` | `inputs/Hydro-Inflow.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816), [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `PEMMDB2.zip` | `inputs/PEMMDB_2.X.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Demand-Profiles.zip` | `inputs/Demand.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `EV-Modelling-Inputs.zip` | `inputs/Electric-Vehicle.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Hybrid-Heat-Pump-Modelling-Inputs.zip` | `inputs/Heat.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Line-data.zip`, `Nodes.zip`, `Hydrogen.zip` | `inputs/Line-data.zip`, `inputs/Nodes.zip`, `inputs/Hydrogen.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| none | `outputs/Supply-Tool.zip`, `outputs/NT+_TimeSeriesDashboard.zip`, `additional/TYNDP_2026_Scenario_Report_Data_Figures.xlsx` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| none | `cba_corrections/Hydrogen.zip`, only for the P2G shares | [#1055](https://github.com/open-energy-transition/open-tyndp/pull/1055) |
| `Investment-Datasets.zip`, `StartingGrid2030.xlsx.zip` | no longer retrieved | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| `tyndp_nuclear_profiles` | no longer retrieved | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `tyndp_vis_plfm` | no longer retrieved | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `countries_centroids` | no longer retrieved | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |

From the 2024 dataset, `retrieve_tyndp` still downloads `Nodes.zip`, `Offshore-hubs.zip`, `20240518-Supply-Tool.xlsm.zip`, the 2024 Scenarios Package and the 2024 market model outputs ([#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051)).

### Rules, wildcards and files

Upstream wildcard and network file changes, such as `{planning_horizons}` to `{horizon}` and `base_s_all___{horizon}.nc` to `solved_{horizon}.nc`, are listed in the upstream [Migration Guide](migration.md#files).

| Before | Now | PR |
|---|---|---|
| none | wildcard `{demand_type}` in `build_tyndp_demand` | [#808](https://github.com/open-energy-transition/open-tyndp/pull/808) |
| `build_pemmdb_and_trajectories` | `build_tyndp_pemmdb_datas` | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| `build_pemmdb_data` | `build_tyndp_pemmdb_data` | [#865](https://github.com/open-energy-transition/open-tyndp/pull/865) |
| `clean_pecd_data` | `clean_tyndp_pecd_data` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887) |
| `prepare_pecd_release` | `prepare_tyndp_pecd_release` | [#1106](https://github.com/open-energy-transition/open-tyndp/pull/1106) |
| `clean_tyndp_output_benchmark` | `clean_tyndp_dashboard_benchmark` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `clean_tyndp_electricity_demand`, `build_electricity_demand_tyndp` | `build_tyndp_demand` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| `build_tyndp_h2_demand`, `build_tyndp_h2_demands` | `build_tyndp_demand` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| `clean_tyndp_h2_imports` | merged into `build_tyndp_h2_imports` | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| none | `retrieve_tyndp_2026`, `build_tyndp_electricity_ntc`, `build_tyndp_wheeling_charges`, `build_tyndp_ev_fleets`, `build_tyndp_nuclear_profiles`, `build_tyndp_synfuels`, `build_tyndp_p2g_shares` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816), [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#935](https://github.com/open-energy-transition/open-tyndp/pull/935), [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039), [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073), [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098), [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
| `build_tyndp_offshore_hubs`, `plot_offshore_network`, `plot_base_offshore_network`, `build_tyndp_transmission_projects` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `build_tyndp_trajectories` | removed | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| `retrieve_tyndp_nuclear_profiles` | removed | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `retrieve_tyndp_vp_data`, `clean_tyndp_vp_data` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `retrieve_countries_centroids` | removed | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| `electricity_demand_{planning_horizons}.csv` | `demand_tyndp_electricity_market_{horizon}.csv` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| none | `demand_tyndp_electricity_prosumer_{horizon}.csv`, `demand_tyndp_electricity_prosumer_btm_{horizon}.csv` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| `h2_demand_tyndp_{planning_horizons}.csv` | `demand_tyndp_h2_z1_{horizon}.csv`, `demand_tyndp_h2_z2_{horizon}.csv` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `profile_pemmdb_hydro.nc` | `pecd_hydro_profile.nc` | [#1106](https://github.com/open-energy-transition/open-tyndp/pull/1106) |
| `tyndp_trajectories.csv`, `offshore_zone_trajectories.csv` | removed | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| none | `tyndp_electricity_ntc_{horizon}.csv`, `wheeling_charges_tyndp.csv` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| none | `ev_fleets_tyndp_{horizon}.csv`, `ev_demand_tyndp_{horizon}.csv` | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| none | `h2_import_profiles_{horizon}.csv`, `p2g_shares_{horizon}.csv` | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939), [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
| none | `synfuel_generators_tyndp_{horizon}.csv`, `synfuel_links_tyndp_{horizon}.csv` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| `nuclear_p_max_pu_{planning_horizons}.csv` (retrieved) | `tyndp_nuclear_profiles/nuclear_p_max_pu_{horizon}.csv` (built) | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `config/test/config.cyears.yaml`, task `tyndp-cyears-test` | `config/test/config.wscenarios.yaml`, task `tyndp-wscenarios-test` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |

### Benchmarking and CBA names

| Before | Now | PR |
|---|---|---|
| benchmark source `market_out` | `dashboard_out` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| benchmark source `vp` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| table `generation_profiles` | removed | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| none | tables `electricity_prosumer_demand`, `ev_charging`, `synfuel_supply` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960), [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039), [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| output folder `benchmarks/tyndp-2024/` | `benchmarks/tyndp-2026/` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `benchmarks_tyndp_output_{scenario}{planning_horizons}.csv` | `benchmarks_tyndp_dashboard_{scenario}{horizon}.csv` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `benchmarks_tyndp_output_crossborder_{scenario}{planning_horizons}.csv`, `benchmarks_tyndp_output_h2_demand_{scenario}{planning_horizons}.csv`, `benchmarks_tyndp_output_elec_demand_{scenario}{planning_horizons}.csv` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| CSV suffix `_cy{snapshot start year}` | `_ws{snapshot start year}` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| PDF suffix `_cy{climate year}` | `_ws{weather scenario index}` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| CBA runs `NT-cy1995`, `NT-cy2008`, `NT-cy2009`, `NT-cyears` | `NT-ws1995`, `NT-ws2008`, `NT-ws2009`, `NT-wscenarios` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `cba.projects` includes `t335` | `t225` | [#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079) |
| `cba.projects` includes `t1106` | `t1106` temporarily removed | [#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079) |

### New features

| Feature | PR |
|---|---|
| Planning horizons 2035 and 2050 | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#1041](https://github.com/open-energy-transition/open-tyndp/pull/1041) |
| Weather scenarios per planning horizon (`wscenarios_tyndp`) | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| PECD v4.2.1 renewable profiles, including CSP | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#1148](https://github.com/open-energy-transition/open-tyndp/pull/1148) |
| SRES and DRES generation with own buses feeding electrolysers | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| Market, prosumer and BTM electricity demand | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| Prosumer buses and distribution grid with wheeling charges | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| Home batteries at the prosumer bus | [#892](https://github.com/open-energy-transition/open-tyndp/pull/892) |
| Per-horizon NTCs and reference grids | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Hourly H2 import profiles in two price bands | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| H2 storage by flexibility class | [#931](https://github.com/open-energy-transition/open-tyndp/pull/931) |
| SMR ramp limits (`smr_ramp_limits_tyndp`) | [#931](https://github.com/open-energy-transition/open-tyndp/pull/931) |
| P2G capacity shares between H2 zones | [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
| New demand types in `build_tyndp_demand` | [#808](https://github.com/open-energy-transition/open-tyndp/pull/808) |
| EV fleets with smart charging and V2G | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| Synthetic fuels: domestic production and import of e-liquids and SNG | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| Nuclear availability profiles from the 2026 Time Series Dashboards | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| New countries `MD`, `TR`, `UA` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| HEV and LEV price scenarios | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066), [#1092](https://github.com/open-energy-transition/open-tyndp/pull/1092) |
| Benchmarking against TYNDP 2026 data | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |

### Removed features

| Before | PR |
|---|---|
| Offshore hubs with offshore electrolysers and hub wind technologies | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Investment-candidate corrections and capacity trajectories | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| Patching of demand with Market Model outputs | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944), [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| DE and GA scenarios | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| ENTSO-E H2 starting grid | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Toggle for H2 zoning (`h2_zones_tyndp`) | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| CBA reference grid corrections (now placeholders) | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Visualisation Platform benchmarking | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| PyPSA-Eur land transport in TYNDP runs (replaced by a flat EU27 oil demand) | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |

### Known limitations

| Limitation |
|---|
| Benchmarking is disabled for HEV and LEV. These scenarios differ from NT+ only in fuel and CO2 prices. |
| The CBA run names `NT-ws1995`, `NT-ws2008`, `NT-ws2009` and `NT-wscenarios` keep the 2024 calendar years. They are not TYNDP 2026 weather scenarios. They change only the snapshots and cutout years, and all use the same 2026 weather scenario profiles. |
| The CBA still uses the TYNDP 2024 project list and bus corrections with 2024 offshore hub buses. |
| **TODO (placeholder, depends on [#1125](https://github.com/open-energy-transition/open-tyndp/pull/1125)):** describe how to use pre-solved Scenario Building networks of releases before v0.10, which use the old network file names, with the CBA. |
