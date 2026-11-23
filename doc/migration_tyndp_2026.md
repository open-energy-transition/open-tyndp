<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Migration guide from Open-TYNDP 2024 to 2026

This guide is for users who ran Open-TYNDP on TYNDP 2024 data, up to release `v0.9`. It explains what changed for the TYNDP 2026 cycle and lists all renamed, added and removed names. Read the upstream [Migration Guide](migration.md) first, as it covers the workflow restructuring of the same release.

!!! danger "Removed and renamed keys are ignored, not rejected"
    Removed or renamed config keys are silently ignored, so their settings have no effect. Old values of choice lists, such as `tyndp_scenario: DE`, fail validation at start-up. List unknown keys outside `plotting`, `conventional` and `benchmarking` with:

    ```console
    $ pixi run -e open-tyndp validate-config config/config.yaml
    ```

## Part A: What changed and why

**Data cycle.** TYNDP 2026 is a new dataset under `data: tyndp_2026:` ([#893](https://github.com/open-energy-transition/open-tyndp/pull/893)). Its version is the date on which Open-TYNDP retrieved the data from the official download page. ENTSO-E does not version these datasets and updates them now and then. Open-TYNDP therefore keeps each retrieval in its own data bucket, so results stay reproducible. The supported versions are listed in `data/tyndp_versions.csv`. PECD also moves to the release of the 2026 cycle. With the default data config, no change is needed. If you pin versions, keep both `data: tyndp:` and `data: tyndp_2026:`.

**Naming conventions.** Bus names now follow the TYNDP 2026 node codes instead of the PyPSA-Eur defaults, for example `DEh2Z1`, `BEO1_OFF` or `{node}RETE`. Carrier names keep the PyPSA-Eur convention. Update analyses that select components by bus name.

**Scenarios.** HEV and LEV replace the DE and GA scenarios ([#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066)). Following the TYNDP 2026 methodology, they add +10 % and -10 % to the NT+ fuel and CO2 prices. Replace `DE` or `GA` with `NT`, `HEV` or `LEV` in `tyndp_scenario` and your scenario files.

**Weather scenarios and full year.** Numbered weather scenarios per horizon replace the 2024 climate years, and `cyear` became `wscenario` ([#878](https://github.com/open-energy-transition/open-tyndp/pull/878)). TYNDP 2026 models the full year, so snapshots cover 8760 h instead of 52 weeks ([#1065](https://github.com/open-energy-transition/open-tyndp/pull/1065)). Set `wscenarios_tyndp` and rename `pecd_renewable_profiles.pre_built.cyears` to `wscenarios`. Set `snapshots.end` to the first day of the next year. A grid-only cutout `europe-grid` without weather data replaces the ERA5 and SARAH-3 cutouts ([#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133)). It only provides the grid for the population layouts. Set `atlite.default_cutout` to `europe-grid` and `data.cutout.source` to `tyndp-archive`, and remove `atlite.default_cutout` from CBA run scenarios.

**Horizons and countries.** TYNDP 2026 provides data for 2030, 2035, 2040 and 2050, so the TYNDP config solves all four. The 2026 node list adds Moldova, Turkey and Ukraine ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). Set the top-level `planning_horizons`, as described in the [upstream guide](migration.md#scenario-block). Add `MD`, `TR` and `UA` to `countries`.

**Offshore.** 2026 offshore nodes are real substations with a real reference-grid topology ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). The 2024 offshore hubs are therefore removed, and offshore buses carry 2026 codes such as `BEO1_OFF`. Offshore wind becomes `offwind-r` and `offwind-oh` plus SRES and DRES variants ([#999](https://github.com/open-energy-transition/open-tyndp/pull/999)). Offshore buses without connections, such as `DKN7_OFF`, are now kept when the new key `electricity.remove_unconnected_components` is `false`, as in the TYNDP config ([#1135](https://github.com/open-energy-transition/open-tyndp/pull/1135)). Remove `sector.offshore_hubs_tyndp` and the old offshore technologies. Update analyses that use `...OH001` buses.

**Generation.** Capacities follow PEMMDB in all horizons without expansion, so trajectories and investment candidates are removed ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051)). SRES generation gets own buses to feed electrolysers, and DRES generation feeds the DRES buses ([#999](https://github.com/open-energy-transition/open-tyndp/pull/999)). SRES generation feeds its electrolysers first, and the surplus goes to the market node ([#1083](https://github.com/open-energy-transition/open-tyndp/pull/1083)). Nuclear profiles come from the 2026 Time Series Dashboards, and hydro inflows use PECD codes. Remove `tyndp_investment_candidates` and update `tyndp_renewable_carriers` to the 2026 names. Rename `electricity.pemmdb_hydro_profiles` to `pecd_hydro_profiles` and use the PECD codes.

**Electricity demand and distribution.** Demand splits into market, prosumer and BTM, with prosumer demand on a `RETE` bus ([#944](https://github.com/open-energy-transition/open-tyndp/pull/944), [#935](https://github.com/open-energy-transition/open-tyndp/pull/935)). The distribution grid uses wheeling charges, and batteries are no longer extendable. Keep `sector.electricity_distribution_grid` enabled, since only it adds the prosumer buses. Add `home battery` to `tyndp_stores` and do not list batteries as extendable stores. Remove `electricity distribution grid` from `sector.transmission_efficiency.enable`, as TYNDP demand includes grid losses. Remove `load.patch_demand_with_mm` and `load.available_years_tyndp`.

**Hydrogen.** The 2026 node list always splits H2 buses into Zone 1 and Zone 2 ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). Imports become hourly generators ([#939](https://github.com/open-energy-transition/open-tyndp/pull/939)). Storage follows flexibility classes and SMR can follow ramp limits ([#931](https://github.com/open-energy-transition/open-tyndp/pull/931)). Electrolysers now also sit at the SRES buses, and DRES buses move to the electricity nodes. P2G shares split the electrolyser capacity between zones ([#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130)). Remove `h2_zones_tyndp`, `h2_reference_grid_source` and `h2_demand_patch_with_mm`. Rename the stores in `tyndp_stores` and update analyses that use ` H2 Z1` or ` H2 Z2` buses.

**Heating and methane.** Hybrid heat pumps with CH4 and H2 boilers now cover the TYNDP heat demand ([#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090)). Their COP profiles come from the TYNDP 2026 data. The exogenous methane demand now comes from the 2026 Supply Tool. Keep `sector.heating` disabled, as a TYNDP run adds the hybrid heat pumps without it. The TYNDP config sets `sector.gas_distribution_grid` to `false`.

**Transport and synthetic fuels.** EV fleets and a flat EU27 oil demand replace the PyPSA-Eur land transport ([#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039)). The oil demand is flat for EU27, because TYNDP reports EU27 only. Keep `sector.transport` enabled, since fixed capacities would otherwise oversupply the system. Remove the `land_transport_*_share` and `bev_charge_efficiency` overrides. Synthetic fuels need `sector.synfuels_tyndp.enable` and `h2_topology_tyndp`, without `industry` or `methanation` ([#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098)).

**Grid.** Reference grids keep one topology in all horizons, and only the NTCs change ([#920](https://github.com/open-energy-transition/open-tyndp/pull/920)). The network is therefore built once from `electricity.tyndp_reference_year`. Hurdle costs now also apply in Scenario Building, on `DC` and `H2 pipeline` links ([#935](https://github.com/open-energy-transition/open-tyndp/pull/935)). Move `cba.hurdle_costs` to the top-level `hurdle_costs`.

**Benchmarking.** Benchmarks now compare against the 2026 Time Series Dashboards and Scenario Report figures ([#960](https://github.com/open-energy-transition/open-tyndp/pull/960)). The `generation_profiles` table served only DE and GA and is removed ([#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066)). New tables cover prosumer demand, EV charging, synthetic fuels and hybrid heating. Rename `market_out` to `dashboard_out` and remove `vp` and `generation_profiles`. Read outputs from `benchmarks/tyndp-2026/`.

**CBA.** Run names changed their prefix from `cy` to `ws` ([#878](https://github.com/open-energy-transition/open-tyndp/pull/878)). In `cba.projects`, replace `t335` with `t225` and remove `t1106`, as both used 2024 offshore hub buses ([#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079)). Pre-solved Scenario Building networks are now retrieved by release, also from releases before `v0.10` with the old network file names ([#1125](https://github.com/open-energy-transition/open-tyndp/pull/1125)). Rename runs such as `NT-cy2009` to `NT-ws2009`. Leave `cba.cba_scenario_input.sb_version` empty to use `data.open_tyndp_prelim.version`.

**Workflow.** Several TYNDP rules and files were renamed, see [Rules, wildcards and files](#rules-wildcards-and-files). One rule builds all TYNDP demand with the new `{demand_type}` wildcard ([#808](https://github.com/open-energy-transition/open-tyndp/pull/808)). The weather scenario test config and its pixi task were renamed. The `all` target now collects the TYNDP outputs of each scenario run ([#1136](https://github.com/open-energy-transition/open-tyndp/pull/1136)). Update custom Snakemake targets, scripts that read resource files and calls of the test task.

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
| `data.countries_centroids` | removed | [#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079) |
| `sector.electricity_distribution_grid_cost_factor` | removed (upstream) | [PyPSA-Eur #2299](https://github.com/PyPSA/pypsa-eur/pull/2299) |
| `electricity.pemmdb_hydro_profiles` | `electricity.pecd_hydro_profiles` | [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `electricity.pecd_renewable_profiles.pre_built.cyears` | `electricity.pecd_renewable_profiles.pre_built.wscenarios` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `cba.hurdle_costs` | `hurdle_costs` | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |

Added keys:

| Now | PR |
|---|---|
| `data.tyndp_2026` | [#893](https://github.com/open-energy-transition/open-tyndp/pull/893) |
| `wscenarios_tyndp` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `electricity.remove_unconnected_components` | [#1135](https://github.com/open-energy-transition/open-tyndp/pull/1135) |
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
| `sector.heating` must be disabled for a TYNDP scenario. |
| `atlite.default_cutout: europe-grid` requires `data.cutout.source: tyndp-archive`, and the reverse. |
| `sector.synfuels_tyndp.enable` requires `sector.h2_topology_tyndp`. |
| `sector.synfuels_tyndp.enable` cannot be combined with `sector.industry` or `sector.methanation`. |

Changed values:

| Key | Before | Now | PR |
|---|---|---|---|
| `sector.biogas_upgrading_cc` (TYNDP config) | `false` | `true` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| `sector.gas_distribution_grid` (TYNDP config) | `true` | `false` | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| `atlite.default_cutout` (TYNDP config) | `europe-2009-sarah3-era5` | `europe-grid` | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
| `atlite.default_cutout` (CBA runs in `config/scenarios.tyndp.yaml`) | `europe-1995-sarah3-era5`, `europe-2008-sarah3-era5`, `europe-2009-sarah3-era5` | removed | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
| `data.cutout.source` (TYNDP data config) | `archive` | `tyndp-archive` | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
| `cba.cba_scenario_input.sb_version` (default) | `latest` | empty, which uses `data.open_tyndp_prelim.version` | [#1125](https://github.com/open-energy-transition/open-tyndp/pull/1125) |

Many values in `config/config.tyndp.yaml` were updated to TYNDP 2026 NT+, among them CO2 prices, CO2 sequestration potential, biomass and biogas values and the shipping oil share. See the [release notes](release_notes.md) for the details.

### Scenarios, horizons and weather scenarios

| Before | Now | PR |
|---|---|---|
| `tyndp_scenario`: `NT`, `DE`, `GA` | `tyndp_scenario`: `NT`, `HEV`, `LEV` | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| climate years `1995`, `2008`, `2009` | weather scenarios in `wscenarios_tyndp`, e.g. `2030: [3, 21, 29]` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| planning horizons 2030, 2040 | 2030, 2035, 2040, 2050 | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#1041](https://github.com/open-energy-transition/open-tyndp/pull/1041) |
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
| none | `H2 Electrolysis Z1`, `H2 Electrolysis Z2`, `H2 Electrolysis SRES Z1`, `H2 Electrolysis SRES Z2`, `H2 Electrolysis DRES Z1`, `H2 Electrolysis DRES Z2` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999), [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
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
| `investment_dataset_carrier` | removed | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
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
| none | `CH4 hybrid heat`, `CH4 hybrid air heat pump`, `CH4 hybrid gas boiler`, `H2 hybrid heat`, `H2 hybrid air heat pump`, `H2 hybrid H2 boiler` | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| none | `e-liquids`, `sng`, `e-liquids import`, `sng import`, `H2 to e-liquids`, `H2 to sng`, `e-liquids to oil`, `sng to gas`, `biogas to gas CC` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |

### Bus names

| Before | Now | PR |
|---|---|---|
| offshore hubs, e.g. `BEOH001`, `DEOH002` | 2026 offshore nodes, e.g. `BEO1_OFF`, `DEKF_OFF` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| H2 buses `{country} H2 Z1`, `{country} H2 Z2`, or `{country} H2` without zoning | 2026 node codes, e.g. `DEh2Z1` (Z1), `DEh2` (Z2), `FRh2N`, plus a `category` tag | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| prosumer bus `{node} low voltage` | `{node}RETE` | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| none | SRES buses with carrier `AC_SRES` | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| DRES buses `{H2 Z2 bus} DRES` | `{node} DRES` | [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |

### Data bundles and versions

| Before | Now | PR |
|---|---|---|
| none | dataset `tyndp_2026`, versions in `data/tyndp_versions.csv` | [#893](https://github.com/open-energy-transition/open-tyndp/pull/893) |
| PECD release of the 2024 cycle | PECD release of the 2026 cycle, versions in `data/tyndp_versions.csv` | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887) |
| `Hydro-Inflows.zip` | `inputs/Hydro-Inflow.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816), [#1059](https://github.com/open-energy-transition/open-tyndp/pull/1059) |
| `PEMMDB2.zip` | `inputs/PEMMDB_2.X.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Demand-Profiles.zip` | `inputs/Demand.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `EV-Modelling-Inputs.zip` | `inputs/Electric-Vehicle.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Hybrid-Heat-Pump-Modelling-Inputs.zip` | `inputs/Heat.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `Line-data.zip`, `Nodes.zip`, `Hydrogen.zip` | `inputs/Line-data.zip`, `inputs/Nodes.zip`, `inputs/Hydrogen.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `20240518-Supply-Tool.xlsm.zip` | `outputs/Supply-Tool.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| 2024 market model output files | `outputs/NT+_TimeSeriesDashboard.zip` | [#816](https://github.com/open-energy-transition/open-tyndp/pull/816) |
| `TYNDP-2024-Scenarios-Package-20250128.zip` (benchmarking) | `additional/TYNDP_2026_Scenario_Report_Data_Figures.xlsx` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| none | `cba_corrections/Hydrogen.zip`, only for the P2G shares | [#1055](https://github.com/open-energy-transition/open-tyndp/pull/1055) |
| cutout `europe-2009-sarah3-era5` | grid-only cutout `europe-grid` from `tyndp-archive` | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
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
| none | `build_tyndp_cop_profiles`, `build_tyndp_hybrid_heating` | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| `build_tyndp_offshore_hubs`, `plot_offshore_network`, `plot_base_offshore_network`, `build_tyndp_transmission_projects` | removed | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `build_tyndp_trajectories` | removed | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| `retrieve_tyndp_nuclear_profiles` | removed | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `retrieve_tyndp_vp_data`, `clean_tyndp_vp_data` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `retrieve_presolved_sb_networks` | merged into `retrieve_presolved_networks` | [#1125](https://github.com/open-energy-transition/open-tyndp/pull/1125) |
| `retrieve_countries_centroids` | removed | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939) |
| `electricity_demand_{planning_horizons}.csv` | `demand_tyndp_electricity_market_{horizon}.csv` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| none | `demand_tyndp_electricity_prosumer_{horizon}.csv`, `demand_tyndp_electricity_prosumer_btm_{horizon}.csv` | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| `h2_demand_tyndp_{planning_horizons}.csv` | `demand_tyndp_h2_z1_{horizon}.csv`, `demand_tyndp_h2_z2_{horizon}.csv` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| `profile_pemmdb_hydro.nc` | `pecd_hydro_profile.nc` | [#1106](https://github.com/open-energy-transition/open-tyndp/pull/1106) |
| `tyndp_trajectories.csv`, `offshore_zone_trajectories.csv` | removed | [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| none | `tyndp_electricity_ntc_{horizon}.csv`, `wheeling_charges_tyndp.csv` | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| none | `ev_fleets_tyndp_{horizon}.csv`, `ev_demand_tyndp_{horizon}.csv` | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |
| none | `h2_import_profiles_{horizon}.csv`, `p2g_shares_{horizon}.csv` | [#939](https://github.com/open-energy-transition/open-tyndp/pull/939), [#1130](https://github.com/open-energy-transition/open-tyndp/pull/1130) |
| none | `cop_profiles_tyndp_{horizon}.nc`, `hybrid_heating_boilers_{horizon}.csv`, `demand_tyndp_thermal_ch4_{horizon}.csv`, `demand_tyndp_thermal_h2_{horizon}.csv` | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| none | `synfuel_generators_tyndp_{horizon}.csv`, `synfuel_links_tyndp_{horizon}.csv` | [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098) |
| `nuclear_p_max_pu_{planning_horizons}.csv` (retrieved) | `tyndp_nuclear_profiles/nuclear_p_max_pu_{horizon}.csv` (built) | [#1073](https://github.com/open-energy-transition/open-tyndp/pull/1073) |
| `config/test/config.cyears.yaml`, task `tyndp-cyears-test` | `config/test/config.wscenarios.yaml`, task `tyndp-wscenarios-test` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |

### Benchmarking and CBA names

| Before | Now | PR |
|---|---|---|
| benchmark source `market_out` | `dashboard_out` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| benchmark source `vp` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| table `generation_profiles` | removed | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| none | tables `electricity_prosumer_demand`, `ev_charging`, `synfuel_supply`, `hybrid_heating` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960), [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039), [#1098](https://github.com/open-energy-transition/open-tyndp/pull/1098), [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| output folder `benchmarks/tyndp-2024/` | `benchmarks/tyndp-2026/` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `benchmarks_tyndp_output_{scenario}{planning_horizons}.csv` | `benchmarks_tyndp_dashboard_{scenario}{horizon}.csv` | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| `benchmarks_tyndp_output_crossborder_{scenario}{planning_horizons}.csv`, `benchmarks_tyndp_output_h2_demand_{scenario}{planning_horizons}.csv`, `benchmarks_tyndp_output_elec_demand_{scenario}{planning_horizons}.csv` | removed | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| CSV suffix `_cy{snapshot start year}` | `_ws{snapshot start year}` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| PDF suffix `_cy{climate year}` | `_ws{weather scenario index}` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| CBA runs `NT-cy1995`, `NT-cy2008`, `NT-cy2009`, `NT-cyears` | `NT-ws1995`, `NT-ws2008`, `NT-ws2009`, `NT-wscenarios` | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| `cba.projects` includes `t335` | `t225` | [#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079) |
| `cba.projects` includes `t1106` | removed | [#1079](https://github.com/open-energy-transition/open-tyndp/pull/1079) |

### New features

| Feature | PR |
|---|---|
| Planning horizons 2035 and 2050 | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887), [#1041](https://github.com/open-energy-transition/open-tyndp/pull/1041) |
| Weather scenarios per planning horizon (`wscenarios_tyndp`) | [#878](https://github.com/open-energy-transition/open-tyndp/pull/878) |
| PECD renewable profiles of the 2026 cycle, including CSP | [#887](https://github.com/open-energy-transition/open-tyndp/pull/887) |
| SRES and DRES generation with own buses feeding electrolysers | [#999](https://github.com/open-energy-transition/open-tyndp/pull/999) |
| SRES generation feeds electrolysers first, the surplus goes to the market node | [#1083](https://github.com/open-energy-transition/open-tyndp/pull/1083) |
| Offshore buses without connections are kept (`electricity.remove_unconnected_components`) | [#1135](https://github.com/open-energy-transition/open-tyndp/pull/1135) |
| Hybrid heat pumps with CH4 and H2 boilers | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| Exogenous methane demand from the 2026 Supply Tool | [#1090](https://github.com/open-energy-transition/open-tyndp/pull/1090) |
| Grid-only cutout `europe-grid` without weather data | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
| Market, prosumer and BTM electricity demand | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944) |
| Prosumer buses and distribution grid with wheeling charges | [#935](https://github.com/open-energy-transition/open-tyndp/pull/935) |
| Home batteries at the prosumer bus | [#892](https://github.com/open-energy-transition/open-tyndp/pull/892) |
| Per-horizon NTCs on a fixed reference grid topology | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
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
| Version-aware retrieval of pre-solved networks, shared by Scenario Building and CBA | [#1125](https://github.com/open-energy-transition/open-tyndp/pull/1125) |
| Benchmarking against TYNDP 2026 data | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |

### Removed features

| Before | PR |
|---|---|
| Offshore hubs with offshore electrolysers and hub wind technologies | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Investment-candidate corrections and capacity trajectories | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920), [#1051](https://github.com/open-energy-transition/open-tyndp/pull/1051) |
| Patching of demand with market model outputs | [#944](https://github.com/open-energy-transition/open-tyndp/pull/944), [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| DE and GA scenarios | [#1066](https://github.com/open-energy-transition/open-tyndp/pull/1066) |
| ENTSO-E H2 starting grid | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Toggle for H2 zoning (`h2_zones_tyndp`) | [#920](https://github.com/open-energy-transition/open-tyndp/pull/920) |
| Visualisation Platform benchmarking | [#960](https://github.com/open-energy-transition/open-tyndp/pull/960) |
| ERA5 and SARAH-3 weather cutouts in TYNDP runs | [#1133](https://github.com/open-energy-transition/open-tyndp/pull/1133) |
| PyPSA-Eur land transport in TYNDP runs (replaced by a flat EU27 oil demand) | [#1039](https://github.com/open-energy-transition/open-tyndp/pull/1039) |

### Known limitations

| Limitation |
|---|
| HEV and LEV differ from NT+ only in fuel and CO2 prices. Their demand adjustments and hydrogen import prices are not yet included. |
| HEV and LEV are not benchmarked yet. Benchmarking for them is planned. |
| The CBA run names `NT-ws1995`, `NT-ws2008`, `NT-ws2009` and `NT-wscenarios` keep the 2024 calendar years. They are not TYNDP 2026 weather scenarios. They change only the snapshots, and all use the same 2026 weather scenario profiles. |
| The CBA still uses the TYNDP 2024 project list and bus corrections with 2024 offshore hub buses. Its reference grid corrections are placeholders. |
| `t1106` is temporarily removed from `cba.projects`. |
| The TYNDP 2024 node list is still used for the CBA projects. |
| The methane demand uses only the NT+ Supply Tool workbook. The HEV and LEV workbooks are not yet included, so `build_tyndp_gas_demand` fails for these scenarios. |
| H2 boilers of hybrid heat pumps use an efficiency of 1/1.08. The TYNDP 2026 Time Series Dashboards imply 1.08. |
