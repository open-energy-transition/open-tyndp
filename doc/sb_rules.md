<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Scenario Building (SB)

The Scenario Building (SB) workflow is implemented in `rules/sb.smk`. Rules are organised into stages: `Retrieve`, `Development`, `Build electricity`, `Build sector`, `Postprocess`, `Benchmark` and `Explore`. In addition, a `Collect` stage provides convenient target rules that run the workflow up to a specific stage, aggregating each stage's output across all configured wildcards.

## Retrieve

### Rule `retrieve_tyndp_pecd`

Downloads the PECD dataset.

### Rule `retrieve_tyndp_vp_data`

Downloads the TYNDP Visualisation Platform data used by the benchmarking framework.

### Rule `retrieve_tyndp_nuclear_profiles`

Downloads per-country nuclear availability profiles.

### Rule `retrieve_presolved_networks`

Downloads pre-solved networks from a previous Open-TYNDP release (*preliminary outcomes*
published on [Zenodo](https://doi.org/10.5281/zenodo.18608105)) and extracts the solved network
for each planning horizon. These can be investigated with PyPSA-Explorer's web interface
using the [`launch_presolved_explorer`](#rule-launch_presolved_explorer) rule without having to re-run the workflow.

**Relevant Settings**

```yaml
data:
    open_tyndp_prelim:
        source:
        version:
```

**Outputs**

- `data/open_tyndp_prelim/{source}/{version}/base_s_all___{planning_horizons}.nc`

### Rule `retrieve_countries_centroids`

Downloads country centroid geometry data by Copyright (c) 2021 Gavin Rehkemper from
<https://cdn.jsdelivr.net/gh/gavinr/world-countries-centroids@v1.0.0/dist/countries.geojson>.

**Relevant Settings**

None.

**Outputs**

- `data/countries_centroids.geojson`

## Development

### Rule `prepare_pecd_release`

::: prepare_pecd_release

## Build electricity

### Rule `clean_tyndp_electricity_demand`

::: clean_tyndp_electricity_demand

### Rule `build_electricity_demand_tyndp`

Extends the upstream [`build_electricity_demand`](preparation.md#electricity_demand) rule with TYNDP-specific load data. Builds
per-country load time series from the TYNDP electricity demand prepared by
[`clean_tyndp_electricity_demand`](#rule-clean_tyndp_electricity_demand).

### Rule `clean_pecd_data`

::: clean_pecd_data

### Rule `build_renewable_profiles_pecd`

::: build_renewable_profiles_pecd

### Rule `build_pemmdb_data`

::: build_pemmdb_data

### Rule `build_tyndp_transmission_projects`

::: build_tyndp_transmission_projects

### Rule `build_tyndp_trajectories`

::: build_tyndp_trajectories

### Rule `clean_tyndp_hydro_inflows`

::: clean_tyndp_hydro_inflows

### Rule `build_tyndp_hydro_profile`

::: build_tyndp_hydro_profile

### Rule `build_electricity_demand_base_tyndp`

Extends the upstream [`build_electricity_demand_base`](preparation.md#rule-build_electricity_demand_base) rule with TYNDP-specific load data. Builds
the electricity demand for base regions from the TYNDP electricity demand prepared by
[`build_electricity_demand_tyndp`](#rule-build_electricity_demand_tyndp).

## Build sector

### Rule `build_tyndp_gas_demand`

::: build_tyndp_gas_demand

### Rule `build_tyndp_h2_demand`

::: build_tyndp_h2_demand

### Rule `build_tyndp_h2_network`

::: build_tyndp_h2_network

### Rule `clean_tyndp_h2_imports`

::: clean_tyndp_h2_imports

### Rule `build_tyndp_h2_imports`

::: build_tyndp_h2_imports

### Rule `clean_tyndp_smr`

::: clean_tyndp_smr

### Rule `clean_tyndp_h2_storages`

::: clean_tyndp_h2_storages

### Rule `build_tyndp_offshore_hubs`

::: build_tyndp_offshore_hubs

### Rule `group_tyndp_conventionals`

::: group_tyndp_conventionals

## Postprocess

### Rule `plot_base_hydrogen_network`

::: plot_base_hydrogen_network

### Rule `plot_base_offshore_network`

::: plot_offshore_network

### Rule `plot_offshore_network`

::: plot_offshore_network

## Benchmark

### Rule `clean_tyndp_output_benchmark`

::: clean_tyndp_output_benchmark

### Rule `clean_tyndp_report_benchmark`

::: clean_tyndp_report_benchmark

### Rule `clean_tyndp_vp_data`

::: clean_tyndp_vp_data

### Rule `build_statistics`

::: build_statistics

### Rule `make_benchmark`

::: make_benchmark

### Rule `plot_benchmark`

::: plot_benchmark

## Collect

Aggregate rules that run the corresponding base rule across all configured wildcards. They do not have dedicated scripts.

### Rule `clean_pecd_datas`

Aggregate [`clean_pecd_data`](#rule-clean_pecd_data) outputs.

### Rule `build_renewable_profiles_pecds`

Aggregate [`build_renewable_profiles_pecd`](#rule-build_renewable_profiles_pecd) outputs.

### Rule `prepare_benchmarks`

Aggregate benchmark inputs before the benchmarking stage.

### Rule `make_benchmarks`

Aggregate [`make_benchmark`](#rule-make_benchmark) outputs.

### Rule `plot_benchmarks`

Aggregate [`plot_benchmark`](#rule-plot_benchmark) outputs.

### Rule `build_pemmdb_and_trajectories`

Aggregate [`build_pemmdb_data`](#rule-build_pemmdb_data) and [`build_tyndp_trajectories`](#rule-build_tyndp_trajectories) outputs.

### Rule `build_tyndp_h2_demands`

Aggregate [`build_tyndp_h2_demand`](#rule-build_tyndp_h2_demand) outputs.

### Rule `build_tyndp_gas_demands`

Aggregate [`build_tyndp_gas_demand`](#rule-build_tyndp_gas_demand) outputs.

## Explore

### Rule `launch_explorer`

::: launch_explorer

### Rule `launch_presolved_explorer`

Mirrors the [`launch_explorer`](#rule-launch_explorer) rule to launch the PyPSA-Explorer web interface with
pre-solved SB networks from previous Open-TYNDP release runs (see [`retrieve_presolved_networks`](#rule-retrieve_presolved_networks)).

### Rule `close_explorers`

Closes all open local instances of launched PyPSA-Explorers and frees up used ports again.
