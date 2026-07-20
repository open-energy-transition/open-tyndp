<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Cost-Benefit Analysis (CBA)

The Cost-Benefit Analysis (CBA) workflow is implemented in `rules/cba.smk`. Rules are organised into stages: `Retrieve`, `Build MSV`, `Build rolling horizon`, `Postprocess` and `Benchmark`. In addition, a `Collect` stage provides convenient target rules that run the workflow up to a specific stage, aggregating each stage's output across all configured wildcards.

## Retrieve

### Rule `retrieve_tyndp_cba_projects`

Downloads the CBA project explorer dataset, containing project definitions.

### Rule `retrieve_tyndp_cba_non_co2_emissions`

Downloads non-CO₂ emission factors (NOx, SO₂, PM, etc.) used to compute B4 indicators.

### Rule `retrieve_cba_guidelines_reference_projects`

Downloads the CBA Implementation Guidelines reference project table
(`table_B1_CBA_Implementations_Guidelines_TYNDP2024.csv`) used to reconcile SB investments
with CBA project definitions.

### Rule `retrieve_presolved_sb_networks`

Downloads pre-solved SB networks from a previous Open-TYNDP release for use as CBA inputs
when `cba.cba_scenario_input.use_presolved` is `true`.

## Build MSV

### Rule `clean_projects` *(checkpoint)*

::: clean_projects

### Rule `clean_tyndp_indicators`

::: clean_tyndp_indicators

### Rule `simplify_sb_network`

::: simplify_sb_network

### Rule `fix_reference_sb_to_cba`

::: fix_reference_sb_to_cba

### Rule `prepare_reference`

::: prepare_reference

### Rule `build_msv_snapshot_weightings`

::: build_msv_snapshot_weightings

## Build rolling horizon

### Rule `prepare_rolling_horizon`

::: prepare_rolling_horizon

### Rule `prepare_project`

::: prepare_project

## Solve

### Rule `solve_cba_msv_extraction`

::: solve_cba_msv_extraction

### Rule `solve_cba_reference_network`

Solves the reference network using the rolling horizon approach. Shares the script with
[`solve_cba_network`](#rule-solve_cba_network).

::: solve_cba_network

### Rule `solve_cba_network`

Solves the CBA network using the rolling horizon approach. Shares the script with [`solve_cba_reference_network`](#rule-solve_cba_reference_network).

::: solve_cba_network

## Postprocess

### Rule `make_indicators`

::: make_indicators

### Rule `combine_indicators`

::: combine_indicators

### Rule `plot_indicators`

::: plot_indicators

## Benchmark

### Rule `plot_cba_benchmark`

Plots per-project indicator benchmarking charts comparing computed indicators against
TYNDP 2024 reference values. Shares the script with [`plot_weather_benchmark`](#rule-plot_weather_benchmark).

::: plot_benchmark_indicators

### Rule `plot_weather_benchmark`

Plots weather ensemble benchmarking charts from indicators aggregated across climate years.
Shares the script with [`plot_cba_benchmark`](#rule-plot_cba_benchmark).

::: plot_benchmark_indicators

### Rule `average_indicators_per_project_and_planning_horizon`

::: average_indicators

### Rule `summarize_indicators_per_project`

::: summarize_indicators

### Rule `summarize_all_indicators`

::: summarize_all

## Collect

Aggregate rules that run the corresponding base rule across all configured wildcards. They do not have dedicated scripts.

### Rule `prepare_references`

Aggregate [`prepare_reference`](#rule-prepare_reference) outputs.

### Rule `collect_cba_scenario`

Collects all per-scenario outputs (indicator plots, benchmark charts) into a single target
for a single climate year run (e.g. `NT-cy2009`).

### Rule `cba`

Top-level target rule. Collects ensemble outputs from all climate year runs in a collection
scenario (e.g. `NT-cyears`) and the per-scenario results from nested runs.
