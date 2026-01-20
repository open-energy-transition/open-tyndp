..
  SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Benchmarking
##########################################

The Open-TYNDP introduces a benchmarking framework for continuous and systematic validation of Open TYNDP model outputs against TYNDP 2024 scenarios. This framework provides flexible and scalable validation across multiple metrics and benchmarking methods.

Introduction
------------

The following metrics from the `TYNDP 2024 Scenarios report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Report_FInal_Version_250128_web.pdf>`_ are considered relevant for benchmarking:

* Exogenous Inputs:

  * Benchmark Final Energy demand by fuel, EU27 (TWh), (Fig 5, p24 and Fig 51, p63)
  * Benchmark Electricity demand per sector, EU27 (TWh), (Fig 6, p25 and Fig 52, p63)
  * Benchmark Methane demand by sector, EU27 (TWh), (Fig 8, p27 and Fig 53, p64)
  * Benchmark Hydrogen demand by sector, EU27 (TWh), (Fig 10, p28 and Fig 54, p64)

* Investment and dispatch modelling outputs:

  * Benchmark of net installed capacity for electricity generation, EU27 (GW), (Fig 25, p39 and Fig 55, p65)
  * Benchmark of electricity generation, EU27 (TWh), (Fig 26, p39 and Fig 56, p65)
  * Benchmark methane supply, EU27 (TWh), (Fig 32, p45 and Fig 57, p66)
  * Benchmark hydrogen supply, EU27 (TWh), (Fig 33, p46 and Fig 58, p67)
  * Benchmark biomass supply, EU27 (TWh), (Fig 59, p67)
  * Benchmark energy imports, EU27 (TWh), (Fig 40, p51 and Fig 60, p68)
  * Hourly generation profile of power generation, Fig 30, p35

The data is published in the `Scenarios package <https://2024-data.entsos-tyndp-scenarios.eu/files/reports/TYNDP-2024-Scenarios-Package-20250128.zip>`_. In addition to the TYNDP 2024 Scenarios Report data, data from the `Visualisation Platform <https://2024.entsos-tyndp-scenarios.eu/visualisation-platform/>`_ is also processed and included in the relevant figures. The Open-TYNDP data is benchmarked against the published TYNDP 2024 Scenarios Report data, while figures present both data sources for comparison purposes. This addition resolves discrepancies observed in the main data source.

The benchmarking is based on the methodology proposed by `Wen et al. (2022) <https://www.sciencedirect.com/science/article/pii/S0306261922011667>`_. This methodology provides a multi-criteria approach to ensure:

- the **diversity** (each indicator has its own added value),
- the **effectiveness** (each indicator provides essential and correct information),
- the **robustness** (against diverse units and orders of magnitude), and
- the **compatibility** (can be used to compare across countries) of the selected set of indicators.

This methodology defines the following indicators:

- **Missing**: Count of carriers / sectors dropped due to missing values
- **sMPE** (Symmetric Mean Percentage Error): Indicates the direction of the deviation between modeled scenarios and TYNDP 2024 outcomes, showing if the output is overall overestimated or underestimated.
- **sMAPE** (Symmetric Mean Absolute Percentage Error): Indicates the absolute magnitude of the deviations, avoiding the cancellation of negative and positive errors.
- **sMdAPE** (Symmetric Median Absolute Percentage Error): Provides skewness information to complement sMAPE.
- **RMSLE** (Root Mean Square Logarithmic Error): Complements the percentage errors since it shows the logarithmic deviation values.
- **Growth error**: Shows the error on the temporal scale. This indicator is ignored for dynamic time series (i.e., hourly generation profiles).

Hourly time series from the TYNDP 2024 will be aggregated to match the temporal resolution of Open-TYNDP.

Summary tables are computed for both the overall and per-carrier results.

Workflow
--------

#. New configuration files `config/benchmarking.default.yaml`.
#. `retrieve_additional_tyndp_data`: Retrieve the TYNDP 2024 Scenarios Report Data Figures package for benchmarking purposes. This rule will be deprecated once the data bundle has been updated (https://github.com/open-energy-transition/open-tyndp/issues/87).
#. (new) `clean_tyndp_report_benchmark`: Read and process the raw TYNDP 2024 Scenarios Report data. The output data structure is a long-format table.
#. (new) `clean_tyndp_vp_data`: Read and process the TYNDP 2024 Visualisation Platform data for benchmarking purposes. The output data structure is a long-format table.
#. (new) `build_statistics`: Compute the benchmark statistics from the optimised network. Run for every planning horizon. The output data structure is a long-format table.
#. (new) `make_benchmark`: Compute accuracy indicators for comparing model results against reference data from TYNDP 2024 Scenarios Report data.
#. (new) `make_benchmarks` to collect `make_benchmark` outputs
#. (new) `plot_benchmark`: Generate visualisation outputs for model validation.
#. (new) `plot_benchmarks` to collect `plot_benchmarks` outputs
#. The full set of files produced for the benchmarking are stored in the `results/validation/` folder. This includes:

   * `results/validation/resources/` for processed inputs information from both Open-TYNDP and TYNDP 2024.
   * `results/validation/csvs_s_{clusters}_{opts}_{sector_opts}_all_years/` for quantitative information for each table
   * `results/validation/graphics_s_{clusters}_{opts}_{sector_opts}_all_years/` for figures of each table
   * `results/validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.csv` as summary table
   * `results/validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.pdf` as summary figure
   * the structure of these outputs can be validated in the artifacts of the GitHub CI (e.g. artifacts section `here <https://github.com/open-energy-transition/open-tyndp/actions/runs/17715799690?pr=73>`_)

.. image:: img/tyndp/benchmarking_workflow.png

Outputs
-------

.. warning::
    Open-TYNDP is under active development and is not yet feature-complete. The current `development status <https://open-tyndp.readthedocs.io/en/latest/index.html#development-status>`__ and the general `Limitations <https://open-tyndp.readthedocs.io/en/latest/limitations.html>`__ are important to understand before using the model. The following outputs are presented for illustrative purposes and do not reflect the quality of the results.

Example of indicators extracted from `power_generation_s_all__all_years.csv` for NT scenario with hourly resolution:

================================  =======  =======  ========  =======  ================  ==================
Carrier                           sMPE     sMAPE    sMdAPE    RMSLE    Growth Error      version
================================  =======  =======  ========  =======  ================  ==================
**Coal + other fossil**           0.77     0.77     0.77      0.97     -0.09             v0.3+g8776987e7
**Hydro (exc. pump storage)**     0.14     0.14     0.14      0.14     0.01              v0.3+g8776987e7
**Hydrogen**                      -0.16    0.16     0.16      0.22     0.03              v0.3+g8776987e7
**Methane**                       -0.28    0.28     0.28      0.28     -0.01             v0.3+g8776987e7
**Nuclear**                       0.32     0.32     0.32      0.32     0                 v0.3+g8776987e7
**Oil**                           -0.54    0.54     0.54      0.67     -0.07             v0.3+g8776987e7
**Solar**                         -0.09    0.09     0.09      0.09     0                 v0.3+g8776987e7
**Wind offshore**                 -0.06    0.06     0.06      0.07     -0.01             v0.3+g8776987e7
**Wind onshore**                  -0.09    0.09     0.09      0.1      -0.01             v0.3+g8776987e7
**Other res**                     —        —        —         —        —                 v0.3+g8776987e7
**Other non-res**                 —        —        —         —        —                 v0.3+g8776987e7
**Biofuels**                      —        —        —         —        —                 v0.3+g8776987e7
**OCGT**                          —        —        —         —        —                 v0.3+g8776987e7
================================  =======  =======  ========  =======  ================  ==================

Example of figure created for the final energy demand for NT scenario in 2030 with hourly resolution:

.. image:: img/tyndp/benchmarking_fed_NT_2030.png

Example of figure including Visualisation Platform data created for the power capacity for NT scenario in 2030 with hourly resolution:

.. image:: img/tyndp/benchmarking_power_capacity_NT_2030_w_vis_pltfm.png

Example of figure created for the generation profiles for DE scenario in 2040 with 45SEG:

.. image:: img/tyndp/benchmarking_gen_profiles_DE_2040.png

Example of indicators extracted from `kpis_eu27_s_all__all_years.csv` for NT scenario with hourly resolution:

===========================  =======  =======  ========  =======  ================  =========  ==================
Metric                       sMPE     sMAPE    sMdAPE    RMSLE    Growth Error      Missing    version
===========================  =======  =======  ========  =======  ================  =========  ==================
biomass_supply               0.12     0.12     0.09      0.17     0                 1          v0.3+g8776987e7
elec_demand                  0        0        0         0        0                 0          v0.3+g8776987e7
energy_imports               0.39     0.39     0.2       0.73     0.01              1          v0.3+g8776987e7
final_energy_demand          0.07     0.09     0.08      0.12     0.01              1          v0.3+g8776987e7
generation_profiles          —        —        —         —        —                 NA         v0.3+g8776987e7
hydrogen_demand              -0.79    0.84     0.45      4.03     -0.03             0          v0.3+g8776987e7
hydrogen_supply              -0.23    0.36     0.22      0.74     -0.02             2          v0.3+g8776987e7
methane_demand               -0.12    0.22     0.16      0.27     -0.01             0          v0.3+g8776987e7
methane_supply               -0.58    0.73     0.15      2.65     0                 3          v0.3+g8776987e7
power_capacity               0        0.08     0.02      0.17     -0.01             5          v0.3+g8776987e7
power_generation             0        0.27     0.19      0.43     -0.01             4          v0.3+g8776987e7
Total (excl. time series)    -0.07    0.38     0.14      1.58     0                 14         v0.3+g8776987e7
===========================  =======  =======  ========  =======  ================  =========  ==================

Example of summary figure created for NT scenario:

.. image:: img/tyndp/benchmarking_overview_NT.png
