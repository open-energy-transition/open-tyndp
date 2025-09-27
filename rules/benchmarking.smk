# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT


rule clean_tyndp_benchmark:
    params:
        benchmarking=config_provider("benchmarking"),
        scenario=config_provider("tyndp_scenario"),
    input:
        scenarios_figures="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
    output:
        benchmarks=RESULTS + "validation/resources/benchmarks_tyndp.csv",
    log:
        logs("clean_tyndp_benchmark.log"),
    benchmark:
        benchmarks("clean_tyndp_benchmark")
    threads: 4
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/clean_tyndp_benchmark.py"


rule build_statistics:
    params:
        benchmarking=config_provider("benchmarking"),
        scenario=config_provider("tyndp_scenario"),
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        loss_factors="data/tyndp_electricity_loss_factors.csv",
    output:
        RESULTS
        + "validation/resources/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    log:
        logs(
            "build_statistics_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "build_statistics_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    threads: 4
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_statistics.py"


rule make_benchmark:
    params:
        benchmarking=config_provider("benchmarking"),
        scenario=config_provider("tyndp_scenario"),
        snapshots=config_provider("snapshots"),
    input:
        results=expand(
            RESULTS
            + "validation/resources/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
            planning_horizons=config_provider("scenario", "planning_horizons"),
            allow_missing=True,
        ),
        benchmarks=RESULTS + "validation/resources/benchmarks_tyndp.csv",
    output:
        benchmarks=directory(
            RESULTS + "validation/csvs_s_{clusters}_{opts}_{sector_opts}_all_years/"
        ),
        kpis=RESULTS
        + "validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.csv",
    threads: 4
    resources:
        mem_mb=8000,
    log:
        logs("make_benchmark_s_{clusters}_{opts}_{sector_opts}_all_years.log"),
    benchmark:
        benchmarks("make_benchmark_s_{clusters}_{opts}_{sector_opts}_all_years")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_benchmark.py"


rule plot_benchmark:
    params:
        benchmarking=config_provider("benchmarking"),
        scenario=config_provider("tyndp_scenario"),
        snapshots=config_provider("snapshots"),
        colors=config_provider("plotting", "tech_colors"),
    input:
        results=expand(
            RESULTS
            + "validation/resources/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
            planning_horizons=config_provider("scenario", "planning_horizons"),
            allow_missing=True,
        ),
        benchmarks=RESULTS + "validation/resources/benchmarks_tyndp.csv",
        kpis=RESULTS
        + "validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.csv",
    output:
        dir=directory(
            RESULTS
            + "validation/graphics_s_{clusters}_{opts}_{sector_opts}_all_years/"
        ),
        kpis=RESULTS
        + "validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.pdf",
    threads: 4
    resources:
        mem_mb=8000,
    log:
        logs("plot_benchmark_s_{clusters}_{opts}_{sector_opts}_all_years.log"),
    benchmark:
        benchmarks("plot_benchmark_s_{clusters}_{opts}_{sector_opts}_all_years")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_benchmark.py"
