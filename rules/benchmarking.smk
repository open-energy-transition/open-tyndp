# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT


rule retrieve_tyndp_benchmark:
    output:
        scenarios_figures="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
    log:
        logs("retrieve_tyndp_benchmark.log"),
    retries: 2
    script:
        "../scripts/retrieve_tyndp_benchmark.py"


rule clean_tyndp_benchmark:
    params:
        benchmarking=config_provider("benchmarking"),
    input:
        scenarios_figures="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
    output:
        benchmarks=resources("benchmarks.csv"),
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


rule build_benchmark:
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
    output:
        RESULTS
        + "benchmarks/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    log:
        logs(
            "build_benchmark_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "build_benchmark_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_benchmark.py"


rule make_benchmark:
    input:
        results=RESULTS
        + "benchmarks/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
        benchmarks=resources("benchmarks.csv"),
    output:
        RESULTS
        + "benchmarks/csvs/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=8000,
    log:
        logs(
            "make_benchmark/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "make_benchmark/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_benchmark.py"


rule plot_benchmark:
    params:
        plotting=config_provider("plotting"),
    input:
        elec_demand=RESULTS
        + "benchmarks/csvs/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
    output:
        RESULTS
        + "benchmarks/graphics/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.pdf",
    threads: 1
    resources:
        mem_mb=8000,
    log:
        logs(
            "plot_benchmark/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "benchmarks/plot_benchmark/{table}_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_benchmark.py"
