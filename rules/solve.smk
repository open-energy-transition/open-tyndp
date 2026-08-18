# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule solve_network:
    input:
        network=resources("networks/composed_{horizon}.nc"),
        offshore_zone_trajectories=branch(
            config_provider("sector", "offshore_hubs_tyndp", "enable"),
            resources("offshore_zone_trajectories.csv"),
        ),
    output:
        network=RESULTS + "networks/solved_{horizon}.nc",
        model=(
            RESULTS + "models/solved_{horizon}.nc"
            if config["solving"]["options"]["store_model"]
            else []
        ),
    log:
        solver=normpath(RESULTS + "logs/solve_network/solver_{horizon}.log"),
        memory=RESULTS + "logs/solve_network/memory_{horizon}.log",
        python=RESULTS + "logs/solve_network/python_{horizon}.log",
    benchmark:
        (RESULTS + "benchmarks/performances/solve_network_{horizon}.log")
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("planning_horizons"),
        sector=config_provider("sector"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential"
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        renewable_carriers_tyndp=config_provider(
            "electricity", "tyndp_renewable_carriers"
        ),
    script:
        "../scripts/solve_network.py"


rule solve_operations_network:
    input:
        network=RESULTS + "networks/solved_{horizon}.nc",
    output:
        network=RESULTS + "networks/operations_{horizon}.nc",
    log:
        solver=normpath(RESULTS + "logs/solve_operations_network/solver_{horizon}.log"),
        memory=RESULTS + "logs/solve_operations_network/memory_{horizon}.log",
        python=RESULTS + "logs/solve_operations_network/python_{horizon}.log",
    benchmark:
        (RESULTS + "benchmarks/performances/solve_operations_network_{horizon}.log")
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("planning_horizons"),
        sector=config_provider("sector"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential"
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    script:
        "../scripts/solve_operations_network.py"
