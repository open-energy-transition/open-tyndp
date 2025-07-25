# SPDX-FileCopyrightText: Open Energy Transition gGmbH and contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule solve_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        carriers_tyndp=config_provider("electricity", "tyndp_renewable_carriers"),
    input:
        network=resources("networks/base_s_{clusters}_elec_{opts}.nc"),
        offshore_zone_trajectories=branch(
            config_provider("sector", "offshore_hubs_tyndp", "enable"),
            resources("offshore_zone_trajectories.csv"),
        ),
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
        config=RESULTS + "configs/config.base_s_{clusters}_elec_{opts}.yaml",
    log:
        solver=normpath(
            RESULTS + "logs/solve_network/base_s_{clusters}_elec_{opts}_solver.log"
        ),
        memory=RESULTS + "logs/solve_network/base_s_{clusters}_elec_{opts}_memory.log",
        python=RESULTS + "logs/solve_network/base_s_{clusters}_elec_{opts}_python.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network/base_s_{clusters}_elec_{opts}")
    threads: solver_threads
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        shadow_config
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule solve_operations_network:
    params:
        options=config_provider("solving", "options"),
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
    output:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}_op.nc",
    log:
        solver=normpath(
            RESULTS
            + "logs/solve_operations_network/base_s_{clusters}_elec_{opts}_op_solver.log"
        ),
        python=RESULTS
        + "logs/solve_operations_network/base_s_{clusters}_elec_{opts}_op_python.log",
    benchmark:
        (RESULTS + "benchmarks/solve_operations_network/base_s_{clusters}_elec_{opts}")
    threads: 4
    resources:
        mem_mb=memory,
        runtime=config_provider("solving", "runtime", default="6h"),
    shadow:
        shadow_config
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_operations_network.py"
