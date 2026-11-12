# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,


rule process_costs:
    input:
        expand_runs(resources("costs_{horizon}_processed.csv"), runs=BASE_RUNS),


rule cluster_networks:
    input:
        expand_runs(
            [resources("networks/clustered.nc"), resources("busmap.csv")],
            runs=BASE_RUNS,
        ),
    message:
        "Collecting clustered network files"


rule compose_networks:
    input:
        expand_runs(resources("networks/composed_{horizon}.nc")),
    message:
        "Collecting composed network files"


rule solve_networks:
    input:
        expand_runs(RESULTS + "networks/solved_{horizon}.nc", horizons=last_horizon),
    message:
        "Collecting solved network files"


rule solve_operations_networks:
    input:
        expand_runs(
            RESULTS + "networks/operations_{horizon}.nc",
            horizons=(
                last_horizon
                if config["foresight"] == "perfect"
                else config_provider("planning_horizons")
            ),
        ),
    message:
        "Collecting operational dispatch network files"


def balance_map_paths(kind, w):
    """
    kind = "static" or "interactive"
    """
    cfg_key = "balance_map" if kind == "static" else "balance_map_interactive"
    ext = "pdf" if kind == "static" else "html"

    if config["foresight"] == "perfect":
        return []

    return expand_runs(
        RESULTS + f"maps/{kind}/balance_map_{{carrier}}_{{horizon}}.{ext}",
        carrier=config_provider("plotting", cfg_key, "bus_carriers")(w),
    )


rule plot_balance_maps:
    input:
        static=lambda w: balance_map_paths("static", w),
        interactive=lambda w: balance_map_paths("interactive", w),
    message:
        "Plotting energy balance maps"


rule plot_balance_maps_static:
    input:
        lambda w: balance_map_paths("static", w),


rule plot_balance_maps_interactive:
    input:
        lambda w: balance_map_paths("interactive", w),


rule plot_power_networks:
    input:
        (
            expand_runs(resources("maps/clustered_network.pdf"), runs=BASE_RUNS)
            if config["foresight"] != "perfect"
            else []
        ),
    message:
        "Plotting clustered power network topology"


rule rulegraphs:
    input:
        expand_runs(resources("dag_rulegraph.pdf"), runs=BASE_RUNS),


rule filegraphs:
    input:
        expand_runs(resources("dag_filegraph.pdf"), runs=BASE_RUNS),
