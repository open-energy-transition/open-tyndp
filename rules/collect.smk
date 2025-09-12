# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


localrules:
    all,
    cluster_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,


rule cluster_networks:
    input:
        expand(
            resources("networks/base_s_{clusters}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_elec_networks:
    input:
        expand(
            resources("networks/base_s_{clusters}_elec_{opts}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule prepare_sector_networks:
    input:
        expand(
            resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "maps/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule plot_balance_maps:
    input:
        lambda w: expand(
            (
                RESULTS
                + "maps/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-balance_map_{carrier}.pdf"
            ),
            **config["scenario"],
            run=config["run"]["name"],
            carrier=config_provider("plotting", "balance_map", "bus_carriers")(w),
        ),


rule plot_power_networks_clustered:
    input:
        expand(
            resources("maps/power-network-s-{clusters}.pdf"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule clean_pecd_datas:
    input:
        lambda w: expand(
            resources("pecd_data_{technology}_{planning_horizons}.csv"),
            **config["scenario"],
            run=config["run"]["name"],
            technology=config_provider(
                "electricity", "pecd_renewable_profiles", "technologies"
            )(w),
        ),


rule build_renewable_profiles_pecds:
    input:
        lambda w: expand(
            resources("profile_pecd_{clusters}_{technology}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
            technology=config_provider(
                "electricity", "pecd_renewable_profiles", "technologies"
            )(w),
        ),


rule prepare_benchmarks:
    input:
        lambda w: expand(
            RESULTS
            + "validation/resources/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        RESULTS + "validation/resources/benchmarks_tyndp.csv",


rule make_benchmarks:
    input:
        lambda w: expand(
            RESULTS
            + "validation/kpis_s_{clusters}_{opts}_{sector_opts}_all_years.csv",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule plot_benchmarks:
    input:
        lambda w: expand(
            RESULTS
            + "validation/graphics_s_{clusters}_{opts}_{sector_opts}_all_years/",
            **config["scenario"],
            run=config["run"]["name"],
        ),
