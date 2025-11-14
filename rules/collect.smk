# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
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


rule process_costs:
    input:
        lambda w: (
            expand(
                resources(
                    f"costs_{config_provider('costs', 'year')(w)}_processed.csv"
                ),
                run=config["run"]["name"],
            )
            if config_provider("foresight")(w) == "overnight"
            else expand(
                resources("costs_{planning_horizons}_processed.csv"),
                **config["scenario"],
                run=config["run"]["name"],
            )
        ),


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
        expand(
            RESULTS
            + "validation/resources/benchmarks_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.csv",
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS + "validation/resources/benchmarks_tyndp.csv",
            run=config["run"]["name"],
        ),
        expand(
            RESULTS + "validation/resources/vp_data_tyndp.csv",
            run=config["run"]["name"],
        ),


rule make_benchmarks:
    input:
        expand(
            RESULTS
            + "validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.csv",
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule plot_benchmarks:
    input:
        expand(
            RESULTS
            + "validation/kpis_eu27_s_{clusters}_{opts}_{sector_opts}_all_years.pdf",
            **config["scenario"],
            run=config["run"]["name"],
        ),


def input_pemmdb_datas(w):
    available_years = config_provider(
        "electricity", "pemmdb_capacities", "available_years"
    )(w)
    return list(
        {
            safe_pyear(year, available_years, verbose=False)
            for year in config_provider("scenario", "planning_horizons")(w)
        }
    )


rule build_pemmdb_and_trajectories:
    input:
        expand(
            resources("pemmdb_capacities_{planning_horizons}.csv"),
            planning_horizons=input_pemmdb_datas,
            run=config["run"]["name"],
        ),
        expand(
            resources("tyndp_trajectories.csv"),
            run=config["run"]["name"],
        ),


rule build_tyndp_h2_demands:
    input:
        expand(
            resources("h2_demand_tyndp_{planning_horizons}.csv"),
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule build_tyndp_gas_demands:
    input:
        expand(
            resources("gas_demand_tyndp_{planning_horizons}.csv"),
            **config["scenario"],

          
rule rulegraphs:
    input:
        expand(
            resources("dag_rulegraph.pdf"),
            run=config["run"]["name"],
        ),


rule filegraphs:
    input:
        expand(
            resources("dag_filegraph.pdf"),
            run=config["run"]["name"],
        ),
