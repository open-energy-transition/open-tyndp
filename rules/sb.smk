# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT


from scripts.sb.build_tyndp_demand import DEMAND_TYPE_MAP
from scripts.sb._helpers import safe_planning_horizon, find_free_port
from scripts._helpers import (
    select_weather_clusters,
    weather_cluster_run,
    weather_cluster_weights,
)
from shutil import unpack_archive, copy2

# Retrieve
##########


if (PECD_DATASET := dataset_version("tyndp_pecd"))["source"] in ARCHIVE_SOURCES:

    rule retrieve_tyndp_pecd:
        input:
            zip_file=storage(
                PECD_DATASET["url"] + f"PECD_{PECD_DATASET['version']}.zip"
            ),
        output:
            dir=directory(PECD_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_pecd.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


if (PRESOLVED_NETWORKS_DATASET := dataset_version("open_tyndp_prelim"))[
    "source"
] in ARCHIVE_SOURCES:

    rule retrieve_presolved_networks:
        input:
            zip_file=storage(PRESOLVED_NETWORKS_DATASET["url"]),
        output:
            network=f"{PRESOLVED_NETWORKS_DATASET['folder']}/base_s_all___{{horizon}}.nc",
        log:
            "logs/retrieve_presolved_networks_{horizon}.log",
        run:
            from pathlib import Path
            from shutil import copyfileobj
            from zipfile import ZipFile

            target_suffix = f"networks/base_s_all___{wildcards.horizon}.nc"
            with ZipFile(input["zip_file"], "r") as zf:
                matches = [m for m in zf.namelist() if m.endswith(target_suffix)]
                if not matches:
                    raise ValueError(
                        f"Could not find '{target_suffix}' in {input['zip_file']}."
                    )
                out_path = Path(output["network"])
                out_path.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(matches[0]) as src, out_path.open("wb") as dst:
                    copyfileobj(src, dst)


# Development
#############
if not "pre-built" in PECD_DATASET["version"]:

    def get_pecd_prebuilt_version(increment_minor=True):
        prebuilt_prefix = f"{PECD_DATASET['version']}+pre-built."
        versions = (
            dataset_version("tyndp_pecd", all_versions=True)
            .query("version.str.contains(@prebuilt_prefix, regex=False)")
            .version.sort_values()
        )

        if versions.empty:
            return "0.1"

        major, minor = versions.iloc[-1].removeprefix(prebuilt_prefix).rsplit(".", 1)

        if increment_minor:
            return f"{major}.{str(int(minor)+1)}"
        else:
            return f"{str(int(major)+1)}.0"

    rule prepare_tyndp_pecd_release:
        input:
            pecd_raw=PECD_DATASET["folder"],
        output:
            pecd_prebuilt=directory(
                f"{PECD_DATASET['folder']}+pre-built.{get_pecd_prebuilt_version(increment_minor= True)}"
            ),
        log:
            "logs/prepare_tyndp_pecd_release.log",
        benchmark:
            benchmarks("performances/prepare_tyndp_pecd_release")
        threads: 4
        resources:
            mem_mb=1000,
        params:
            wscenarios=config_provider(
                "electricity",
                "pecd_renewable_profiles",
                "pre_built",
                "wscenarios",
            ),
            available_years=config_provider(
                "electricity", "pecd_renewable_profiles", "available_years"
            ),
        script:
            scripts("sb/prepare_tyndp_pecd_release.py")


# Build electricity
###################


def weather_scenario(w: Wildcards) -> int:
    """Return the single TYNDP 2026 weather scenario (climate year column index) of the run, also for data of an earlier available year."""
    horizons = config_provider("planning_horizons")(w)
    selected = select_weather_clusters(
        config_provider("weather_clusters_tyndp")(w), horizons
    ).weather_scenario.tolist()
    if len(selected) != 1:
        raise ValueError(
            f"The run with planning horizons {horizons} selects the weather scenarios "
            f"{selected}, but its weather-dependent jobs need exactly one."
        )
    return selected[0]


rule build_tyndp_demand:
    input:
        demand=rules.retrieve_tyndp_2026.output.demand_profiles,
    output:
        demand=resources("demand_tyndp_{demand_type}_{horizon}.csv"),
    log:
        logs("build_tyndp_demand_{demand_type}_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_tyndp_demand_{demand_type}_{horizon}")
    wildcard_constraints:
        # Limited to the demand types the workflow knows how to process
        demand_type="|".join(DEMAND_TYPE_MAP),
    threads: 1
    resources:
        mem_mb=4000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        wscenario=weather_scenario,
    script:
        scripts("sb/build_tyndp_demand.py")


def get_pecd_prebuilt(w):
    if "pre-built" in PECD_DATASET["version"]:
        return rules.retrieve_tyndp_pecd.output.dir
    else:
        return rules.prepare_tyndp_pecd_release.output.pecd_prebuilt


rule clean_tyndp_pecd_data:
    input:
        pecd_prebuilt=get_pecd_prebuilt,
        buses_tyndp=from_parent(rules.build_tyndp_network.output.substations),
        pemmdb_capacities=resources("pemmdb_capacities_{horizon}.csv"),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        pecd_data_clean=resources("pecd_data_{technology}_{horizon}.csv"),
    log:
        logs("clean_tyndp_pecd_data_{technology}_{horizon}.log"),
    benchmark:
        benchmarks("performances/clean_tyndp_pecd_data_{technology}_{horizon}")
    threads: 4
    resources:
        mem_mb=4000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        fill_gaps_method=config_provider(
            "electricity", "pecd_renewable_profiles", "fill_gaps_method"
        ),
        available_years=config_provider(
            "electricity", "pecd_renewable_profiles", "available_years"
        ),
        wscenario=weather_scenario,
    script:
        scripts("sb/clean_tyndp_pecd_data.py")


def input_data_pecd(w):
    available_years = config_provider(
        "electricity", "pecd_renewable_profiles", "available_years"
    )(w)
    planning_horizons = config_provider("planning_horizons")(w)
    safe_planning_horizons = set(
        safe_planning_horizon(year, available_years, "PECD", verbose=False)
        for year in planning_horizons
    )
    return {
        f"pecd_data_{planning_horizon}": resources(
            "pecd_data_{technology}_" + str(planning_horizon) + ".csv"
        )
        for planning_horizon in safe_planning_horizons
    }


rule build_renewable_profiles_pecd:
    input:
        unpack(input_data_pecd),
    output:
        profile=resources("pecd_profile_{technology}.nc"),
    log:
        logs("build_renewable_profile_pecd_{technology}.log"),
    benchmark:
        benchmarks("performances/build_renewable_profile_pecd_{technology}")
    wildcard_constraints:
        technology="(?!hydro).*",  # Any technology other than hydro
    threads: 1
    resources:
        mem_mb=4000,
    params:
        planning_horizons=config_provider("planning_horizons"),
        available_years=config_provider(
            "electricity", "pecd_renewable_profiles", "available_years"
        ),
    script:
        scripts("sb/build_renewable_profiles_pecd.py")


pemmdb_techs = branch(
    config_provider("electricity", "pemmdb_capacities", "enable"),
    config_provider("electricity", "pemmdb_capacities", "technologies"),
)


rule build_tyndp_pemmdb_data:
    input:
        pemmdb_dir=rules.retrieve_tyndp_2026.output.pemmdb,
        carrier_mapping="data/tyndp_technology_map.csv",
        buses_tyndp=from_parent(rules.build_tyndp_network.output.substations),
    output:
        pemmdb_capacities=resources("pemmdb_capacities_{horizon}.csv"),
        pemmdb_profiles=resources("pemmdb_profiles_{horizon}.nc"),
    log:
        logs("build_tyndp_pemmdb_data_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_tyndp_pemmdb_data_{horizon}")
    threads: config_provider("electricity", "pemmdb_capacities", "nprocesses")
    resources:
        mem_mb=16000,
    params:
        pemmdb_techs=pemmdb_techs,
        wscenario=weather_scenario,
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        available_years=config_provider(
            "electricity", "pemmdb_capacities", "available_years"
        ),
    script:
        scripts("sb/build_tyndp_pemmdb_data.py")


rule build_tyndp_nuclear_profiles:
    input:
        dashboard_dir=rules.retrieve_tyndp_2026.output.market_outputs,
    output:
        resources("tyndp_nuclear_profiles/nuclear_p_max_pu_{horizon}.csv"),
    log:
        logs("build_tyndp_nuclear_profiles_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_tyndp_nuclear_profiles_{horizon}")
    threads: 4
    resources:
        mem_mb=2000,
    params:
        wscenario=weather_scenario,
        snapshots=config_provider("snapshots"),
    script:
        scripts("sb/build_tyndp_nuclear_profiles.py")


rule clean_tyndp_hydro_inflows:
    input:
        hydro_inflows_dir=rules.retrieve_tyndp_2026.output.hydro_inflows,
        busmap=from_parent(resources("busmap_cluster_network.csv")),
    output:
        hydro_inflows_tyndp=resources("hydro_inflows_tyndp_{tech}_{horizon}.csv"),
    log:
        logs("clean_tyndp_hydro_inflows_{tech}_{horizon}.log"),
    benchmark:
        benchmarks("performances/clean_tyndp_hydro_inflows_{tech}_{horizon}")
    retries: 2
    threads: 4
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        available_years=config_provider(
            "electricity", "pecd_hydro_profiles", "available_years"
        ),
        wscenario=weather_scenario,
    script:
        scripts("sb/clean_tyndp_hydro_inflows.py")


def input_data_hydro_tyndp(w):
    available_years = config_provider(
        "electricity", "pecd_hydro_profiles", "available_years"
    )(w)
    planning_horizons = config_provider("planning_horizons")(w)
    safe_planning_horizons = set(
        safe_planning_horizon(
            year,
            available_years,
            "PECD hydro",
            verbose=False,
        )
        for year in planning_horizons
    )
    technologies = config_provider(
        "electricity", "pecd_hydro_profiles", "technologies"
    )(w)
    return {
        f"hydro_inflow_tyndp_{tech}_{planning_horizon}": resources(
            f"hydro_inflows_tyndp_{tech}_{str(planning_horizon)}.csv"
        )
        for planning_horizon in safe_planning_horizons
        for tech in technologies
    }


rule build_tyndp_hydro_profile:
    input:
        unpack(input_data_hydro_tyndp),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        profile=resources("pecd_hydro_profile.nc"),
    log:
        logs("build_tyndp_hydro_profile.log"),
    benchmark:
        benchmarks("performances/build_tyndp_hydro_profile")
    resources:
        mem_mb=5000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        planning_horizons=config_provider("planning_horizons"),
        available_years=config_provider(
            "electricity", "pecd_hydro_profiles", "available_years"
        ),
        technologies=config_provider(
            "electricity", "pecd_hydro_profiles", "technologies"
        ),
    script:
        scripts("sb/build_tyndp_hydro_profile.py")


use rule build_electricity_demand_base as build_electricity_demand_base_tyndp with:
    input:
        unpack(input_elec_demand_base),
        raster=[],
        gb_excel=[],
        gb_geojson=[],
        nuts3=[],
        load=resources("demand_tyndp_electricity_market_{horizon}.csv"),
    output:
        resources("electricity_demand_simplified_{horizon}.nc"),
    log:
        logs("build_electricity_demand_simplified_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_electricity_demand_simplified_{horizon}")


use rule cluster_electricity_demand as cluster_electricity_demand_tyndp with:
    input:
        load=resources("electricity_demand_simplified_{horizon}.nc"),
        busmap=from_parent(resources("busmap_cluster_network.csv")),
    output:
        resources("electricity_demand_clustered_{horizon}.nc"),
    log:
        logs("cluster_electricity_demand_{horizon}.log"),
    benchmark:
        benchmarks("performances/cluster_electricity_demand_{horizon}")


# Build sector
##############


rule build_tyndp_gas_demand:
    input:
        supply_tool=rules.retrieve_tyndp.output.supply_tool,
    output:
        gas_demand=resources("gas_demand_tyndp_{horizon}.csv"),
    log:
        logs("build_tyndp_gas_demand_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_tyndp_gas_demand_{horizon}")
    conda:
        "../envs/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    params:
        planning_horizons=config_provider("planning_horizons"),
    script:
        scripts("sb/build_tyndp_gas_demand.py")


rule build_tyndp_wheeling_charges:
    input:
        wheeling_charges=rules.retrieve_tyndp_2026.output.wheeling_charges,
        demand=rules.retrieve_tyndp_2026.output.demand_profiles,
    output:
        wheeling_charges=resources("wheeling_charges_tyndp.csv"),
    log:
        logs("build_tyndp_wheeling_charges.log"),
    benchmark:
        benchmarks("performances/build_tyndp_wheeling_charges")
    conda:
        "../envs/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("sb/build_tyndp_wheeling_charges.py")


rule build_tyndp_ev_fleets:
    input:
        ev_market=resources("demand_tyndp_ev_market_{horizon}.csv"),
        ev_prosumer=resources("demand_tyndp_ev_prosumer_{horizon}.csv"),
        ev_modelling=rules.retrieve_tyndp_2026.output.ev_modelling,
    output:
        ev_demand=resources("ev_demand_tyndp_{horizon}.csv"),
        ev_fleets=resources("ev_fleets_tyndp_{horizon}.csv"),
        avail_profile=resources("avail_profile_tyndp_{horizon}.csv"),
        dsm_profile=resources("dsm_profile_tyndp_{horizon}.csv"),
    log:
        logs("build_tyndp_ev_fleets_{horizon}.log"),
    benchmark:
        benchmarks("performances/build_tyndp_ev_fleets_{horizon}")
    conda:
        "../envs/environment.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("sb/build_tyndp_ev_fleets.py")


if config["sector"]["h2_topology_tyndp"]:

    rule build_tyndp_h2_network:
        input:
            h2_reference_grid=rules.retrieve_tyndp_2026.output.h2_reference_grid_entsos,
        output:
            h2_grid_prepped=resources("h2_reference_grid_tyndp_{horizon}.csv"),
        log:
            logs("build_tyndp_h2_network_{horizon}.log"),
        benchmark:
            benchmarks("performances/build_tyndp_h2_network_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("sb/build_tyndp_h2_network.py")

    rule build_tyndp_h2_imports:
        input:
            import_potentials_raw=rules.retrieve_tyndp_2026.output.h2_imports,
            import_profiles_raw=rules.retrieve_tyndp_2026.output.h2_import_profiles,
        output:
            import_potentials=resources("h2_import_potentials_{horizon}.csv"),
            import_profiles=resources("h2_import_profiles_{horizon}.csv"),
        log:
            logs("build_tyndp_h2_imports_{horizon}.log"),
        benchmark:
            benchmarks("performances/build_tyndp_h2_imports_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=4000,
        params:
            snapshots=config_provider("snapshots"),
            drop_leap_day=config_provider("enable", "drop_leap_day"),
        script:
            scripts("sb/build_tyndp_h2_imports.py")

    rule clean_tyndp_smr:
        input:
            smr=rules.retrieve_tyndp_2026.output.smr,
        output:
            smr_prepped=resources("smr_data_prepped_{horizon}.csv"),
        log:
            logs("clean_tyndp_smr_{horizon}.log"),
        benchmark:
            benchmarks("performances/clean_tyndp_smr_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("sb/clean_tyndp_smr.py")

    rule clean_tyndp_h2_storages:
        input:
            h2_storages=rules.retrieve_tyndp_2026.output.h2_storages,
        output:
            h2_storages_prepped=resources("h2_storages_prepped_{horizon}.csv"),
        log:
            logs("clean_tyndp_h2_storages_{horizon}.log"),
        benchmark:
            benchmarks("performances/clean_tyndp_h2_storages_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("sb/clean_tyndp_h2_storages.py")


if (
    config["sector"]["h2_topology_tyndp"]
    and config["sector"]["synfuels_tyndp"]["enable"]
):

    rule build_tyndp_synfuels:
        input:
            synfuel_lines=rules.retrieve_tyndp_2026.output.synfuel_lines,
            synfuel_generators=rules.retrieve_tyndp_2026.output.synfuel_generators,
        output:
            synfuel_links=resources("synfuel_links_tyndp_{horizon}.csv"),
            synfuel_generators=resources("synfuel_generators_tyndp_{horizon}.csv"),
        log:
            logs("build_tyndp_synfuels_{horizon}.log"),
        benchmark:
            benchmarks("performances/build_tyndp_synfuels_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=1000,
        params:
            synfuels=config_provider("sector", "synfuels_tyndp", "carriers"),
        script:
            scripts("sb/build_tyndp_synfuels.py")


rule group_tyndp_conventionals:
    input:
        pemmdb_capacities=resources("pemmdb_capacities_{horizon}.csv"),
        pemmdb_profiles=resources("pemmdb_profiles_{horizon}.nc"),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        pemmdb_capacities_grouped=resources("pemmdb_capacities_{horizon}_grouped.csv"),
        pemmdb_profiles_grouped=resources("pemmdb_profiles_{horizon}_grouped.nc"),
    log:
        logs("group_tyndp_conventionals_{horizon}.log"),
    benchmark:
        benchmarks("performances/group_tyndp_conventionals_{horizon}")
    conda:
        "../envs/environment.yaml"
    threads: 1
    resources:
        mem_mb=2000,
    params:
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
    script:
        scripts("sb/group_tyndp_conventionals.py")


# Postprocess
#############

if config["foresight"] != "perfect":

    rule plot_base_hydrogen_network:
        input:
            network=resources("networks/composed_{horizon}.nc"),
            regions_onshore=from_parent(resources("onshore_regions.geojson")),
        output:
            map=resources("maps/base_h2_network_{horizon}.pdf"),
        log:
            RESULTS + "logs/plot_base_hydrogen_network_{horizon}.log",
        benchmark:
            benchmarks("performances/plot_base_hydrogen_network_{horizon}")
        conda:
            "../envs/environment.yaml"
        threads: 1
        resources:
            mem_mb=4000,
        params:
            plotting=config_provider("plotting"),
        script:
            scripts("sb/plot_base_hydrogen_network.py")


# Benchmark
###########

if config["benchmarking"]["enable"]:

    def benchmark_inputs(w: Wildcards, fn: str) -> list[str]:
        """Return the benchmark file `fn` per planning horizon of the run, weighted over the weather-cluster runs for a parent."""
        weighted = "weighted/" if w.run in run_parents.values() else ""
        return expand(
            RESULTS + "benchmarks/tyndp-2026/resources/" + weighted + fn,
            run=w.run,
            horizon=config_provider("planning_horizons")(w),
            scenario=config_provider("tyndp_scenario")(w),
        )

    def weather_weights(w: Wildcards) -> pd.Series:
        """Return the renormalized weights of the weather scenarios of the parent run for the `horizon` wildcard."""
        return weather_cluster_weights(
            config_provider("weather_clusters_tyndp")(w), int(w.horizon)
        )

    def price_weighting(w: Wildcards) -> str | None:
        """Return the variant of the weighted prices benchmarked for a parent, None for a weather-cluster run."""
        if w.run in run_parents.values():
            return "cluster_weighted"
        return None

    def weather_label(w: Wildcards) -> str:
        """Return the weather scenario of a weather-cluster run, or `weighted` for a parent."""
        if w.run in run_parents.values():
            return "weighted"
        return str(weather_scenario(w))

    rule weight_benchmark:
        input:
            lambda w: expand(
                RESULTS + "benchmarks/tyndp-2026/resources/{name}{horizon}.csv",
                run=[
                    weather_cluster_run(w.run, int(w.horizon), ws)
                    for ws in weather_weights(w).index
                ],
                name=w.name,
                horizon=w.horizon,
            ),
        output:
            RESULTS + "benchmarks/tyndp-2026/resources/weighted/{name}{horizon}.csv",
        log:
            logs("weight_benchmark_{name}{horizon}.log"),
        benchmark:
            benchmarks("performances/weight_benchmark_{name}{horizon}")
        threads: 1
        resources:
            mem_mb=2000,
        params:
            weights=lambda w: weather_weights(w).tolist(),
        script:
            scripts("sb/weight_benchmark.py")

    rule clean_tyndp_dashboard_benchmark:
        input:
            dashboard_dir=rules.retrieve_tyndp_2026.output.market_outputs,
            carrier_mapping="data/tyndp_technology_map.csv",
        output:
            benchmarks=RESULTS
            + "benchmarks/tyndp-2026/resources/benchmarks_tyndp_dashboard_{scenario}{horizon}.csv",
        log:
            logs("clean_tyndp_dashboard_benchmark_{scenario}{horizon}.log"),
        benchmark:
            benchmarks(
                "performances/clean_tyndp_dashboard_benchmark_{scenario}{horizon}"
            )
        threads: 4
        resources:
            mem_mb=8000,
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            countries=config_provider("countries"),
            wscenario=weather_scenario,
        script:
            scripts("sb/clean_tyndp_dashboard_benchmark.py")

    rule clean_tyndp_report_benchmark:
        input:
            scenarios_figures=rules.retrieve_tyndp_2026.output.benchmark,
            carrier_mapping="data/tyndp_technology_map.csv",
        output:
            benchmarks=RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_tyndp.csv",
        log:
            logs("clean_tyndp_report_benchmark.log"),
        benchmark:
            benchmarks("performances/clean_tyndp_report_benchmark")
        threads: 4
        resources:
            mem_mb=8000,
        params:
            benchmarking=config_provider("benchmarking"),
        script:
            scripts("sb/clean_tyndp_report_benchmark.py")

    rule build_statistics:
        input:
            network=RESULTS + "networks/solved_{horizon}.nc",
            carrier_mapping="data/tyndp_technology_map.csv",
        output:
            RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_{horizon}.csv",
        log:
            python=logs("build_statistics_{horizon}.log"),
        benchmark:
            benchmarks("performances/build_statistics_{horizon}")
        threads: 1
        resources:
            mem_mb=8000,
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            load_shedding=config_provider(
                "solving", "options", "load_shedding", "carriers"
            ),
            low_voltage=config_provider("sector", "electricity_distribution_grid"),
            group_tyndp_conventionals=config_provider(
                "electricity", "group_tyndp_conventionals"
            ),
            synfuels=config_provider("sector", "synfuels_tyndp", "carriers"),
        script:
            scripts("sb/build_statistics.py")

    rule make_benchmark:
        input:
            results=lambda w: benchmark_inputs(w, "benchmarks_{horizon}.csv"),
            benchmarks=from_parent(
                RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_tyndp.csv"
            ),
            dashboard_data=lambda w: (
                benchmark_inputs(
                    w, "benchmarks_tyndp_dashboard_{scenario}{horizon}.csv"
                )
                if config_provider("tyndp_scenario")(w)
                == "NT"  # Only NT has dashboard files for now
                else []
            ),
        output:
            benchmarks=directory(RESULTS + "benchmarks/tyndp-2026/csvs_all_years/"),
            kpis_by_bus=RESULTS + "benchmarks/tyndp-2026/kpis_all_years_by_bus.csv",
            kpis_by_country=RESULTS
            + "benchmarks/tyndp-2026/kpis_all_years_by_country.csv",
        log:
            logs("make_benchmark_all_years.log"),
        benchmark:
            benchmarks("performances/make_benchmark_all_years")
        threads: 4
        resources:
            mem_mb=8000,
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
            price_weighting=price_weighting,
        script:
            scripts("sb/make_benchmark.py")

    rule plot_benchmark:
        input:
            results=lambda w: benchmark_inputs(w, "benchmarks_{horizon}.csv"),
            dashboard_data=lambda w: (
                benchmark_inputs(
                    w, "benchmarks_tyndp_dashboard_{scenario}{horizon}.csv"
                )
                if config_provider("tyndp_scenario")(w)
                == "NT"  # Only NT has dashboard files for now
                else []
            ),
            benchmarks=from_parent(
                RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_tyndp.csv"
            ),
            kpis_by_bus=RESULTS + "benchmarks/tyndp-2026/kpis_all_years_by_bus.csv",
            kpis_by_country=RESULTS
            + "benchmarks/tyndp-2026/kpis_all_years_by_country.csv",
        output:
            dir=directory(RESULTS + "benchmarks/tyndp-2026/graphics_all_years/"),
            kpis_by_bus=RESULTS + "benchmarks/tyndp-2026/kpis_all_years_by_bus.pdf",
            kpis_by_country=RESULTS
            + "benchmarks/tyndp-2026/kpis_all_years_by_country.pdf",
        log:
            logs("plot_benchmark_all_years.log"),
        benchmark:
            benchmarks("performances/plot_benchmark_all_years")
        threads: 4
        resources:
            mem_mb=8000,
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
            tech_colors=config_provider("plotting", "tech_colors"),
            bench_colors=config_provider("plotting", "benchmarking", "colors"),
            wscenario=weather_label,
            price_weighting=price_weighting,
        script:
            scripts("sb/plot_benchmark.py")


# Collect
#########


rule clean_pecd_datas:
    input:
        lambda w: expand_runs(
            resources("pecd_data_{technology}_{horizon}.csv"),
            technology=config_provider(
                "electricity", "pecd_renewable_profiles", "technologies"
            )(w),
        ),


rule build_renewable_profiles_pecds:
    input:
        lambda w: expand(
            resources("pecd_profile_{technology}.nc"),
            run=RUNS,
            technology=config_provider(
                "electricity", "pecd_renewable_profiles", "technologies"
            )(w),
        ),


rule prepare_benchmarks:
    input:
        expand_runs(
            RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_{horizon}.csv"
        ),
        expand(
            RESULTS + "benchmarks/tyndp-2026/resources/benchmarks_tyndp.csv",
            run=BASE_RUNS,
        ),


rule make_benchmarks:
    input:
        kpis_by_bus=expand(
            rules.make_benchmark.output.kpis_by_bus,
            run=RUNS + WEIGHTED_RUNS,
        ),
        kpis_by_country=expand(
            rules.make_benchmark.output.kpis_by_country,
            run=RUNS + WEIGHTED_RUNS,
        ),


rule plot_benchmarks:
    input:
        kpis_by_bus=expand(
            rules.plot_benchmark.output.kpis_by_bus,
            run=RUNS + WEIGHTED_RUNS,
        ),
        kpis_by_country=expand(
            rules.plot_benchmark.output.kpis_by_country,
            run=RUNS + WEIGHTED_RUNS,
        ),


def pemmdb_horizons(w: Wildcards) -> list[int]:
    available_years = config_provider(
        "electricity", "pemmdb_capacities", "available_years"
    )(w)
    return list(
        {
            safe_planning_horizon(year, available_years, verbose=False)
            for year in config_provider("planning_horizons")(w)
        }
    )


rule build_tyndp_pemmdb_datas:
    input:
        expand_runs(
            rules.build_tyndp_pemmdb_data.output.pemmdb_capacities,
            horizons=pemmdb_horizons,
        ),


rule build_tyndp_gas_demands:
    input:
        expand_runs(resources("gas_demand_tyndp_{horizon}.csv"), runs=BASE_RUNS),


# Explore
###########


rule launch_explorer:
    input:
        lambda w: expand_runs(
            RESULTS + "networks/solved_{horizon}.nc",
            runs=[c for c in child_runs(w.get("run")) if c in RUNS] or RUNS,
        ),
    output:
        RESULTS + "logs/explorer_launched.log",
    params:
        port=find_free_port(start_port=8050, max_attempts=50),
        launch_msg="Launching PyPSA-Explorer...",
    run:
        import platform
        import subprocess
        import sys

        output_log = str(output[0])
        input_files = list(input)
        # Define command line executable
        cmd = [
            sys.executable,
            "scripts/sb/launch_explorer.py",
            output_log,
            str(params.port),
        ] + input_files
        print(params.launch_msg)
        # Open logfile before Popen so the log exists when the subprocess validates its path
        popen_kwargs = {
            "stdout": open(output_log, "w"),
            "stderr": subprocess.STDOUT,
        }
        # Use creationflags for Windows and start_new_session for Linux/Unix
        if platform.system() == "Windows":
            popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
        else:
            popen_kwargs["start_new_session"] = True
        process = subprocess.Popen(cmd, **popen_kwargs)
        print(f"Explorer subprocess started with PID: {process.pid}")
        print(f"PyPSA-Explorer is running at http://127.0.0.1:{params.port}.")
        print(
            f"Your browser should open automatically. If not, click the link above."
        )


rule close_explorers:
    run:
        import psutil

        print("Closing all explorer instances...")
        killed_count = 0
        for proc in psutil.process_iter(["pid", "name", "cmdline"]):
            try:
                cmdline = proc.info.get("cmdline", [])
                if cmdline and "launch_explorer.py" in " ".join(cmdline):
                    proc.kill()
                    print(f"Killed explorer process (PID: {proc.info['pid']})")
                    killed_count += 1
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        if killed_count == 0:
            print("No explorer processes found running.")
        else:
            print(f"Closed {killed_count} explorer instance(s).")


if (PRESOLVED_NETWORKS_DATASET := dataset_version("open_tyndp_prelim"))[
    "source"
] in ARCHIVE_SOURCES:

    use rule launch_explorer as launch_presolved_explorer with:
        input:
            expand(
                f"{PRESOLVED_NETWORKS_DATASET['folder']}/base_s_all___{{horizon}}.nc",
                horizon=config["planning_horizons"],
            ),
        output:
            "logs/presolved_explorer_launched.log",
        params:
            launch_msg=f"Launching PyPSA-Explorer with presolved networks for release v{PRESOLVED_NETWORKS_DATASET['version']}...",
