# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT


from scripts._helpers import safe_pyear
from shutil import unpack_archive, rmtree, copy2


# Retrieve
##########

if (CUTOUT_ADDITIONAL_DATASET := dataset_version("cutout_additional"))["source"] in [
    "archive"
]:

    rule retrieve_additional_cutout:
        input:
            storage(CUTOUT_ADDITIONAL_DATASET["url"] + "{cutout}.nc"),
        output:
            CUTOUT_DATASET["folder"] + "/{cutout}.nc",
        log:
            "logs/retrieve_additional_cutout/{cutout}.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            copy2(input[0], output[0])

    ruleorder: retrieve_additional_cutout > retrieve_cutout


if (PECD_DATASET := dataset_version("tyndp_pecd"))["source"] in ["archive"]:

    rule retrieve_tyndp_pecd:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(
                PECD_DATASET["url"] + f"PECD_{PECD_DATASET["version"]}.zip"
            ),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(PECD_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_pecd.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


if (TYNDP_HYDRO_INFLOWS_DATASET := dataset_version("tyndp_hydro_inflows"))[
    "source"
] in ["archive"]:

    rule retrieve_tyndp_hydro_inflows:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(TYNDP_HYDRO_INFLOWS_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(TYNDP_HYDRO_INFLOWS_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_hydro_inflows.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


if (PEMMDB_DATASET := dataset_version("tyndp_pemmdb"))["source"] in ["archive"]:

    rule retrieve_tyndp_pemmdb_data:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(PEMMDB_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(PEMMDB_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_pemmdb_data.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


if (SUPPLY_TOOL_DATASET := dataset_version("tyndp_supply_tool"))["source"] in [
    "archive"
]:

    rule retrieve_tyndp_supply_tool:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(SUPPLY_TOOL_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(SUPPLY_TOOL_DATASET["folder"]),
            file=f"{SUPPLY_TOOL_DATASET["folder"]}/20240518-Supply-Tool.xlsm",
        log:
            "logs/retrieve_tyndp_supply_tool.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")

            # Remove __MACOSX directory if it exists
            macosx_dir = f"{output["dir"]}/__MACOSX"
            rmtree(macosx_dir, ignore_errors=True)


if (BENCHMARK_DATASET := dataset_version("tyndp_benchmark"))["source"] in ["archive"]:

    rule retrieve_tyndp_benchmark:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            file=storage(BENCHMARK_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(BENCHMARK_DATASET["folder"]),
            file=f"{BENCHMARK_DATASET["folder"]}/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
        log:
            "logs/retrieve_tyndp_benchmark.log",
        run:
            copy2(input["file"], output["file"])


if (VIS_PLFM_DATASET := dataset_version("tyndp_vis_plfm"))["source"] in ["archive"]:

    rule retrieve_tyndp_vp_data:
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(VIS_PLFM_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(VIS_PLFM_DATASET["folder"]),
            elec_demand=f"{VIS_PLFM_DATASET["folder"]}/250117_TYNDP2024Scenarios_Electricity_Demand.xlsx",
            elec_flex=f"{VIS_PLFM_DATASET["folder"]}/250117_TYNDP2024Scenarios_Electricity_Flexibility.xlsx",
            elec_supply=f"{VIS_PLFM_DATASET["folder"]}/250117_TYNDP2024Scenarios_Electricity_SupplyMix.xlsx",
        log:
            "logs/retrieve_tyndp_vp_data.log",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


# Versioning not implemented as the dataset is used only for plotting
# License - MIT - Copyright (c) 2021 Gavin Rehkemper
# Website: https://github.com/gavinr/world-countries-centroids
rule retrieve_countries_centroids:
    output:
        "data/countries_centroids.geojson",
    log:
        "logs/retrieve_countries_centroids.log",
    run:
        from scripts._helpers import progress_retrieve

        progress_retrieve(
            "https://cdn.jsdelivr.net/gh/gavinr/world-countries-centroids@v1.0.0/dist/countries.geojson",
            output[0],
            disable=True,
        )


# Development
#############
if not "pre-built" in PECD_DATASET["version"]:

    def get_pecd_prebuilt_version(increment_minor=True):
        prebuilt_prefix = f"{PECD_DATASET["version"]}+pre-built."
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

    rule prepare_pecd_release:
        params:
            cyears=config_provider(
                "electricity", "pecd_renewable_profiles", "pre_built", "cyears"
            ),
            available_pyears=config_provider(
                "electricity", "pecd_renewable_profiles", "available_years"
            ),
        input:
            pecd_raw=PECD_DATASET["folder"],
        output:
            pecd_prebuilt=directory(
                f"{PECD_DATASET["folder"]}+pre-built.{get_pecd_prebuilt_version(increment_minor= True)}"
            ),
        log:
            "logs/prepare_pecd_release.log",
        benchmark:
            "benchmarks/prepare_pecd_release"
        threads: 4
        resources:
            mem_mb=1000,
        script:
            "../scripts/sb/prepare_pecd_release.py"


# Build electricity
###################

if config["load"]["source"] == "tyndp":

    rule clean_tyndp_electricity_demand:
        params:
            planning_horizons=config_provider("scenario", "planning_horizons"),
            snapshots=config_provider("snapshots"),
            scenario=config_provider("tyndp_scenario"),
            available_years=config_provider("load", "available_years_tyndp"),
        input:
            electricity_demand=rules.retrieve_tyndp.output.demand_profiles,
        output:
            electricity_demand_prepped=resources("electricity_demand_raw_tyndp.csv"),
        log:
            logs("clean_tyndp_electricity_demand.log"),
        benchmark:
            benchmarks("clean_tyndp_electricity_demand")
        threads: 4
        resources:
            mem_mb=4000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/clean_tyndp_electricity_demand.py"


use rule build_electricity_demand as build_electricity_demand_tyndp with:
    input:
        unpack(input_elec_demand),
        reported=rules.clean_tyndp_electricity_demand.output.electricity_demand_prepped,
    output:
        resources("electricity_demand_{planning_horizons}.csv"),
    log:
        logs("build_electricity_demand_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_electricity_demand_{planning_horizons}")


def get_pecd_prebuilt(w):
    if "pre-built" in PECD_DATASET["version"]:
        return rules.retrieve_tyndp_pecd.output.dir
    else:
        return rules.prepare_pecd_release.output.pecd_prebuilt


rule clean_pecd_data:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        fill_gaps_method=config_provider(
            "electricity", "pecd_renewable_profiles", "fill_gaps_method"
        ),
        available_years=config_provider(
            "electricity", "pecd_renewable_profiles", "available_years"
        ),
        prebuilt_years=config_provider(
            "electricity", "pecd_renewable_profiles", "pre_built", "cyears"
        ),
    input:
        pecd_prebuilt=get_pecd_prebuilt,
        offshore_buses=rules.retrieve_tyndp.output.offshore_nodes,
        onshore_buses=resources("busmap_base_s_all.csv"),
    output:
        pecd_data_clean=resources("pecd_data_{technology}_{planning_horizons}.csv"),
    log:
        logs("clean_pecd_data_{technology}_{planning_horizons}.log"),
    benchmark:
        benchmarks("clean_pecd_data_{technology}_{planning_horizons}")
    threads: 4
    resources:
        mem_mb=4000,
    script:
        "../scripts/sb/clean_pecd_data.py"


def input_data_pecd(w):
    available_years = config_provider(
        "electricity", "pecd_renewable_profiles", "available_years"
    )(w)
    planning_horizons = config_provider("scenario", "planning_horizons")(w)
    safe_pyears = set(
        safe_pyear(year, available_years, "PECD", verbose=False)
        for year in planning_horizons
    )
    return {
        f"pecd_data_{pyear}": resources("pecd_data_{technology}_" + str(pyear) + ".csv")
        for pyear in safe_pyears
    }


rule build_renewable_profiles_pecd:
    params:
        planning_horizons=config_provider("scenario", "planning_horizons"),
        available_years=config_provider(
            "electricity", "pecd_renewable_profiles", "available_years"
        ),
    input:
        unpack(input_data_pecd),
    output:
        profile=resources("profile_pecd_{clusters}_{technology}.nc"),
    log:
        logs("build_renewable_profile_pecd_{clusters}_{technology}.log"),
    benchmark:
        benchmarks("build_renewable_profile_pecd_{clusters}_{technology}")
    threads: 1
    resources:
        mem_mb=4000,
    wildcard_constraints:
        technology="(?!hydro).*",  # Any technology other than hydro
    script:
        "../scripts/sb/build_renewable_profiles_pecd.py"


pemmdb_techs = branch(
    config_provider("electricity", "pemmdb_capacities", "enable"),
    config_provider("electricity", "pemmdb_capacities", "technologies"),
)


rule build_pemmdb_data:
    params:
        pemmdb_techs=pemmdb_techs,
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        available_years=config_provider(
            "electricity", "pemmdb_capacities", "available_years"
        ),
        tyndp_scenario=config_provider("tyndp_scenario"),
    input:
        pemmdb_dir=rules.retrieve_tyndp_pemmdb_data.output.dir,
        carrier_mapping="data/tyndp_technology_map.csv",
        busmap=resources("busmap_base_s_all.csv"),
    output:
        pemmdb_capacities=resources("pemmdb_capacities_{planning_horizons}.csv"),
        pemmdb_profiles=resources("pemmdb_profiles_{planning_horizons}.nc"),
    log:
        logs("build_pemmdb_data_{planning_horizons}.log"),
    threads: config_provider("electricity", "pemmdb_capacities", "nprocesses")
    resources:
        mem_mb=16000,
    benchmark:
        benchmarks("build_pemmdb_data_{planning_horizons}")
    script:
        "../scripts/sb/build_pemmdb_data.py"


rule build_tyndp_transmission_projects:
    input:
        buses=rules.build_tyndp_network.output.substations_geojson,
        invest_grid=rules.retrieve_tyndp.output.invest_grid,
    output:
        resources("tyndp/new_links.csv"),
    log:
        logs("build_tyndp_transmission_projects.log"),
    benchmark:
        benchmarks("build_tyndp_transmission_projects")
    resources:
        mem_mb=1000,
    threads: 1
    script:
        "../scripts/sb/build_tyndp_transmission_projects.py"


rule build_tyndp_trajectories:
    params:
        tyndp_scenario=config_provider("tyndp_scenario"),
    input:
        trajectories=rules.retrieve_tyndp.output.trajectories,
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        tyndp_trajectories=resources("tyndp_trajectories.csv"),
    log:
        logs("build_tyndp_trajectories.log"),
    threads: 4
    benchmark:
        benchmarks("build_tyndp_trajectories")
    script:
        "../scripts/sb/build_tyndp_trajectories.py"


rule clean_tyndp_hydro_inflows:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        available_years=config_provider(
            "electricity", "pemmdb_hydro_profiles", "available_years"
        ),
    input:
        hydro_inflows_dir=rules.retrieve_tyndp_hydro_inflows.output.dir,
        busmap=resources("busmap_base_s_all.csv"),
    output:
        hydro_inflows_tyndp=resources(
            "hydro_inflows_tyndp_{tech}_{planning_horizons}.csv"
        ),
    log:
        logs("clean_tyndp_hydro_inflows_{tech}_{planning_horizons}.log"),
    threads: 4
    retries: 2
    benchmark:
        benchmarks("clean_tyndp_hydro_inflows_{tech}_{planning_horizons}")
    script:
        "../scripts/sb/clean_tyndp_hydro_inflows.py"


def input_data_hydro_tyndp(w):
    available_years = config_provider(
        "electricity", "pemmdb_hydro_profiles", "available_years"
    )(w)
    planning_horizons = config_provider("scenario", "planning_horizons")(w)
    safe_pyears = set(
        safe_pyear(
            year,
            available_years,
            "PEMMDB hydro",
            verbose=False,
        )
        for year in planning_horizons
    )
    technologies = config_provider(
        "electricity", "pemmdb_hydro_profiles", "technologies"
    )(w)
    return {
        f"hydro_inflow_tyndp_{tech}_{pyear}": resources(
            f"hydro_inflows_tyndp_{tech}_{str(pyear)}.csv"
        )
        for pyear in safe_pyears
        for tech in technologies
    }


rule build_tyndp_hydro_profile:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        available_years=config_provider(
            "electricity", "pemmdb_hydro_profiles", "available_years"
        ),
        technologies=config_provider(
            "electricity", "pemmdb_hydro_profiles", "technologies"
        ),
    input:
        unpack(input_data_hydro_tyndp),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        profile=resources("profile_pemmdb_hydro.nc"),
    log:
        logs("build_tyndp_hydro_profile.log"),
    benchmark:
        benchmarks("build_tyndp_hydro_profile")
    resources:
        mem_mb=5000,
    script:
        "../scripts/sb/build_tyndp_hydro_profile.py"


use rule build_electricity_demand_base as build_electricity_demand_base_tyndp with:
    input:
        unpack(input_elec_demand_base),
        load=resources("electricity_demand_{planning_horizons}.csv"),
    output:
        resources("electricity_demand_base_s_{planning_horizons}.nc"),
    log:
        logs("build_electricity_demand_base_s_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_electricity_demand_base_s_{planning_horizons}")


# Build sector
##############


rule build_tyndp_gas_demand:
    params:
        scenario=config_provider("tyndp_scenario"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
    input:
        supply_tool=rules.retrieve_tyndp_supply_tool.output.file,
    output:
        gas_demand=resources("gas_demand_tyndp_{planning_horizons}.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_tyndp_gas_demand_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_tyndp_gas_demand_{planning_horizons}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/sb/build_tyndp_gas_demand.py"


rule build_tyndp_h2_demand:
    params:
        snapshots=config_provider("snapshots"),
        scenario=config_provider("tyndp_scenario"),
    input:
        h2_demand=rules.retrieve_tyndp.output.demand_profiles,
    output:
        h2_demand=resources("h2_demand_tyndp_{planning_horizons}.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_tyndp_h2_demand_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_tyndp_h2_demand_{planning_horizons}")
    script:
        "../scripts/sb/build_tyndp_h2_demand.py"


if config["sector"]["h2_topology_tyndp"]:

    rule build_tyndp_h2_network:
        params:
            snapshots=config_provider("snapshots"),
            scenario=config_provider("tyndp_scenario"),
        input:
            h2_reference_grid=rules.retrieve_tyndp.output.h2_reference_grid,
        output:
            h2_grid_prepped=resources("h2_reference_grid_tyndp_{planning_horizons}.csv"),
            interzonal_prepped=resources("h2_interzonal_tyndp_{planning_horizons}.csv"),
        log:
            logs("build_tyndp_h2_network_{planning_horizons}.log"),
        benchmark:
            benchmarks("build_tyndp_h2_network_{planning_horizons}")
        threads: 1
        resources:
            mem_mb=4000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/build_tyndp_h2_network.py"

    rule clean_tyndp_h2_imports:
        input:
            import_potentials_raw=rules.retrieve_tyndp.output.h2_imports,
            countries_centroids=rules.retrieve_countries_centroids.output,
        output:
            import_potentials_prepped=resources("h2_import_potentials_prepped.csv"),
        log:
            logs("clean_tyndp_h2_imports.log"),
        benchmark:
            benchmarks("clean_tyndp_h2_imports")
        threads: 1
        resources:
            mem_mb=4000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/clean_tyndp_h2_imports.py"

    rule build_tyndp_h2_imports:
        params:
            scenario=config_provider("tyndp_scenario"),
        input:
            import_potentials_prepped=rules.clean_tyndp_h2_imports.output.import_potentials_prepped,
        output:
            import_potentials_filtered=resources(
                "h2_import_potentials_{planning_horizons}.csv"
            ),
        log:
            logs("build_tyndp_h2_imports_{planning_horizons}.log"),
        benchmark:
            benchmarks("build_tyndp_h2_imports_{planning_horizons}")
        threads: 1
        resources:
            mem_mb=4000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/build_tyndp_h2_imports.py"


if config["sector"]["offshore_hubs_tyndp"]["enable"]:

    rule build_tyndp_offshore_hubs:
        params:
            planning_horizons=config_provider("scenario", "planning_horizons"),
            scenario=config_provider("tyndp_scenario"),
            countries=config_provider("countries"),
            offshore_hubs_tyndp=config_provider("sector", "offshore_hubs_tyndp"),
            extendable_carriers=config_provider("electricity", "extendable_carriers"),
            h2_zones_tyndp=config_provider("sector", "h2_zones_tyndp"),
        input:
            nodes=rules.retrieve_tyndp.output.offshore_nodes,
            grid=rules.retrieve_tyndp.output.offshore_grid,
            electrolysers=rules.retrieve_tyndp.output.offshore_electrolysers,
            generators=rules.retrieve_tyndp.output.offshore_generators,
        output:
            offshore_buses=resources("offshore_buses.csv"),
            offshore_grid=resources("offshore_grid.csv"),
            offshore_electrolysers=resources("offshore_electrolysers.csv"),
            offshore_generators=resources("offshore_generators.csv"),
            offshore_zone_trajectories=resources("offshore_zone_trajectories.csv"),
        log:
            logs("build_tyndp_offshore_hubs.log"),
        benchmark:
            benchmarks("build_tyndp_offshore_hubs")
        threads: 1
        resources:
            mem_mb=4000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/build_tyndp_offshore_hubs.py"


rule group_tyndp_conventionals:
    params:
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
    input:
        pemmdb_capacities=resources("pemmdb_capacities_{planning_horizon}.csv"),
        pemmdb_profiles=resources("pemmdb_profiles_{planning_horizon}.nc"),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        pemmdb_capacities_grouped=resources(
            "pemmdb_capacities_{planning_horizon}_grouped.csv"
        ),
        pemmdb_profiles_grouped=resources(
            "pemmdb_profiles_{planning_horizon}_grouped.nc"
        ),
    log:
        logs("group_tyndp_conventionals_{planning_horizon}.log"),
    benchmark:
        benchmarks("group_tyndp_conventionals_{planning_horizon}")
    threads: 1
    resources:
        mem_mb=2000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/sb/group_tyndp_conventionals.py"


# Postprocess
#############

if config["foresight"] != "perfect":

    rule plot_base_hydrogen_network:
        params:
            plotting=config_provider("plotting"),
        input:
            network=resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
            ),
            regions_onshore=resources("regions_onshore.geojson"),
        output:
            map=resources(
                "maps/base_h2_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}.pdf"
            ),
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            benchmarks(
                "plot_base_hydrogen_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
            )
        log:
            RESULTS
            + "logs/plot_base_hydrogen_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log",
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/plot_base_hydrogen_network.py"

    rule plot_base_offshore_network:
        params:
            plotting=config_provider("plotting"),
            expanded=False,
        input:
            network=resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
            ),
            regions_offshore=resources("regions_offshore.geojson"),
        output:
            map=resources(
                "maps/base_offshore_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}.pdf"
            ),
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            benchmarks(
                "plot_base_offshore_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}"
            )
        log:
            RESULTS
            + "logs/plot_base_offshore_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}.log",
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/sb/plot_offshore_network.py"

    use rule plot_base_offshore_network as plot_offshore_network with:
        params:
            expanded=True,
        input:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        output:
            map=RESULTS
            + "maps/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-offshore_network_{carrier}.pdf",
        benchmark:
            benchmarks(
                "plot_offshore_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}"
            )
        log:
            RESULTS
            + "logs/plot_offshore_network_{clusters}_{opts}_{sector_opts}_{planning_horizons}_{carrier}.log",


# Benchmarking
##############

if config["benchmarking"]["enable"]:

    rule clean_tyndp_output_benchmark:
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
        input:
            tyndp_output_file="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Outputs/MMStandardOutputFile_NT/MMStandardOutputFile_{scenario}{planning_horizons}_Plexos_CY2009_2.5_v40.xlsx",
        output:
            benchmarks=RESULTS
            + "validation/resources/benchmarks_tyndp_output_{scenario}{planning_horizons}.csv",
        log:
            logs("clean_tyndp_output_benchmark_{scenario}{planning_horizons}.log"),
        benchmark:
            benchmarks("clean_tyndp_output_benchmark_{scenario}{planning_horizons}")
        threads: 4
        resources:
            mem_mb=8000,
        script:
            "../scripts/sb/clean_tyndp_output_benchmark.py"

    rule clean_tyndp_report_benchmark:
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
        input:
            scenarios_figures=rules.retrieve_tyndp_benchmark.output.file,
        output:
            benchmarks=RESULTS + "validation/resources/benchmarks_tyndp.csv",
        log:
            logs("clean_tyndp_report_benchmark.log"),
        benchmark:
            benchmarks("clean_tyndp_report_benchmark")
        threads: 4
        resources:
            mem_mb=8000,
        script:
            "../scripts/sb/clean_tyndp_report_benchmark.py"

    rule clean_tyndp_vp_data:
        params:
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
            unit_conversion=config_provider("benchmarking", "unit_conversion"),
        input:
            elec_demand=rules.retrieve_tyndp_vp_data.output.elec_demand,
            elec_supplymix=rules.retrieve_tyndp_vp_data.output.elec_supply,
            elec_flex=rules.retrieve_tyndp_vp_data.output.elec_flex,
        output:
            RESULTS + "validation/resources/vp_data_tyndp.csv",
        log:
            logs("clean_tyndp_vp_data.log"),
        benchmark:
            benchmarks("clean_tyndp_vp_data")
        threads: 4
        resources:
            mem_mb=8000,
        script:
            "../scripts/sb/clean_tyndp_vp_data.py"

    rule build_statistics:
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            tyndp_renewable_carriers=config_provider(
                "electricity", "tyndp_renewable_carriers"
            ),
        input:
            network=RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
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
            mem_mb=8000,
        script:
            "../scripts/sb/build_statistics.py"

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
                RESULTS
                + "validation/csvs_s_{clusters}_{opts}_{sector_opts}_all_years/"
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
        script:
            "../scripts/sb/make_benchmark.py"

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
            vp_data=RESULTS + "validation/resources/vp_data_tyndp.csv",
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
        script:
            "../scripts/sb/plot_benchmark.py"


# Collect
#########


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
            rules.build_pemmdb_data.output.pemmdb_capacities,
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
            run=config["run"]["name"],
        ),
