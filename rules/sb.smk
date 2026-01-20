# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT


import zipfile
from shutil import move
from scripts._helpers import safe_pyear

# Retrieve
##########

if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_additional_cutout:
        input:
            storage("https://storage.googleapis.com/open-tyndp-data-store/{cutout}.nc"),
        output:
            CDIR.joinpath("{cutout}.nc").as_posix(),
        log:
            Path("logs").joinpath(CDIR, "retrieve_cutout_{cutout}.log").as_posix(),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])

    ruleorder: retrieve_additional_cutout > retrieve_cutout


if config["enable"]["retrieve"]:

    rule retrieve_tyndp_bundle:
        input:
            zip_file=storage(
                "https://zenodo.org/records/14230568/files/TYNDP_2024_data_bundle.zip"
            ),
        output:
            dir=directory("data/tyndp_2024_bundle"),
            elec_reference_grid="data/tyndp_2024_bundle/Line data/ReferenceGrid_Electricity.xlsx",
            buses="data/tyndp_2024_bundle/Nodes/LIST OF NODES.xlsx",
            h2_reference_grid="data/tyndp_2024_bundle/Line data/ReferenceGrid_Hydrogen.xlsx",
            electricity_demand=directory("data/tyndp_2024_bundle/Demand Profiles"),
            h2_imports="data/tyndp_2024_bundle/Hydrogen/H2 IMPORTS GENERATORS PROPERTIES.xlsx",
            offshore_nodes="data/tyndp_2024_bundle/Offshore hubs/NODE.xlsx",
            offshore_grid="data/tyndp_2024_bundle/Offshore hubs/GRID.xlsx",
            offshore_electrolysers="data/tyndp_2024_bundle/Offshore hubs/ELECTROLYSER.xlsx",
            offshore_generators="data/tyndp_2024_bundle/Offshore hubs/GENERATOR.xlsx",
            trajectories="data/tyndp_2024_bundle/Investment Datasets/TRAJECTORY.xlsx",
        run:
            with zipfile.ZipFile(input.zip_file, "r") as zip_ref:
                zip_ref.extractall(output.dir)

    rule retrieve_tyndp_pecd_data_raw:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_version}.zip",
            source="PECD raw",
        input:
            "data/tyndp_2024_bundle",
        output:
            dir=directory("data/tyndp_2024_bundle/PECD/PECD_{pecd_version}"),
        log:
            "logs/retrieve_tyndp_pecd_data_raw_{pecd_version}.log",
        retries: 2
        wildcard_constraints:
            pecd_version="(?!.*pre-built).*",  # Cannot be pre-built version
        script:
            "../scripts/sb/retrieve_additional_tyndp_data.py"

    use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_hydro_inflows with:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/Hydro_Inflows.zip",
            source="Hydro Inflows",
        output:
            dir=directory("data/tyndp_2024_bundle/Hydro Inflows"),
        log:
            "logs/retrieve_tyndp_hydro_inflows.log",

    use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_pemmdb_data with:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/PEMMDB.zip",
            source="PEMMDB",
        output:
            dir=directory("data/tyndp_2024_bundle/PEMMDB2"),
        log:
            "logs/retrieve_tyndp_pemmdb_data.log",

    use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_supply_tool with:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/20240518-Supply-Tool.xlsm.zip",
            source="Supply Tool",
        output:
            dir=directory("data/tyndp_2024_bundle/Supply Tool"),
            file="data/tyndp_2024_bundle/Supply Tool/20240518-Supply-Tool.xlsm",
        log:
            "logs/retrieve_tyndp_supply_tool.log",

    if config["electricity"]["pecd_renewable_profiles"]["pre_built"]["retrieve"]:

        use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_pecd_data_prebuilt with:
            params:
                url="https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_prebuilt_version}.zip",
                source="PECD prebuilt",
            output:
                dir=directory(
                    "data/tyndp_2024_bundle/PECD/PECD_{pecd_prebuilt_version}"
                ),
            log:
                "logs/retrieve_tyndp_pecd_data_raw_{pecd_prebuilt_version}.log",

    use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_benchmark with:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
            source="Benchmarks",
        output:
            dir=directory("data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package"),
            file="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package/TYNDP_2024-Scenario-Report-Data-Figures_240522.xlsx",
        log:
            "logs/retrieve_tyndp_benchmark.log",

    use rule retrieve_tyndp_pecd_data_raw as retrieve_tyndp_vp_data with:
        params:
            # TODO Integrate into Zenodo tyndp data bundle
            url="https://storage.googleapis.com/open-tyndp-data-store/250117-TYNDP-2024-Visualisation-Platform.zip",
            source="Visualisation Platform",
        output:
            dir=directory("data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform"),
            elec_demand="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_Demand.xlsx",
            elec_flex="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_Flexibility.xlsx",
            elec_supply="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_SupplyMix.xlsx",
        log:
            "logs/retrieve_tyndp_vp_data.log",

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

if not config["electricity"]["pecd_renewable_profiles"]["pre_built"]["retrieve"]:

    def pecd_version(w):
        version = config_provider("electricity", "pecd_renewable_profiles", "version")(
            w
        )
        return {"pecd_raw": f"data/tyndp_2024_bundle/PECD/PECD_{version}"}

    rule prepare_pecd_release:
        params:
            cyears=config_provider(
                "electricity", "pecd_renewable_profiles", "pre_built", "cyears"
            ),
            available_pyears=config_provider(
                "electricity", "pecd_renewable_profiles", "available_years"
            ),
        input:
            unpack(pecd_version),
        output:
            pecd_prebuilt=directory(
                "data/tyndp_2024_bundle/PECD/PECD_{pecd_prebuilt_version}"
            ),
        log:
            "logs/prepare_pecd_release_{pecd_prebuilt_version}.log",
        benchmark:
            "benchmarks/prepare_pecd_release_{pecd_prebuilt_version}"
        threads: 4
        resources:
            mem_mb=1000,
        script:
            "../scripts/sb/prepare_pecd_release.py"


# Build electricity
###################


use rule build_electricity_demand as build_electricity_demand_tyndp with:
    input:
        unpack(input_elec_demand),
        reported=resources("electricity_demand_raw_tyndp.csv"),
    output:
        resources("electricity_demand_{planning_horizons}.csv"),
    log:
        logs("build_electricity_demand_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_electricity_demand_{planning_horizons}")


def pecd_prebuilt_version(w):
    pre_built_version = config_provider(
        "electricity", "pecd_renewable_profiles", "pre_built", "pecd_prebuilt_version"
    )(w)
    pecd_raw_version = config_provider(
        "electricity", "pecd_renewable_profiles", "version"
    )(w)
    return {
        "pecd_prebuilt": f"data/tyndp_2024_bundle/PECD/PECD_{pecd_raw_version}+pre-built.{pre_built_version}"
    }


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
        unpack(pecd_prebuilt_version),
        offshore_buses="data/tyndp_2024_bundle/Offshore hubs/NODE.xlsx",
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
        pemmdb_dir="data/tyndp_2024_bundle/PEMMDB2",
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


rule build_tyndp_trajectories:
    params:
        tyndp_scenario=config_provider("tyndp_scenario"),
    input:
        trajectories="data/tyndp_2024_bundle/Investment Datasets/TRAJECTORY.xlsx",
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
        hydro_inflows_dir="data/tyndp_2024_bundle/Hydro Inflows",
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


if config["load"]["source"] == "tyndp":

    rule clean_tyndp_electricity_demand:
        params:
            planning_horizons=config_provider("scenario", "planning_horizons"),
            snapshots=config_provider("snapshots"),
            scenario=config_provider("tyndp_scenario"),
            available_years=config_provider("load", "available_years_tyndp"),
        input:
            electricity_demand="data/tyndp_2024_bundle/Demand Profiles",
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


# Build sector
##############


rule build_tyndp_gas_demand:
    params:
        scenario=config_provider("tyndp_scenario"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
    input:
        supply_tool="data/tyndp_2024_bundle/Supply Tool/20240518-Supply-Tool.xlsm",
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
        h2_demand="data/tyndp_2024_bundle/Demand Profiles",
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
            tyndp_reference_grid="data/tyndp_2024_bundle/Line data/ReferenceGrid_Hydrogen.xlsx",
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
            import_potentials_raw="data/tyndp_2024_bundle/Hydrogen/H2 IMPORTS GENERATORS PROPERTIES.xlsx",
            countries_centroids="data/countries_centroids.geojson",
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
            import_potentials_prepped=resources("h2_import_potentials_prepped.csv"),
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
            nodes="data/tyndp_2024_bundle/Offshore hubs/NODE.xlsx",
            grid="data/tyndp_2024_bundle/Offshore hubs/GRID.xlsx",
            electrolysers="data/tyndp_2024_bundle/Offshore hubs/ELECTROLYSER.xlsx",
            generators="data/tyndp_2024_bundle/Offshore hubs/GENERATOR.xlsx",
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

    rule plot_base_network:
        params:
            plotting=config_provider("plotting"),
        input:
            network=resources("networks/base.nc"),
            regions_onshore=resources("regions_onshore.geojson"),
        output:
            map=resources("maps/power-network.pdf"),
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            benchmarks("plot_base_network/base")
        script:
            "../scripts/plot_base_network.py"

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

    rule clean_plexos_benchmark:
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
        input:
            plexos_file="data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Outputs/MMStandardOutputFile_NT/MMStandardOutputFile_{scenario}{planning_horizons}_Plexos_CY2009_2.5_v40.xlsx",
        output:
            benchmarks=RESULTS
            + "validation/resources/benchmarks_plexos_{scenario}{planning_horizons}.csv",
        log:
            logs("clean_plexos_benchmark_{scenario}{planning_horizons}.log"),
        benchmark:
            benchmarks("clean_plexos_benchmark_{scenario}{planning_horizons}")
        threads: 4
        resources:
            mem_mb=8000,
        script:
            "../scripts/sb/clean_plexos_benchmark.py"

    rule clean_tyndp_benchmark:
        params:
            benchmarking=config_provider("benchmarking"),
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
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
            mem_mb=8000,
        script:
            "../scripts/sb/clean_tyndp_benchmark.py"

    rule clean_tyndp_vp_data:
        params:
            scenario=config_provider("tyndp_scenario"),
            snapshots=config_provider("snapshots"),
            unit_conversion=config_provider("benchmarking", "unit_conversion"),
        input:
            elec_demand="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_Demand.xlsx",
            elec_supplymix="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_SupplyMix.xlsx",
            elec_flex="data/tyndp_2024_bundle/TYNDP-2024-Visualisation-Platform/250117_TYNDP2024Scenarios_Electricity_Flexibility.xlsx",
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
            run=config["run"]["name"],
        ),
