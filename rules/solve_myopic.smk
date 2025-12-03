# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule add_existing_baseyear:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        carriers=config_provider("electricity", "renewable_carriers"),
        costs=config_provider("costs"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        load_source=config_provider("load", "source"),
        conventional=config_provider("conventional"),
        conventional_carriers=config_provider("electricity", "conventional_carriers"),
    input:
        unpack(input_conventional),
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        busmap_s=resources("busmap_base_s.csv"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        costs=lambda w: resources(
            f"costs_{config_provider("scenario", "planning_horizons",0)(w)}_processed.csv"
        ),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv"
        ),
        heating_efficiencies=resources("heating_efficiencies.csv"),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
    wildcard_constraints:
        # TODO: The first planning_horizon needs to be aligned across scenarios
        # snakemake does not support passing functions to wildcard_constraints
        # reference: https://github.com/snakemake/snakemake/issues/2703
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=3000,
    log:
        logs(
            "add_existing_baseyear_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_existing_baseyear/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/add_existing_baseyear.py"


def input_profile_tech_brownfield(w):
    return {
        f"profile_{tech}": resources("profile_{clusters}_" + tech + ".nc")
        for tech in (set(config_provider("electricity", "renewable_carriers")(w)))
        if tech != "hydro"
    }


def input_profile_tech_brownfield_pecd(w):
    return {
        f"profile_{tech}": resources("profile_pecd_{clusters}_" + tech + ".nc")
        for tech in pecd_techs(w)
    }


rule add_brownfield:
    params:
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        H2_retrofit_capacity_per_CH4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        threshold_capacity=config_provider("existing_capacities", "threshold_capacity"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        carriers=config_provider("electricity", "renewable_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        tes=config_provider("sector", "tes"),
        dynamic_ptes_capacity=config_provider(
            "sector", "district_heating", "ptes", "dynamic_capacity"
        ),
        offshore_hubs_tyndp=config_provider("sector", "offshore_hubs_tyndp", "enable"),
        carriers_tyndp=config_provider("electricity", "tyndp_renewable_carriers"),
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
        hydrogen_fuel_cell=config_provider("sector", "hydrogen_fuel_cell"),
        hydrogen_turbine=config_provider("sector", "hydrogen_turbine"),
        group_tyndp_conventionals=config_provider(
            "electricity", "group_tyndp_conventionals"
        ),
    input:
        unpack(input_profile_tech_brownfield),
        unpack(input_profile_tech_brownfield_pecd),
        simplify_busmap=resources("busmap_base_s.csv"),
        cluster_busmap=resources("busmap_base_s_{clusters}.csv"),
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs=resources("costs_{planning_horizons}_processed.csv"),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        carrier_mapping="data/tyndp_technology_map.csv",
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
    threads: 4
    resources:
        mem_mb=10000,
    log:
        logs(
            "add_brownfield_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_brownfield/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/add_brownfield.py"


ruleorder: add_existing_baseyear > add_brownfield


rule solve_sector_network_myopic:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        renewable_carriers_tyndp=config_provider(
            "electricity", "tyndp_renewable_carriers"
        ),
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
        costs=resources("costs_{planning_horizons}_processed.csv"),
        offshore_zone_trajectories=branch(
            config_provider("sector", "offshore_hubs_tyndp", "enable"),
            resources("offshore_zone_trajectories.csv"),
        ),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        shadow_config
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/solve_network.py"


rule prepare_uncertainty_scenarios:
    message:
        "Preparing different scenarios for uncertainty analysis."
    params:
        uncertainty_scenarios=config["uncertainty_scenarios"],
    input:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
    output:
        networks=[
            resources(s)
            for s in [
                uncertainty_scenario
                + "/networks"
                + "/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
                for uncertainty_scenario in config["uncertainty_scenarios"]
            ]
        ],
    log:
        logs(
            "prepare_uncertainty_scenarios_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    script:
        "../scripts/prepare_uncertainty_scenarios.py"


rule solve_uncertainty_scenarios_myopic:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        renewable_carriers_tyndp=config_provider(
            "electricity", "tyndp_renewable_carriers"
        ),
    input:
        network=resources(
            "{uncertainty_scenario}/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        offshore_zone_trajectories=branch(
            config_provider("sector", "offshore_hubs_tyndp", "enable"),
            resources("offshore_zone_trajectories.csv"),
        ),
        line_limits=resources(
            "uncertainty_scenarios/{clusters}_lluk_{sector_opts}_{planning_horizons}.csv"
        ),
    output:
        network=RESULTS
        + "{uncertainty_scenario}/networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "{uncertainty_scenario}/configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        shadow_config
    log:
        solver=RESULTS
        + "{uncertainty_scenario}/logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "{uncertainty_scenario}/logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "{uncertainty_scenario}/logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: 1
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "{uncertainty_scenario}/benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/solve_network.py"


rule extract_uncertainty_results:
    message:
        "Extracting results from uncertainty analysis scenarios for subsequent processing."
    params:
        uncertainty_scenarios=config["uncertainty_scenarios"],
    input:
        networks=[
            RESULTS
            + uncertainty_scenario
            + "/networks/base_s_{clusters}__{sector_opts}_{planning_horizons}.nc"
            for uncertainty_scenario in config["uncertainty_scenarios"]
        ],
    output:
        line_limits=resources(
            "uncertainty_scenarios/{clusters}_lluk_{sector_opts}_{planning_horizons}.csv"
        ),
    log:
        logs(
            "extract_uncertainty_results_base_s_{clusters}__{sector_opts}_{planning_horizons}.log"
        ),
    script:
        "../scripts/extract_uncertainty_results.py"


rule solve_sector_network_myopic_line_limited:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        renewable_carriers_tyndp=config_provider(
            "electricity", "tyndp_renewable_carriers"
        ),
    input:
        network=resources(
            "networks/base_s_{clusters}__{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
        offshore_zone_trajectories=branch(
            config_provider("sector", "offshore_hubs_tyndp", "enable"),
            resources("offshore_zone_trajectories.csv"),
        ),
        line_limits=resources(
            "uncertainty_scenarios/{clusters}_lluk_{sector_opts}_{planning_horizons}.csv"
        ),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        shadow_config
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network_myopic/base_s_{clusters}_lluk_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/solve_network.py"


# Need to set the ruleorder to have the `opts=lluk` wildcard version be processed
# by the correct rule, not the `rule solve_sector_network_myopic`.
ruleorder: solve_sector_network_myopic_line_limited > solve_sector_network_myopic
