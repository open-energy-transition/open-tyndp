# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

"""
Production implementation of streamlined PyPSA-EUR workflow compose rules.

This file implements the streamlined workflow structure:
base → simplified → clustered → composed → solved

All configuration is now driven by config sections rather than wildcards.
"""

from scripts._helpers import safe_pyear


def get_tyndp_compose_inputs(w):
    """Open-TYNDP additions to the compose_network inputs."""
    cfg = get_config(w)
    horizon = int(w.horizon)
    scenario = cfg["tyndp_scenario"]
    sector = cfg["sector"]
    elec = cfg["electricity"]
    h2_tyndp = sector["h2_topology_tyndp"]
    mm_years = [2030, 2040]  # only years with Market Model output data
    mm_dir = RESULTS + "benchmarks/tyndp-2024/resources/"
    mm_available = scenario == "NT" and horizon in mm_years

    inputs = {
        "carrier_mapping": "data/tyndp_technology_map.csv",
        "load_tyndp": (
            resources(f"electricity_demand_clustered_{horizon}.nc")
            if cfg["load"]["source"] == "tyndp"
            else []
        ),
        "h2_grid_tyndp": (
            resources(f"h2_reference_grid_tyndp_{horizon}.csv") if h2_tyndp else []
        ),
        "interzonal_prepped": (
            resources(f"h2_interzonal_tyndp_{horizon}.csv") if h2_tyndp else []
        ),
        "buses_h2": (
            resources("tyndp/build/geojson/buses_h2.geojson") if h2_tyndp else []
        ),
        "h2_imports_tyndp": (
            resources(f"h2_import_potentials_{horizon}.csv") if h2_tyndp else []
        ),
        "tyndp_smr": resources(f"smr_data_prepped_{horizon}.csv") if h2_tyndp else [],
        "tyndp_h2_storages": (
            resources(f"h2_storages_prepped_{horizon}.csv") if h2_tyndp else []
        ),
        "profile_pemmdb_hydro": (
            resources("pemmdb_hydro_profile.nc")
            if elec["pemmdb_hydro_profiles"]["enable"]
            else []
        ),
        "gas_demand": resources(f"gas_demand_tyndp_{horizon}.csv") if scenario else [],
        "tyndp_trajectories": (
            resources("tyndp_trajectories.csv")
            if elec["tyndp_renewable_carriers"]
            or "uranium" in elec["tyndp_conventional_carriers"]
            else []
        ),
        "tyndp_projects": (
            resources(f"tyndp/new_links_{horizon}.csv")
            if horizon in (cfg["tyndp_investment_candidates"]["elec_projects"] or {})
            else []
        ),
        "tyndp_projects_fix": (
            resources(f"cba/reference_sb_to_cba_{horizon}.csv")
            if cfg["tyndp_investment_candidates"]["patch_sb_with_annexe"]
            else []
        ),
        "tyndp_nuclear_profiles": (
            rules.retrieve_tyndp_nuclear_profiles.output[
                f"nuclear_p_max_pu_{safe_pyear(horizon, available_years=mm_years, verbose=False)}"
            ]
            if scenario and cfg["conventional"]["tyndp_availability_profiles"]
            else []
        ),
        "h2_demand": (
            f"{mm_dir}benchmarks_tyndp_output_h2_demand_{scenario}{horizon}.csv"
            if mm_available and sector["h2_demand_patch_with_mm"]
            else (resources(f"h2_demand_tyndp_{horizon}.csv") if scenario else [])
        ),
        "elec_demand_mm": (
            f"{mm_dir}benchmarks_tyndp_output_elec_demand_{scenario}{horizon}.csv"
            if mm_available and cfg["load"]["patch_demand_with_mm"]
            else []
        ),
    }

    if elec["pecd_renewable_profiles"]["enable"]:
        for tech in elec["pecd_renewable_profiles"]["technologies"]:
            inputs[f"profile_pecd_{tech}"] = resources(f"pecd_profile_{tech}.nc")

    if elec["pemmdb_capacities"]["enable"]:
        pemmdb_year = safe_pyear(
            horizon, elec["pemmdb_capacities"]["available_years"], verbose=False
        )
        grouped = "_grouped" if elec["group_tyndp_conventionals"] else ""
        inputs["pemmdb_capacities"] = resources(
            f"pemmdb_capacities_{pemmdb_year}{grouped}.csv"
        )
        inputs["pemmdb_profiles"] = resources(
            f"pemmdb_profiles_{pemmdb_year}{grouped}.nc"
        )

    if sector["offshore_hubs_tyndp"]["enable"]:
        for f in (
            "offshore_buses",
            "offshore_grid",
            "offshore_electrolysers",
            "offshore_generators",
        ):
            inputs[f] = resources(f"{f}.csv")
        inputs["tyndp_offshore_fix"] = (
            f"{mm_dir}benchmarks_tyndp_output_crossborder_{scenario}{horizon}.csv"
            if sector["offshore_hubs_tyndp"]["patch_crossborder_with_mm"]
            and horizon in mm_years
            else []
        )

    return inputs


def get_compose_inputs(w):
    """Determine inputs for compose rule based on foresight and horizon."""
    cfg = get_config(w)
    foresight = cfg["foresight"]
    horizon = int(w.horizon)
    sector_enabled = cfg["sector"]["enabled"]
    horizons = cfg["planning_horizons"]

    # Electricity-only inputs (always included)
    inputs = dict(
        **input_profile_tech(w),
        **input_class_regions(w),
        **input_conventional(w),
        tech_costs=resources(f"costs_{horizon}_processed.csv"),
        powerplants=resources("powerplants.csv"),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
        unit_commitment="data/unit_commitment.csv",
        fuel_price=(
            resources("monthly_fuel_price.csv")
            if cfg["conventional"]["dynamic_fuel_price"]
            else []
        ),
        co2_price=resources("co2_price.csv"),
        eurostat=(
            resources("eurostat_energy_balances.csv")
            if cfg["co2_budget"]["relative"]
            else []
        ),
        co2=(
            rules.retrieve_ghg_emissions.output["csv"]
            if cfg["co2_budget"]["relative"]
            else []
        ),
        load=(
            resources("electricity_demand.nc")
            if cfg["load"]["source"] != "tyndp"
            else []
        ),
        snapshot_weightings=resources("snapshot_weightings.csv"),
        network=resources("networks/clustered.nc"),
        solar_rooftop_potentials=(
            resources("solar_rooftop_potentials.csv")
            if "solar" in cfg["electricity"]["renewable_carriers"]
            else []
        ),
    )

    # Sector-specific inputs (only when sector coupling is enabled)
    if sector_enabled:
        sector_inputs = dict(
            **input_heat_source_power(w),
            clustered_gas_network=(
                resources("gas_network_clustered.csv")
                if cfg["sector"]["gas_network"] or cfg["sector"]["H2_retrofit"]
                else []
            ),
            gas_input_nodes_simplified=(
                resources("gas_input_locations_simplified.csv")
                if cfg["sector"]["gas_network"]
                or (
                    cfg["sector"]["imports"]["enable"]
                    and (
                        "gas" in cfg["sector"]["imports"]["carriers"]
                        or (
                            "H2" in cfg["sector"]["imports"]["carriers"]
                            and not cfg["sector"]["h2_topology_tyndp"]
                        )
                    )
                )
                else []
            ),
            pop_weighted_energy_totals=resources("pop_weighted_energy_totals.csv"),
            pop_weighted_heat_totals=resources("pop_weighted_heat_totals.csv"),
            shipping_demand=resources("shipping_demand.csv"),
            transport_demand=resources("transport_demand.csv"),
            transport_data=resources("transport_data.csv"),
            avail_profile=resources("avail_profile.csv"),
            dsm_profile=resources("dsm_profile.csv"),
            heat_dsm_profile=resources("residential_heat_dsm_profile.csv"),
            co2_totals_name=resources("co2_totals.csv"),
            biomass_potentials=resources("biomass_potentials_{horizon}.csv"),
            h2_cavern=resources("salt_cavern_potentials.csv"),
            clustered_pop_layout=resources("pop_layout.csv"),
            industrial_demand=resources("industrial_energy_demand_{horizon}.csv"),
            hourly_heat_demand_total=resources("hourly_heat_demand_total.nc"),
            industrial_production=resources("industrial_production_{horizon}.csv"),
            district_heat_share=resources("district_heat_share_{horizon}.csv"),
            heating_efficiencies=resources("heating_efficiencies.csv"),
            existing_heating_distribution=resources(
                "existing_heating_distribution_{horizon}.csv"
            ),
            temp_soil_total=resources("temp_soil_total.nc"),
            temp_air_total=resources("temp_air_total.nc"),
            cop_profiles=resources("cop_profiles_{horizon}.nc"),
            direct_heat_source_utilisation_profiles=resources(
                "direct_heat_source_utilisation_profiles_{horizon}.nc"
            ),
            retro_cost=(
                resources("retro_cost.csv")
                if cfg["sector"]["retrofitting"]["retro_endogen"]
                else []
            ),
            floor_area=(
                resources("floor_area.csv")
                if cfg["sector"]["retrofitting"]["retro_endogen"]
                else []
            ),
            biomass_transport_costs=(
                resources("biomass_transport_costs.csv")
                if cfg["sector"]["biomass_transport"]
                or cfg["sector"]["biomass_spatial"]
                else []
            ),
            sequestration_potential=(
                resources("co2_sequestration_potential.csv")
                if cfg["sector"]["regional_co2_sequestration_potential"]["enable"]
                else []
            ),
            ptes_e_max_pu_profiles=(
                resources("ptes_e_max_pu_profiles_{horizon}.nc")
                if cfg["sector"]["district_heating"]["ptes"]["dynamic_capacity"]
                else []
            ),
            ptes_direct_utilisation_profiles=(
                resources("ptes_direct_utilisation_profiles_{horizon}.nc")
                if cfg["sector"]["district_heating"]["ptes"]["supplemental_heating"][
                    "enable"
                ]
                else []
            ),
            solar_thermal_total=(
                resources("solar_thermal_total.nc")
                if cfg["sector"]["solar_thermal"]
                else []
            ),
            egs_potentials=(
                resources("egs_potentials.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            egs_overlap=(
                resources("egs_overlap.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            egs_capacity_factors=(
                resources("egs_capacity_factors.csv")
                if cfg["sector"]["enhanced_geothermal"]["enable"]
                else []
            ),
            ates_potentials=(
                resources("ates_potentials_{horizon}.csv")
                if cfg["sector"]["district_heating"]["ates"]["enable"]
                else []
            ),
        )
        inputs.update(sector_inputs)

    inputs.update(get_tyndp_compose_inputs(w))

    # Add brownfield inputs for non-first horizons
    if foresight == "overnight" and len(horizons) > 1:
        raise ValueError(
            "Overnight optimization can only be run for a single planning horizon."
        )

    if horizon != horizons[0]:
        # Not first horizon - need previous network
        prev_horizon = horizons[horizons.index(horizon) - 1]

        if foresight == "myopic":
            # Myopic uses solved network from previous horizon
            inputs["network_previous"] = RESULTS + f"networks/solved_{prev_horizon}.nc"
        elif foresight == "perfect":
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = resources(
                f"networks/composed_{prev_horizon}.nc"
            )
        else:
            raise ValueError(f"Invalid foresight type: {foresight}")

    return inputs


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        resources("networks/composed_{horizon}.nc"),
    log:
        logs("compose_network_{horizon}.log"),
    benchmark:
        benchmarks("performances/compose_network_{horizon}")
    threads: 1
    resources:
        mem_mb=8000,
    params:
        foresight=config_provider("foresight"),
        electricity=config_provider("electricity"),
        sector=config_provider("sector"),
        clustering=config_provider("clustering"),
        clustering_temporal=config_provider("clustering", "temporal"),
        existing_capacities=config_provider("existing_capacities"),
        pypsa_eur=config_provider("pypsa_eur"),
        renewable=config_provider("renewable"),
        conventional=config_provider("conventional"),
        costs=config_provider("costs"),
        emission_prices=config_provider("costs", "emission_prices"),
        load=config_provider("load"),
        lines=config_provider("lines"),
        links=config_provider("links"),
        transmission_losses=config_provider("solving", "options", "transmission_losses"),
        industry=config_provider("industry"),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        horizons=config_provider("planning_horizons"),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        conventional_carriers=config_provider("electricity", "conventional_carriers"),
        fuel_carriers=config_provider("existing_capacities", "conventional_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        h2_retrofit=config_provider("sector", "H2_retrofit"),
        h2_retrofit_capacity_per_ch4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        capacity_threshold=config_provider("existing_capacities", "threshold_capacity"),
        temperature_limited_stores=config_provider(
            "sector", "district_heating", "temperature_limited_stores"
        ),
        tes=config_provider("sector", "tes"),
        dynamic_ptes_capacity=config_provider(
            "sector", "district_heating", "ptes", "dynamic_capacity"
        ),
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        co2_budget=config_provider("co2_budget"),
        adjustments=config_provider("adjustments"),
        tyndp_scenario=config_provider("tyndp_scenario"),
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
        tyndp_renewable_carriers=config_provider(
            "electricity", "tyndp_renewable_carriers"
        ),
        tyndp_stores=config_provider("electricity", "tyndp_stores"),
        group_tyndp_conventionals=config_provider(
            "electricity", "group_tyndp_conventionals"
        ),
        h2_topology_tyndp=config_provider("sector", "h2_topology_tyndp"),
        offshore_hubs_tyndp=config_provider("sector", "offshore_hubs_tyndp", "enable"),
        hydrogen_fuel_cell=config_provider("sector", "hydrogen_fuel_cell"),
        hydrogen_turbine=config_provider("sector", "hydrogen_turbine"),
        patch_load_mm=config_provider("load", "patch_demand_with_mm"),
        uniform_renewable_profiles=config_provider(
            "existing_capacities", "uniform_renewable_profiles"
        ),
    message:
        "Composing network for horizon {wildcards.horizon}"
    script:
        "../scripts/compose_network.py"
