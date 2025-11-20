# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean TYNDP Scenario Building offshore hubs data to be used in the PyPSA-Eur workflow. Depending on the scenario, different planning years (`pyear`) are available. DE and GA are defined for 2030, 2040 and 2050. NT scenario is only defined for 2030 and 2040. All the planning years are read at once.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd

from scripts._helpers import SCENARIO_DICT, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"


def load_offshore_hubs(fn: str):
    """
    Load and process offshore hub coordinates from Excel file.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore hub data.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the offshore hub data.

        The GeoDataFrame uses the coordinate reference system defined by `GEO_CRS`.
    """
    column_names = {
        "OFFSHORE_NODE": "Bus",
        "OFFSHORE_NODE_TYPE": "type",
        "LAT": "y",
        "LON": "x",
    }

    nodes = (
        pd.read_excel(fn, sheet_name="NODE")
        .rename(columns=column_names)
        .assign(
            location=lambda x: x.Bus,
            country=lambda x: x.location.str[:2],
        )
    )

    nodes = gpd.GeoDataFrame(
        nodes, geometry=gpd.points_from_xy(nodes.x, nodes.y), crs=GEO_CRS
    )

    return nodes


def expand_all_scenario(df: pd.DataFrame, scenarios: list):
    all_mask = df["scenario"] == "All"
    all_rows = (
        df[all_mask]
        .drop(columns="scenario")
        .merge(pd.DataFrame({"scenario": list(set(scenarios))}), how="cross")
    )
    return pd.concat([df[~all_mask], all_rows], ignore_index=True)


def load_offshore_grid(
    fn: str,
    scenario: str,
    planning_horizons: list[int],
    countries: list[str],
    max_capacity: dict[str, int],
):
    """
    Load offshore grid (electricity and hydrogen) and format data.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore grid data.
    scenario : str
        Scenario identifier to filter the grid data. Must be one of the scenario
        codes: "DE" (Distributed Energy), "GA" (Global Ambition), or
        "NT" (National Trends).
    planning_horizons : list[int]
        List of planning years to include in the cost data filtering.
    countries : list[str]
        List of country codes used to clean data.
    max_capacity : dict[str, int]
        Maximum transmission capacity between two offshore hubs per carrier

    Returns
    -------
    pd.DataFrame
        DataFrame containing the merged offshore grid data.
    """
    column_names = {
        "FROM": "bus0",
        "TO": "bus1",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "MARKET": "carrier",
        "CAPACITY": "p_nom_min",
        "CAPEX": "capex",
        "OPEX": "opex",
    }

    # Load reference grid
    grid = (
        pd.read_excel(fn, sheet_name="Reference grid")
        .rename(columns=column_names)
        .replace(
            {
                "carrier": {"E": "DC_OH", "H2": "H2 pipeline OH"},
                "scenario": SCENARIO_DICT,
            }
        )
    )
    grid = expand_all_scenario(grid, SCENARIO_DICT.values()).query(
        "scenario == @scenario"
    )

    # Load costs data
    grid_costs = (
        pd.read_excel(
            fn,
            sheet_name="COST",
        )
        .rename(columns=column_names)
        .replace({"scenario": SCENARIO_DICT})
        .query("pyear in @planning_horizons and scenario == @scenario")
    )
    grid_costs[["capex", "opex"]] = grid_costs[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW
    grid_costs["carrier"] = grid_costs["carrier"].replace(
        {"E": "DC_OH", "H2": "H2 pipeline OH"}
    )

    # Merge information
    grid = grid.merge(
        grid_costs, how="outer", on=["bus0", "bus1", "pyear", "scenario", "carrier"]
    ).assign(
        p_min_pu=0,
        p_max_pu=1,
    )

    suffix = "H2 Z1" if options["h2_zones_tyndp"] else "H2"
    # Filter out radial nodes and Convert to explicit hydrogen buses
    grid = grid.query("~bus0.str.contains('OR') and ~bus1.str.contains('OR')").assign(
        bus0=lambda df: np.where(
            df.carrier == "H2 pipeline OH",
            np.where(
                df.bus0.str.contains("OH"),
                df.bus0 + " H2",
                df.bus0.str[:2] + f" {suffix}",
            ),
            df.bus0,
        ),
        bus1=lambda df: np.where(
            df.carrier == "H2 pipeline OH",
            np.where(
                df.bus1.str.contains("OH"),
                df.bus1 + " H2",
                df.bus1.str[:2] + f" {suffix}",
            ),
            df.bus1,
        ),
    )

    # Handle missing data
    # TODO Validate assumption
    grid["p_nom_extendable"] = ~grid[["capex", "opex"]].isna().any(axis=1)
    grid[["capex", "opex"]] = grid[["capex", "opex"]].fillna(0)
    grid["p_nom_min"] = grid["p_nom_min"].fillna(0)

    # Add maximum transmission capacities
    grid["p_nom_max"] = np.where(
        grid.bus0.str.contains("OH") & grid.bus1.str.contains("OH"),
        np.where(
            (grid.carrier == "DC_OH"),
            max_capacity["DC_OH"],
            max_capacity["H2 pipeline OH"],
        )
        * 1e3,
        grid.get("p_nom_max", np.inf),
    )  # GW > MW

    # Rename UK in GB
    grid[["bus0", "bus1"]] = grid[["bus0", "bus1"]].replace("UK", "GB", regex=True)

    # Filter selected countries and nodes
    grid = grid.assign(
        country0=lambda x: x.bus0.str[:2],
        country1=lambda x: x.bus1.str[:2],
    ).query("country0 in @countries and country1 in @countries")

    return grid


def load_offshore_electrolysers(
    fn: str,
    scenario: str,
    planning_horizons: list[int],
    countries: list[str],
    options: dict,
):
    """
    Load offshore electrolysers data and format data.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore electrolyser data.
    scenario : str
        Scenario identifier to filter the grid data. Must be one of the scenario
        codes: "DE" (Distributed Energy), "GA" (Global Ambition), or
        "NT" (National Trends).
    planning_horizons : list[int]
        List of planning years to include in the cost data filtering.
    countries : list[str]
        List of country codes used to clean data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the formatted offshore electrolyser data.
    """
    column_names = {
        "NODE": "bus0",
        "OFFSHORE_NODE_TYPE": "type",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "CAPEX": "capex",
        "OPEX": "opex",
    }

    # Load electrolysers data
    electrolysers = (
        pd.read_excel(
            fn,
            sheet_name="COST",
        )
        .rename(columns=column_names)
        .query("pyear in @planning_horizons")
        .replace({"scenario": SCENARIO_DICT})
        .query("scenario == @scenario")
        .assign(country=lambda x: x.bus0.str[:2], bus1=lambda x: x.bus0 + " H2")
        .drop(columns="OFFSHORE_NODE")
    )

    suffix = "H2 Z1" if options["h2_zones_tyndp"] else "H2"
    mask = electrolysers["type"] == "Radial"
    electrolysers.loc[mask, "bus1"] = electrolysers.loc[mask, "country"] + f" {suffix}"

    electrolysers[["capex", "opex"]] = electrolysers[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW

    # rename UK in GB
    electrolysers[["bus0", "bus1", "country"]] = electrolysers[
        ["bus0", "bus1", "country"]
    ].replace("UK", "GB", regex=True)

    # filter selected countries
    electrolysers = electrolysers.query("country in @countries")

    return electrolysers


def collect_from_layer(generators_e, generators_l, nodes):
    """
    Combine existing capacities with potentials and resolve bus allocations.

    This function merges generator data from two sources: existing capacities (EXISTING sheet)
    and potential capacities (LAYER sheet). It handles reallocation of radial wind farms by
    correcting inconsistencies between EXISTING and LAYER data, particularly for offshore
    radial connections that need to be mapped to hub connections.

    Parameters
    ----------
    generators_e : pd.DataFrame
       Existing generator capacities.
    generators_l : pd.DataFrame
       Layer potential capacities. Contains candidate generators without explicit bus assignments.
    nodes : pd.DataFrame
       Node definitions. Used to deduce bus assignments for candidates lacking explicit bus information.

    Returns
    -------
    pd.DataFrame
       Combined generator dataframe with corrected bus allocations and merged existing
       and potential capacities.
    """
    # Identify reallocations of radial wind farms using LAYER_POTENTIAL
    idx = ["location", "bus", "type", "pyear", "scenario", "carrier"]
    generators_el = generators_e.merge(
        generators_l,
        how="outer",
        on=["location", "type", "pyear", "scenario", "carrier"],
        suffixes=("", "_l"),
    )

    radial_inconsistent = generators_el.query(
        "carrier.str.contains('-r') "  # only radial connection
        "and p_nom_min != p_nom_min_l "  # when EXISTING and LAYER_POTENTIAL values are inconsistent
        "and ~(p_nom_min.isna() and p_nom_min_l == 0)"  # treat missing values and zeros as equivalent
    )[idx + ["p_nom_min"]]

    # Fix EXISTING technologies by reallocating radial to hubs
    corrections_radial = radial_inconsistent.assign(
        bus=lambda x: x.location, carrier=lambda x: x.carrier.str.replace("-r", "-oh")
    )
    generators_e_fixed = (
        pd.concat(
            [
                generators_e.set_index(idx).drop(
                    radial_inconsistent.set_index(idx).index
                ),
                corrections_radial.set_index(idx),
            ]
        )
        .groupby(level=list(range(len(idx))))
        .sum()
        .reset_index()
        .assign(carrier_mapped=lambda x: x.carrier.str.replace("h2", "dc", regex=True))
    )
    generators_l_fixed = (
        generators_l.drop(columns="p_nom_min")
        .query("p_nom_max != 0")
        .rename(columns={"carrier": "carrier_mapped"})
    )

    # Combine existing capacities with potentials
    # Set identical potentials for both hydrogen- and electricity-generating farms
    generators = (
        generators_e_fixed.merge(
            generators_l_fixed,
            how="outer",
            on=["location", "type", "pyear", "scenario", "carrier_mapped"],
        )
        .assign(
            p_nom_min=lambda df: df.p_nom_min.fillna(0),
            carrier=lambda df: df.carrier.fillna(df.carrier_mapped),
        )
        .drop(columns="carrier_mapped")
    )

    # Fill missing buses
    generators = (
        generators.merge(nodes[["location", "HOME_NODE"]], how="left", on="location")
        .assign(
            bus=lambda df: df.bus.fillna(
                df.HOME_NODE.where(df.carrier.str.contains("-r"), df.location)
            )
        )
        .drop(columns="HOME_NODE")
    )

    # Remove duplicates introduced by LAYER
    generators = generators.sort_values(
        by="p_nom_min", ascending=False
    ).drop_duplicates(subset=idx)

    return generators


def load_offshore_generators(
    fn: str,
    nodes: pd.DataFrame,
    scenario: str,
    planning_horizons: list[int],
    countries: list[str],
    extendable_carriers: dict[str, list[str]],
):
    """
    Load offshore generators data and format data.

    The `EXISTING` sheet is assumed to contain the existing capacities collected prior to any reallocations intended to align with the PEMMDB. This sheet appears to be excluded from the modelling exercise, except for hydrogen-generating capacities.

    The `LAYER_POTENTIAL` sheet is viewed as containing the reallocated existing capacities (excluding hydrogen-generating specific information) and the theoretical potentials per technology. Existing capacities are specified for both electricity- and hydrogen-generating offshore wind farms. Technology shares from `EXISTING` will be used to supplement the data.

    The `ZONE_POTENTIAL` sheet is considered as the source for achievable potentials for each node across all planning horizons. It establishes a nodal constraint on top of the theoretical potentials outlined by `LAYER_POTENTIAL`.

    **Existing capacities** will be read from the `LAYER_POTENTIAL` sheet, utilizing technology shares specified in `EXISTING` for hydrogen-generating capacities. A discrepancy of 526 MW for `DEOH002` in 2045 (across all scenarios) is noted when comparing existing capacities with `ZONE_POTENTIAL`. It remains uncertain which of the two values is correct: 5828.55 MW from `LAYER_POTENTIAL` or 6354.55 MW from `ZONE_POTENTIAL`. Currently, the value of 5828.55 MW is used. Additionally, the existing capacity of 3768.25 MW in 2040 exceeds the potential of 3242.25MW shown in the `LAYER_POTENTIAL` sheet.  This value has been adjusted to match the maximum potential value across all scenarios.

    **Potentials** will be obtained from both the `LAYER_POTENTIAL` and the `ZONE_POTENTIAL` sheets. `LAYER_POTENTIAL` will establish a technology level constraint, while `ZONE_POTENTIAL` will restrict expansion across all technologies at each node. The same 526 MW discrepancy for `DEOH002` in 2045 and 2050 (across all planning horizons and scenarios) has been noted and needs to be addressed to ensure that existing capacities do not exceed their potential. Currently, the value `ZONE_POTENTIAL` value is corrected at 5828.55 MW.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore generators data.
    nodes : pd.DataFrame
        DataFrame containing node information.
    scenario : str
        Scenario identifier to filter the grid data. Must be one of the scenario
        codes: "DE" (Distributed Energy), "GA" (Global Ambition), or
        "NT" (National Trends).
    planning_horizons : list[int]
        List of planning years to include in the cost data filtering.
    countries : list[str]
        List of country codes used to clean data.
    extendable_carriers : dict[str, list[str]]
        Dictionary mapping components to list of extendable carriers.

    Returns
    -------
    generators : pd.DataFrame
        DataFrame containing the formatted offshore generators data

    zone_trajectories : pd.DataFrame
        DataFrame containing the zone potentials trajectories
    """
    column_names = {
        "NODE": "bus",
        "OFFSHORE_NODE": "location",
        "OFFSHORE_NODE_TYPE": "type",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "TECHNOLOGY": "carrier",
        "TECH1": "carrier",
        "CAPEX": "capex",
        "OPEX": "opex",
        "MW": "p_nom_min",
        "EXISTING_MW": "p_nom_min",
        "MAX_MW": "p_nom_max",
    }

    column_del = [
        "TECH2",
        "TECH3",
        "TECH4",
        "TECH5",
        "TECH6",
        "MARGIN_MW",
        "LAYER",
    ]

    # Load data
    def load_generators(sheet_name, tech_switch=None):
        generators = pd.read_excel(
            fn,
            sheet_name=sheet_name,
        )
        if tech_switch:
            generators = generators.dropna(subset=tech_switch).assign(
                TECH1=lambda df: df[tech_switch]
            )
        generators = (
            generators.rename(columns=column_names)
            .replace({"scenario": SCENARIO_DICT})
            .query("pyear in @planning_horizons and scenario == @scenario")
            .assign(
                carrier=lambda x: "offwind-"
                + x.carrier.str.lower().replace("_", "-", regex=True)
            )
            .drop(columns=column_del, errors="ignore")
        )
        return generators

    generators_e = load_generators("EXISTING")
    generators_l_e = load_generators("LAYER_POTENTIAL")
    generators_l_h2 = load_generators("LAYER_POTENTIAL", tech_switch="TECH2").drop(
        columns="p_nom_min"
    )
    generators_l = pd.concat([generators_l_e, generators_l_h2])
    generators_z = load_generators("ZONE_POTENTIAL").drop(
        columns=["carrier", "p_nom_min"]
    )
    generators_c = load_generators("COST")
    generators_c[["capex", "opex"]] = generators_c[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW

    # Collect existing capacities and potentials in LAYER_POTENTIAL using H2 tech shares from EXISTING
    generators = collect_from_layer(generators_e, generators_l, nodes)

    # Collect potentials trajectories in ZONE_POTENTIAL
    zone_trajectories = generators_z

    # Resolve known DEOH002 data discrepancy
    # This is a temporary fix for a 526 MW discrepancy between LAYER_POTENTIAL
    # and ZONE_POTENTIAL data sources.
    # TODO: Remove this once upstream TYNDP data is corrected
    DEOH002_DISCREPANCY_MW = 526
    idx_l = generators.query(
        "location=='DEOH002' and pyear == 2040 and carrier=='offwind-ac-fb-r'"
    ).index
    generators.loc[idx_l, "p_nom_min"] = (
        generators.loc[idx_l, "p_nom_min"] - DEOH002_DISCREPANCY_MW
    )
    idx_z = zone_trajectories.query(
        "location=='DEOH002' and pyear in [2045, 2050]"
    ).index
    zone_trajectories.loc[idx_z, "p_nom_max"] = (
        zone_trajectories.loc[idx_z, "p_nom_max"] - DEOH002_DISCREPANCY_MW
    )

    # Collect cost assumptions
    generators = generators.merge(
        generators_c,
        how="left",
        on=["bus", "location", "pyear", "scenario", "type", "carrier"],
    )

    # Convert to explicit hydrogen buses
    mask = generators["carrier"].str.contains("h2")
    generators.loc[mask, "bus"] = generators.loc[mask, "bus"] + " H2"

    # Validate that all required cost assumptions are defined
    if generators[["capex", "opex"]].isna().any().any():
        raise ValueError("Missing generator cost data in input dataset.")
    generators["p_nom_extendable"] = generators["carrier"].isin(
        extendable_carriers["Generator"]
    )

    # Rename UK in GB
    generators[["bus", "location"]] = generators[["bus", "location"]].replace(
        "UK", "GB", regex=True
    )
    zone_trajectories["location"] = zone_trajectories["location"].replace(
        "UK", "GB", regex=True
    )

    # Filter selected countries
    generators = generators.assign(country=lambda x: x.location.str[:2]).query(
        "country in @countries"
    )
    zone_trajectories = zone_trajectories.assign(
        country=lambda x: x.location.str[:2]
    ).query("country in @countries")

    return generators, zone_trajectories


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_offshore_hubs", configfiles="config/config.tyndp-test.yaml"
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    planning_horizons = snakemake.params["planning_horizons"]
    countries = snakemake.params["countries"]
    extendable_carriers = snakemake.params["extendable_carriers"]
    options = snakemake.params.sector

    nodes = load_offshore_hubs(snakemake.input.nodes)

    grid = load_offshore_grid(
        snakemake.input.grid,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
        snakemake.params["offshore_hubs_tyndp"]["max_capacity"],
    )

    electrolysers = load_offshore_electrolysers(
        snakemake.input.electrolysers,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
        options,
    )

    generators, zone_trajectories = load_offshore_generators(
        snakemake.input.generators,
        nodes,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
        extendable_carriers,
    )

    # Convert country codes and retain only specified countries and offshore wind hubs
    nodes[["Bus", "location", "country"]] = nodes[
        ["Bus", "location", "country"]
    ].replace("UK", "GB", regex=True)
    nodes = nodes.query("type != 'Radial' and country in @countries").drop(
        columns="HOME_NODE"
    )

    # Save data
    nodes.to_csv(snakemake.output.offshore_buses, index=False)
    grid.to_csv(snakemake.output.offshore_grid, index=False)
    electrolysers.to_csv(snakemake.output.offshore_electrolysers, index=False)
    generators.to_csv(snakemake.output.offshore_generators, index=False)
    zone_trajectories.to_csv(snakemake.output.offshore_zone_trajectories, index=False)
