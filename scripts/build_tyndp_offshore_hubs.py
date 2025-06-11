# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean TYNDP Scenario Building offshore hubs data to be used in the PyPSA-Eur workflow. Depending on the scenario, different planning years (`pyear`) are available. DE and GA are defined for 2030, 2040 and 2050. NT scenario is only defined for 2030 and 2040. All the planning years are read at once.
"""

import logging

import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from shapely.geometry import Point

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"


def load_offshore_hubs(fn: str, countries: list[str]):
    """
    Load offshore hubs coordinates and format data.

    Offshore Hubs (OH) nodes are situated offshore, while Offshore Radial (OR) nodes are located in the homeland market node.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore hub data.
    countries : list[str]
        List of country codes used to clean data.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the offshore hub data.

        The GeoDataFrame uses the coordinate reference system defined by `GEO_CRS`.
    """
    column_dict = {
        "OFFSHORE_NODE": "Bus",
        "OFFSHORE_NODE_TYPE": "type",
        "HOME_NODE": "location",
        "LAT": "y",
        "LON": "x",
    }

    nodes = pd.read_excel(
        fn,
        sheet_name="NODE",
    ).rename(columns=column_dict)

    mask = nodes["Bus"].str.contains("OH")
    nodes.loc[mask, "location"] = nodes.loc[mask, "Bus"]
    nodes.loc[:, "country"] = nodes.loc[:, "location"].str[:2]

    nodes["geometry"] = nodes.apply(lambda row: Point(row["x"], row["y"]), axis=1)
    nodes = gpd.GeoDataFrame(nodes, geometry="geometry", crs=GEO_CRS)

    # rename UK in GB
    nodes[["Bus", "location", "country"]] = nodes[
        ["Bus", "location", "country"]
    ].replace("UK", "GB", regex=True)

    # filter selected countries
    nodes = nodes.query("country in @countries")

    return nodes


def expand_all_scenario(df: pd.DataFrame, scenarios: list):
    all_mask = df["scenario"] == "All"
    all_rows = (
        df[all_mask]
        .drop(columns="scenario")
        .merge(pd.DataFrame({"scenario": scenarios}), how="cross")
    )
    return pd.concat([df[~all_mask], all_rows], ignore_index=True)


def load_offshore_grid(
    fn: str,
    nodes: pd.DataFrame,
    scenario: str,
    planning_horizons: list[int],
    countries: list[str],
):
    """
    Load offshore grid (electricity and hydrogen) and format data.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore grid data.
    nodes : pd.DataFrame
        DataFrame containing node information (currently not used in function body
        but may be needed for validation or future functionality).
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
        DataFrame containing the merged offshore grid data.
    """
    column_dict = {
        "FROM": "bus0",
        "TO": "bus1",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "MARKET": "carrier",
        "CAPACITY": "p_nom",
        "CAPEX": "capex",
        "OPEX": "opex",
    }

    scenario_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
    }

    # Load reference grid
    grid = (
        pd.read_excel(
            fn,
            sheet_name="Reference grid",
        )
        .rename(columns=column_dict)
        .assign(
            p_min_pu=0,
            p_max_pu=1,
        )
    )
    grid["carrier"] = grid["carrier"].replace("E", "DC")
    grid = expand_all_scenario(grid, scenario_dict.values()).query(
        "scenario == @scenario"
    )

    # Load costs data
    grid_costs = (
        pd.read_excel(
            fn,
            sheet_name="COST",
        )
        .rename(columns=column_dict)
        .query("pyear in @planning_horizons")
        .replace({"scenario": scenario_dict})
    )
    grid_costs[["capex", "opex"]] = grid_costs[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW
    grid_costs["carrier"] = grid_costs["carrier"].replace("E", "DC")

    # Merge information
    grid = grid.merge(
        grid_costs, how="left", on=["bus0", "bus1", "pyear", "scenario", "carrier"]
    )

    # Assume non-extendable when missing data
    # TODO Validate assumption
    grid["p_nom_extendable"] = ~grid[["capex", "opex"]].isna().any(axis=1)
    grid[["capex", "opex"]] = grid[["capex", "opex"]].fillna(0)

    # Rename UK in GB
    grid[["bus0", "bus1"]] = grid[["bus0", "bus1"]].replace("UK", "GB", regex=True)

    # Filter selected countries
    grid = grid.assign(
        country0=lambda x: x.bus0.str[:2],
        country1=lambda x: x.bus1.str[:2],
    ).query("country0 in @countries and country1 in @countries")

    return grid


def load_offshore_electrolysers(
    fn: str, scenario: str, planning_horizons: list[int], countries: list[str]
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
    column_dict = {
        "NODE": "location",
        "OFFSHORE_NODE": "bus0",
        "OFFSHORE_NODE_TYPE": "type",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "CAPEX": "capex",
        "OPEX": "opex",
    }

    scenario_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
    }

    # Load electrolysers data
    electrolysers = (
        pd.read_excel(
            fn,
            sheet_name="COST",
        )
        .rename(columns=column_dict)
        .query("pyear in @planning_horizons")
        .replace({"scenario": scenario_dict})
        .query("scenario == @scenario")
        .assign(bus1=lambda x: x.bus0 + " H2")
    )

    electrolysers[["capex", "opex"]] = electrolysers[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW

    # rename UK in GB
    electrolysers[["bus0", "bus1", "location"]] = electrolysers[
        ["bus0", "bus1", "location"]
    ].replace("UK", "GB", regex=True)

    # filter selected countries
    electrolysers = electrolysers.assign(country=lambda x: x.bus0.str[:2]).query(
        "country in @countries"
    )

    return electrolysers


def get_shares(df, values, dropna=True):
    df_share = (
        df.pivot_table(
            index=["bus0", "type", "pyear", "scenario"],
            values=values,
            columns="carrier",
        )
        .pipe(lambda df: df.div(df.sum(axis=1), axis=0))
        .melt(ignore_index=False, value_name="tech_share")
        .reset_index()
    )
    if dropna:
        df_share = df_share.dropna(subset="tech_share")
    return df_share


def collect_from_layer(generators_e, generators_l):
    # Identify reallocations of radial wind farms using LAYER_POTENTIAL
    idx = ["location", "bus0", "type", "pyear", "scenario", "carrier"]
    generators_el = generators_e.merge(
        generators_l,
        how="outer",
        on=["bus0", "type", "pyear", "scenario", "carrier"],
        suffixes=("", "_l"),
    )

    radial_inconsistent = generators_el.query(
        "carrier.str.contains('-r') "  # only radial connection
        "and p_nom_min != p_nom_min_l "  # when EXISTING and LAYER_POTENTIAL values are inconsistent
        "and ~(p_nom_min.isna() and p_nom_min_l == 0)"  # treat missing values and zeros as equivalent
    )

    # Fix EXISTING technologies by reallocating radial to hubs
    generators_e_fixed = generators_e.copy().set_index(idx)
    corrections = radial_inconsistent.assign(
        location=lambda x: x.bus0, carrier=lambda x: x.carrier.str.replace("-r", "-oh")
    ).set_index(idx)
    generators_e_fixed = (
        pd.concat(
            [
                generators_e_fixed.drop(radial_inconsistent.set_index(idx).index),
                corrections[generators_e_fixed.columns],
            ]
        )
        .groupby(level=list(range(len(idx))))
        .sum()
    )

    # Calculate technology shares, considering PEMMDB related reallocations
    tech_shares = get_shares(generators_e_fixed, values="p_nom_min")

    # Apply H2 to DC carrier mapping and get relevant shares
    h2_shares = (
        tech_shares.query("carrier.str.contains('h2')")
        .assign(
            carrier_mapped=lambda x: x.carrier.str.replace(
                "h2", "dc", regex=True
            ).str.replace("-r", "-oh", regex=True)
        )
        .drop(columns=["tech_share", "carrier"])
        .merge(tech_shares, how="left")
    )

    # Apply shares to data to calculate existing capacities
    generators = (
        generators_l.merge(
            h2_shares,
            how="outer",
            left_on=["bus0", "type", "pyear", "scenario", "carrier"],
            right_on=["bus0", "type", "pyear", "scenario", "carrier_mapped"],
            suffixes=("_x", ""),
        )
        .assign(
            tech_share=lambda x: x.tech_share.fillna(1),
            carrier=lambda x: x.carrier.fillna(x.carrier_x),
            p_nom_min=lambda x: x.p_nom_min * x.tech_share,
        )
        .drop(columns=["carrier_x", "tech_share", "carrier_mapped"])
    )

    return generators


def load_offshore_generators(
    fn: str, scenario: str, planning_horizons: list[int], countries: list[str]
):
    """
    Load offshore generators data and format data.

    The `EXISTING` sheet is assumed to contain the collected existing capacities collected prior to any reallocations intended to align with the PEMMDB. This sheet appears to be excluded from the modelling exercise, except for hydrogen-generating capacities.

    The `LAYER_POTENTIAL` sheet is viewed as containing the reallocated existing capacities (excluding hydrogen-generating specific information) and the theoretical potentials per technology. Existing capacities are specified for both electricity- and hydrogen-generating offshore wind farms. Technology shares from `EXISTING` will be used to supplement the data.

    The `ZONE_POTENTIAL` sheet is considered as the source for achievable potentials for each node across all planning horizons. It establishes a nodal constraint on top of the theoretical potentials outlined by `LAYER_POTENTIAL`.

    **Existing capacities** will be read from the `LAYER_POTENTIAL` sheet, utilizing technology shares specified in `EXISTING` for hydrogen-generating capacities. A discrepancy of 526 MW for `DEOH002` in 2045 (across all scenarios) is noted when comparing existing capacities with `ZONE_POTENTIAL`. It remains uncertain which of the two values is correct: 5828.55 MW from `LAYER_POTENTIAL` or 6354.55 MW from `ZONE_POTENTIAL`. Currently, the value of 5828.55 MW is used.

    **Potentials** will be obtained from both the `LAYER_POTENTIAL` and the `ZONE_POTENTIAL` sheets. `LAYER_POTENTIAL` will establish a technology level constraint, while `ZONE_POTENTIAL` will restrict expansion across all technologies at each node. The same 526 MW discrepancy in `DEOH002` (across all planning horizons and scenarios) has been noted and needs to be addressed to ensure that existing capacities do not exceed their potential. Currently, the value `ZONE_POTENTIAL` value is corrected at 5828.55 MW.

    Parameters
    ----------
    fn : str
        Path to the Excel file containing offshore generators data.
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
    generators : pd.DataFrame
        DataFrame containing the formatted offshore generators data

    zone_trajectories : pd.DataFrame
        DataFrame containing the zone potentials trajectories
    """
    column_names = {
        "NODE": "location",
        "OFFSHORE_NODE": "bus0",
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

    scenario_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
    }

    # Load data
    def load_generators(sheet_name):
        generators = (
            pd.read_excel(
                fn,
                sheet_name=sheet_name,
            )
            .rename(columns=column_names)
            .query("pyear in @planning_horizons")
            .replace({"scenario": scenario_dict})
            .query("scenario == @scenario")
            .assign(
                carrier=lambda x: "offwind-"
                + x.carrier.str.lower().replace("_", "-", regex=True)
            )
            .drop(columns=column_del, errors="ignore")
        )
        return generators

    generators_e = load_generators("EXISTING")
    generators_l = load_generators("LAYER_POTENTIAL")
    generators_z = load_generators("ZONE_POTENTIAL").drop(
        columns=["carrier", "p_nom_min"]
    )
    generators_c = load_generators("COST")
    generators_c[["capex", "opex"]] = generators_c[["capex", "opex"]].mul(
        1e3
    )  # kEUR/MW to EUR/MW

    # Collect existing capacities and potentials in LAYER_POTENTIAL using H2 tech shares from EXISTING
    generators = collect_from_layer(generators_e, generators_l)

    # Collect potentials trajectories in ZONE_POTENTIAL
    zone_trajectories = generators_z

    # Resolve discrepancy in DEOH002
    idx = zone_trajectories.query("bus0=='DEOH002' and pyear in [2045, 2050]").index
    zone_trajectories.loc[idx, "p_nom_max"] = (
        zone_trajectories.loc[idx, "p_nom_max"] - 526
    )

    # Collect cost assumptions
    generators = generators.merge(
        generators_c,
        how="left",
        on=["bus0", "pyear", "scenario", "type", "carrier"],
    )

    # Ensure that all cost assumptions are present
    idx = generators[["capex", "opex"]].isna().any(axis=1)
    if not (generators.loc[idx].p_nom_min == 0).all():
        raise RuntimeError("Missing generator cost data in input dataset.")
    generators = generators.dropna(subset=["capex", "opex"])
    generators.loc[:, "p_nom_extendable"] = True

    # Rename UK in GB
    generators[["bus0", "location"]] = generators[["bus0", "location"]].replace(
        "UK", "GB", regex=True
    )
    zone_trajectories["bus0"] = zone_trajectories["bus0"].replace(
        "UK", "GB", regex=True
    )

    # Filter selected countries
    generators = generators.assign(country=lambda x: x.bus0.str[:2]).query(
        "country in @countries"
    )
    zone_trajectories = zone_trajectories.assign(
        country=lambda x: x.bus0.str[:2]
    ).query("country in @countries")

    return generators, zone_trajectories


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_offshore_hubs", configfiles="config/test/config.tyndp.yaml"
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    planning_horizons = snakemake.params["planning_horizons"]
    countries = snakemake.params["countries"]

    nodes = load_offshore_hubs(snakemake.input.nodes, countries)

    grid = load_offshore_grid(
        snakemake.input.grid,
        nodes,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
    )

    electrolysers = load_offshore_electrolysers(
        snakemake.input.electrolysers,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
    )

    generators, zone_trajectories = load_offshore_generators(
        snakemake.input.generators,
        snakemake.params["scenario"],
        planning_horizons,
        countries,
    )

    # Save data
    nodes.to_csv(snakemake.output.offshore_buses, index=False)
    grid.to_csv(snakemake.output.offshore_grid, index=False)
    electrolysers.to_csv(snakemake.output.offshore_electrolysers, index=False)
    generators.to_csv(snakemake.output.offshore_generators, index=False)
    zone_trajectories.to_csv(snakemake.output.offshore_zone_trajectories, index=False)
