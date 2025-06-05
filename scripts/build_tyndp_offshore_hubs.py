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


def load_offshore_hubs(fn: str):
    """
    Load offshore hubs coordinates and format data.

    Offshore Hubs (OH) nodes are situated offshore, while Offshore Radial (OR) nodes are located in the homeland market node.

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

    nodes["geometry"] = nodes.apply(lambda row: Point(row["x"], row["y"]), axis=1)
    nodes = gpd.GeoDataFrame(nodes, geometry="geometry", crs=GEO_CRS)

    # rename UK in GB
    nodes[["Bus", "location"]] = nodes[["Bus", "location"]].replace(
        "UK", "GB", regex=True
    )

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
    fn: str, nodes: pd.DataFrame, scenario: str, planning_horizons: list[int]
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
    )  # EUR/kW to EUR/MW
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

    return grid


def load_offshore_electrolysers(fn: str, scenario: str, planning_horizons: list[int]):
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
    )  # EUR/kW to EUR/MW

    # rename UK in GB
    electrolysers[["bus0", "bus1", "location"]] = electrolysers[
        ["bus0", "bus1", "location"]
    ].replace("UK", "GB", regex=True)

    return electrolysers


def load_offshore_generators(fn: str, scenario: str, planning_horizons: list[int]):
    """
    Load offshore generators data and format data.

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

    Returns
    -------
    pd.DataFrame
        DataFrame containing the formatted offshore generators data.
    """
    column_dict = {  # TODO
        "NODE": "location",
        "OFFSHORE_NODE": "bus0",
        "OFFSHORE_NODE_TYPE": "type",
        "YEAR": "pyear",
        "SCENARIO": "scenario",
        "TECHNOLOGY": "carrier",
        "CAPEX": "capex",
        "OPEX": "opex",
        "MW": "p_nom_min",
    }

    scenario_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
    }

    # Load generators data
    generators = (
        pd.read_excel(
            fn,
            sheet_name="EXISTING",
        )
        .rename(columns=column_dict)
        .query("pyear in @planning_horizons")
        .replace({"scenario": scenario_dict})
        .query("scenario == @scenario")
        .assign(
            carrier=lambda x: "offwind-"
            + x.carrier.str.lower().replace("_", "-", regex=True)
        )
    )

    # Load costs data
    generators_costs = (
        pd.read_excel(
            fn,
            sheet_name="COST",
        )
        .rename(columns=column_dict)
        .query("pyear in @planning_horizons")
        .replace({"scenario": scenario_dict})
        .query("scenario == @scenario")
        .assign(
            carrier=lambda x: "offwind-"
            + x.carrier.str.lower().replace("_", "-", regex=True)
        )
    )

    generators_costs[["capex", "opex"]] = generators_costs[["capex", "opex"]].mul(
        1e3
    )  # EUR/kW to EUR/MW

    # Merge information
    generators = generators_costs.merge(
        generators,
        how="outer",
        on=["bus0", "location", "pyear", "scenario", "type", "carrier"],
    ).assign(p_nom_min=lambda x: x["p_nom_min"].fillna(0))

    # Assume non-extendable when missing data
    # TODO Validate assumption
    generators["p_nom_extendable"] = ~generators[["capex", "opex"]].isna().any(axis=1)
    generators[["capex", "opex"]] = generators[["capex", "opex"]].fillna(0)

    # Rename UK in GB
    generators[["bus0", "location"]] = generators[["bus0", "location"]].replace(
        "UK", "GB", regex=True
    )

    return generators


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

    nodes = load_offshore_hubs(snakemake.input.nodes)

    grid = load_offshore_grid(
        snakemake.input.grid, nodes, snakemake.params["scenario"], planning_horizons
    )

    electrolysers = load_offshore_electrolysers(
        snakemake.input.electrolysers, snakemake.params["scenario"], planning_horizons
    )

    generators = load_offshore_generators(
        snakemake.input.generators, snakemake.params["scenario"], planning_horizons
    )

    # Save data
    nodes.to_csv(snakemake.output.offshore_buses, index=False)
    grid.to_csv(snakemake.output.offshore_grid, index=False)
    electrolysers.to_csv(snakemake.output.offshore_electrolysers, index=False)
    generators.to_csv(snakemake.output.offshore_generators, index=False)
