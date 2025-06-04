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
    )
    grid_costs[["capex", "opex"]] = grid_costs[["capex", "opex"]].mul(
        1e3
    )  # EUR/kW to EUR/MW
    grid_costs["carrier"] = grid_costs["carrier"].replace("E", "DC")
    grid_costs["scenario"] = grid_costs["scenario"].replace(scenario_dict)

    # Rename UK in GB
    grid[["bus0", "bus1"]] = grid[["bus0", "bus1"]].replace("UK", "GB", regex=True)
    grid_costs[["bus0", "bus1"]] = grid_costs[["bus0", "bus1"]].replace(
        "UK", "GB", regex=True
    )

    # Merge information
    grid = grid.merge(
        grid_costs, how="left", on=["bus0", "bus1", "pyear", "scenario", "carrier"]
    )

    # Assume non-extendable when missing data
    # TODO Validate assumption
    grid["p_nom_extendable"] = ~grid.isna().any(axis=1)
    grid[["capex", "opex"]] = grid[["capex", "opex"]].fillna(0)

    return grid


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

    # Save data
    nodes.to_csv(snakemake.output.offshore_buses, index=False)
    grid.to_csv(snakemake.output.offshore_grid, index=False)
