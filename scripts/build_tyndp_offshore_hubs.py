# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean TYNDP Scenario Building offshore hubs data to be used in the PyPSA-Eur workflow. Depending on the scenario, different planning years (`pyear`) are available. DE and GA are defined for 2030, 2040 and 2050. NT scenario is only defined for 2030 and 2040. All the planning years are read at once.
"""

import logging

import geopandas as gpd
import numpy as np
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
    )  # kEUR/MW to EUR/MW

    # rename UK in GB
    electrolysers[["bus0", "bus1", "location"]] = electrolysers[
        ["bus0", "bus1", "location"]
    ].replace("UK", "GB", regex=True)

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


def collect_generators_capacities(generators_e, generators_l):
    # Determine technology shares
    generators_e_share_raw = get_shares(generators_e, values="p_nom_min")

    # Add missing H2 shares
    generators_e_share = generators_e_share_raw.query(
        "carrier.str.contains('h2')"
    ).assign(carrier_rfc=lambda x: x.carrier.str.replace("h2", "dc", regex=True))

    generators_e_share_rfc = (
        generators_e_share.drop(columns=["tech_share", "carrier"])
        .assign(carrier=lambda x: x.carrier_rfc)
        .merge(generators_e_share_raw, how="left")
    )
    generators_e_share = pd.concat([generators_e_share, generators_e_share_rfc])

    # Compute existing capacities using shares
    generators = (
        generators_l.merge(
            generators_e_share,
            how="outer",
            left_on=["bus0", "type", "pyear", "scenario", "carrier"],
            right_on=["bus0", "type", "pyear", "scenario", "carrier_rfc"],
            suffixes=("_x", ""),
        )
        .assign(
            tech_share=lambda x: x.tech_share.fillna(1),
            p_nom_min=lambda x: x.p_nom_min * x.tech_share,
            carrier=lambda x: x.carrier.fillna(x.carrier_x),
        )
        .drop(columns=["carrier_x", "tech_share", "carrier_rfc", "p_nom_max"])
    )

    return generators


def collect_generators_potentials(generators, generators_l, generators_z):
    # Get technology shares
    generators_l_share = get_shares(generators_l, values="p_nom_max", dropna=False)

    # Compute potentials using shares
    generators_z_tech = (
        generators_z.merge(
            generators_l_share,
            how="left",
            on=["bus0", "type", "pyear", "scenario"],
        )
        .assign(
            tech_share=lambda x: x.tech_share.fillna(0),
            p_nom_max=lambda x: x.p_nom_max * x.tech_share,
        )
        .drop(columns=["tech_share"])
        .query("p_nom_max != 0")
    )

    # Integrate potentials in data
    generators = (
        generators.merge(
            generators_z_tech,
            how="outer",
            on=["bus0", "type", "pyear", "scenario", "carrier"],
        )
        .replace(0, np.nan)
        .dropna(subset=["p_nom_min", "p_nom_max"], how="all")
    )

    # Copy DC potentials for H2 wind farms
    idx = generators.carrier.str.contains("h2")
    generators_h2 = (
        generators.loc[idx]
        .drop(columns="p_nom_max")
        .assign(
            carrier_ori=lambda x: x.carrier,
            carrier=lambda x: x.carrier.str.replace("h2", "dc", regex=True),
        )
        .merge(
            generators.drop(columns="p_nom_min"),
            how="left",
            on=["bus0", "type", "pyear", "scenario", "carrier"],
        )
        .drop(columns="carrier")
        .rename(columns={"carrier_ori": "carrier"})
    )

    generators = pd.concat([generators.loc[~idx], generators_h2]).fillna(0)

    return generators


def load_offshore_generators(fn: str, scenario: str, planning_horizons: list[int]):
    """
    Load offshore generators data and format data.

    The `COST` sheet provides techno-economic assumptions for offshore generators. It is assumed that only the specified generators can be expanded.

    The `EXISTING` sheet is assumed to contain the collected existing capacities collected prior to any reallocations intended to align with the PEMMDB. This sheet appears to be excluded from the modelling exercise.

    The `LAYER_POTENTIAL` sheet is viewed as containing the reallocated existing capacities and the theoretical potentials per technology. Existing capacities are specified for both electricity- and hydrogen-generating offshore wind farms. Technology shares from `EXISTING` will be used to supplement the data.

    The `ZONE_POTENTIAL` sheet is considered as the source for achievable potentials at each planning horizon. Technology shares from `LAYER_POTENTIAL` will be used to supplement the data.

    **Existing capacities** will be read from the `LAYER_POTENTIAL` sheet, utilizing technology shares specified in `EXISTING`. A discrepancy of 526 MW for `DEOH002` in 2045 (across all scenarios) is noted when comparing existing capacities with `ZONE_POTENTIAL`. It remains uncertain which of the two values is correct: 5828.55 MW from `LAYER_POTENTIAL` or 6354.55 MW from `ZONE_POTENTIAL`.

    **Potentials** will be read from the `ZONE_POTENTIAL` sheet, utilizing technology shares specified in `LAYER_POTENTIAL`. The same 526 MW discrepancy in `DEOH002` (across all planning horizons and scenarios) has been identified and needs to be addressed to ensure that existing capacities do not exceed their potential.

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

    # Collect existing capacities in LAYER_POTENTIAL using H2 tech shares from EXISTING
    generators = collect_generators_capacities(generators_e, generators_l)

    # Collect potentials in ZONE_POTENTIAL using tech shares from LAYER_POTENTIAL
    generators = collect_generators_potentials(generators, generators_l, generators_z)

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
