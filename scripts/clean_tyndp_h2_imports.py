# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to load and clean TYNDP H2 import data.
"""

import logging

import geopandas as gpd
import pandas as pd
from _helpers import (
    configure_logging,
    make_index,
    set_scenario_config,
)
from shapely.geometry import Point

logger = logging.getLogger(__name__)


def match_centroids(df, countries_centroids):
    """
    Matches coordinates of country centroids to bus0 countries.
    Manually matches coordinates next to Faroe Island ("FO") to Ammonia import node.

    Parameters
    ----------
    df : pd:DataFrame
        Dataframe containing import data with bus0 as import nodes
    countries_centroids : gpd.GeoDataFrame
        GeoDataFrame containing country centroid information as geometry

    Returns
    -------
    pd.DataFrame
        The function returns the input Dataframe df with matched coordinates inside new columns bus0_x and bus0_y
    """
    import_nodes = df.bus0.unique()

    # Match coordinates
    coordinates = (
        countries_centroids.replace(
            {"FO": "Ammonia"}
        )  # manually match coordinates next to Faroe Islands with Ammonia imports
        .query("ISO in @import_nodes")
        .set_index("ISO")
        .geometry
    )
    if "Ammonia" in import_nodes:
        # manually match coordinates next to Faroe Islands with Ammonia imports
        coordinates.loc["Ammonia"] = Point(
            coordinates.loc["Ammonia"].x - 1, coordinates.loc["Ammonia"].y + 1
        )
    if not coordinates.empty:
        logger.info(
            f"Found coordinates for import nodes: {', '.join(coordinates.index.values)}."
        )
    else:
        logger.warning(
            f"Can't match centroid coordinates as none of the import nodes are defined as ISO countries: {', '.join(import_nodes)}."
        )
    return df.assign(
        bus0_x=df.bus0.map(coordinates.x), bus0_y=df.bus0.map(coordinates.y)
    )


def load_import_data(fn, countries_centroids):
    """
    Load and clean TYNDP H2 import potentials, maximum capacity, offer quantity and marginal cost for pipeline and shipping
    Returns the cleaned data as dataframe.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP H2 imports data.

    Returns
    -------
    pd.DataFrame
        The function returns cleaned TYNDP H2 import potentials, maximum capacity, offer quantity and marginal cost.
    """

    column_dict = {
        "YEAR": "Year",
        "SCENARIO": "Scenario",
        "CORRIDOR": "Corridor",
        "NODE FROM": "bus0",
        "NODE TO": "bus1",
        "MAX CAPACITY [MW]": "p_nom",
        "OFFER QUANTITY [MW]": "offer_quantity",
        "OFFER PRICE [€/MWh]": "marginal_cost",
        "MAX ENERGY YEAR [GWh]": "e_sum_max",
    }

    replace_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
        "Lh2": "LH2",
    }

    # Read data, rename and convert to MWh
    imports = pd.read_excel(fn).rename(columns=column_dict).replace(replace_dict)
    imports.loc[:, "e_sum_max"] *= 1e3  # convert from GWh to MWh

    # Convert marginal_cost to numeric, replacing non-numeric values with NaN, then fill NaN with 0.0
    imports.loc[:, "marginal_cost"] = pd.to_numeric(
        imports.marginal_cost, errors="coerce"
    ).fillna(0.0)

    # Extract Band information into new column
    imports.loc[:, "Band"] = (
        imports.Corridor.str.split("-", expand=True).iloc[:, -1].str.lower()
    )

    # Match countries centroids for defining coordinates of import nodes
    imports = match_centroids(imports, countries_centroids)
    imports.index = imports.apply(make_index, axis=1, prefix="H2 import") + " - " + imports.Band

    return imports


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_tyndp_h2_imports")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load countries centroids for defining coordinates of import nodes later on
    countries_centroids = gpd.read_file(snakemake.input.countries_centroids)

    # Load and prep import potentials
    import_potentials = load_import_data(
        snakemake.input.import_potentials_raw, countries_centroids
    )

    # Save prepped H2 import potentials
    import_potentials.to_csv(snakemake.output.import_potentials_prepped)
