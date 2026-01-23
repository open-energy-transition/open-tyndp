# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT


"""
Build TYNDP transmission projects from the Grid Investment Dataset.

This script processes the TYNDP 2024 Grid Investment Dataset to extract
transmission projects for inclusion in the reference grid.
Projects are filtered by commissioning year and category, then formatted
as network links with geometry attributes derived from bus locations.

Inputs
------
- ``data/tyndp_2024_bundle/Investment Datasets/GRID.xlsx``: Grid investment dataset
  containing electricity and hydrogen transmission projects with capacities and
  commissioning years.
- ``resources/tyndp/substations.geojson``: Bus with geometry.

Outputs
-------
- ``resources/tyndp/new_links.csv``: Processed transmission projects with bus
  connections, capacities (p_nom), lengths, and geometry attributes ready for
  network integration.
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import (
    configure_logging,
    extract_grid_data_tyndp,
    set_scenario_config,
)
from scripts.build_tyndp_network import (
    MAP_GRID_TYNDP,
    add_links_missing_attributes,
)

logger = logging.getLogger(__name__)


def read_invest_projects(
    fn_invest: str,
    buses: gpd.GeoDataFrame,
    carrier: str = "Electricity",
    year: int = 2035,
    op: str = "==",
    category: str = "Real",
) -> gpd.GeoDataFrame:
    """
    Read grid investment dataset for Electricity or Hydrogen.

    For Electricity, only 'Real' projects are considered by default (excluding 'Concept' projects).

    Parameters
    ----------
    fn_invest : str
        Path to the investment dataset file.
    buses : gpd.GeoDataFrame
        GeoDataFrame of TYNDP buses with geometry.
    carrier : str, optional
        Energy carrier to extract: 'Electricity' or 'Hydrogen'. Default is 'Electricity'.
    year : int, optional
        Commissioning year to filter projects. Default is 2035.
    op : str, optional
        Comparison operator for year filtering (e.g., '==', '<=', '>='). Default is '=='.
    category : str, optional
        Project category filter for Electricity. Default is 'Real'.

    Returns
    -------
    gpd.GeoDataFrame
        Processed links with geometry.
    """
    if carrier not in ["Electricity", "Hydrogen"]:
        raise ValueError(
            f"Carrier must be 'Electricity' or 'Hydrogen', got '{carrier}'."
        )
    border_condition = (
        f"BORDER.str.contains('{category}') & " if carrier == "Electricity" else ""
    )

    projects = (
        pd.read_excel(fn_invest, sheet_name=carrier)
        .query(f"{border_condition}YEAR {op} {year}")
        .rename(
            columns={
                "DIRECT CAPACITY INCREASE (MW)": "Summary Direction 1",
                "INDIRECT CAPACITY INCREASE (MW)": "Summary Direction 2",
                "FROM NODE": "bus0",
                "TO NODE": "bus1",
            }
        )
        .query("bus0 in @buses.index & bus1 in @buses.index")
        .groupby(["bus0", "bus1"])
        .sum()[["Summary Direction 1", "Summary Direction 2"]]
        .reset_index()
    )

    links = extract_grid_data_tyndp(
        links=projects,
        replace_dict=MAP_GRID_TYNDP,
        expand_from_index=False,
        idx_connector="-",
        idx_suffix="-DC",
        idx_separator="",
    )

    links = add_links_missing_attributes(links, buses)

    return links


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_transmission_projects")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    buses = gpd.read_file(snakemake.input.buses).set_index("bus_id")
    new_links = read_invest_projects(snakemake.input.invest_grid, buses)

    new_links.to_csv(snakemake.output[0], quotechar="'")
