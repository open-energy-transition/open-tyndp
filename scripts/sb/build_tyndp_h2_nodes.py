# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds the TYNDP 2026 hydrogen node list.

Combines the ``H2_Demand``, ``H2_Offshore``, ``H2_Imports`` and
``H2_Bottlenecks`` sheets of the TYNDP node list into a single hydrogen bus
list. Node identities are kept exactly as given in the raw data (no
renaming, splitting or merging), except for the "UK" -> "GB" country-code
convention used throughout Open-TYNDP. This is a deliberate departure from
the TYNDP 2024 workflow, which built one synthetic country-level hydrogen
bus per country (optionally split into "Z1"/"Z2" zones via the
``h2_zones_tyndp`` option): TYNDP 2026 already provides the real,
zone-specific hydrogen nodes (e.g. ``FRh2N``/``FRh2S``/``FRh2SW``), so no
synthetic splitting is needed anymore.

Node coordinates are not given directly in the TYNDP node list (unlike
electricity, hydrogen nodes have no matching bidding-zone shape), so they are
approximated by the representative point of the node's country, derived from
the electricity bidding-zone shapes. This is a coarse approximation for
zone-split countries and import/bottleneck nodes; it is only used for
plotting and downstream distance-based calculations, not for the pipeline
topology itself.

Inputs
------
- `rules.retrieve_tyndp_2026.output.nodes`: TYNDP 2026 node list Excel file
  ("LIST OF NODES.xlsx"), sheets `H2_Demand`, `H2_Offshore`, `H2_Imports`,
  `H2_Bottlenecks`.
- `resources/bidding_zones.geojson`: Electricity bidding zone shapes, used to
  approximate hydrogen node coordinates at country level.

Outputs
-------
- `resources/tyndp/build/buses_h2.csv` / `.../geojson/buses_h2.geojson`:
  Hydrogen buses with approximate coordinates and country attribution.
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_tyndp_network import BUSES_COLUMNS, GEO_CRS, build_shapes

logger = logging.getLogger(__name__)

H2_NODE_SHEETS = ["H2_Demand", "H2_Offshore", "H2_Imports", "H2_Bottlenecks"]

# Raw H2_Imports entries that are already bare (non-EU) country codes rather
# than TYNDP node codes, and thus need no further country extraction.
BARE_COUNTRY_CODES = {"DZ", "MA", "NO", "Y_NO", "TN", "TR", "IL", "UA"}


def extract_country(bus_id: str) -> str:
    """
    Approximate the country of a raw TYNDP hydrogen node code.

    Every raw code (H2_Demand, H2_Offshore) starts with its two-letter
    country code, except for a handful of import/bottleneck node naming
    conventions (H2_Imports, H2_Bottlenecks) handled explicitly below.

    Parameters
    ----------
    bus_id : str
        Raw hydrogen node code (e.g. "DEh2Z1", "Ammonia_BE", "IB_ITh2").

    Returns
    -------
    str
        Two-letter country code.
    """
    if bus_id in BARE_COUNTRY_CODES:
        return bus_id[-2:]
    if bus_id.startswith("Ammonia_"):
        return bus_id.removeprefix("Ammonia_")[:2]
    if bus_id.startswith("IB_"):
        return bus_id[3:5]
    return bus_id[:2]


def build_country_shapes(bidding_shapes: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Derive a representative point per country from electricity bidding zone shapes.

    Parameters
    ----------
    bidding_shapes : gpd.GeoDataFrame
        Bidding zone shapes as returned by `build_tyndp_network.build_shapes`.

    Returns
    -------
    gpd.GeoDataFrame
        Country shapes with a representative point per country.
    """
    return bidding_shapes.dissolve(by="country")[["geometry", "x", "y"]]


def build_buses_h2(
    nodes_fn: str,
    bidding_shapes: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
) -> gpd.GeoDataFrame:
    """
    Build the hydrogen node list with attributes, incl. approximate coordinates.

    Parameters
    ----------
    nodes_fn : str
        Path to the TYNDP node list Excel file ("LIST OF NODES.xlsx").
    bidding_shapes : gpd.GeoDataFrame
        Electricity bidding zone shapes, used to approximate hydrogen node
        coordinates at country level.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.

    Returns
    -------
    gpd.GeoDataFrame
        Hydrogen buses as used in Open-TYNDP.
    """
    country_shapes = build_country_shapes(bidding_shapes)

    nodes = pd.concat(
        [pd.read_excel(nodes_fn, sheet_name=sheet) for sheet in H2_NODE_SHEETS],
        ignore_index=True,
    ).replace("UK", "GB", regex=True)

    nodes["country"] = nodes["NODE"].map(extract_country)
    missing_shape = set(nodes["country"]) - set(country_shapes.index)
    if missing_shape:
        logger.warning(
            "No bidding-zone shape for countries, dropping coordinates for "
            f"hydrogen nodes in: {', '.join(sorted(missing_shape))}"
        )

    buses = (
        nodes.merge(
            country_shapes[["x", "y"]], how="left", left_on="country", right_index=True
        )
        .rename({"NODE": "bus_id"}, axis=1)
        .assign(
            station_id=lambda df: df["bus_id"],
            voltage=None,
            dc="f",
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
            geometry=lambda df: gpd.points_from_xy(df["x"], df["y"]),
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )

    return gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_h2_nodes")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    bidding_shapes = build_shapes(snakemake.input.bidding_shapes, countries)
    buses_h2 = build_buses_h2(snakemake.input.nodes, bidding_shapes)

    buses_h2.to_csv(snakemake.output["buses_h2"], quotechar="'")
    buses_h2.to_file(snakemake.output["buses_h2_geojson"])
