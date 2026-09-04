# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString

from scripts._helpers import (
    configure_logging,
    extract_grid_data_tyndp,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
BUSES_COLUMNS = [
    "station_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
    # Marks which raw node-list sheet a bus came from: "onshore"/"offshore"
    # for electricity buses (see `ELEC_SHEET_CATEGORIES`), or
    # "offshore"/"import"/"bottleneck"/"Z1"/"Z2" for hydrogen buses (see
    # `H2_SHEET_CATEGORIES`).
    "category",
]
LINES_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "circuits",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
LINKS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "station_id",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]

# Countries bordering the modelled area that never have their own node in the
# TYNDP electricity node list, but still appear as a link endpoint in the
# reference grid.
KNOWN_EXCEPTIONS = {
    "DZ00",  # Algeria
    "EG00",  # Egypt
    "FR15",  # Corsica (French side of the FR15-ITCO interconnector)
    "IS00",  # Iceland
    "IL00",  # Israel
    "LY00",  # Libya
    "MA00",  # Morocco
    "PS00",  # Palestine
    "TN00",  # Tunisia
}

ELEC_NODE_SHEETS = ["Electricity", "Electricity_Offshore"]

# Category assigned to electricity nodes from each raw node-list sheet.
ELEC_SHEET_CATEGORIES = {
    "Electricity": "onshore",
    "Electricity_Offshore": "offshore",
}

H2_NODE_SHEETS = ["H2_Demand", "H2_Offshore", "H2_Imports", "H2_Bottlenecks"]

# Category assigned to hydrogen nodes from each raw node-list sheet, except
# "H2_Demand" whose nodes are categorised as "Z1"/"Z2" instead (see
# `categorize_h2_demand_zone`).
H2_SHEET_CATEGORIES = {
    "H2_Offshore": "offshore",
    "H2_Imports": "import",
    "H2_Bottlenecks": "bottleneck",
}

# Raw H2_Imports entries that are already country codes rather
# than TYNDP node codes, and thus need no further country extraction.
BARE_COUNTRY_CODES = {"DZ", "MA", "NO", "Y_NO", "TN", "TR", "IL", "UA"}


def format_bz_names(s: str) -> str:
    """
    Standardize bidding zone name formats to Open-TYNDP conventions.

    Parameters
    ----------
    s : str
        Raw bidding zone name string to format.

    Returns
    -------
    str
        Formatted bidding zone name with standardized region codes.
    """
    s = s.replace("UK-N", "UKNI").replace("UK", "GB")
    return s


def extract_shape_by_bbox(
    gdf: gpd.GeoDataFrame,
    country: str,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
    region_id: str,
):
    """
    Extracts a shape from a country's GeoDataFrame based on latitude and longitude bounds.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame containing country geometries.
    country : str
        The country code or name to filter.
    min_lon : float
        Minimum longitude bound for extraction.
    max_lon : float
        Maximum longitude bound for extraction.
    min_lat : float
        Minimum latitude bound for extraction.
    max_lat : float
        Maximum latitude bound for extraction.
    region_id : str
        String to assign an ID to the extracted region.

    Returns
    -------
    gpd.GeoDataFrame
        Updated GeoDataFrame with the extracted shape separated.
    """
    country_gdf = gdf.explode().query(f"country == '{country}'").reset_index(drop=True)

    extracted_region = country_gdf.cx[min_lon:max_lon, min_lat:max_lat].assign(
        id=region_id
    )

    remaining_country = (
        country_gdf.drop(extracted_region.index).dissolve(by="country").reset_index()
    )

    return pd.concat(
        [
            gdf.query(f"country != '{country}'"),
            remaining_country,
            extracted_region.dissolve(by="country").reset_index(),
        ]
    ).reset_index(drop=True)


def build_shapes(
    bz_fn: str,
    countries: list[str],
    geo_crs: str = GEO_CRS,
):
    """
    Process bidding zones from the shape file and calculate representative point.

    Parameters
    ----------
    bz_fn : str
        Path to bidding zone shape file.
    countries : list[str]
        List of countries to consider.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.

    Returns
    -------
    gpd.GeoDataFrame
        Bidding zone shapes with a representative point per zone.
    """
    bidding_zones = gpd.read_file(bz_fn)

    bidding_shapes = bidding_zones.assign(
        bz_id=lambda df: df["zone_name"].apply(format_bz_names),
        node=lambda df: (
            df.geometry.to_crs(DISTANCE_CRS).representative_point().to_crs(geo_crs)
        ),
        x=lambda df: df["node"].x,
        y=lambda df: df["node"].y,
    ).set_index("bz_id")

    return bidding_shapes


def build_buses(
    buses_fn: str,
    countries: list[str],
    bidding_shapes: gpd.GeoDataFrame,
    offshore_bus_locations_fn: str,
    geo_crs: str = GEO_CRS,
):
    """
    Build the electricity node list with attributes, incl. country and coordinates.

    Node identities are kept exactly as given in the TYNDP node list (no
    renaming, splitting or merging), except for the "UK" -> "GB" country-code
    convention used throughout Open-TYNDP.

    Each bus is additionally tagged with a ``category``: "onshore" or
    "offshore", depending on whether it came from the "Electricity" or
    "Electricity_Offshore" node-list sheet (see `ELEC_SHEET_CATEGORIES`).

    Parameters
    ----------
    buses_fn : str
        Path to the TYNDP node list Excel file ("LIST OF NODES.xlsx").
    countries : list[str]
        List of countries to consider.
    bidding_shapes : gpd.GeoDataFrame
        A GeoDataFrame including bidding zone geometry, representative point and id.
    offshore_bus_locations_fn : str
        Path to a CSV of manually-guessed ``x``/``y`` coordinates (see
        ``data/tyndp_offshore_bus_location.csv``), keyed by ``bus_id``, used
        to fill in coordinates for offshore/virtual nodes that have no
        matching bidding-zone shape.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.

    Returns
    -------
    gpd.GeoDataFrame
        Electricity buses as used in Open-TYNDP.
    """
    nodes = pd.concat(
        [
            pd.read_excel(buses_fn, sheet_name=sheet).assign(sheet=sheet)
            for sheet in ELEC_NODE_SHEETS
        ],
        ignore_index=True,
    )

    buses = (
        nodes.replace("UK", "GB", regex=True)
        .merge(
            bidding_shapes[["country", "node", "x", "y"]],
            how="left",
            left_on="NODE",
            right_index=True,
        )
        .rename({"NODE": "bus_id", "node": "geometry"}, axis=1)
        .assign(
            station_id=lambda df: df["bus_id"],
            voltage=380,  # TODO Improve assumption
            dc=None,
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
            # Bidding-zone shapes don't cover offshore/virtual/sub-zone node
            # codes (e.g. "BEO1_OFF", "PL00E"); fall back to the node code's
            # country prefix
            country=lambda df: df["country"].fillna(df["bus_id"].map(extract_country)),
            category=lambda df: df["sheet"].map(ELEC_SHEET_CATEGORIES),
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs)

    # Assume the same coordinates for all LU buses
    if "LU" in countries:
        buses.loc["LUB1"] = buses.loc["LUB1"].fillna(buses.loc["LUG1"])
        buses.loc["LUF1"] = buses.loc["LUF1"].fillna(buses.loc["LUG1"])
        buses.loc["LUV1"] = buses.loc["LUV1"].fillna(buses.loc["LUG1"])

    # Fill in manually-guessed coordinates for offshore/virtual nodes that
    # have no matching bidding-zone shape (see data/tyndp_offshore_bus_location.csv)
    offshore_locations = pd.read_csv(offshore_bus_locations_fn, index_col="bus_id")
    missing = buses.index[
        buses["geometry"].isna() & buses.index.isin(offshore_locations.index)
    ]
    if not missing.empty:
        buses.loc[missing, "x"] = offshore_locations.loc[missing, "x"]
        buses.loc[missing, "y"] = offshore_locations.loc[missing, "y"]
        buses.loc[missing, "geometry"] = gpd.points_from_xy(
            offshore_locations.loc[missing, "x"], offshore_locations.loc[missing, "y"]
        )

    return buses


def build_country_shapes(bidding_shapes: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Derive a representative point per country from electricity bidding zone shapes.

    Parameters
    ----------
    bidding_shapes : gpd.GeoDataFrame
        Bidding zone shapes as returned by `build_shapes`.

    Returns
    -------
    gpd.GeoDataFrame
        Country shapes with a representative point per country.
    """
    return bidding_shapes.dissolve(by="country")[["geometry", "x", "y"]]


def extract_country(bus_id: str) -> str:
    """
    Approximate the country of a raw TYNDP node code.

    Every raw code starts with its two-letter country code, except for a
    handful of hydrogen import/bottleneck node naming conventions (H2_Imports,
    H2_Bottlenecks) handled explicitly below; those never match electricity
    node codes, so this also works unchanged for electricity nodes without a
    matching bidding-zone shape (e.g. offshore/virtual/sub-zone codes).

    Parameters
    ----------
    bus_id : str
        Raw TYNDP node code (e.g. "DEh2Z1", "Ammonia_BE", "IB_ITh2", "BEO1_OFF").

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


def build_buses_h2(
    nodes_fn: str,
    bidding_shapes: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
) -> gpd.GeoDataFrame:
    """
    Build the hydrogen node list with attributes, incl. approximate coordinates.

    Combines the ``H2_Demand``, ``H2_Offshore``, ``H2_Imports`` and
    ``H2_Bottlenecks`` sheets of the TYNDP node list into a single hydrogen
    bus list. Node identities are kept exactly as given in the raw data (no
    renaming, splitting or merging), except for the "UK" -> "GB" country-code
    convention used throughout Open-TYNDP.

    Each bus is additionally tagged with a ``category``: "offshore", "import"
    or "bottleneck" for nodes from the corresponding sheet, and "Z1"/"Z2" for
    ``H2_Demand`` nodes depending on whether their node code names the "Z1"
    zone (e.g. "DEh2Z1") or not (implicitly "Z2").

    Node coordinates are not given directly in the TYNDP node list (unlike
    electricity, hydrogen nodes have no matching bidding-zone shape), so they
    are approximated by the representative point of the node's country,
    derived from the electricity bidding-zone shapes. This is a coarse
    approximation for zone-split countries and import/bottleneck nodes; it is
    only used for plotting and downstream distance-based calculations, not
    for the pipeline topology itself.

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
        [
            pd.read_excel(nodes_fn, sheet_name=sheet).assign(sheet=sheet)
            for sheet in H2_NODE_SHEETS
        ],
        ignore_index=True,
    ).replace("UK", "GB", regex=True)

    is_demand = nodes["sheet"] == "H2_Demand"
    nodes["category"] = nodes["sheet"].map(H2_SHEET_CATEGORIES)
    nodes.loc[is_demand, "category"] = nodes.loc[is_demand, "NODE"].apply(
        lambda node: "Z1" if "Z1" in node else "Z2"
    )

    nodes["country"] = nodes["NODE"].map(extract_country)
    missing_shape = set(nodes["country"]) - set(country_shapes.index)
    if missing_shape:
        logger.warning(
            "No bidding-zone shape for countries, dropping coordinates for "
            f"hydrogen nodes in: {', '.join(sorted(missing_shape))}"
        )

    buses_h2 = (
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

    return gpd.GeoDataFrame(buses_h2, geometry="geometry", crs=geo_crs)


def add_links_missing_attributes(
    links: pd.DataFrame,
    buses: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
    distance_crs: str = DISTANCE_CRS,
):
    """
    Add geometry attributes to links based on connected bus locations.

    Parameters
    ----------
    links : pd.DataFrame
        DataFrame of links with bus0 and bus1 columns.
    buses : gpd.GeoDataFrame
        GeoDataFrame of electrical buses including country and coordinates.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
    distance_crs : str, optional
        Coordinate reference system for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame of links with added geometry columns.
    """
    links = links.merge(
        buses["geometry"], how="left", left_on="bus0", right_index=True
    ).merge(
        buses["geometry"],
        how="left",
        left_on="bus1",
        right_index=True,
        suffixes=("0", "1"),
    )

    unknown_buses = set(
        links["bus0"][links[["bus0", "geometry0"]].isna().any(axis=1)]
    ).union(set(links["bus1"][links[["bus1", "geometry1"]].isna().any(axis=1)]))
    if unknown_buses - KNOWN_EXCEPTIONS:
        logger.warning(
            f"Dropping links connected to unknown buses: "
            f"{', '.join(sorted(unknown_buses - KNOWN_EXCEPTIONS))}"
        )
    links = links.dropna()  # TODO Remove this when all nodes are known

    links["geometry"] = gpd.GeoSeries(
        [
            LineString([p0, p1])
            for p0, p1 in zip(links["geometry0"], links["geometry1"])
        ],
        index=links.index,
    )
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=geo_crs)

    links = (
        links.assign(
            link_id=lambda df: df["bus0"] + "-" + df["bus1"] + "-DC",
            voltage=380,  # TODO Improve assumption
            length=lambda df: df.geometry.to_crs(distance_crs).length,
            underground="t",
            under_construction="f",
            tags=lambda df: df["bus0"] + " -> " + df["bus1"],
        )
        .groupby(by="link_id")
        .agg(
            {
                **{col: "first" for col in LINKS_COLUMNS if col != "p_nom"},
                "p_nom": "sum",
            }
        )[LINKS_COLUMNS]
    )
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=geo_crs)

    return links


def build_links(
    grid_fn,
    buses: gpd.GeoDataFrame,
    reference_year: int,
):
    """
    Process reference grid information to produce link data. p_nom are NTC values.

    Parameters
    ----------
    grid_fn : str | Path
        Path to the TYNDP electricity reference grid Excel file.
    buses : gpd.GeoDataFrame
        A GeoDataFrame of electrical buses including country and coordinates.
    reference_year : int
        Planning horizon sheet to read from the reference grid workbook. The
        border topology is identical across all available planning horizons,
        only NTC capacities differ; `build_tyndp_network` builds the static
        (unwildcarded) base network topology from a single, fixed year, given
        via `electricity: tyndp_reference_year`, while `build_tyndp_electricity_ntc`
        re-reads per-horizon NTC for `tyndp_scenario` runs.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame including NTC from the reference grid.
    """
    links = pd.read_excel(grid_fn, sheet_name=f"Year_{reference_year}")
    links["Border"] = links["Border"].replace("UK", "GB", regex=True)
    links = extract_grid_data_tyndp(links=links, idx_connector="->")

    # Add missing attributes
    links = add_links_missing_attributes(links, buses)

    return links


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_network", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    # Build node coordinates
    bidding_shapes = build_shapes(snakemake.input.bidding_shapes, countries)
    buses = build_buses(
        snakemake.input.buses,
        countries,
        bidding_shapes,
        snakemake.input.offshore_bus_locations,
    )
    buses_h2 = build_buses_h2(snakemake.input.buses, bidding_shapes)

    # Build links
    links = build_links(
        snakemake.input.elec_reference_grid,
        buses,
        reference_year=snakemake.params.reference_year,
    )

    # Build placeholder lines, converters and transformers as empty DataFrames
    lines = gpd.GeoDataFrame(columns=LINES_COLUMNS, geometry="geometry").set_index(
        pd.Index([], name="line_id")
    )
    converters = gpd.GeoDataFrame(
        columns=CONVERTERS_COLUMNS, geometry="geometry"
    ).set_index(pd.Index([], name="converter_id"))
    transformers = gpd.GeoDataFrame(
        columns=TRANSFORMERS_COLUMNS, geometry="geometry"
    ).set_index(pd.Index([], name="transformer_id"))

    # Export to csv for base_network
    buses.to_csv(snakemake.output["substations"], quotechar="'")
    buses_h2.to_csv(snakemake.output["substations_h2"], quotechar="'")
    lines.to_csv(snakemake.output["lines"], quotechar="'")
    links.to_csv(snakemake.output["links"], quotechar="'")
    converters.to_csv(snakemake.output["converters"], quotechar="'")
    transformers.to_csv(snakemake.output["transformers"], quotechar="'")

    # Export to GeoJSON for quick validations
    buses.to_file(snakemake.output["substations_geojson"])
    buses_h2.to_file(snakemake.output["substations_h2_geojson"])
    lines.to_file(snakemake.output["lines_geojson"])
    links.to_file(snakemake.output["links_geojson"])
    converters.to_file(snakemake.output["converters_geojson"])
    transformers.to_file(snakemake.output["transformers_geojson"])
