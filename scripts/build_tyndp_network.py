# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString, Point

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

# Poland organizes its lines in three sections,
# PL00 for demand/generation, -E for exporting lines and -I for importing lines
MAP_GRID_TYNDP = {
    "PL00E": "PL00",
    "PL00I": "PL00",
    "UK": "GB",
}

AC_VIRTUAL_NODES_IT = {
    "ITCO": "FR15",
    "ITVI": "ITSI",
}

IBFI_COORD = (63.0, 25.0)


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
    s = s.replace("FR-C", "FR15").replace("UK-N", "UKNI").replace("UK", "GB")
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
    distance_crs: str = DISTANCE_CRS,
):
    """
    Process bidding zones from the shape file and calculate representative point. Deduce the country shapes and their representative point.

    Parameters
    ----------
    bz_fn : str
        Path to bidding zone shape file.
    countries : list[str]
        List of countries to consider.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
    distance_crs : str, optional
        Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
    tuple
        A tuple of (bidding_shapes, country_shapes) GeoDataFrames.
    """
    bidding_zones = gpd.read_file(bz_fn)

    # Bidding zone shapes
    bidding_shapes = bidding_zones.assign(
        bz_id=lambda df: df["zone_name"].apply(format_bz_names),
        node=lambda df: (
            df.geometry.to_crs(distance_crs).representative_point().to_crs(geo_crs)
        ),
        x=lambda df: df["node"].x,
        y=lambda df: df["node"].y,
    ).set_index("bz_id")

    # Country shapes
    country_shapes = bidding_shapes.dissolve(by="country")[["geometry"]].assign(
        node=lambda df: (
            df.geometry.to_crs(distance_crs).representative_point().to_crs(geo_crs)
        ),
        x=lambda df: df["node"].x,
        y=lambda df: df["node"].y,
    )

    # Correct DK, IT, GR, SE and GB coordinates
    if "DK" in countries:
        country_shapes.loc["DK", ["node", "x", "y"]] = bidding_shapes.loc[
            "DKE1", ["node", "x", "y"]
        ]
    if "IT" in countries:
        country_shapes.loc["IT", ["node", "x", "y"]] = bidding_shapes.loc[
            "ITCA", ["node", "x", "y"]
        ]
    if "GR" in countries:
        country_shapes.loc["GR", ["node", "x", "y"]] = bidding_shapes.loc[
            "GR00", ["node", "x", "y"]
        ]
    if "SE" in countries:
        country_shapes.loc["SE", ["node", "x", "y"]] = bidding_shapes.loc[
            "SE01", ["node", "x", "y"]
        ]
    if "GB" in countries:
        country_shapes.loc["GB", ["node", "x", "y"]] = bidding_shapes.loc[
            "GB00", ["node", "x", "y"]
        ]

    return bidding_shapes, country_shapes


def _add_virtual_node(
    target_gdf: gpd.GeoDataFrame,
    new_bus: str,
    ref_bus: str,
    source_gdf: gpd.GeoDataFrame | None = None,
    **overrides,
) -> None:
    """
    Add a virtual node to target Dataframe as a copy of an existing reference bus.

    The virtual node inherits every attribute of the reference bus and takes its
    own name as ``station_id`` and ``tags``. Any remaining attribute is set
    through ``overrides``.

    Parameters
    ----------
    target_gdf : gpd.GeoDataFrame
        Bus GeoDataFrame the virtual node is appended to, modified in place.
    new_bus : str
        Name of the virtual node.
    ref_bus : str
        Name of the reference bus whose attributes are copied.
    source_gdf : gpd.GeoDataFrame, optional
        Bus GeoDataFrame holding the reference bus. Defaults to None in which case ``target_gdf`` is used.
    **overrides, optional
        Attribute values overriding those inherited from the reference bus.
    """
    source_gdf = target_gdf if source_gdf is None else source_gdf
    target_gdf.loc[new_bus] = (
        source_gdf.loc[[ref_bus]]
        .assign(station_id=new_bus, tags=new_bus, **overrides)
        .loc[ref_bus]
    )


def build_buses(
    buses_fn: str,
    countries: list[str],
    bidding_shapes: gpd.GeoDataFrame,
    country_shapes: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
):
    """
    Extend the node list for both electricity and hydrogen with attributes, incl. country and coordinates.

    Parameters
    ----------
    buses_fn : str
        Path to bidding zone shape file.
    countries : list[str]
        List of countries to consider.
    bidding_shapes : gpd.GeoDataFrame
        A GeoDataFrame including bidding zone geometry, representative point and id.
    country_shapes : gpd.GeoDataFrame
        A GeoDataFrame including country geometry and representative point.
    geo_crs : str, optional
        Coordinate reference system for geographic calculations. Defaults to GEO_CRS.

    Returns
    -------
    tuple
        A tuple of (buses, buses_h2) GeoDataFrames.
    """
    buses = (
        pd.read_excel(buses_fn)
        .replace("UK", "GB", regex=True)
        .merge(
            bidding_shapes[["country", "node", "x", "y"]],
            how="outer",
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
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs)

    # Assume the same coordinates for all LU buses
    if "LU" in countries:
        buses.loc["LUB1"] = buses.loc["LUB1"].fillna(buses.loc["LUG1"])
        buses.loc["LUF1"] = buses.loc["LUF1"].fillna(buses.loc["LUG1"])
        buses.loc["LUV1"] = buses.loc["LUV1"].fillna(buses.loc["LUG1"])

    # Manually add Italian virtual nodes  # TODO Refine assumptions
    if "IT" in countries:
        for node, location in AC_VIRTUAL_NODES_IT.items():
            _add_virtual_node(
                target_gdf=buses, new_bus=node, ref_bus=location, country="IT"
            )

    buses_h2 = (
        country_shapes[["node", "x", "y"]]
        .reset_index()
        .rename({"node": "geometry"}, axis=1)
        .assign(
            bus_id=lambda df: df[["country"]] + " H2",
            station_id=lambda df: df["bus_id"],
            voltage=None,
            dc="f",
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )
    buses_h2 = gpd.GeoDataFrame(buses_h2, geometry="geometry", crs=geo_crs)

    # Manually add IBIT and IBFI nodes  # TODO Refine assumptions
    if "IT" in countries:
        _add_virtual_node(
            target_gdf=buses_h2,
            new_bus="IBIT H2",
            ref_bus="ITN1",
            source_gdf=buses,
            voltage=None,
            dc="f",
        )
    if "FI" in countries:
        ibfi_lat, ibfi_long = IBFI_COORD
        _add_virtual_node(
            target_gdf=buses_h2,
            new_bus="IBFI H2",
            ref_bus="FI H2",
            x=ibfi_long,
            y=ibfi_lat,
            geometry=Point(ibfi_long, ibfi_lat),
        )

    return buses, buses_h2


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
    known_exceptions = {
        "DEKF",  # Connection from DE to the Kriegers Flak offshore wind farm
        "DKKF",  # Connection from DK to the Kriegers Flak offshore wind farm
        "DZ00",  # Algeria
        "EG00",  # Egypt
        "IS00",  # Iceland
        "IL00",  # Israel
        "LY00",  # Libya
        "MA00",  # Morocco
        "MD00",  # Moldova
        "PS00",  # Palestine
        "TN00",  # Tunisia
        "TR00",  # Turkey
        "UA00",  # Ukraine
        "UA01",  # Ukraine
    }
    if unknown_buses - known_exceptions:
        logger.warning(
            f"Dropping links connected to unknown buses: "
            f"{', '.join(sorted(unknown_buses - known_exceptions))}"
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
):
    """
    Process reference grid information to produce link data. p_nom are NTC values.

    Parameters
    ----------
    grid_fn : str | Path
        Path to bidding zone shape file.
    buses : gpd.GeoDataFrame
        A GeoDataFrame of electrical buses including country and coordinates.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame including NTC from the reference grid.
    """
    links = pd.read_excel(grid_fn)
    links = extract_grid_data_tyndp(
        links=links,
        idx_prefix="Transmission line",
        replace_dict=MAP_GRID_TYNDP,
        idx_connector="->",
    )

    # Add missing attributes
    links = add_links_missing_attributes(links, buses)

    return links


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    # Build node coordinates
    bidding_shapes, country_shapes = build_shapes(
        snakemake.input.bidding_shapes, countries
    )
    buses, buses_h2 = build_buses(
        snakemake.input.buses, countries, bidding_shapes, country_shapes
    )

    # Build links
    links = build_links(snakemake.input.elec_reference_grid, buses)

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
