# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create interactive energy balance maps for the defined carriers using `n.explore()`.
"""

import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import pydeck as pdk
import pypsa
from pypsa.plot.maps.interactive import PydeckPlotter
from pypsa.statistics import get_transmission_carriers
from shapely.geometry import box

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import sanitize_carriers
from scripts.build_tyndp_network import IBFI_COORD

VALID_MAP_STYLES = PydeckPlotter.VALID_MAP_STYLES


def scalar_to_rgba(
    value: float,
    *,
    norm: mcolors.Normalize,
    cmap: mcolors.Colormap,
    alpha: float = 1.0,
) -> list[int]:
    """
    Map a scalar float value to an RGBA color encoded as 8-bit integers.

    Parameters
    ----------
    value : float
        Scalar input to map through the normalization and colormap.
    norm : matplotlib.colors.Normalize
        Normalization defining vmin and vmax used for scaling.
    cmap : matplotlib.colors.Colormap
        Colormap used to convert normalized values to RGBA colors.
    alpha : float, optional (default = 1.0)
        Opacity in the range [0, 1]. Overrides the colormap's alpha.

    Returns
    -------
    List[int]
        A list ``[R, G, B, A]`` where each channel is an integer in the 0–255 range.
    """

    # Clamp to normalization bounds
    p = max(norm.vmin, min(norm.vmax, value))

    # Convert to RGBA floats (0–1)
    r, g, b, _ = cmap(norm(p))

    # Clamp and apply alpha
    a = alpha if 0.0 <= alpha <= 1.0 else 1.0

    # Convert to 8-bit integers
    return [
        int(r * 255),
        int(g * 255),
        int(b * 255),
        int(a * 255),
    ]


def dissolve_h2_regions_tyndp(regions: gpd.GeoDataFrame, buses_h2_fn: str):
    """
    Dissolve hydrogen regions to align with the TYNDP topology.

    Two zones require dedicated treatments, as defined in ``build_tyndp_network``:
    IBIT and IBFI. The IBIT zone uses the ITN1 shape. The Finland shape is divided
    in two: the northern part is assigned to FI H2 and the southern part to IBFI H2.

    Parameters
    ----------
    regions : gpd.GeoDataFrame
        Regions shape of the electrical topology.
    buses_h2_fn : str
        File path to the H2 buses used.

    Returns
    -------
    gpd.GeoDataFrame
        Regions dissolved to the hydrogen topology.
    """
    buses_h2 = gpd.read_file(buses_h2_fn).set_index("bus_id")
    buses_h2_real = buses_h2[~buses_h2.index.str.startswith("IB")]
    country_to_bus = buses_h2_real.reset_index()[["bus_id", "country"]].set_index(
        "country"
    )["bus_id"]
    regions["country"] = regions.index.str[:2]
    regions["bus_id"] = regions["country"].map(country_to_bus)
    if "ITN1" in regions.index:
        regions.loc["ITN1", "bus_id"] = "IBIT H2"
    regions = regions.dissolve("bus_id")[["geometry"]]

    if "FI H2" in regions.index:
        ibfi_lat, ibfi_long = IBFI_COORD
        fi_lat = buses_h2.loc["FI H2"].y
        split_lat = (ibfi_lat + fi_lat) / 2
        fi_geom = regions.loc["FI H2", "geometry"]

        south = fi_geom.intersection(box(minx=-180, miny=-90, maxx=180, maxy=split_lat))
        north = fi_geom.intersection(box(minx=-180, miny=split_lat, maxx=180, maxy=90))

        fi_splitted = gpd.GeoDataFrame(
            {
                "geometry": [south, north],
                "bus_id": ["IBFI H2", "FI H2"],
            },
            crs=regions.crs,
        ).set_index("bus_id")
        regions = gpd.GeoDataFrame(
            pd.concat([regions[regions.index != "FI H2"], fi_splitted]),
            crs=regions.crs,
        )

    return regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_balance_map_interactive",
            clusters=50,
            opts="",
            sector_opts="",
            planning_horizons="2050",
            carrier="H2",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Interactive map settings
    settings = snakemake.params.settings
    unit_conversion = settings["unit_conversion"]
    cmap = settings["cmap"]
    region_alpha = settings["region_alpha"]
    region_unit = settings["region_unit"]
    branch_color = settings["branch_color"]
    arrow_size_factor = settings["arrow_size_factor"]
    bus_size_max = settings["bus_size_max"]
    branch_width_max = settings["branch_width_max"]
    map_style = settings.get("map_style")
    map_style = VALID_MAP_STYLES.get(map_style, "road")
    tooltip = settings["tooltip"]

    # Import
    n = pypsa.Network(snakemake.input.network)
    sanitize_carriers(n, snakemake.config)
    pypsa.options.params.statistics.round = 8
    pypsa.options.params.statistics.drop_zero = True
    pypsa.options.params.statistics.nice_names = False

    carrier = snakemake.wildcards.carrier
    carrier = carrier.replace("_", " ")
    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    regions.geometry = regions.geometry.simplify(
        0.05
    )  # reduces file size

    if carrier == "H2" and snakemake.params.h2_topology_tyndp:
        regions = dissolve_h2_regions_tyndp(regions, snakemake.input.buses_h2)

    # Fill missing carrier colors
    missing_color = "#808080"
    b_missing = n.carriers.query("color == '' or color.isnull()").index
    n.carriers.loc[b_missing, "color"] = missing_color

    transmission_carriers = get_transmission_carriers(n, bus_carrier=carrier).rename(
        {"name": "carrier"}
    )
    components = transmission_carriers.unique("component")
    carriers = transmission_carriers.unique("carrier")

    ### Pie charts
    eb = n.statistics.energy_balance(
        bus_carrier=carrier,
        groupby=["bus", "carrier"],
    )

    # Only carriers that are also in the energy balance
    carriers_in_eb = carriers[carriers.isin(eb.index.get_level_values("carrier"))]

    eb.loc[components] = eb.loc[components].drop(index=carriers_in_eb, level="carrier")
    eb = eb.dropna()
    bus_size = eb.groupby(level=["bus", "carrier"]).sum()

    # line and links widths according to optimal capacity
    flow = n.statistics.transmission(groupby=False, bus_carrier=carrier, at_port=0)
    if not flow.empty:
        flow_reversed_mask = flow.index.get_level_values(1).str.contains("reversed")
        flow_reversed = flow[flow_reversed_mask].rename(
            lambda x: x.replace("-reversed", "")
        )
        flow = flow[~flow_reversed_mask].subtract(flow_reversed, fill_value=0)

    # only line first index
    line_flow = (
        flow.loc[flow.index.get_level_values(0).str.contains("Line")]
        .copy()
        .droplevel(0)
    )
    link_flow = (
        flow.loc[flow.index.get_level_values(0).str.contains("Link")]
        .copy()
        .droplevel(0)
    )

    branch_components = ["Link"]
    if carrier == "AC":
        branch_components = ["Line", "Link"]

    ### Prices
    buses = n.buses.query("carrier in @carrier").index
    demand = (
        n.statistics.energy_balance(
            bus_carrier=carrier, aggregate_time=False, groupby=["bus", "carrier"]
        )
        .clip(lower=0)
        .groupby("bus")
        .sum()
        .reindex(buses)
        .rename(n.buses.location)
        .T
    )

    weights = n.snapshot_weightings.generators
    price = (
        weights
        @ n.buses_t.marginal_price.reindex(buses, axis=1).rename(
            n.buses.location, axis=1
        )
        / weights.sum()
    )

    if carrier == "co2 stored" and "CO2Limit" in n.global_constraints.index:
        co2_price = n.global_constraints.loc["CO2Limit", "mu"]
        price = price - co2_price

    # if only one price is available, use this price for all regions
    if price.size == 1:
        regions["price"] = price.values[0]
        shift = round(abs(price.values[0]) / 20, 0)
    else:
        regions["price"] = price.reindex(regions.index).fillna(0)
        shift = 0

    vmin, vmax = regions.price.min() - shift, regions.price.quantile(0.98) + shift
    if settings["vmin"] is not None:
        vmin = settings["vmin"]
    if settings["vmax"] is not None:
        vmax = settings["vmax"]

    # Map colors
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap)

    regions["color"] = regions["price"].apply(
        scalar_to_rgba,
        norm=norm,
        cmap=cmap,
        alpha=region_alpha,
    )

    # Create tooltips
    regions["tooltip_html"] = (
        "<b>"
        + regions.index
        + "</b><br>"
        + "<b>Weighted price:</b> "
        + regions["price"].round(2).astype(str)
        + " "
        + region_unit
    )
    # regions["tooltip_html"] = regions["price"].round(2).astype(str)
    # Create layer
    regions_layer = pdk.Layer(
        "GeoJsonLayer",
        regions,
        stroked=True,
        filled=True,
        get_fill_color="color",
        get_line_color=[255, 255, 255, 255],
        line_width_min_pixels=1,
        pickable=True,
        auto_highlight=True,
    )

    map = n.explore(
        branch_components=branch_components,
        bus_size=bus_size.div(unit_conversion).round(3),
        bus_split_circle=True,
        line_width=line_flow.div(unit_conversion).round(3),
        line_flow=line_flow.div(unit_conversion).round(3),
        line_color="rosybrown",
        link_width=link_flow.div(unit_conversion).round(3),
        link_flow=link_flow.div(unit_conversion).round(3),
        link_color=branch_color,
        arrow_size_factor=arrow_size_factor,
        tooltip=tooltip,
        auto_scale=True,
        branch_width_max=branch_width_max,
        bus_size_max=bus_size_max,
        map_style=map_style,
    )

    map.layers.insert(0, regions_layer)

    map.to_html(snakemake.output[0], offline=False)
