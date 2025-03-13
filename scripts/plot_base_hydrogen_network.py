import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pypsa
import os
import yaml
import matplotlib.animation as animation
import geopandas as gpd
import cartopy.crs as ccrs
from pypsa.plot import add_legend_lines, add_legend_patches, add_legend_circles
from plot_hydrogen_network import load_projection
plt.style.use(["ggplot"])

from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

def plot_h2_map_base(network):
    n = network.copy()

    linewidth_factor = 4e3

    n.links.drop(
        n.links.index[~n.links.carrier.str.contains("H2 pipeline")], inplace=True
    )
    h2_pipes = n.links[n.links.carrier == "H2 pipeline"].p_nom
    link_widths_total = h2_pipes / linewidth_factor
    if link_widths_total.notnull().empty:
        logger.info("No base H2 pipeline capacities to plot.")
        return
    link_widths_total = link_widths_total.reindex(n.links.index).fillna(0.0)

    proj = load_projection(dict(name="EqualEarth"))
    n.buses.drop(
        n.buses.index[~n.buses.carrier.str.contains("H2")], inplace=True
    )

    logger.info("Plotting base H2 pipeline capacities.")
    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})
    color_h2_pipe = "#499a9c"
    color_h2_node = "#ff29d9"

    n.plot(
        geomap=True,
        bus_sizes=0.1,
        bus_colors=color_h2_node,
        link_colors=color_h2_pipe,
        link_widths=link_widths_total,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / 4e3
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.01, 1.07),
        frameon=False,
        ncol=2,
        labelspacing=0.8,
        handletextpad=1,
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.25, 1.13),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes=[0.2],
        labels=["H2 Node"],
        srid=n.srid,
        patch_kw=dict(facecolor=color_h2_node),
        legend_kw=legend_kw,
    )

    colors = [color_h2_pipe]
    labels = ["H2 Pipeline"]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=1,
        frameon=False,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)

    ax.set_facecolor("white")

    plt.savefig(snakemake.output.map, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_base_hydrogen_network",
            opts="",
            clusters="100",
            sector_opts="",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    map_opts = snakemake.params.plotting["map"]

    if map_opts["boundaries"] is None:
        regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
        map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

    proj = load_projection(snakemake.params.plotting)

    plot_h2_map_base(n)