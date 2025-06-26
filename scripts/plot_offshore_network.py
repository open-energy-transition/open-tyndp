# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Plot offshore transmission network.
"""

import logging

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from _helpers import configure_logging, set_scenario_config
from plot_power_network import load_projection
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

plt.style.use(["ggplot"])


logger = logging.getLogger(__name__)


def plot_offshore_map(network, map_opts, map_fn, carrier="DC", expanded=False):
    """
    Plots the offshore network hydrogen and electricity capacities and offshore-hubs buses.
    If expanded is enabled, the optimal capacities are plotted instead.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network for plotting the offshore grid. Can be either presolving or post solving.
    map_opts : dict
        Map options for plotting.
    map_fn : str
        Path to save the final map plot to.
    carrier : str, optional
        Carrier to plot
    expanded : bool, optional
        Whether to plot expanded capacities. Defaults to plotting only base network (p_nom).

    Returns
    -------
    None
        Saves the map plot as figure.
    """
    n = network.copy()

    linewidth_factor = 1e4

    mask = {"DC": "Offshore DC", "H2": "Offshore H2 pipeline"}

    n.links.drop(
        n.links.index[~(n.links.index.str.contains(mask[carrier]))],
        inplace=True,
    )

    p_nom = "p_nom_opt" if expanded else "p_nom"
    # transmission capacities
    links = n.links[n.links.index.str.contains(mask[carrier])][p_nom]

    # set link widths
    link_widths = links / linewidth_factor
    if link_widths.notnull().empty:
        logger.info(f"No offshore capacities for {carrier}, skipping plot.")
        return
    link_widths = link_widths.reindex(n.links.index).fillna(0.0)

    # keep relevant buses
    bus_carriers = [carrier] + (["AC"] if carrier == "DC" else [])
    n.buses.drop(
        n.buses.index[
            (~n.buses.carrier.isin(bus_carriers))
            | (n.buses.index.str.contains("Z1"))
            | (n.buses.index.str.contains("DRES"))
        ],
        inplace=True,
    )
    n_oh = n.copy()
    n_oh.buses.drop(
        n_oh.buses.index[~n_oh.buses.index.str.contains("OH")], inplace=True
    )

    # plot transmission network
    logger.info("Plotting offshore transmission network.")
    proj = load_projection(dict(name="EqualEarth"))
    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})
    color_h2 = "#f081dc"
    color_dc = "darkseagreen"
    color = color_dc if carrier == "DC" else color_h2
    color_oh_nodes = "#ff29d9"
    color_hm_nodes = "darkgray"

    n.plot(
        geomap=True,
        bus_sizes=0.05,
        bus_colors=color_hm_nodes,
        link_colors=color,
        link_widths=link_widths,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n_oh.plot(
        geomap=True,
        bus_sizes=0.05,
        bus_colors=color_oh_nodes,
        branch_components=[],
        ax=ax,
        **map_opts,
    )

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.32, 1.13),
        frameon=False,
        ncol=1,
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
        bbox_to_anchor=(0.55, 1.13),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes=[0.1],
        labels=["Home market"],
        srid=n.srid,
        patch_kw=dict(facecolor=color_hm_nodes),
        legend_kw=legend_kw,
    )

    legend_kw["bbox_to_anchor"] = (0.55, 1.08)

    add_legend_circles(
        ax,
        sizes=[0.1],
        labels=["Offshore hubs"],
        srid=n.srid,
        patch_kw=dict(facecolor=color_oh_nodes),
        legend_kw=legend_kw,
    )

    label = "DC link" if carrier == "DC" else "H2 pipeline"

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=1,
        frameon=False,
    )

    add_legend_patches(ax, color, label, legend_kw=legend_kw)

    ax.set_facecolor("white")

    plt.savefig(map_fn, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_offshore_network",
            opts="",
            clusters="all",
            sector_opts="",
            planning_horizons=2050,
            carrier="DC",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    map_opts = snakemake.params.plotting["map"]

    if map_opts["boundaries"] is None:
        regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
        map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

    proj = load_projection(snakemake.params.plotting)
    map_fn = snakemake.output.map

    plot_offshore_map(
        n,
        map_opts,
        map_fn,
        carrier=snakemake.wildcards.carrier,
        expanded=snakemake.params.expanded,
    )
