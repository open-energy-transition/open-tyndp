# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Plot offshore transmission network with existing capacities. If `expanded` is enabled, the optimal capacities are plotted instead.
"""

import logging
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from pypsa.plot.maps.static import add_legend_circles, add_legend_lines, add_legend_patches

from scripts._helpers import configure_logging, set_scenario_config
from scripts.plot_power_network import load_projection

plt.style.use(["ggplot"])


logger = logging.getLogger(__name__)


def plot_offshore_map(
    network,
    map_opts,
    map_fn,
    planning_horizons,
    carrier="DC_OH",
    p_nom="p_nom",
    legend=True,
    hubs_only=False,
):
    """
    Plots the offshore network hydrogen or electricity capacities and offshore-hubs buses.
    If `p_nom` parameter is set as `p_nom_opt`, optimal capacities are plotted instead.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network for plotting the offshore grid. Can be either presolving or post solving.
    map_opts : dict
        Map options for plotting.
    map_fn : str
        Path to save the final map plot to.
    planning_horizons : int
        The planning horizon year
    carrier : str, optional
        Carrier to plot
    p_nom : str | float, optional
        Nominal power parameter for determining link thickness. If str, must be "p_nom" or "p_nom_opt".
        If float, uses fixed value for all links. Defaults to plotting only base network (p_nom).
    legend : bool, optional
        Whether to display a legend on the plot. Defaults to display the legend.
    hubs_only : bool, optional
        Whether to only plot the offshore hubs. Defaults to plot both home market nodes and offshore hubs.

    Returns
    -------
    None
        Saves the map plot as figure to the provided map_fn path.
    """
    n = network.copy()

    lw_factor = 1e4 if carrier == "DC_OH" else 5e3
    link_lower_threshold = 1e2  # MW below which not drawn

    n.links.drop(
        n.links.index[n.links.carrier != carrier],
        inplace=True,
    )

    # transmission capacities
    if isinstance(p_nom, str):
        links = (
            n.links[n.links.carrier == carrier][p_nom]
            .rename(index=lambda x: re.sub(r"-\d{4}$", f"-{planning_horizons}", x))
            .groupby(level=0)
            .sum()
        )
        # set link widths
        links[links < link_lower_threshold] = 0.0
        link_widths = links / lw_factor
        if link_widths.notnull().empty:
            logger.info(f"No offshore capacities for {carrier}, skipping plot.")
            return
        link_widths = link_widths.reindex(n.links.index).fillna(0.0)
    elif isinstance(p_nom, float) or isinstance(p_nom, int):
        link_widths = p_nom
    else:
        raise ValueError("Parameter 'p_nom' must be either str or float.")

    # keep relevant buses
    bus_carriers = [carrier.replace("DC", "AC")] + (
        ["AC"] if carrier == "DC_OH" else ["H2", "H2_OH"]
    )
    n.buses.drop(
        n.buses.index[~n.buses.carrier.isin(bus_carriers)],
        inplace=True,
    )
    n_oh = n.copy()
    n_oh.buses.drop(
        n_oh.buses.index[~n_oh.buses.carrier.str.contains("OH")], inplace=True
    )

    if hubs_only:
        n.buses = n_oh.buses

    # plot transmission network
    logger.info("Plotting offshore transmission network.")
    proj = load_projection(dict(name="EqualEarth"))
    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})
    color_h2 = "#f081dc"
    color_dc = "darkseagreen"
    color = color_dc if carrier == "DC_OH" else color_h2
    color_oh_nodes = "#ff29d9"
    color_hm_nodes = "darkgray"

    n.plot.map(
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

    if legend:
        sizes = [30, 10]
        labels = [f"{s} GW" for s in sizes]
        scale = 1e3 / lw_factor
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

        label = "DC link" if carrier == "DC_OH" else "H2 pipeline"

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
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_offshore_network",
            opts="",
            clusters="all",
            sector_opts="",
            planning_horizons=2050,
            carrier="DC_OH",
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

    p_nom = "p_nom_opt" if snakemake.params.expanded else "p_nom"

    plot_offshore_map(
        n,
        map_opts,
        map_fn,
        snakemake.wildcards.planning_horizons,
        carrier=snakemake.wildcards.carrier,
        p_nom=p_nom,
    )
