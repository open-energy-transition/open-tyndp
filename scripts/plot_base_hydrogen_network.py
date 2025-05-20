# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates map of base hydrogen network, added before the optimization.
"""

import logging

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pypsa
from _helpers import configure_logging, set_scenario_config
from plot_hydrogen_network import load_projection
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

plt.style.use(["ggplot"])

logger = logging.getLogger(__name__)


def group_import_corridors(df):
    """
    Group pipes which connect same buses and return overall capacity.
    """
    df = df.copy()

    # there are pipes for each investment period rename to AC buses name for plotting
    df["index_orig"] = df.index
    df.rename(index=lambda x: x.split(" - ")[0], inplace=True)
    return df.groupby(level=0).agg(
        {"p_nom": "sum", "p_nom_opt": "sum", "index_orig": "first"}
    )


def plot_h2_map_base(
    network, map_opts, map_fn, expanded=False, regions_for_storage=None
):
    """
    Plots the base hydrogen network pipelines capacities, hydrogen buses and import potentials.
    If expanded is enabled, the optimal capacities are plotted instead.
    If regions are given, hydrogen storage capacities are plotted for those regions with aggregated H2 tank storage
    and underground H2 cavern capacities.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network for plotting the hydrogen grid. Can be either presolving or post solving.
    map_opts : dict
        Map options for plotting.
    map_fn : str
        Path to save the final map plot to.
    expanded : bool, optional
        Whether to plot expanded capacities. Defaults to plotting only base network (p_nom).
    regions_for_storage : gpd.GeoDataframe, optional
        Geodataframe of regions to use for plotting hydrogen storage capacities. Index needs to match storage locations.
        If none is given, no hydrogen storage capacities are plotted.

    Returns
    -------
    None
        Saves the map plot as figure.
    """
    n = network.copy()

    linewidth_factor = 4e3

    n.links.drop(
        n.links.index[
            ~(
                n.links.carrier.str.contains("H2 pipeline")
                | n.links.carrier.str.contains("H2 import")
            )
        ],
        inplace=True,
    )

    p_nom = "p_nom_opt" if expanded else "p_nom"
    # capacity of pipes and imports
    h2_pipes = n.links[n.links.carrier == "H2 pipeline"][p_nom]
    h2_imports = n.links[n.links.carrier.str.contains("H2 import")]

    # group high and low import corridors together
    h2_imports = group_import_corridors(h2_imports)[p_nom]
    n.links.rename(index=lambda x: x.split(" - ")[0], inplace=True)
    # group links by summing up p_nom values and taking the first value of the rest of the columns
    other_cols = dict.fromkeys(n.links.columns.drop(["p_nom_opt", "p_nom"]), "first")
    n.links = n.links.groupby(level=0).agg(
        {"p_nom_opt": "sum", "p_nom": "sum", **other_cols}
    )

    # set link widths
    link_widths_pipes = h2_pipes / linewidth_factor
    link_widths_imports = h2_imports / linewidth_factor
    if link_widths_pipes.notnull().empty:
        logger.info("No base H2 pipeline capacities to plot.")
        return
    link_widths_pipes = link_widths_pipes.reindex(n.links.index).fillna(0.0)
    link_widths_imports = link_widths_imports.reindex(n.links.index).fillna(0.0)

    # drop non H2 buses
    n.buses.drop(n.buses.index[~n.buses.carrier.str.contains("H2")], inplace=True)

    # optionally add hydrogen storage capacities onto the map
    if regions_for_storage is not None:
        h2_storage = n.stores.query("carrier.str.contains('H2')")
        regions_for_storage["H2"] = (
            h2_storage.rename(index=h2_storage.bus.map(n.buses.location))
            .e_nom_opt.groupby(level=0)
            .sum()
            .div(1e6)
        )  # TWh
        regions_for_storage["H2"] = regions_for_storage["H2"].where(
            regions_for_storage["H2"] > 0.1
        )

    # plot H2 pipeline capacities and imports
    logger.info("Plotting base H2 pipeline and import capacities.")
    proj = load_projection(dict(name="EqualEarth"))
    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})
    color_h2_pipe = "#499a9c"
    color_h2_imports = "#FFA500"
    color_h2_node = "#ff29d9"

    n.plot(
        geomap=True,
        bus_sizes=0.1,
        bus_colors=color_h2_node,
        link_colors=color_h2_pipe,
        link_widths=link_widths_pipes,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    if regions_for_storage is not None:
        regions_for_storage = regions_for_storage.to_crs(proj.proj4_init)
        regions_for_storage.plot(
            ax=ax,
            column="H2",
            cmap="Blues",
            linewidths=0,
            legend=True,
            vmax=6,
            vmin=0,
            legend_kwds={
                "label": "Hydrogen Storage [TWh]",
                "shrink": 0.7,
                "extend": "max",
            },
        )

    if not h2_imports.empty:
        n.plot(
            geomap=True,
            bus_sizes=0,
            link_colors=color_h2_imports,
            link_widths=link_widths_imports,
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
        sizes=[0.2],
        labels=["H2 Node"],
        srid=n.srid,
        patch_kw=dict(facecolor=color_h2_node),
        legend_kw=legend_kw,
    )

    colors = (
        [color_h2_pipe, color_h2_imports] if not h2_imports.empty else [color_h2_pipe]
    )
    labels = ["H2 Pipeline", "H2 import"] if not h2_imports.empty else ["H2 Pipeline"]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=1,
        frameon=False,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)

    ax.set_facecolor("white")

    plt.savefig(map_fn, bbox_inches="tight")
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

    if n.buses.country.isin(["MA", "DZ"]).any():
        map_opts["boundaries"] = list(np.add(map_opts["boundaries"], [0, 0, -6, 0]))
    proj = load_projection(snakemake.params.plotting)
    map_fn = snakemake.output.map

    plot_h2_map_base(n, map_opts, map_fn)
