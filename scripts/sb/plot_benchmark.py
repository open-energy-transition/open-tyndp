# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Plot benchmark figures.
"""

import logging
import multiprocessing as mp
import textwrap
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    convert_units,
    get_snapshots,
    get_version,
    set_scenario_config,
)
from scripts.sb.make_benchmark import (
    SOURCES_MAP,
    get_bus_col_name,
    load_data,
    match_temporal_resolution,
)

logger = logging.getLogger(__name__)
plt.style.use("bmh")

FIGURE_WIDTH_DEFAULT = 12
FIGURE_HEIGHT_DEFAULT = 8


def add_metadata(
    ax: plt.Axes,
    fig: plt.Figure,
    model_col: str = "",
    rfc_source: str = "",
    rfc_cols: list[str] = [],
):
    # Version
    version = get_version()
    fig.canvas.draw()
    bbox_fig = fig.get_tightbbox(fig.canvas.get_renderer())
    fig_width_inches, fig_height_inches = fig.get_size_inches()
    x0_fig = (
        bbox_fig.x0 / fig_width_inches
    )  # Convert bbox coordinates from inches to figure coordinates
    x1_fig = (
        bbox_fig.x1 / fig_width_inches
    )  # Convert bbox coordinates from inches to figure coordinates
    y0_fig = bbox_fig.y0 / fig_height_inches

    ax.text(
        x1_fig,
        y0_fig - 0.05,
        f"version: {version}",
        transform=fig.transFigure,
        ha="right",
        va="bottom",
        fontsize=8,
        alpha=0.7,
    )

    # Reference source
    if model_col != "" and rfc_source != "":
        additional_sources = (
            "" if len(rfc_cols) <= 1 else " Other sources shown for comparison."
        )
        ax.text(
            x0_fig,
            y0_fig - 0.05,
            f"Model outputs ({model_col}) benchmarked against {rfc_source}.{additional_sources}",
            transform=fig.transFigure,
            ha="left",
            va="bottom",
            fontsize=8,
            alpha=0.7,
        )


def _plot_scenario_comparison(
    df: pd.DataFrame,
    table: str,
    year: int,
    bus: str,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_cols: list[str],
    rfc_source: str,
    source_unit: str,
    bench_colors: dict,
):
    table_title = table.replace("_", " ").title()
    if table == "power_generation":
        table_title += " (pre-curtailment)"
    idx = [model_col, rfc_source] + [
        c for c in rfc_cols if c in df.columns and c != rfc_source
    ]

    tyndp_str = "TYNDP 2024 Scenarios Report"
    if "TYNDP 2024 Vis Pltfm" in idx and tyndp_str in idx:
        tyndp_str_ext = "TYNDP 2024 Scenarios Report"
        idx = [tyndp_str_ext if i == tyndp_str else i for i in idx]
        df = df.rename(columns={tyndp_str: tyndp_str_ext})

    # Wrap long x-axis labels
    df = df.set_index("carrier")
    df.index = [textwrap.fill(label, width=30) for label in df.index]

    fig_width = FIGURE_WIDTH_DEFAULT + max((df.shape[0] - FIGURE_WIDTH_DEFAULT) * 3, 0)
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_DEFAULT, FIGURE_HEIGHT_DEFAULT))

    bar_colors = [bench_colors.get(col, "grey") for col in idx]
    df[idx].plot.bar(
        ax=ax,
        color=bar_colors,
        width=0.7,
        xlabel="",
        ylabel=f"{table_title} [{source_unit}]",
        title=f"{table_title} - {bus} - Scenario {scenario} - CY {cyear} - Year {year}",
    )
    ax.tick_params(axis="x", labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")
    legend = ax.legend(facecolor="white", frameon=True, edgecolor="black")
    for txt in legend.get_texts():
        if txt.get_text() in [model_col, rfc_source]:
            txt.set_fontweight("bold")

    for c in ax.containers:
        max_val = df[idx].replace(np.inf, np.nan).fillna(0).max().max()
        rotation = 90 if max_val > 100 and fig_width > FIGURE_WIDTH_DEFAULT else 0
        ax.bar_label(
            c,
            fmt=lambda x: f"{x:.1f}",
            padding=3,
            fontsize=8,
            rotation=rotation,
        )
    if rotation != 0:
        ax.set_ylim(top=ax.get_ylim()[1] * 1.05)

    add_metadata(ax, fig, model_col=model_col, rfc_source=rfc_source, rfc_cols=rfc_cols)

    output_filename = Path(
        output_dir, f"benchmark_{table}_{bus.replace(' ', '_')}_cy{cyear}_{year}.pdf"
    )
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def _plot_time_series(
    df: pd.DataFrame,
    table: str,
    year: int,
    bus: str,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_col: str,
    source_unit: str,
    tech_colors: dict,
):
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_DEFAULT, FIGURE_HEIGHT_DEFAULT))

    # Remove rows where either value is NaN
    df_clean = df.dropna(subset=[model_col, rfc_col])

    if df_clean.empty:
        logger.warning(f"No valid data points for time series plot: {table} {year}")
        plt.close(fig)
        return

    # Create scatter plot
    df.carrier = df.carrier.astype("category")
    carriers = df.carrier.cat.categories

    for i, carrier in enumerate(carriers):
        carrier_data = df[df.carrier == carrier]
        ax.scatter(
            carrier_data[model_col],
            carrier_data[rfc_col],
            c=tech_colors.get(carrier, "grey"),
            label=carrier,
            edgecolors="grey",
            s=20,
            alpha=0.7,
        )

    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    # Add 1:1 diagonal line for reference
    min_val = min(df_clean[rfc_col].min(), df_clean[model_col].min())
    max_val = max(df_clean[rfc_col].max(), df_clean[model_col].max())
    ax.plot(
        [min_val, max_val],
        [min_val, max_val],
        "r--",
        linewidth=2,
    )
    ax.set_xlim(min_val)
    ax.set_ylim(min_val)

    correlation = df_clean[model_col].corr(df_clean[rfc_col])
    table_title = table.replace("_", " ").title()
    ax.set_title(
        f"{table_title} - EU27 - Scenario {scenario} - CY {cyear} - Year {year}"
    )
    ax.set_xlabel(f"{rfc_col} [{source_unit}]")
    ax.set_ylabel(f"{model_col} [{source_unit}]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)

    # Add text box with statistics
    n_points = len(df_clean)
    missing_series = len(df[df[model_col].isna()].carrier.unique())
    textstr = (
        f"n = {n_points}\nR² = {correlation**2:.3f}\nmissing series = {missing_series}"
    )
    props = dict(boxstyle="round", facecolor="white")
    ax.text(
        0.05,
        0.95,
        textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    add_metadata(ax, fig)

    output_filename = Path(
        output_dir, f"benchmark_{table}_{bus.replace(' ', '_')}_cy{cyear}_{year}.pdf"
    )
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def _plot_prices(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_source: str,
    source_unit: str,
    bench_colors: dict,
    eps: float = 1e-6,
):
    fig, ax = plt.subplots(
        figsize=(FIGURE_WIDTH_DEFAULT * 1.7, FIGURE_HEIGHT_DEFAULT * 0.7)
    )
    table_title = (
        table.replace("_", " ").replace("excl shed", "excl. load shedding").title()
    )
    bar_colors = [bench_colors.get(col, "grey") for col in [model_col, rfc_source]]
    df.index = df.index.get_level_values("spatial")
    df[[model_col, rfc_source]].plot.bar(
        title=f"{table_title} - Scenario {scenario} - CY {cyear} - Year {year}",
        ylabel=source_unit,
        color=bar_colors,
        ax=ax,
    )

    ax.grid(axis="y", linestyle="--")

    errors = (
        (df[model_col] - df[rfc_source]).abs()
        / ((df[model_col].abs() + df[rfc_source].abs()) / 2 + eps)
        * 100
    )

    ax2 = ax.twinx()
    ax2.grid(False)
    ax2.plot(
        ax.get_xticks(),
        errors.values,
        color="red",
        marker="o",
        linestyle="",
        label=f"sMAPE [%] (Overall sMAPE {errors.mean():.1f}% / sMdAPE {errors.median():.1f}%)",
    )
    ax2.set_ylim(0, max(ax2.get_ylim()[1], 100))
    ax2.set_ylabel("%")

    means = df.mean()
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    labels = [f"{label} (Avg. {means[label]:.1f} {source_unit})" for label in labels]
    ax.legend(
        handles + handles2,
        labels + labels2,
        loc="upper left",
        frameon=True,
        facecolor="white",
    )

    add_metadata(ax, fig, model_col=model_col, rfc_source=rfc_source)

    output_filename = Path(output_dir, f"benchmark_{table}_cy{cyear}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def _plot_flows(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_source: str,
    source_unit: str,
    bench_colors: dict,
):
    fig, ax = plt.subplots(
        figsize=(
            FIGURE_WIDTH_DEFAULT,
            FIGURE_HEIGHT_DEFAULT * np.ceil(df.shape[0] / 45),
        )
    )
    table_title = (
        table.replace("_", " ")
        .replace("crossborder", "cross-border exchanges for")
        .title()
    )
    bar_colors = [bench_colors.get(col, "grey") for col in [model_col, rfc_source]]
    df.index = df.index.get_level_values("spatial")

    # remove flows between identical locations
    df = df[df.index.str.split("->").map(lambda x: x[0] != x[1])]

    df[[model_col, rfc_source]].sort_index(ascending=False).plot.barh(
        title=f"{table_title} - Scenario {scenario} - CY {cyear} - Year {year}",
        xlabel=source_unit,
        color=bar_colors,
        ax=ax,
    )

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis="x", which="major", linestyle="-")
    ax.grid(axis="x", which="minor", linestyle="--", alpha=0.7)
    ax.grid(axis="y", linestyle="-", alpha=0.7)

    ax.legend(
        frameon=True,
        facecolor="white",
    )

    add_metadata(ax, fig, model_col=model_col, rfc_source=rfc_source)

    output_filename = Path(output_dir, f"benchmark_{table}_cy{cyear}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def plot_benchmark(
    table: str,
    bus: str,
    benchmarks_raw: pd.DataFrame,
    output_dir: str,
    scenario: str,
    snapshots: dict[str, str],
    options: dict,
    tech_colors: dict,
    bench_colors: dict,
    model_col: str = "Open-TYNDP",
    bus_col_name: str = "bus",
):
    """
    Create benchmark comparison figures and export one file per year.

    Parameters
    ----------
    table : str
        Benchmark table to plot.
    bus : str
        Bus of the current figure.
    benchmarks_raw: pd.DataFrame
        Combined DataFrame containing both model and reference data.
    output_dir: str
        Output directory.
    scenario: str
        Scenario name.
    snapshots : dict[str, str]
        Dictionary defining the temporal range with 'start' and 'end' keys.
    options : dict
        Full benchmarking configuration containing table units and conversions.
    tech_colors : dict
        Dictionary mapping technology/carrier names to colors for time series.
    bench_colors : dict
        Dictionary mapping data source names to colors for scenario comparisons.
    model_col : str, default "Open-TYNDP"
        Column name for model values.
    bus_col_name : str, default "bus"
        Bus column name.
    """

    # Parameters
    opt = options["tables"][table]
    table_type = opt["table_type"]
    source_unit = "TWh" if "crossborder" in table else opt["unit"]
    rfc_cols = [SOURCES_MAP.get(s, s) for s in opt["rfc_sources"]]
    rfc_source = rfc_cols[0]
    cyear = get_snapshots(snapshots)[0].year
    bus_col_name = get_bus_col_name(bus_col_name, table)

    # Filter data and Convert back to source unit
    logger.debug(
        f"Making benchmark for {table} at {bus} using {rfc_cols} and {model_col}"
    )

    filter_by_bus = "price" not in table and "crossborder" not in table
    condition_str = f" and {bus_col_name}==@bus" if filter_by_bus else ""
    benchmarks = (
        benchmarks_raw.query(f"table==@table{condition_str}")
        .dropna(how="all", axis=1)
        .assign(spatial=lambda df: df[bus_col_name])
        .drop(columns=["bus", "country", "border", "corridor"], errors="ignore")
    )
    op = "sum" if "price" not in table else "mean"
    benchmarks = (
        benchmarks.groupby([c for c in benchmarks.columns if c != "value"])
        .value.agg(op)
        .reset_index()
        .assign(unit=opt["unit"])
    )

    if benchmarks.empty:
        logger.warning(
            f"No data available for table '{table}' and bus {bus} in Open-TYNDP or TYNDP 2024 datasets"
        )
        return

    if "price" not in table:
        benchmarks.loc[benchmarks.table.str.contains("crossborder"), "unit"] = "TWh"
        benchmarks = convert_units(benchmarks, invert=True)

    available_columns = [
        c for c in benchmarks.columns if c not in ["value", "source", "unit"]
    ]
    bench_wide = benchmarks.pivot_table(
        index=available_columns, values="value", columns="source", dropna=False
    )

    # Check if at least two sources are available to compare
    if len(bench_wide.columns) < 2:
        logger.info(
            f"Skipping table {table} at {bus}, need at least two sources to compare."
        )
        return

    for year in bench_wide.index.get_level_values("year").unique():
        bench_year = bench_wide.query("year==@year").copy()

        if table_type == "scenario_comparison":
            _plot_scenario_comparison(
                df=bench_year.reset_index(),
                table=table,
                year=year,
                bus=bus,
                output_dir=output_dir,
                scenario=scenario,
                cyear=cyear,
                model_col=model_col,
                rfc_cols=[c for c in rfc_cols if c in bench_year.columns],
                rfc_source=rfc_source,
                source_unit=source_unit,
                bench_colors=bench_colors,
            )
        elif table_type == "time_series":
            rfc_col_str = [c for c in rfc_cols if c in bench_year.columns][0]
            bench_agg = match_temporal_resolution(
                bench_year, snapshots, model_col, rfc_col_str
            ).reset_index()
            _plot_time_series(
                df=bench_agg,
                table=table,
                year=year,
                bus=bus,
                output_dir=output_dir,
                scenario=scenario,
                cyear=cyear,
                model_col=model_col,
                rfc_col=rfc_col_str,
                source_unit=source_unit,
                tech_colors=tech_colors,
            )
        elif table_type == "prices":
            _plot_prices(
                df=bench_year,
                table=table,
                year=year,
                output_dir=output_dir,
                scenario=scenario,
                cyear=cyear,
                model_col=model_col,
                rfc_source=rfc_source,
                source_unit=source_unit,
                bench_colors=bench_colors,
            )
        elif table_type == "flows":
            _plot_flows(
                df=bench_year,
                table=table,
                year=year,
                output_dir=output_dir,
                scenario=scenario,
                cyear=cyear,
                model_col=model_col,
                rfc_source=rfc_source,
                source_unit=source_unit,
                bench_colors=bench_colors,
            )
        else:
            raise ValueError(f"Unknown table type {table_type}.")


def orchestrate_benchmark(
    bus_col_name: str,
    benchmarks_raw: pd.DataFrame,
    output_dir: str,
    scenario: str,
    snapshots: dict[str, str],
    options: dict,
    tech_colors: dict,
    bench_colors: dict,
    threads: int,
):
    logger.info(f"Producing benchmark figures by {bus_col_name}")

    output_dir_bus_col = Path(output_dir, f"by_{bus_col_name}")
    output_dir_bus_col.mkdir(parents=True, exist_ok=True)

    table_bus_col_pairs = [
        (table, bus_col)
        for table in options["tables"]
        for bus_col in (
            [""]
            if "price" in table or "crossborder" in table
            else benchmarks_raw.query("table == @table")[bus_col_name].unique()
        )
    ]

    tqdm_kwargs = {
        "ascii": False,
        "unit": " figure",
        "total": len(table_bus_col_pairs),
        "desc": "Producing benchmark figures",
    }

    func = partial(
        plot_benchmark,
        benchmarks_raw=benchmarks_raw,
        output_dir=output_dir_bus_col,
        scenario=scenario,
        snapshots=snapshots,
        options=options,
        tech_colors=tech_colors,
        bench_colors=bench_colors,
        bus_col_name=bus_col_name,
    )

    with mp.Pool(processes=threads) as pool:
        list(tqdm(pool.starmap(func, table_bus_col_pairs), **tqdm_kwargs))


def plot_overview(
    indicators: pd.DataFrame,
    fn: str,
    scenario: str,
    snapshots: dict[str, str],
    metric: str = "sMAPE",
    bus_col_name: str = "bus",
):
    """
    Plot benchmark overview figure.

    Parameters
    ----------
    indicators : pd.DataFrame
        Indicators DataFrame.
    fn : str
        Output filename.
    scenario : str
        Scenario name.
    snapshots : dict[str, str]
        Dictionary defining the temporal range with 'start' and 'end' keys.
    metric : str, default "sMAPE"
        Metric to plot.
    bus_col_name : str, default "bus"
        Bus column name.
    """
    fig, ax = plt.subplots(figsize=(12, FIGURE_HEIGHT_DEFAULT))
    cyear = get_snapshots(snapshots)[0].year

    # Keep relevant indicators and rows
    df_clean = indicators[[metric, "Missing carriers"]].dropna()
    df_clean.index = df_clean.index.str.replace("_", " ").str.title()

    # Wrap long x-axis labels
    df_clean.index = [textwrap.fill(label, width=30) for label in df_clean.index]

    # Create bar plot with metric
    df_clean.plot.bar(
        ax=ax,
        y=metric,
        width=0.7,
        xlabel="",
        ylabel=metric,
        title=f"Comparison of Open-TYNDP and TYNDP 2024 outputs by {bus_col_name}, CY {cyear} and {scenario} scenario\n{metric} accuracy indicator (a lower error is better)",
        legend=True,
        ylim=[0, max(df_clean[metric].max() + 0.1, 1)],
    )

    # Add missing carriers information
    df_clean.plot(
        ax=ax,
        y="Missing carriers",
        secondary_y=True,
        mark_right=True,
        legend=True,
        marker=".",
        color="red",
        markersize=10,
        linestyle="None",
        ylabel="Missing carriers",
        ylim=0,
    )

    ax.tick_params(axis="x", labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")

    # Combine legends from both axes
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax.right_ax.get_legend_handles_labels()
    ax.legend(
        h1 + h2,
        l1 + l2,
        facecolor="white",
        frameon=True,
        edgecolor="black",
        loc="upper left",
    )

    add_metadata(ax, fig)

    fig.savefig(fn, bbox_inches="tight")

    plt.close(fig)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_benchmark",
            opts="",
            clusters="all",
            sector_opts="",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    tech_colors = snakemake.params["tech_colors"]
    bench_colors = snakemake.params["bench_colors"]
    scenario = snakemake.params["scenario"]
    snapshots = snakemake.params.snapshots
    benchmarks_fn = snakemake.input.benchmarks
    vp_data_fn = snakemake.input.vp_data
    mm_data_fn = snakemake.input.mm_data
    results_fn = snakemake.input.results
    output_dir = Path(snakemake.output.dir)
    threads = snakemake.threads

    # Load data
    benchmarks_raw = load_data(
        benchmarks_fn=benchmarks_fn,
        results_fn=results_fn,
        scenario="TYNDP " + scenario,
        vp_data_fn=vp_data_fn,
        mm_data_fn=mm_data_fn,
    )

    # Produce benchmark figures
    for bus_col_name, kpis_in, kpis_out, enabled in [
        (
            "bus",
            snakemake.input.kpis_by_bus,
            snakemake.output.kpis_by_bus,
            options["spatial"]["by_bus"],
        ),
        (
            "country",
            snakemake.input.kpis_by_country,
            snakemake.output.kpis_by_country,
            options["spatial"]["by_country"],
        ),
    ]:
        if enabled:
            orchestrate_benchmark(
                bus_col_name=bus_col_name,
                benchmarks_raw=benchmarks_raw,
                output_dir=output_dir,
                scenario=scenario,
                snapshots=snapshots,
                options=options,
                tech_colors=tech_colors,
                bench_colors=bench_colors,
                threads=threads,
            )
            indicators = pd.read_csv(kpis_in, index_col=0)
            plot_overview(
                indicators, kpis_out, scenario, snapshots, bus_col_name=bus_col_name
            )
        else:
            Path(kpis_out).touch()

    logger.info("Benchmark plotting completed successfully")
