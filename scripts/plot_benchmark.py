# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Plot benchmark figures.
"""

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    get_version,
    set_scenario_config,
)
from scripts.clean_tyndp_benchmark import _convert_units
from scripts.make_benchmark import load_data, match_temporal_resolution

logger = logging.getLogger(__name__)
plt.style.use("bmh")


def add_version(ax: plt.Axes, fig: plt.Figure):
    version = get_version()
    fig.canvas.draw()
    bbox_fig = fig.get_tightbbox(fig.canvas.get_renderer())
    fig_width_inches, fig_height_inches = fig.get_size_inches()
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


def _plot_scenario_comparison(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_col: list[str],
    source_unit: str,
):
    fig, ax = plt.subplots(figsize=(12, 8))

    table_title = table.replace("_", " ").title()
    idx = [model_col] + [c for c in rfc_col if c in df.columns]
    df.set_index("carrier")[idx].plot.bar(
        ax=ax,
        color=["#1f77b4", "#ff7f0e", "#aeff39"],
        width=0.7,
        xlabel="",
        ylabel=f"{table_title} [{source_unit}]",
        title=f"{table_title} - EU27 - Scenario {scenario} - CY {cyear} - Year {year}",
    )
    ax.tick_params(axis="x", labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.legend(facecolor="white", frameon=True, edgecolor="black")

    for c in ax.containers:
        ax.bar_label(c, fmt="%.0f", padding=3, fontsize=8)

    add_version(ax, fig)

    output_filename = Path(output_dir, f"benchmark_{table}_eu27_cy{cyear}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def _plot_time_series(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
    scenario: str,
    cyear: int,
    model_col: str,
    rfc_col: str,
    source_unit: str,
    colors: dict,
):
    fig, ax = plt.subplots(figsize=(12, 8))

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
            c=colors.get(carrier, "grey"),
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

    add_version(ax, fig)

    output_filename = Path(output_dir, f"benchmark_{table}_eu27_cy{cyear}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def plot_benchmark(
    table: str,
    benchmarks_raw: pd.DataFrame,
    output_dir: str,
    scenario: str,
    cyear: int,
    options: dict,
    colors: dict,
    model_col: str = "Open-TYNDP",
    rfc_col: list[str] = ["TYNDP 2024", "TYNDP 2024 VP"],
):
    """
    Create benchmark comparison figures and export one file per year.

    Parameters
    ----------
    table : str
        Benchmark table to plot.
    benchmarks_raw: pd.DataFrame
        Combined DataFrame containing both model and reference data.
    output_dir: str
        Output directory.
    scenario: str
        Scenario name.
    cyear : int
        Climate year.
    options : dict
        Full benchmarking configuration containing table units and conversions.
    colors : dict,
        Dictionary of colors to be used for each technology.
    model_col : str, default "Open-TYNDP"
        Column name for model values.
    rfc_col : list[str], default ["TYNDP 2024", "TYNDP 2024 Vis Plfm"]
        Column names for reference values.
    """

    # Parameters
    opt = options["tables"][table]
    table_type = opt["table_type"]
    source_unit = opt["unit"]
    unit_conversion = options["unit_conversion"]
    unit_conversion = {
        k: 1 / v for k, v in unit_conversion.items()
    }  # Inverse conversion factor to revert unit

    # Filter data and Convert back to source unit
    logger.info(f"Making benchmark for {table} using {rfc_col} and {model_col}")
    benchmarks_raw = benchmarks_raw.query("table==@table").dropna(how="all", axis=1)
    benchmarks = _convert_units(benchmarks_raw, source_unit, unit_conversion)

    available_columns = [
        c for c in benchmarks.columns if c not in ["value", "source", "unit"]
    ]
    bench_wide = benchmarks.pivot_table(
        index=available_columns, values="value", columns="source", dropna=False
    )

    # Check if at least two sources are available to compare
    if len(bench_wide.columns) < 2:
        logging.info(f"Skipping table {table}, need at least two sources to compare.")
        return

    for year in bench_wide.index.get_level_values("year").unique():
        bench_year = bench_wide.query("year==@year").copy()

        if table_type == "scenario_comparison":
            _plot_scenario_comparison(
                bench_year.reset_index(),
                table,
                year,
                output_dir,
                scenario,
                cyear,
                model_col,
                rfc_col,
                source_unit,
            )
        elif table_type == "time_series":
            rfc_col_str = [c for c in rfc_col if c in bench_year.columns][0]
            bench_agg = match_temporal_resolution(
                bench_year, model_col, rfc_col_str
            ).reset_index()
            _plot_time_series(
                bench_agg,
                table,
                year,
                output_dir,
                scenario,
                cyear,
                model_col,
                rfc_col_str,
                source_unit,
                colors,
            )


def plot_overview(
    indicators: pd.DataFrame, fn: str, scenario: str, metric: str = "sMAPE"
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
    metric : str, default "sMAPE"
        Metric to plot.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Keep relevant indicators and rows
    df_clean = indicators[[metric, "Missing"]].dropna()
    df_clean.index = df_clean.index.str.replace("_", " ").str.title()

    # Create bar plot with metric
    df_clean.plot.bar(
        ax=ax,
        y=metric,
        width=0.7,
        xlabel="",
        ylabel=metric,
        title=f"Comparison of Open-TYNDP and TYNDP 2024 results for EU27, CY {cyear} and {scenario} scenario\n{metric} accuracy indicator (a lower error is better)",
        legend=True,
    )

    # Add missing carriers information
    df_clean.plot(
        ax=ax,
        y="Missing",
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

    add_version(ax, fig)

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
    colors = snakemake.params["colors"]
    scenario = snakemake.params["scenario"]
    cyear = get_snapshots(snakemake.params.snapshots)[0].year
    benchmarks_fn = snakemake.input.benchmarks
    data_vp_fn = snakemake.input.data_vp
    results_fn = snakemake.input.results
    output_dir = Path(snakemake.output.dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    kpis_in = snakemake.input.kpis
    kpis_out = snakemake.output.kpis

    # Load data
    benchmarks_raw = load_data(
        benchmarks_fn, results_fn, "TYNDP " + scenario, data_vp_fn
    )

    # Produce benchmark figures
    logger.info("Producing benchmark figures")

    tqdm_kwargs = {
        "ascii": False,
        "unit": " figure",
        "total": len(options["tables"]),
        "desc": "Producing benchmark figures",
    }

    func = partial(
        plot_benchmark,
        benchmarks_raw=benchmarks_raw,
        output_dir=output_dir,
        scenario=scenario,
        cyear=cyear,
        options=options,
        colors=colors,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        results = list(tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs))

    # Plot overview
    indicators = pd.read_csv(kpis_in, index_col=0)
    plot_overview(indicators, kpis_out, scenario)

    logger.info("Benchmark plotting completed successfully")
