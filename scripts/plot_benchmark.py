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

from scripts._helpers import configure_logging, set_scenario_config
from scripts.clean_tyndp_benchmark import _convert_units
from scripts.make_benchmark import load_data, match_temporal_resolution

logger = logging.getLogger(__name__)
plt.style.use("bmh")


def _plot_scenario_comparison(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
    model_col: str,
    rfc_col: str,
    source_unit: str,
):
    fig, ax = plt.subplots(figsize=(12, 8))

    df.set_index("carrier")[[model_col, rfc_col]].plot.bar(
        ax=ax,
        color=["#1f77b4", "#ff7f0e"],
        width=0.7,
        rot=45,
        xlabel="",
    )
    table_title = table.replace("_", " ").title()
    ax.set_title(f"{table_title} - Year {year}")
    ax.set_ylabel(f"{table_title} [{source_unit}]")

    for c in ax.containers:
        ax.bar_label(c, fmt="%.0f", padding=3, fontsize=8)

    output_filename = Path(output_dir, f"benchmark_{table}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def _plot_time_series(
    df: pd.DataFrame,
    table: str,
    year: int,
    output_dir: str,
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
    ax.set_title(f"{table_title} - Year {year}")
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

    output_filename = Path(output_dir, f"benchmark_{table}_{year}.pdf")
    fig.savefig(output_filename, bbox_inches="tight")

    plt.close(fig)


def plot_benchmark(
    table: str,
    benchmarks_raw: pd.DataFrame,
    output_dir: str,
    options: dict,
    colors: dict,
    model_col: str = "Open-TYNDP",
    rfc_col: str = "TYNDP 2024",
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
    options : dict
        Full benchmarking configuration containing table units and conversions.
    colors : dict,
        Dictionary of colors to be used for each technology.
    model_col : str, default "Open-TYNDP"
        Column name for model values.
    rfc_col : str, default "TYNDP 2024"
        Column name for reference values.
    """

    # Parameters
    opt = options["tables"][table]
    table_type = opt["table_type"]
    source_unit = opt["unit"]
    unit_conversion = options["processing"]["unit_conversion"]
    unit_conversion = {
        k: 1 / v for k, v in unit_conversion.items()
    }  # Inverse conversion factor to revert unit

    # Filter data and Convert back to source unit
    logger.info(f"Making benchmark for {table} using {rfc_col} and {model_col}")
    benchmarks_raw = benchmarks_raw.query("table==@table").copy()
    benchmarks = _convert_units(benchmarks_raw, source_unit, unit_conversion)

    available_columns = [
        c for c in benchmarks.columns if c not in ["value", "source", "unit"]
    ]
    bench_wide = benchmarks.pivot_table(
        index=available_columns, values="value", columns="source", dropna=False
    )

    for year in bench_wide.index.get_level_values("year").unique():
        bench_year = bench_wide.query("year==@year").copy()

        if table_type == "scenario_comparison":
            _plot_scenario_comparison(
                bench_year.reset_index(),
                table,
                year,
                output_dir,
                model_col,
                rfc_col,
                source_unit,
            )
        elif table_type == "time_series":
            bench_agg = match_temporal_resolution(
                bench_year, model_col, rfc_col
            ).reset_index()
            _plot_time_series(
                bench_agg,
                table,
                year,
                output_dir,
                model_col,
                rfc_col,
                source_unit,
                colors,
            )


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
    scenario = "TYNDP " + snakemake.params["scenario"]  # noqa: F841
    benchmarks_fn = snakemake.input.benchmarks
    results_fn = snakemake.input.results
    output_dir = Path(snakemake.output[0])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    benchmarks_raw = load_data(benchmarks_fn, results_fn, scenario)

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
        options=options,
        colors=colors,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        results = list(tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs))

    logger.info("Benchmark plotting completed successfully")
