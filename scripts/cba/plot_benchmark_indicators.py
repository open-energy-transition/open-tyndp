# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Create benchmark plots for CBA indicators.

This script reads a indicators CSV file (single project or collected) that
includes model rows and TYNDP benchmark rows, and generates plots comparing
model values (dot) to TYNDP min/mean/max bands.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

MODEL_SUBINDEX_PREFERENCE = {
    "B2_societal_cost_variation": ["central", "mean", ""],
}


def select_model_values(df: pd.DataFrame, indicator: str) -> pd.Series:
    model = df[(df["source"] == "model") & (df["indicator"] == indicator)].copy()
    if model.empty:
        return pd.Series(dtype=float)

    if indicator.startswith("B4"):
        candidates = ["mean", ""]
    else:
        candidates = MODEL_SUBINDEX_PREFERENCE.get(indicator, ["", "mean", "central"])

    for subindex in candidates:
        subset = model[model["subindex"] == subindex]
        if not subset.empty:
            return subset.groupby("project_id")["value"].mean()

    return model.groupby("project_id")["value"].mean()


def plot_benchmark_indicator(df: pd.DataFrame, indicator: str, output_dir: Path) -> None:
    benchmark = df[
        (df["source"] == "2024 tyndp")
        & (df["indicator"] == indicator)
        & (df["subindex"].isin(["min", "mean", "max"]))
    ].copy()

    model_values = select_model_values(df, indicator)
    if benchmark.empty or model_values.empty:
        logger.info("Skipping %s benchmark plot (missing data)", indicator)
        return

    benchmark = benchmark[benchmark["project_id"].isin(model_values.index)]
    if benchmark.empty:
        logger.info("Skipping %s benchmark plot (no matching projects)", indicator)
        return

    pivot = benchmark.pivot_table(
        index="project_id", columns="subindex", values="value", aggfunc="mean"
    ).reindex(model_values.index)
    for col in ["min", "mean", "max"]:
        if col not in pivot.columns:
            logger.info("Skipping %s benchmark plot (missing %s)", indicator, col)
            return
    pivot = pivot.dropna(subset=["min", "mean", "max"])

    if pivot.empty:
        logger.info("Skipping %s benchmark plot (no complete min/mean/max)", indicator)
        return

    projects = pivot.index.tolist()
    x = range(len(projects))
    mean = pivot["mean"].to_numpy()
    ymin = pivot["min"].to_numpy()
    ymax = pivot["max"].to_numpy()

    plt.figure(figsize=(10, 4))
    plt.errorbar(
        x,
        mean,
        yerr=[mean - ymin, ymax - mean],
        fmt="o",
        color="gray",
        ecolor="lightgray",
        capsize=3,
        label="2024 TYNDP mean +/- min/max",
    )
    plt.scatter(
        x,
        model_values.loc[projects].to_numpy(),
        color="tab:blue",
        zorder=3,
        label="Model",
    )

    units = (
        df.loc[(df["indicator"] == indicator) & (df["source"] == "model"), "units"]
        .dropna()
        .unique()
    )
    unit_label = units[0] if len(units) else ""
    plt.xticks(list(x), projects, rotation=45, ha="right")
    plt.ylabel(unit_label)
    plt.title(indicator)
    plt.tight_layout()

    output_path = output_dir / f"benchmark_{indicator}.png"
    plt.savefig(output_path, dpi=150)
    plt.close()


def create_plots(indicators_file, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(indicators_file)
    if df.empty:
        logger.warning("No indicators data to plot")
        return

    if "source" not in df.columns:
        df["source"] = "model"

    benchmark_indicators = [
        "B1_total_system_cost_change",
        "co2_variation",
        "B2_societal_cost_variation",
        "B3_res_generation_change_mwh",
        "B3_res_capacity_change_mw",
        "B4a_nox_t_per_year",
        "B4b_nh3_t_per_year",
        "B4c_sox_t_per_year",
        "B4d_pm25_t_per_year",
        "B4e_pm10_t_per_year",
        "B4f_nmvoc_t_per_year",
    ]

    for indicator in benchmark_indicators:
        plot_benchmark_indicator(df, indicator, output_dir)

    logger.info("Benchmark plots saved to %s", output_dir)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("plot_benchmark_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    create_plots(snakemake.input.indicators, snakemake.output.plot_dir)
