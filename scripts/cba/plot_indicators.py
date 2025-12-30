# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Create plots for CBA indicators.

This script reads the collected indicators CSV file and generates various
plots to visualize the cost-benefit analysis results, including the B1
indicator (Total System Cost difference).
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

COLORS = {
    "capex": "#2980b9",
    "opex": "#e67e22",
    "beneficial": "#27ae60",
    "not_beneficial": "#c0392b",
}


def load_and_merge_data(indicators_path, projects_path):
    """Load indicators and merge with project metadata."""
    indicators = pd.read_csv(indicators_path)
    projects = pd.read_csv(projects_path)

    border_counts = projects.groupby("project_id").size().rename("border_count")
    projects_agg = (
        projects.groupby("project_id")
        .agg(
            {
                "project_name": "first",
                "is_crossborder": "first",
                "p_nom 0->1": "sum",
                "p_nom 1->0": "sum",
            }
        )
        .reset_index()
    )
    projects_agg = projects_agg.merge(border_counts, on="project_id")
    projects_agg["total_capacity_MW"] = (
        projects_agg["p_nom 0->1"] + projects_agg["p_nom 1->0"]
    )

    merged = indicators.merge(projects_agg, on="project_id", how="left")
    merged["B1_billion_EUR"] = merged["B1_total_system_cost_change"] / 1e9
    merged["capex_change_billion"] = merged["capex_change"] / 1e9
    merged["opex_change_billion"] = merged["opex_change"] / 1e9

    return merged, projects["project_id"].nunique()


def plot_b1_top_projects(df, output_dir, method, n_top=20):
    """Waterfall chart showing B1 CAPEX/OPEX breakdown for top N projects."""
    df_top = df.nlargest(n_top, "B1_billion_EUR", keep="first")
    df_sorted = df_top.sort_values("B1_billion_EUR", ascending=True).reset_index(
        drop=True
    )

    fig, ax = plt.subplots(figsize=(16, 10))
    y_pos = np.arange(len(df_sorted))
    bar_height = 0.35

    ax.barh(
        y_pos - bar_height / 2,
        df_sorted["capex_change_billion"],
        height=bar_height,
        color=COLORS["capex"],
        alpha=0.9,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.barh(
        y_pos + bar_height / 2,
        df_sorted["opex_change_billion"],
        height=bar_height,
        color=COLORS["opex"],
        alpha=0.9,
        edgecolor="white",
        linewidth=0.5,
    )

    for i, (_, row) in enumerate(df_sorted.iterrows()):
        color = (
            COLORS["beneficial"]
            if row["B1_billion_EUR"] >= 0
            else COLORS["not_beneficial"]
        )
        ax.plot(
            row["B1_billion_EUR"],
            i,
            "D",
            color=color,
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            zorder=5,
        )

    labels = []
    for _, row in df_sorted.iterrows():
        multi = "*" if row["border_count"] > 1 else ""
        labels.append(
            f"{row['project_id']}: {row['project_name']}{multi} "
            f"[{row['p_nom 0->1']:.0f}/{row['p_nom 1->0']:.0f} MW]"
        )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)

    ax.axvline(x=0, color="black", linewidth=1.2, zorder=3)
    ax.set_xlabel("Cost Change (Billion EUR)", fontsize=12, fontweight="bold")
    ax.set_title(
        f"{method.upper()} B1 Indicator: Top {n_top} Projects by Benefit\n(Positive = Beneficial)",
        fontsize=14,
        fontweight="bold",
    )

    max_val = df_sorted["B1_billion_EUR"].abs().max()
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        total = row["B1_billion_EUR"]
        color = COLORS["beneficial"] if total >= 0 else COLORS["not_beneficial"]
        ax.text(
            max_val * 1.1,
            i,
            f"{total:+.2f}B",
            va="center",
            ha="left",
            fontsize=9,
            color=color,
            fontweight="medium",
        )

    ax.grid(axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)

    legend_elements = [
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=COLORS["capex"], alpha=0.9, label="CAPEX Change"
        ),
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=COLORS["opex"], alpha=0.9, label="OPEX Change"
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=COLORS["beneficial"],
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            label="B1 (Beneficial)",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=COLORS["not_beneficial"],
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            label="B1 (Not Beneficial)",
        ),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=10)
    ax.set_xlim(ax.get_xlim()[0], max_val * 1.4)

    plt.tight_layout()
    fig.savefig(
        output_dir / f"b1_top_{n_top}_projects.png", dpi=150, bbox_inches="tight"
    )
    plt.close(fig)
    logger.info(f"Saved B1 top {n_top} projects plot")


def plot_b1_summary(df, output_dir, method, total_projects):
    """Summary plot with B1 histogram."""
    beneficial = df[df["is_beneficial"]]
    not_beneficial = df[~df["is_beneficial"]]

    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(
        df["B1_billion_EUR"].min() - 0.1, df["B1_billion_EUR"].max() + 0.1, 30
    )
    ax.hist(
        beneficial["B1_billion_EUR"],
        bins=bins,
        color=COLORS["beneficial"],
        alpha=0.7,
        label=f"Beneficial ({len(beneficial)})",
        edgecolor="white",
    )
    ax.hist(
        not_beneficial["B1_billion_EUR"],
        bins=bins,
        color=COLORS["not_beneficial"],
        alpha=0.7,
        label=f"Not Beneficial ({len(not_beneficial)})",
        edgecolor="white",
    )

    ax.axvline(x=0, color="black", linewidth=1.5, linestyle="--")
    ax.set_xlabel("B1 Indicator (Billion EUR)", fontsize=11)
    ax.set_ylabel("Number of Projects", fontsize=11)
    ax.set_title(
        f"{method.upper()} B1: Distribution of Project Benefits",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    # Force integer y-axis ticks (number of projects must be integers)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.tight_layout()
    fig.savefig(output_dir / "b1_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved B1 summary plot")


def plot_b1_capex_vs_opex(df, output_dir, method):
    """Scatter plot of B1 CAPEX vs OPEX changes."""
    fig, ax = plt.subplots(figsize=(10, 8))

    beneficial = df[df["is_beneficial"]]
    not_beneficial = df[~df["is_beneficial"]]

    ax.scatter(
        beneficial["capex_change_billion"],
        beneficial["opex_change_billion"],
        c=COLORS["beneficial"],
        s=80,
        alpha=0.7,
        label=f"Beneficial ({len(beneficial)})",
        edgecolors="white",
        linewidth=0.5,
    )
    ax.scatter(
        not_beneficial["capex_change_billion"],
        not_beneficial["opex_change_billion"],
        c=COLORS["not_beneficial"],
        s=80,
        alpha=0.7,
        label=f"Not Beneficial ({len(not_beneficial)})",
        edgecolors="white",
        linewidth=0.5,
    )

    ax.axhline(y=0, color="gray", linewidth=0.8, linestyle="--")
    ax.axvline(x=0, color="gray", linewidth=0.8, linestyle="--")

    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, [-x for x in lims], "k:", alpha=0.5, label="CAPEX + OPEX = 0")

    ax.set_xlabel("CAPEX Change (Billion EUR)", fontsize=11)
    ax.set_ylabel("OPEX Change (Billion EUR)", fontsize=11)
    ax.set_title(
        f"{method.upper()} B1: CAPEX vs OPEX Changes\n"
        "(Projects above diagonal are beneficial)",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / "b1_capex_vs_opex.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved B1 CAPEX vs OPEX plot")


# Future indicator plot functions:
# def plot_b2_...
# def plot_b3_...


def create_plots(indicators_path, projects_path, output_dir):
    """Create all CBA indicator plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df, total_projects = load_and_merge_data(indicators_path, projects_path)
    if df.empty:
        logger.warning("No indicators data to plot")
        return

    method = df["cba_method"].iloc[0].lower()
    logger.info(f"Creating {method.upper()} plots for {len(df)} projects")

    # B1 indicator plots
    plot_b1_top_projects(df, output_dir, method, n_top=min(20, len(df)))
    plot_b1_summary(df, output_dir, method, total_projects)
    plot_b1_capex_vs_opex(df, output_dir, method)

    # Future indicators (add when implemented in make_indicators.py):
    # plot_b2_...(df, output_dir, method)
    # plot_b3_...(df, output_dir, method)

    logger.info(f"All plots saved to {output_dir}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_indicators",
            cba_method="toot",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    create_plots(
        snakemake.input.indicators,
        snakemake.input.transmission_projects,
        snakemake.output.plot_dir,
    )
