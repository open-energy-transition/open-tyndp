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

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def plot_b1(df, output_dir):
    """
    Plot the B1 indicator (Total System Cost Change) for all projects.

    The B1 indicator represents the change in total system cost:
    - PINT: positive B1 means beneficial (project reduces costs)
    - TOOT: positive B1 means beneficial (removing project increases costs)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with indicator data
    output_dir : Path
        Path to output directory for plots

    Notes
    -----
    TODO: Implement B1 visualization
    - Bar chart showing B1 values per project
    - Color coding for beneficial (green) vs non-beneficial (red)
    - Separate plots for CAPEX and OPEX changes
    - Comparison between PINT and TOOT methods if both present
    """
    pass


def create_plots(indicators_file, output_dir):
    """
    Create all indicator plots.

    Parameters
    ----------
    indicators_file : str or Path
        Path to collected indicators CSV file
    output_dir : str or Path
        Path to output directory for plots
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read indicators
    df = pd.read_csv(indicators_file)

    if df.empty:
        logger.warning("No indicators data to plot")
        return

    logger.info(f"Creating plots for {len(df)} projects")

    # Generate B1 plots
    plot_b1(df, output_dir)

    logger.info(f"Plots saved to {output_dir}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("plot_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Create plots from collected indicators
    create_plots(snakemake.input.indicators, snakemake.output.plot_dir)
