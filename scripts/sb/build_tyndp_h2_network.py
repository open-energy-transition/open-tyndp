# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script loads and cleans the TYNDP 2026 H2 reference grid for a given wildcard
planning horizon. The reference grid gives one sheet per planning horizon
(``Year_{planning_horizons}``); the border topology is identical across all
available planning horizons, only the NTC capacities differ.

Note this replaces the TYNDP 2024 workflow, which combined a single base-year
reference grid with separate "investment candidate" projects
(`build_tyndp_transmission_projects.py`, sourced from the now-discontinued
Investment Datasets) to approximate later planning horizons. TYNDP 2026
instead gives the final grid directly for each planning horizon, so no
project layering is needed anymore.
"""

import logging

import pandas as pd

from scripts._helpers import (
    configure_logging,
    extract_grid_data_tyndp,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_h2_grid(fn_grid: str, pyear: int) -> pd.DataFrame:
    """
    Load and clean the TYNDP 2026 H2 reference grid for a given planning horizon.

    Parameters
    ----------
    fn_grid : str
        Path to the TYNDP 2026 H2 reference grid Excel file
        ("ReferenceGrid_Hydrogen.xlsx").
    pyear : int
        Planning horizon used to select the corresponding sheet.

    Returns
    -------
    pd.DataFrame
        Cleaned TYNDP H2 reference grid.
    """
    h2_grid_raw = pd.read_excel(fn_grid, sheet_name=f"Year_{pyear}")
    h2_grid_raw["Border"] = h2_grid_raw["Border"].replace("UK", "GB", regex=True)
    h2_grid = extract_grid_data_tyndp(
        links=h2_grid_raw, idx_prefix="H2 pipeline", idx_connector="->"
    )

    return h2_grid


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_network",
            planning_horizons=2030,
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)

    # Load and clean H2 reference grid
    h2_grid = load_h2_grid(fn_grid=snakemake.input.h2_reference_grid, pyear=pyear)

    # Save cleaned H2 grid
    h2_grid.to_csv(snakemake.output.h2_grid_prepped)
