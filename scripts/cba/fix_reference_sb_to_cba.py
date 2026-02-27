# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Build a ``Link`` dataframe of capacity corrections to align SB transmission projects with the CBA reference grid. A positive ``p_nom`` indicates capacity to add (potentially as a new link), while a negative ``p_nom`` indicates capacity to reduce on an existing link.

This only works for 2040, as we are assuming the 2030 reference grid does not need corrections.
"""

import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd

from scripts._helpers import (
    configure_logging,
    extract_grid_data_tyndp,
    set_scenario_config,
)
from scripts.build_tyndp_network import add_links_missing_attributes
from scripts.sb.build_tyndp_transmission_projects import read_invest_file


def read_guidelines(path: Path) -> pd.DataFrame:
    query = "In_ref_grid_2030=='no' and In_ref_grid_2040=='yes'"
    df = (
        pd.read_csv(path)
        .query(query)
        .rename(
            columns={
                "A_B_MW": "Summary Direction 1",
                "B_A_MW": "Summary Direction 2",
            }
        )
        .replace(
            {
                "ITS1-GR00": "GR00-ITS1",
                "MT00-ITSI": "ITSI-MT00",
                "UKNI-UK00": "UK00-UKNI",
            }
        )
        .replace({"UK": "GB"}, regex=True)
        .set_index("Border")[["Summary Direction 1", "Summary Direction 2"]]
        .dropna(how="all")
        .astype({"Summary Direction 1": float, "Summary Direction 2": float})
        .groupby("Border")
        .sum()
    )
    return df[~df.index.str.contains("internal|Offshore|OBZ")]


def build_links(projects_fix: pd.DataFrame, buses_fn: str) -> pd.DataFrame:
    links_fix = extract_grid_data_tyndp(
        projects_fix.reset_index(),
        expand_from_index=True,
        idx_connector="-",
        idx_suffix="-DC",
        idx_separator="",
    )
    buses = gpd.read_file(buses_fn).set_index("bus_id")
    links_fix = add_links_missing_attributes(links_fix, buses)

    links_fix["carrier"] = "DC"
    links_fix["length"] /= 1e3
    links_fix[["underground", "under_construction"]] = (
        links_fix[["underground", "under_construction"]] == "t"
    )

    return links_fix


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("fix_reference_sb_to_cba", planning_horizons="2040")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    invest_path = Path(snakemake.input.invest_grid)
    guidelines_path = Path(snakemake.input.guidelines)
    build_years = snakemake.params.build_years

    if not build_years:
        pd.DataFrame().to_csv(snakemake.output.corrections, index=True)
        sys.exit(0)

    # Read references
    df_invest = (
        read_invest_file(invest_path, years=build_years)
        .assign(Border=lambda df: df.bus0 + "-" + df.bus1)
        .drop(columns=["bus0", "bus1"])
        .set_index("Border")
    )
    df_gl_projects = read_guidelines(guidelines_path)

    # Build fixes
    projects_fix = df_gl_projects.subtract(df_invest, fill_value=0)
    projects_fix = projects_fix[projects_fix.sum(axis=1) != 0]

    # Add EU-GB border specific treatment (see Implementation Guidelines for TYNDP 2024, Appendix B.2, p119)
    current = 6850  # MW (2030)
    desired = 8725  # MW (2035)
    projects_fix.loc["FR00-GB00"] = [
        desired - current,
        desired - current,
    ]

    # Transform fixes to links definitions
    links_fix = build_links(projects_fix, snakemake.input.buses)
    links_fix.to_csv(snakemake.output.corrections, quotechar="'")
