# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Build a table with values to use to correct the CBA reference grid.

This only works for 2040, as we are assuming the 2030 reference grid does not need corrections.
"""

import sys
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
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


def build_table(sb_invest: pd.DataFrame, cba_guidelines: pd.DataFrame) -> pd.DataFrame:
    sb = sb_invest.rename(lambda x: "SB - " + x, axis=1)
    cba = cba_guidelines.rename(lambda x: "CBA Projects - " + x, axis=1)
    merged = pd.concat([sb, cba], axis=1)
    correction = cba_guidelines.subtract(sb_invest, fill_value=0).rename(
        lambda x: "Correction - " + x, axis=1
    )
    df = pd.concat([merged, correction], axis=1)
    df = df[df.filter(like="Correction", axis=1).sum(axis=1) != 0]
    return df


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
        pd.DataFrame().to_csv(snakemake.output.corrections_csv, index=True)
        sys.exit(0)

    df_invest = (
        read_invest_file(invest_path, years=build_years)
        .assign(Border=lambda df: df.bus0 + "-" + df.bus1)
        .drop(columns=["bus0", "bus1"])
        .set_index("Border")
    )
    df_gl_projects = read_guidelines(guidelines_path)
    output = build_table(df_invest, df_gl_projects)

    # Add EU-GB border specific treatment (see Implementation Guidelines for TYNDP 2024, Appendix B.2, p119)
    current = 6850  # MW
    desired = 8725  # MW
    output.loc["FR00-GB00"] = [
        current,
        current,
        desired,
        desired,
        desired - current,
        desired - current,
    ]

    output.to_csv(snakemake.output.corrections_csv, index=True)
