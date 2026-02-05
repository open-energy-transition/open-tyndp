# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Build a table with values to use to correct the CBA reference grid.

This only works for 2040, as we are assuming the 2030 reference grid does not need corrections.
"""

from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config


def read_invest(path: Path, year: int, op: str) -> pd.DataFrame:
    return (
        pd.read_excel(path, sheet_name="Electricity")
        .query(f"BORDER.str.contains('Real') & YEAR {op} {year}")
        .rename(
            columns={
                "DIRECT CAPACITY INCREASE (MW)": "Summary Direction 1",
                "INDIRECT CAPACITY INCREASE (MW)": "Summary Direction 2",
                "FROM NODE": "bus0",
                "TO NODE": "bus1",
            }
        )
        .groupby(["bus0", "bus1"])
        .sum()[["Summary Direction 1", "Summary Direction 2"]]
        .reset_index()
        .assign(Border=lambda df: df.bus0 + "-" + df.bus1)
        .drop(columns=["bus0", "bus1"])
        .set_index("Border")
    )


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
    return pd.concat([merged, correction], axis=1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("fix_reference_sb_to_cba")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    invest_path = Path(snakemake.input.invest_grid)
    guidelines_path = Path(snakemake.input.guidelines)

    df_invest_strict = read_invest(invest_path, year=2035, op="==")
    df_gl_projects = read_guidelines(guidelines_path)
    output = build_table(df_invest_strict, df_gl_projects)
    output.to_csv(snakemake.output.corrections_csv, index=True)
