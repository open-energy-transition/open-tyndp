# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP 2026 Shared RES (SRES) electrolyser profiles for a given planning
horizon.

Inputs
------

- `SRES/{planning_horizons}/SRES_<node>h2_<year>_SumGenH2.csv`: hourly
  electricity input to the SRES electrolyser in MW_e, one column per weather
  scenario (``WS061`` to ``WS090``).

Outputs
-------

- `resources/sres_profiles_{planning_horizons}.csv`: hourly SRES electrolyser
  electricity input in MW_e for the selected weather scenario, indexed by
  ``snapshot``, with one column per e-market node (e.g. ``DEh2`` becomes
  ``DE00``).
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import (
    check_snapshot_year,
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_sres_profiles", planning_horizons=2040)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    year = sns[0].year
    check_snapshot_year(year, snakemake.params.drop_leap_day, data_type="SRES")
    wscenario = f"WS{snakemake.params.wscenario:03d}"
    planning_horizon = snakemake.wildcards.planning_horizons

    profiles = {}
    for fn in Path(snakemake.input.sres_dir, planning_horizon).glob("*_SumGenH2.csv"):
        sres = pd.read_csv(fn, skiprows=10, usecols=["Date", "Hour", wscenario])
        datetime_idx = pd.to_datetime(
            f"{year}."
            + sres["Date"].str.cat((sres["Hour"] - 1).astype(str), sep=" "),
            format="%Y.%d.%m. %H",
        )
        node = fn.name.split("_")[1].replace("h2", "00")
        profiles[node] = sres.set_index(datetime_idx).loc[sns, wscenario]

    pd.DataFrame(profiles).rename_axis("snapshot").to_csv(
        snakemake.output.sres_profiles
    )