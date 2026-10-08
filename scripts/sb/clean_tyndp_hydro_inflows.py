# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available hydro inflow data from TYNDP data bundle for a given

* climate year,
* planning horizon,
* hydro technology.

Input data for TYNDP 2024 comes from PEMMDB v2.5.

Outputs
-------
Cleaned csv file with hourly hydro inflow time series in MW per region.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    get_wscenario,
    safe_planning_horizon,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

HYDRO_TECH_CODES = {
    "Run of River": "HRR",
    "Pondage": "HPI",
    "Reservoir": "HRI",
    "PS Open": "HOL",
    "PS Closed": "HCL",
}


def read_hydro_inflows_file(
    node: str,
    hydro_inflows_dir: str,
    wscenario: int,
    planning_horizon: int,
    hydro_tech: str,
    sns: pd.DatetimeIndex,
    date_index: dict,
) -> pd.Series:
    fn = Path(
        hydro_inflows_dir,
        "Hydro Inflows",
        str(planning_horizon),
        f"{node.replace('GB', 'UK')}_Hydro_Inflows_{HYDRO_TECH_CODES[hydro_tech]}_{planning_horizon}.csv",
    )

    if not os.path.isfile(fn):
        return None

    inflow_tech = pd.read_csv(fn, index_col=0)

    # infer resolution of data for each technology
    tech_res = "w" if inflow_tech.index.name == "WEEK" else "d"

    inflow_tech = (
        inflow_tech[f"WS{wscenario:03d}"]
        .set_axis(date_index[tech_res])
        .reindex(sns)  # filter for hourly subset of snapshots only
        .ffill()  # upsample to hourly data
        .div(  # calculate hourly inflow in MW
            # input value was either in MWh/week or in MWh/day
            24 * 7 if tech_res == "w" else 24
        )
        .rename(node)
    )

    return inflow_tech


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_hydro_inflows",
            clusters="all",
            planning_horizons=2030,
            tech="Run_of_River",
            run="NT",
            configfiles="config/config.tyndp.yaml",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    year = sns[0].year
    date_index = {
        "w": pd.date_range(
            start=f"{year}-01-01",
            periods=53,  # 53 weeks
            freq="7D",
        ),
        "d": get_snapshots(
            {
                "start": f"{year}-01-01",
                "end": f"{year + 1}-01-01",
                "inclusive": "left",
            },
            drop_leap_day=True,
            freq="D",
        ),
    }

    # Planning year
    planning_horizon = safe_planning_horizon(
        snakemake.wildcards.planning_horizons,
        available_years=snakemake.params.available_years,
        source="Hydro inflows",
    )

    # Weather scenario
    wscenario = get_wscenario(snakemake.params.wscenarios, planning_horizon)

    # Parameters
    onshore_buses = pd.read_csv(snakemake.input.busmap, index_col=0)
    nodes = onshore_buses.index
    hydro_inflows_dir = snakemake.input.hydro_inflows_dir
    hydro_tech = str(snakemake.wildcards.tech).replace("_", " ")

    # Load and prep inflow data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading TYNDP hydro inflows data",
    }

    func = partial(
        read_hydro_inflows_file,
        hydro_inflows_dir=hydro_inflows_dir,
        wscenario=wscenario,
        planning_horizon=planning_horizon,
        hydro_tech=hydro_tech,
        sns=sns,
        date_index=date_index,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        inflows = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    inflows_df = (
        # start with empty dataframe so workflow will not crash if no inflows are found (e.g. for PS Closed)
        pd.concat([pd.DataFrame(index=sns), *inflows], axis=1)
        .reindex(
            nodes,
            axis=1,
        )  # include missing node data with empty columns
        .fillna(0.0)  # fill missing data with zero values
    )

    inflows_df.to_csv(snakemake.output.hydro_inflows_tyndp)
