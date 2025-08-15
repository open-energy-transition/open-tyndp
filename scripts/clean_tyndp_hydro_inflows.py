# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available hydro inflow data from TYNDP data bundle for a given

* climate year,
* planning horizon,
* hydro technology.

Input data for TYNDP 2024 comes from PEMMDB v2.4.

Outputs
-------
Cleaned csv file with hourly hydro inflow time series in MW per region.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def read_hydro_inflows_file(
    node: str,
    hydro_inflows_dir: str,
    cyear: str,
    pyear: str,
    hydro_tech: str,
    sns: pd.DatetimeIndex,
    date_index: dict,
) -> pd.Series:
    fn = Path(
        hydro_inflows_dir,
        pyear,
        f"PEMMDB_{node.replace('GB', 'UK')}_Hydro_Inflows_{pyear}.xlsx",
    )

    if not os.path.isfile(fn):
        return None

    inflow_tech = pd.read_excel(
        fn,
        skiprows=1,
        usecols=lambda name: name == "Day"
        or name == "Week"
        or name == "ShortName"
        or name == "Variable"
        or name == int(cyear),
        sheet_name=f"{hydro_tech} - Year Dependent",
    )

    # infer resolution of data for each technology
    tech_res = "w" if "Week" in inflow_tech.columns else "d"

    inflow_tech = (
        inflow_tech.query("ShortName == 'INFLOW'")
        .assign(datetime=date_index[tech_res])
        .set_index("datetime")
        .reindex(sns)  # filter for hourly subset of snapshots only
        .ffill()  # upsample to hourly data
        .assign(
            **{
                node: lambda df: np.where(  # calculate hourly inflow in MWh/h
                    # input value was either in GWh/week or in GWh/day
                    df.Variable.str.contains("week"),
                    df[int(cyear)] / (24 * 7 * 1e-3),
                    df[int(cyear)] / (24 * 1e-3),
                )
            }
        )[node]
    )

    return inflow_tech


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_hydro_inflows",
            clusters="all",
            planning_horizons=2030,
            tech="Run of River",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    date_index = {
        "w": pd.date_range(
            start=f"{cyear}-01-01",
            periods=53,  # 53 weeks
            freq="7D",
        ),
        "d": pd.date_range(
            start=f"{cyear}-01-01",
            periods=366,  # 366 days (incl. first day of next year)
            freq="D",
        ),
    }
    if int(cyear) < 1982 or int(cyear) > 2019:
        logger.warning(
            f"Snapshot year {cyear} doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear = str(snakemake.wildcards.planning_horizons)

    # Parameters
    onshore_buses = pd.read_csv(snakemake.input.busmap, index_col=0)
    nodes = onshore_buses.index
    hydro_inflows_dir = snakemake.input.hydro_inflows_dir
    hydro_tech = str(snakemake.wildcards.tech)

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
        cyear=cyear,
        pyear=pyear,
        hydro_tech=hydro_tech,
        sns=sns,
        date_index=date_index,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        inflows = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    inflows_df = (
        pd.concat(inflows, axis=1)
        .reindex(
            nodes,
            axis=1,
        )  # include missing node data with empty columns
        .fillna(0.0)  # fill missing data with zero values
    )

    inflows_df.to_csv(snakemake.output.hydro_inflows_tyndp)
