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
):
    fn = Path(hydro_inflows_dir, pyear, f"PEMMDB_{node}_Hydro_Inflows_{pyear}.xlsx")

    if not os.path.isfile(fn):
        return None

    tech_res = {
        "Run of River": "d",
        "Pondage": "d",
        "Reservoir": "w",
        "PS Open": "w",
        "PS Closed": "w",
    }

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

    sns_year = sns[0].year
    date_index = {
        "w": pd.date_range(
            start=f"{sns_year}-01-01",
            periods=53,  # 53 weeks
            freq="7D",
        ),
        "d": pd.date_range(
            start=f"{sns_year}-01-01",
            periods=366,  # 366 days (incl. first day of next year)
            freq="D",
        ),
    }

    inflow_tech = (
        inflow_tech.query("ShortName == 'INFLOW'")
        .assign(
            datetime=date_index[tech_res[hydro_tech]],
            p_nom=lambda df: np.where(  # calculate hourly inflow in GWh/h
                df.Variable.str.contains("week"),
                df[int(cyear)].div(24 * 7),  # the value will be either in GWh/week
                df[int(cyear)].div(24),  # the value will be either in GWh/day
            ),
        )
        .set_index("datetime")
        .reindex(sns)  # filter for snapshots only
        .ffill()
        .div(1e3)  # convert from GW to MW
        .p_nom.rename(node)
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
    if int(cyear) < 1982 or int(cyear) > 2019:
        logger.warning(
            "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear = str(snakemake.wildcards.planning_horizons)
    hydro_tech = str(snakemake.wildcards.tech)

    onshore_buses = pd.read_csv(snakemake.input.onshore_buses, index_col=0)

    nodes = onshore_buses.index.str.replace("GB", "UK", regex=True)
    hydro_inflows_dir = snakemake.input.hydro_inflows_dir

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
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        inflows = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    inflows_df = (
        pd.concat(inflows, axis=1)
        .reindex(
            nodes, axis=1, fill_value=0.0
        )  # include missing node data with empty columns
        .rename(
            columns=lambda x: x.replace("UK", "GB")
        )  # replace UK with GB for naming convention
        .fillna(0.0)
    )

    inflows_df.to_csv(snakemake.output.hydro_inflows_tyndp)
