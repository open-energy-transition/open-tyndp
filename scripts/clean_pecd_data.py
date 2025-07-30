# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PECD capacity factor generation time series based on PECD weather data.
The script is executed for a given technology and planning horizon. Technologies can be one of:

   * CSP
   * LFSolarPVUtility,
   * LFSolarPVRooftop,
   * Wind_Offshore,
   * Wind_Onshore.

Outputs
-------
Cleaned csv file with capacity factor generation time series and regions as columns.
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
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def read_pecd_file(
    node: str,
    dir_pecd: str,
    cyear: str,
    pyear: str,
    technology: str,
    sns: pd.DatetimeIndex,
):
    fn = Path(dir_pecd, pyear, f"PECD_{technology}_{pyear}_{node}_edition 2023.2.csv")

    # PECD only differentiates between utility and rooftop PV for some nodes
    if not os.path.isfile(fn) and "LFSolarPV" in technology:
        fn = Path(str(fn).replace(technology, "LFSolarPV"))
    if not os.path.isfile(fn):
        logger.warning(f"Missing data for {technology} in {node} in {pyear}.")
        return None

    # Malta CSP data file has an extra header row that must be skipped
    if node == "MT00" and technology == "CSP_noStorage" and pyear == "2040":
        skiprows = 11
    else:
        skiprows = 10

    pecd_bus = pd.read_csv(
        fn,
        skiprows=skiprows,  # first rows contain only file metadata
        usecols=lambda name: name == "Date"
        or name == "Hour"
        or name == str(cyear)
        or name == str(float(cyear)),
    ).rename(columns={str(float(cyear)): str(cyear)})

    datetime_str = f"{cyear}." + pecd_bus["Date"].str.cat(
        (pecd_bus["Hour"] - 1).astype(str), sep=" "
    )
    cf_pecd = (
        pecd_bus.set_index(pd.to_datetime(datetime_str, format="%Y.%d.%m. %H"))
        .drop(columns=["Date", "Hour"])
        .rename(columns={str(cyear): node})
        .loc[sns]  # filter for snapshots only
    )

    return cf_pecd


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_pecd_data",
            clusters="all",
            technology="Wind_Offshore",
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    if int(cyear) < 1982 or int(cyear) > 2019:
        # TODO: Note that because of this fallback, the snapshots of the profiles will not always match with the model snapshots
        logger.warning(
            "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear = str(snakemake.wildcards.planning_horizons)

    # Technology as in PECD terminology
    pecd_tech = snakemake.wildcards.technology

    offshore_buses = pd.read_excel(snakemake.input.offshore_buses, index_col=0)
    onshore_buses = pd.read_csv(snakemake.input.onshore_buses, index_col=0)

    nodes = (
        offshore_buses.index
        if pecd_tech == "Wind_Offshore"
        else onshore_buses.index.str.replace("GB", "UK", regex=True)
    )
    dir_pecd = snakemake.input.dir_pecd

    # Load and prep pecd data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PECD capacity factor data",
    }

    func = partial(
        read_pecd_file,
        dir_pecd=dir_pecd,
        cyear=cyear,
        pyear=pyear,
        technology=pecd_tech,
        sns=sns,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        pecd = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    if all(data is None for data in pecd):
        raise ValueError(
            f"No PECD data found for {pecd_tech} in {pyear}. Please specify a technology covered within the TYNDP PECD data."
        )
    pecd_df = pd.concat(pecd, axis=1)
    fill_na = (
        pd.Series(0.0, index=pecd_df.index)
        if snakemake.params.fill_gaps_method == "zero"
        else pecd_df.agg(snakemake.params.fill_gaps_method, axis=1)
    )
    pecd_df = (
        pecd_df.reindex(nodes, axis=1)  # include missing node data with empty columns
        .rename(
            columns=lambda x: x.replace("UK", "GB")
        )  # replace UK with GB for naming convention
        .where(
            lambda df: df.notna(), fill_na, axis=0
        )  # fill missing node data with configured aggregation method
    )

    pecd_df.to_csv(snakemake.output.pecd_data_clean)
