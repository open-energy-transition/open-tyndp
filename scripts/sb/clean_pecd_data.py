# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PECD capacity factor generation time series based on PECD weather data.
The script is executed for a given technology and planning horizon. Technologies can be one of:

   * LFSolarPVUtility,
   * LFSolarPVRooftop,
   * Wind_Offshore,
   * Wind_Onshore,
   * CSP_noStorage,
   * CSP_withStorage_7h_dispatched,
   * CSP_withStorage_7h_preDispatch (note: includes cf > 1 for when thermal storage can be used).

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
    safe_pyear,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def read_pecd_file(
    node: str,
    dir_pecd: str,
    cyear: int,
    cyear_i: int,
    pyear: int,
    technology: str,
    sns: pd.DatetimeIndex,
):
    fn = Path(
        dir_pecd,
        str(pyear),
        f"PECD_{technology}_{pyear}_{node.replace('GB', 'UK')}_edition 2023.2.csv",
    )

    # PECD only differentiates between utility and rooftop PV for some nodes
    if not os.path.isfile(fn) and "LFSolarPV" in technology:
        fn = Path(str(fn).replace(technology, "LFSolarPV"))
    if not os.path.isfile(fn):
        logger.warning(f"Missing data for {technology} in {node} in {pyear}.")
        return None

    pecd_bus = pd.read_csv(fn)

    datetime_str = f"{cyear_i}." + pecd_bus["Date"].str.cat(
        (pecd_bus["Hour"] - 1).astype(str), sep=" "
    )
    cf_pecd = (
        pecd_bus.set_index(pd.to_datetime(datetime_str, format="%Y.%d.%m. %H"))
        .drop(columns=["Date", "Hour"])
        .loc[sns, [str(cyear)]]  # filter for snapshots and climate year only
        .rename(columns={str(cyear): node})
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

    # Parameters
    ############

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    # define climate year to use for the Datetime Index later on
    cyear_i = cyear
    prebuilt_years = snakemake.params.prebuilt_years

    if int(cyear) not in prebuilt_years:
        # TODO: Note that because of this fallback, the snapshots of the profiles will not always match with the model snapshots
        fallback_year = (
            2009 if 2009 in prebuilt_years else (prebuilt_years[-1])
        )  # use 2009 as default fallback if one of the filtered cyears
        logger.warning(
            f"Snapshot year doesn't match available TYNDP data. Falling back to {fallback_year}."
        )
        cyear = fallback_year

    # Planning year (falls back to latest available pyear if not in list of available years)
    pyear = safe_pyear(
        snakemake.wildcards.planning_horizons,
        available_years=snakemake.params.available_years,
        source="PECD",
    )

    # Technology as in PECD terminology
    pecd_tech = snakemake.wildcards.technology

    offshore_buses = pd.read_excel(snakemake.input.offshore_buses, index_col=0)
    onshore_buses = pd.read_csv(snakemake.input.onshore_buses, index_col=0)

    nodes = (
        offshore_buses.index.str.replace(
            "UK", "GB", regex=True
        )  # replace UK with GB for naming convention
        if pecd_tech == "Wind_Offshore"
        else onshore_buses.index
    )
    dir_pecd = snakemake.input.pecd_prebuilt

    # Load and prep pecd data
    #########################

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
        cyear_i=cyear_i,
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
        pecd_df.reindex(
            nodes, axis=1
        ).where(  # include missing node data with empty columns
            lambda df: df.notna(), fill_na, axis=0
        )  # fill missing node data with configured aggregation method
    )

    pecd_df.to_csv(snakemake.output.pecd_data_clean)
