# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PECD capacity factor generation time series based on PECD weather data.
The script is executed for a given technology, and planning horizon. Technologies can be one of:

   * CSP_noStorage,
   * CSP_withStorage,
   * LFSolarPV,
   * LFSolarPVRooftop,
   * LFSolarPVUtility,
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


def read_pecd_file(node: list, fn_pecd: str, cyear: str, pyear: str, technology: str):
    fn = Path(fn_pecd, pyear, f"PECD_{technology}_{pyear}_{node}_edition 2023.2.csv")

    if not os.path.isfile(fn):
        return None

    pecd_bus = pd.read_csv(
        fn,
        skiprows=10,
        usecols=lambda name: name == "Date" or name == "Hour" or name == str(cyear),
    )
    datetime_str = f"{cyear}." + pecd_bus["Date"].str.cat(
        (pecd_bus["Hour"] - 1).astype(str), sep=" "
    )
    cf_pecd = (
        pecd_bus.set_index(pd.to_datetime(datetime_str, format="%Y.%d.%m. %H"))
        .drop(columns=["Date", "Hour"])
        .rename(columns={str(cyear): node})
    )

    return cf_pecd


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_pecd_data",
            clusters="all",
            technology="offwind-ac",
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Climate year from snapshots
    cyear = get_snapshots(snakemake.params.snapshots)[0].year
    if int(cyear) < 1982 or int(cyear) > 2019:
        # TODO: Note that because of this fallback, the snapshots of the profiles will not always match with the model snapshots
        logger.warning(
            "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear = str(snakemake.wildcards.planning_horizons)

    # Technology as in PECD terminology
    # TODO: find solution for solar profiles being differentiated between Utility and Rooftop for Italy
    pecd_tech_dict = {
        "offwind-ac": "Wind_Offshore",
        "offwind-dc": "Wind_Offshore",
        "offwind-float": "Wind_Offshore",
        "offwind": "Wind_Offshore",
        "onwind": "Wind_Onshore",
        "solar": "LFSolarPV",
        "solar-hsat": "LFSolarPV",
    }
    pecd_tech = pecd_tech_dict[snakemake.wildcards.technology]

    offshore_buses = pd.read_excel(snakemake.input.offshore_buses, index_col=0)
    onshore_buses = pd.read_csv(snakemake.input.onshore_buses, index_col=0)

    nodes = (
        offshore_buses.index if pecd_tech == "Wind_Offshore" else onshore_buses.index
    )
    fn_pecd = snakemake.input.fn_pecd

    # Load and prep electricity demand
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PECD capacity factor data",
    }

    func = partial(
        read_pecd_file,
        fn_pecd=fn_pecd,
        cyear=cyear,
        pyear=pyear,
        technology=pecd_tech,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        demand = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    pecd_df = (
        pd.concat(demand, axis=1)
        .reindex(nodes, axis=1)
        .fillna(0.0)  # include missing node data with empty columns
        .rename(
            columns=lambda x: x.replace("UK", "GB")
        )  # replace UK with GB for naming convention
    )

    pecd_df.to_csv(snakemake.output.pecd_data_clean)
