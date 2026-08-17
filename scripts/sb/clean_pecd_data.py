# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PECD capacity factor generation time series based on PECD weather data.
The script is executed for a given technology and planning horizon. Technologies can be one of:

   * Solar PV Utility,
   * Solar PV Rooftop,
   * Wind_Offshore,
   * Wind_Onshore,
   * CSP_noStorage_0h_dispatched,

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
    weather_year: int,
    pyear: int,
    technology: str,
    sns: pd.DatetimeIndex,
):

    if "Solar" in technology:
        fn = Path(dir_pecd, str(pyear), f"{technology} {node.replace('UK', 'GB')}.csv")
    else:
        fn = Path(
            dir_pecd,
            str(pyear),
            f"{node.replace('UK', 'GB')}_CapacityFactors_{technology}_{pyear}.csv",
        )

    if not os.path.isfile(fn):
        logger.warning(f"Missing data for {technology} in {node} in {pyear}.")
        return None

    pecd_bus = pd.read_csv(fn)
    year = sns[0].year
    datetime_idx = pd.to_datetime(
        f"{year}."
        + pecd_bus["Date"].str.cat((pecd_bus["Hour"] - 1).astype(str), sep=" "),
        format="%Y.%d.%m. %H",
    )

    cf_pecd = (
        pecd_bus.set_index(datetime_idx)
        .drop(columns=["Date", "Hour"])
        .loc[sns, [weather_year]]  # filter for snapshots and weather scenario only
        .rename(columns={weather_year: node})
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
    weather_year = f"WS{snakemake.params.weather_year:03d}"

    # Planning year (falls back to latest available pyear if not in list of available years)
    pyear = safe_pyear(
        snakemake.wildcards.planning_horizons,
        available_years=snakemake.params.available_years,
        source="PECD",
    )

    # Technology as in PECD terminology
    pecd_tech = snakemake.wildcards.technology

    df_nodes = pd.read_excel(snakemake.input.nodes, sheet_name=None)
    onshore_buses = df_nodes["Electricity"]["NODE"].tolist()
    offshore_buses = onshore_buses + df_nodes["Electricity_Offshore"]["NODE"].tolist()
    busmap = pd.read_csv(snakemake.input.busmap).name.tolist()
    nodes = offshore_buses if pecd_tech == "Wind_Offshore" else onshore_buses
    nodes = [x.replace("UK", "GB") for x in nodes]

    # Nodes present in the TYNDP 2026 node list but absent from the rest of the workflow,
    # which still relies on the TYNDP 2024 node set. Dropped to keep PECD consistent with it.
    # TODO Remove once the TYNDP 2026 nodes are integrated
    # excluded nodes - "MD00", "NOS1", "NOS2", "NOS3", "TR00", "UA00", "PL00E", "PL00I"
    nodes = [x for x in nodes if x in busmap]
    dir_pecd = snakemake.input.pecd_input

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
        weather_year=weather_year,
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
