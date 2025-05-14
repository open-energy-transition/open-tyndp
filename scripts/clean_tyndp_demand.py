# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean TYNDP Scenario Building demand data to be used in the PyPSA-Eur workflow. The `snapshot` year is used as climatic year (`cyear`). For DE and GA, it must be one of the following years: 1995, 2008 or 2009. For NT, it must be between 1982 and 2019. If the `snapshot` is not one of these years, then the demand is set to 2009 electricity demand (2009 being considered as the most representative of the three years).

Depending on the scenario, different planning years (`pyear`) are available. DE and GA are defined for 2030, 2040 and 2050. NT scenario is only defined for 2030 and 2040. All the planning years are read at once.
"""

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import pandas as pd
from _helpers import configure_logging, get_snapshots, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)


def load_elec_demand(fn: str, scenario: str, pyear: int, cyear: int):
    """
    Load electricity demand files into dictionary of dataframes. Filter for specific climatic year and format data.
    """
    if scenario == "NT":
        if pyear == 2050:
            logger.warning(
                "2050 electricity demand data are not defined for NT in 2024 TYNDP cycle. Falling back to 2040."
            )
            pyear = 2040
        demand_fn = Path(
            fn,
            scenario,
            "Electricity demand profiles",
            f"{pyear}_National Trends.xlsx",
        )

        if int(cyear) < 1982 or int(cyear) > 2019:
            logger.warning(
                "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
            )
            cyear = 2009
        data = pd.read_excel(
            demand_fn,
            skiprows=7,
            index_col=0,
            usecols=lambda name: name == "Date" or name == int(cyear),
            sheet_name=None,
        )

    elif scenario in ["DE", "GA"]:
        demand_fn = Path(
            fn,
            scenario,
            str(pyear),
            f"ELECTRICITY_MARKET {scenario} {pyear}.xlsx",
        )

        if int(cyear) not in [1995, 2008, 2009]:
            logger.warning(
                "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
            )
            cyear = 2009
        data = pd.read_excel(
            demand_fn,
            skiprows=11,
            index_col=0,
            usecols=lambda name: name == "Date" or name == int(cyear),
            sheet_name=None,
        )

        # Fix inconsistencies in input data
        if pyear in [2040, 2050]:
            data["PL00"].index = data["AL00"].index
            data["UK00"] = pd.read_excel(
                demand_fn,
                skiprows=11,
                index_col=0,
                usecols=lambda name: name == "Date" or name == int(cyear) - 1,
                sheet_name="UK00",
            )

    demand = pd.concat(data, axis=1).droplevel(1, axis=1)

    # need to reindex load time series to target year
    demand.index = demand.index.map(lambda t: t.replace(year=pyear))

    # rename UK in GB
    demand.columns = demand.columns.str.replace("UK", "GB")

    return demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_demand", configfiles="config/test/config.tyndp.yaml"
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params.get("scenario", "DE")
    cyear = get_snapshots(snakemake.params.snapshots)[0].year
    planning_horizons = snakemake.params["planning_horizons"]

    # Load and prep electricity demand
    tqdm_kwargs = {
        "ascii": False,
        "unit": " pyear",
        "total": len(planning_horizons),
        "desc": "Loading TYNDP demand data",
    }

    func = partial(
        load_elec_demand,
        snakemake.input.electricity_demand,
        scenario,
        cyear=cyear,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        demand = list(tqdm(pool.imap(func, planning_horizons), **tqdm_kwargs))

    demand = pd.concat(demand)

    # Save prepped electricity demand
    demand.to_csv(snakemake.output.electricity_demand_prepped)
