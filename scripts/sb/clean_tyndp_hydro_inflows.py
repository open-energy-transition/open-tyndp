# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available hydro inflow data from TYNDP data bundle for a given

* weather scenario,
* planning horizon,
* hydro technology.

Input data for TYNDP 2026 comes from PECD v4.2.

Outputs
-------
Cleaned csv file with hourly hydro inflow time series in MW per region.
If no inflow data exists for a technology, the output contains zero inflows for all nodes.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config
from scripts.sb._helpers import (
    check_snapshot_year,
    get_wscenario,
    safe_planning_horizon,
)

logger = logging.getLogger(__name__)


def read_hydro_inflows_file(
    node: str,
    hydro_inflows_dir: str,
    wscenario: int,
    planning_horizon: int,
    hydro_tech: str,
    sns: pd.DatetimeIndex,
    date_index: dict,
) -> pd.Series | None:
    fn = Path(
        hydro_inflows_dir,
        "Hydro Inflows",
        str(planning_horizon),
        f"{node.replace('GB', 'UK')}_Hydro_Inflows_{hydro_tech}_{planning_horizon}.csv",
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
            tech="HRR",
            configfiles="config/test/config.tyndp.yaml",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    year = sns[0].year
    check_snapshot_year(year, snakemake.params.drop_leap_day, data_type="hydro inflow")
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
        source="PECD hydro inflow",
    )

    # Weather scenario
    wscenario = get_wscenario(snakemake.params.wscenarios, planning_horizon)

    # Parameters
    onshore_buses = pd.read_csv(snakemake.input.busmap, index_col=0)
    nodes = onshore_buses.index
    hydro_inflows_dir = snakemake.input.hydro_inflows_dir
    hydro_tech = snakemake.wildcards.tech

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
        pd.concat([pd.DataFrame(index=sns), *inflows], axis=1)
        .reindex(
            nodes,
            axis=1,
        )  # include missing node data with empty columns
        .fillna(0.0)  # fill missing data with zero values
    )

    inflows_df.to_csv(snakemake.output.hydro_inflows_tyndp)
