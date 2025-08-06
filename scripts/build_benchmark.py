# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script computes the benchmarks from the optimised network.
"""

import logging
import multiprocessing as mp
from functools import partial

import pandas as pd
import pypsa
from tqdm import tqdm

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def compute_benchmark(n: pypsa.Network, table: str) -> pd.DataFrame:
    """
    Compute benchmark metrics from optimised network.

    Parameters
    ----------
    n : pypsa.Network
        Optimised network.
    table : str
        Benchmark metric to compute.

    Returns
    -------
    pd.DataFrame
        Benchmark data in long format.
    """
    elec_bus_carrier = ["AC", "AC_OH", "low voltage"]
    supply_comps = ["Generator", "Link"]
    demand_comps = ["Link", "Load"]

    if table == "final_energy_demand":
        exclude_carriers = [
            "DC",
            "DC_OH",
            "electricity distribution grid",
            "H2 pipeline",
            "battery charger",
            "home battery charger",
        ]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier=elec_bus_carrier + ["gas", "H2"],
                groupby=["bus_carrier", "carrier"],
                nice_names=False,
                aggregate_across_components=True,
            )
            .loc[lambda x: ~x.index.get_level_values("carrier").isin(exclude_carriers)]
            .groupby(level="bus_carrier")
            .sum()
        )
    elif table == "elec_demand":
        df = n.statistics.withdrawal(
            comps=demand_comps,
            bus_carrier=elec_bus_carrier,
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        ).drop(
            index=[
                "DC",
                "DC_OH",
                "electricity distribution grid",
                "battery charger",
                "home battery charger",
            ]
        )
    elif table == "methane_demand":
        df = n.statistics.withdrawal(
            comps=demand_comps,
            bus_carrier="gas",
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        )
    elif table == "hydrogen_demand":
        df = n.statistics.withdrawal(
            comps=demand_comps,
            bus_carrier="H2",
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        )
    elif table == "power_capacity":
        df = (
            n.statistics.optimal_capacity(
                bus_carrier=elec_bus_carrier,
                groupby="carrier",
                nice_names=False,
                aggregate_across_components=True,
            )
            .loc[lambda x: x > 0]
            .drop(index=["battery discharger"])
        )
    elif table == "power_generation":
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier=elec_bus_carrier,
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        ).drop(
            index=[
                "DC",
                "DC_OH",
                "electricity distribution grid",
                "battery discharger",
                "home battery discharger",
            ]
        )
    elif table == "methane_supply":
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier="gas",
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        )
    elif table == "hydrogen_supply":
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier="H2",
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        )
    elif table == "biomass_supply":
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier="solid biomass",
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        )
    elif table == "energy_imports":
        # TODO cannot extract gas imports
        # TODO no biomass import is assumed
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier=["H2", "oil", "coal", "lignite"],
            groupby="carrier",
            nice_names=False,
            aggregate_across_components=True,
        ).reindex(
            ["coal", "lignite", "oil refining", "H2 import LH2", "H2 import Pipeline"]
        )
    elif table == "generation_profiles":
        df = (
            n.statistics.supply(
                bus_carrier=elec_bus_carrier,
                groupby="carrier",
                nice_names=False,
                aggregate_across_components=True,
                aggregate_time=False,
            )
            .drop(
                index=[
                    "DC",
                    "DC_OH",
                    "electricity distribution grid",
                ]
            )
            .melt(ignore_index=False)
            .reset_index()
            .set_index(["snapshot", "carrier"])
        )
    else:
        logger.warning(f"Unknown benchmark table: {table}")
        df = pd.DataFrame()

    df = (
        df.reset_index()
        .rename(columns={"bus_carrier": "carrier", 0: "value", "objective": "value"})
        .assign(table=table)
    )

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_benchmark",
            configfiles="config/config.tyndp.yaml",
            opts="",
            clusters="all",
            sector_opts="",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]

    # Read networks
    logger.info("Reading network")
    n = pypsa.Network(snakemake.input[0])

    logger.info("Building benchmark from network")
    tqdm_kwargs = {
        "ascii": False,
        "unit": " benchmark",
        "total": len(options["tables"]),
        "desc": "Computing benchmark",
    }

    func = partial(
        compute_benchmark,
        n,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        benchmarks = list(
            tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs)
        )

    # Combine all benchmark data
    benchmarks_combined = pd.concat(benchmarks, ignore_index=True).assign(
        year=snakemake.wildcards.planning_horizons,
        scenario=snakemake.params["scenario"],
        source="Open-TYNDP",
    )
    if benchmarks_combined.empty:
        logger.warning("No benchmark data was successfully processed")

    # Save data
    benchmarks_combined.to_csv(snakemake.output[0], index=False)
