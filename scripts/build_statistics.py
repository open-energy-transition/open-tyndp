# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script computes the benchmark statistics from the optimised network.
"""

import logging
import multiprocessing as mp
from functools import partial

import country_converter as coco
import numpy as np
import pandas as pd
import pypsa
from tqdm import tqdm

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def get_loss_factors(fn: str, n: pypsa.Network, planning_horizons: str) -> pd.Series:
    """
    Load and prepare loss factors

    Parameters
    ----------
    fn : str
        Path to the file containing loss factors.
    n : pypsa.Network
        Network to use.
    planning_horizons : str
        Planning horizon for which to read the data.

    Returns
    -------
    pd.Series
        Loss factors data.
    """
    # Read data
    pyear = np.clip(5 * (planning_horizons // 10), 2030, 2050)
    loss_factors = pd.read_csv(fn, index_col=0)[pyear]

    # Create index map
    idx_map = n.buses.query("Bus.str.contains('low voltage')").country
    loss_factors = idx_map.map(loss_factors).dropna()

    return loss_factors


def compute_benchmark(
    n: pypsa.Network,
    table: str,
    options: dict,
    eu27: list[str],
    loss_factors: pd.Series = pd.Series(),
) -> pd.DataFrame:
    """
    Compute benchmark metrics from optimised network.

    Parameters
    ----------
    n : pypsa.Network
        Optimised network.
    table : str
        Benchmark metric to compute.
    options : dict
        Full benchmarking configuration.
    eu27 : list[str]
        List of member state of European Union (EU27).
    loss_factors : pd.Series, optional
        Series containing loss factors indexed by country.

    Returns
    -------
    pd.DataFrame
        Benchmark data in long format.
    """
    opt = options["tables"][table]
    map = opt.get("mapping", {})
    elec_bus_carrier = ["AC", "AC_OH", "low voltage"]
    supply_comps = ["Generator", "Link"]
    demand_comps = ["Link", "Load"]
    eu27_idx = n.buses[n.buses.country.isin(eu27)].index

    if table == "final_energy_demand":
        # TODO Clarify what renewables encompass
        grouper = ["bus_carrier", "carrier"]
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
                bus_carrier=elec_bus_carrier + ["gas", "H2", "coal", "oil"],
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .loc[lambda x: ~x.index.get_level_values("carrier").isin(exclude_carriers)]
            .groupby(level="bus_carrier")
            .sum()
        )
    elif table == "elec_demand":
        grouper = ["carrier"]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .loc[pd.IndexSlice[:, ["electricity"]]]
            .reindex(eu27_idx, level="bus")
            .dropna()
        )
        if not loss_factors.empty:
            df /= 1 + loss_factors.reindex(df.index, level="bus")
        df = df.groupby(by=grouper).sum()
    elif table == "methane_demand":
        # TODO Energy and non-energy industrial demand are mixed
        grouper = ["carrier"]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier="gas",
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
        )
    elif table == "hydrogen_demand":
        # TODO Energy and non-energy industrial demand are mixed
        # TODO Aviation has no H2 demand
        grouper = ["carrier"]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier="H2",
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .drop(index="H2 pipeline")
        )
    elif table == "power_capacity":
        grouper = ["carrier"]
        df = (
            n.statistics.optimal_capacity(
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .loc[lambda x: x > 0]
            .drop(index=["electricity distribution grid"], errors="ignore")
        )
    elif table == "power_generation":
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .drop(
                index=[
                    "DC",
                    "DC_OH",
                    "electricity distribution grid",
                    "battery discharger",
                    "home battery discharger",
                ]
            )
        )
    elif table == "methane_supply":
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier="gas",
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
        )
    elif table == "hydrogen_supply":
        # TODO Clarify difference between low carbon and renewable imports
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier="H2",
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .drop(index=["H2 pipeline", "H2 pipeline OH"])
        )
    elif table == "biomass_supply":
        # TODO Clarify how to deal with unsustainable sources
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier="solid biomass",
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
        )
    elif table == "energy_imports":
        # TODO Cannot extract gas imports
        # TODO No biomass import is assumed
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier=["H2", "oil", "coal", "lignite"],
                groupby=["bus"] + grouper,
                nice_names=False,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .reindex(
                [
                    "coal",
                    "lignite",
                    "oil refining",
                    "H2 import LH2",
                    "H2 import Pipeline",
                ]
            )
        )
    elif table == "generation_profiles":
        if n.snapshots.year[0] == 2009:
            grouper = ["carrier"]
            df = (
                n.statistics.supply(
                    bus_carrier=elec_bus_carrier,
                    groupby=["bus"] + grouper,
                    nice_names=False,
                    aggregate_across_components=True,
                    aggregate_time=False,
                )
                .reindex(eu27_idx, level="bus")
                .groupby(by=grouper)
                .sum()
                .drop(
                    index=[
                        "DC",
                        "DC_OH",
                        "electricity distribution grid",
                        "H2 Electrolysis",
                        "battery charger",
                        "home battery charger",
                        "methanolisation",
                        "electricity",
                    ]
                )
                .melt(ignore_index=False)
                .reset_index()
                .set_index(["snapshot", "carrier"])["value"]
            )
        else:
            logger.warning(f"Unknown climate year for table: {table}")
            df = pd.DataFrame(columns=["carrier"])
    else:
        logger.warning(f"Unknown benchmark table: {table}")
        df = pd.DataFrame(columns=["carrier"])

    df = (
        df.reset_index()
        .rename(columns={"bus_carrier": "carrier", 0: "value", "objective": "value"})
        .assign(carrier=lambda x: x["carrier"].map(map).fillna(x["carrier"]))
    )
    grouper = [c for c in ["carrier", "snapshot"] if c in df.columns]
    df = df.groupby(by=grouper).sum().reset_index().assign(table=table)

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_statistics",
            opts="",
            clusters="all",
            sector_opts="",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()
    planning_horizons = str(snakemake.wildcards.planning_horizons)
    loss_factors_fn = snakemake.input.loss_factors

    # Read network
    logger.info("Reading network")
    n = pypsa.Network(snakemake.input.network)

    # Load loss factors
    loss_factors = get_loss_factors(loss_factors_fn, n, planning_horizons)

    logger.info("Building benchmark from network")
    tqdm_kwargs = {
        "ascii": False,
        "unit": " benchmark",
        "total": len(options["tables"]),
        "desc": "Computing benchmark",
    }

    func = partial(
        compute_benchmark, n, options=options, eu27=eu27, loss_factors=loss_factors
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        benchmarks = list(
            tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs)
        )

    # Combine all benchmark data
    benchmarks_combined = pd.concat(benchmarks, ignore_index=True).assign(
        year=snakemake.wildcards.planning_horizons,
        scenario="TYNDP " + snakemake.params["scenario"],
        source="Open-TYNDP",
    )
    if benchmarks_combined.empty:
        logger.warning("No benchmark data was successfully processed")

    # Save data
    benchmarks_combined.to_csv(snakemake.output[0], index=False)
