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
import pandas as pd
import pypsa
from tqdm import tqdm

from scripts._helpers import (
    ENERGY_UNITS,
    POWER_UNITS,
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

pypsa.options.params.statistics.nice_names = False


def remove_last_day(sws: pd.Series, nhours: int = 24):
    """
    Remove the last day from snapshots to ensure exactly 52 weeks of data.

    Parameters
    ----------
    sws : pd.Series
        Snapshot weightings.
    nhours : int, default 24
        Number of hours to consider.

    Returns
    -------
    tuple[pd.DatetimeIndex, pd.Series]
        Modified snapshots and snapshot weightings with the last day removed.
    """
    sws = sws.copy()

    remaining_hours = sws.iloc[::-1].cumsum() - nhours
    sws[remaining_hours < 0] = 0
    last_i = remaining_hours[remaining_hours >= 0].index[0]
    sws.loc[last_i] = remaining_hours.loc[last_i]

    return sws


def compute_benchmark(
    n: pypsa.Network,
    table: str,
    options: dict,
    eu27: list[str],
    tyndp_renewable_carriers: list[str],
) -> pd.DataFrame:
    """
    Compute benchmark metrics from optimized network.

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
    tyndp_renewable_carriers : list[str]
        List of renewable carriers in TYNDP 2024.

    Returns
    -------
    pd.DataFrame
        Benchmark data in long format.
    """
    opt = options["tables"][table]
    mapping = opt.get("mapping", {})
    elec_bus_carrier = ["AC", "AC_OH", "low voltage"]
    supply_comps = ["Generator", "Link"]
    demand_comps = ["Link", "Load"]
    eu27_idx = n.buses[n.buses.country.isin(eu27)].index

    # Optionally remove the last day of the year to have exactly 52 weeks
    if options["remove_last_day"]:
        sws = remove_last_day(n.snapshot_weightings.generators)
    else:
        sws = n.snapshot_weightings.generators

    if table == "final_energy_demand":
        # TODO Clarify what renewables encompass
        grouper = ["bus_carrier"]
        df_countries = (
            n.statistics.withdrawal(
                comps="Load",
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                aggregate_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .reindex(eu27_idx, level="bus")
            .groupby(level="bus_carrier")
            .sum()
        )

        # Add EU level demands
        df_eu = (
            n.statistics.withdrawal(
                comps="Load",
                groupby=grouper,
                aggregate_across_components=True,
                aggregate_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[lambda s: ~s.index.isin(df_countries.index)]
            .drop(index=["solid biomass"], errors="ignore")
        )

        # Add Biomass to Liquid
        df_btl = (
            n.statistics.withdrawal(
                comps="Link",
                carrier=["biomass to liquid", "biomass to liquid CC"],
                groupby="carrier",
                aggregate_across_components=True,
                aggregate_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .rename_axis("bus_carrier")
        )

        df = pd.concat([df_countries, df_eu, df_btl])
    elif table == "elec_demand":
        grouper = ["carrier"]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                aggregate_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[pd.IndexSlice[:, ["electricity"]]]
            .reindex(eu27_idx, level="bus")
            .dropna()
        )
        df = df.groupby(by=grouper).sum()
    elif table == "methane_demand":
        # TODO Energy and non-energy industrial demand are mixed
        grouper = ["carrier"]
        if "EU gas" in n.buses.index:
            reindex_idx = pd.Index(["EU gas"])
        else:
            reindex_idx = eu27_idx
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier="gas",
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
            )
            .reindex(reindex_idx, level="bus")
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
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .drop(index="H2 pipeline", errors="ignore")
        )
    elif table == "power_capacity":
        grouper = ["carrier"]
        df = (
            n.statistics.optimal_capacity(
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .loc[lambda x: x > 0]
            .drop(index=["electricity distribution grid"], errors="ignore")
        )

        # Add H2 offwind capacities in MW_e
        off_car = [c for c in tyndp_renewable_carriers if c.startswith("offwind-h2")]  # noqa: F841
        df_offwind_h2 = (
            n.generators.query("carrier.isin(@off_car)")
            .assign(p_nom_opt=lambda df: df.p_nom_opt / df.efficiency_dc_to_h2)
            .groupby(by=["bus"] + grouper)
            .p_nom_opt.sum()
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
        )

        df = pd.concat([df, df_offwind_h2])
    elif table == "power_generation":
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps + ["StorageUnit"],
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
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
                    "PHS",
                    "hydro-phs",  # TODO: double check when correct hydro techs are added
                    "hydro-phs-pure",
                ],
                errors="ignore",
            )
        )
    elif table == "methane_supply":
        grouper = ["carrier"]
        df = (
            n.statistics.supply(
                comps=supply_comps,
                bus_carrier="gas",
                groupby=["bus"] + grouper,
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
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
        )
    elif table == "energy_imports":
        # TODO Account for domestic production of gas, solid and liquid fossil fuels
        # TODO No biomass import is assumed
        grouper = ["carrier"]
        df_countries = (
            n.statistics.supply(
                comps="Link",
                bus_carrier=["H2"],
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus")
            .groupby(by=grouper)
            .sum()
            .drop(
                index=[
                    "H2 Electrolysis",
                    "H2 pipeline",
                    "H2 pipeline OH",
                    "SMR",
                    "SMR CC",
                ],
                errors="ignore",
            )
        )

        # Add EU level demands
        df_eu = n.statistics.supply(
            comps="Generator",
            bus_carrier=["oil primary", "coal", "lignite", "gas"],
            groupby=grouper,
            aggregate_across_components=True,
        ).loc[lambda s: ~s.index.isin(df_countries.index)]

        df = pd.concat([df_countries, df_eu])
    elif table == "generation_profiles":
        if n.snapshots.year[0] == 2009:
            grouper = ["carrier"]
            df = (
                n.statistics.supply(
                    bus_carrier=elec_bus_carrier,
                    groupby=["bus"] + grouper,
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
                    ],
                    errors="ignore",
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
        .rename(
            columns={
                "bus_carrier": "carrier",
                0: "value",
                "objective": "value",
            }
        )
        .assign(carrier=lambda x: x["carrier"].map(mapping).fillna(x["carrier"]))
    )
    grouper = [c for c in ["carrier", "snapshot"] if c in df.columns]
    df = (
        df.groupby(by=grouper)
        .sum()
        .reset_index()
        .assign(
            table=table,
            unit=lambda x: "MWh"
            if opt["unit"] in ENERGY_UNITS
            else "MW"
            if opt["unit"] in POWER_UNITS
            else opt["unit"],
        )
    )

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
    tyndp_renewable_carriers = snakemake.params["tyndp_renewable_carriers"]
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    # Read network
    logger.info("Reading network")
    n = pypsa.Network(snakemake.input.network)

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
        options=options,
        eu27=eu27,
        tyndp_renewable_carriers=tyndp_renewable_carriers,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        benchmarks = list(
            tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs)
        )

    # Combine all benchmark data
    benchmarks_combined = pd.concat(benchmarks, ignore_index=True).assign(
        year=planning_horizons,
        scenario="TYNDP " + snakemake.params["scenario"],
        source="Open-TYNDP",
    )
    if benchmarks_combined.empty:
        logger.warning("No benchmark data was successfully processed")

    # Save data
    benchmarks_combined.to_csv(snakemake.output[0], index=False)
