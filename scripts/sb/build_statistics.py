# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script computes the benchmark statistics from the optimised network.
"""

import logging
import multiprocessing as mp
import re
from functools import partial

import country_converter as coco
import numpy as np
import pandas as pd
import pypsa
from tqdm import tqdm

from scripts._helpers import (
    ENERGY_UNITS,
    POWER_UNITS,
    PRICE_UNITS,
    configure_logging,
    normalize_direction,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

pypsa.options.params.statistics.nice_names = False


def add_benchmarking_mappings(carrier_mapping_fn: str, tables: dict) -> None:
    """
    Load benchmarking mappings from the carrier mapping file and add them into
    each table's config dict under the key ``mapping``.

    Parameters
    ----------
    carrier_mapping_fn : str
        Path to csv file with carrier mapping.
    tables : dict
        Dictionary of statistics tables. Each entry with a ``mapping_col`` key
        will have a ``mapping`` key added containing the loaded carrier mapping.

    Returns
    -------
    None
        Modifies the tables dictionary in place by adding the mapping for each table.
    """
    tech_map = pd.read_csv(carrier_mapping_fn)
    for table, table_opts in tables.items():
        col = table_opts.get("mapping_col")
        if col is None:
            continue
        if col not in tech_map.columns:
            logger.warning(
                f"No existing mapping for table {table} in 'tyndp_technology_mapping.csv'."
            )
            continue
        table_opts["mapping"] = (
            tech_map[["open_tyndp_index", "open_tyndp_carrier", col]]
            .assign(
                open_tyndp=lambda x: x[["open_tyndp_index", "open_tyndp_carrier"]]
                .bfill(axis=1)
                .iloc[:, 0]
            )
            .dropna(subset=["open_tyndp", col])
            .set_index("open_tyndp")[col]
            .to_dict()
        )


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
    planning_horizons: int,
    load_shedding: dict[str, int],
    low_voltage: bool,
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
    planning_horizons : int
        The current planning horizon year.
    load_shedding : dict[str, int]
        Value of lost load per carrier.
    low_voltage: bool
        If True, include low voltage buses in the electricity bus carrier list.

    Returns
    -------
    pd.DataFrame
        Benchmark data in long format.
    """
    opt = options["tables"][table]
    mapping = opt.get("mapping", {})
    elec_bus_carrier = ["AC", "AC_OH"] + (["low voltage"] if low_voltage else [])
    supply_comps = ["Generator", "Link"]
    demand_comps = ["Link", "Load"]
    eu27_idx = n.buses[n.buses.country.isin(eu27)].index

    # Optionally remove the last day of the year to have exactly 52 weeks
    if options["remove_last_day"]:
        sws = remove_last_day(n.snapshot_weightings.generators)
    else:
        sws = n.snapshot_weightings.generators

    if table == "final_energy_demand":
        grouper = ["bus_carrier"]
        df_countries = (
            n.statistics.withdrawal(
                comps="Load",
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                groupby_time=False,
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
                groupby_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[lambda s: ~s.index.isin(df_countries.index)]
        )

        # Biogas not upgraded to biomethane is part of the FED in Open-TYNDP
        biogas_not_upgraded = (
            options["tables"]["biomass_supply"]["biogas_not_upgraded"][
                planning_horizons
            ]
            * 1e6
        )
        df_eu.loc["solid biomass"] -= biogas_not_upgraded

        df = pd.concat([df_countries, df_eu])
    elif table == "elec_demand":
        grouper = ["carrier"]
        df = (
            n.statistics.withdrawal(
                comps=demand_comps,
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                groupby_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[pd.IndexSlice[:, ["electricity"]]]
            .reindex(eu27_idx, level="bus")
            .dropna()
        )
        df = df.groupby(by=grouper).sum()
    elif table == "methane_demand":
        grouper = ["carrier"]
        df_countries = (
            n.statistics.withdrawal(
                comps="Link",
                bus_carrier="gas",
                groupby=["bus1"] + grouper,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus1")
            .groupby(by=grouper)
            .sum()
        )

        df_eu = (
            n.statistics.withdrawal(
                comps="Link",
                bus_carrier="gas",
                groupby=["bus1"] + grouper,
                aggregate_across_components=True,
            )
            .loc[
                lambda s: (
                    ~s.index.get_level_values("carrier").isin(df_countries.index)
                )
                & (s.index.get_level_values("bus1").str.startswith("EU"))
            ]
            .groupby(by=grouper)
            .sum()
        )

        df = pd.concat([df_countries, df_eu])
    elif table == "hydrogen_demand":
        grouper = ["carrier"]
        exclusions = [
            "H2 pipeline",
            "H2 cavern-storage charger",
            "H2 tank-storage charger",
        ]
        df = n.statistics.withdrawal(
            comps=demand_comps,
            bus_carrier="H2",
            groupby=["bus"] + grouper,
            aggregate_across_components=True,
        ).loc[lambda df: ~df.index.get_level_values("carrier").isin(exclusions)]
    elif table == "power_capacity":
        grouper = ["carrier"]
        exclusions = ["electricity distribution grid", "DC", "DC_OH", "load"]
        df = (
            n.statistics.optimal_capacity(
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
            )
            .loc[lambda x: x > 0]
            .reset_index()
            .loc[lambda df: ~df.carrier.isin(exclusions)]
            .assign(bus=lambda df: df.bus.str.split(" ").str[0])
            .groupby(["bus"] + grouper)
            .sum()
            .iloc[:, 0]
        )

        # Add H2 offwind capacities in MW_e
        off_car = [c for c in tyndp_renewable_carriers if c.startswith("offwind-h2")]  # noqa: F841
        df_offwind_h2 = (
            n.generators.query("carrier.isin(@off_car)")
            .assign(
                p_nom_opt=lambda df: df.p_nom_opt / df.efficiency_dc_to_h2,
                bus=lambda df: df.bus.str.split(" ").str[0],
            )
            .groupby(by=["bus"] + grouper)
            .p_nom_opt.sum()
        )

        df = pd.concat([df, df_offwind_h2])
    elif table == "power_generation":
        grouper = ["carrier"]
        exclusions = [
            "DC",
            "DC_OH",
            "electricity distribution grid",
            "battery discharger",
            "battery charger",
            "home battery discharger",
            "home battery charger",
            "PHS",
            "hydro-phs-turbine",
            "hydro-phs-pump",
            "hydro-phs-pure-turbine",
            "hydro-phs-pure-pump",
            "H2 Electrolysis",
        ]
        df = (
            n.statistics.supply(
                comps=supply_comps + ["StorageUnit"],
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                groupby_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[lambda df: ~df.index.get_level_values("carrier").isin(exclusions)]
        )

        # TYNDP 2024 report available generation for renewables (pre-curtailment)
        # and add H2 offwind capacities in MWh_e
        # TODO Review once solar thermals are integrated
        res_carriers = n.carriers.filter(regex="offwind.*|solar.*|onwind", axis=0).index
        res_idx = n.generators[n.generators.carrier.isin(res_carriers)].index
        eff_dc_to_b0 = n.generators.loc[res_idx, "efficiency_dc_to_b0"].fillna(1)

        res_gen = (
            (sws @ (n.generators_t.p_max_pu[res_idx] * n.generators.p_nom_opt[res_idx]))
            .div(eff_dc_to_b0)
            .groupby([n.generators.bus, n.generators.carrier])
            .sum()
        )
        df = res_gen.combine_first(df)

        df = (
            df.rename(index=lambda x: x.split(" ")[0], level=0)
            .groupby(["bus"] + grouper)
            .sum()
        )

        # add curtailment to power generation statistics
        curtailment_exclusions = [
            "other-res-mix",
            "hydro-reservoir",
            "hydro-pondage",
            "hydro-ror",
            "load",
            "dsr",
        ]
        df_curtailment = (
            n.statistics.curtailment(
                bus_carrier=elec_bus_carrier,
                groupby=["bus"] + grouper,
                aggregate_across_components=True,
                groupby_time=False,
            )
            .mul(sws, axis=1)
            .sum(axis=1)
            .loc[
                lambda df: ~df.index.get_level_values("carrier").isin(
                    curtailment_exclusions
                )
            ]
            .rename(index=lambda x: x.removesuffix(" low voltage"), level="bus")
            .rename(index=lambda _: "dumped energy", level="carrier")
            .groupby(["bus"] + grouper)
            .sum()
        )
        df = pd.concat([df, df_curtailment])

    elif table == "methane_supply":
        grouper = ["carrier"]
        df_countries = (
            n.statistics.supply(
                comps="Link",
                bus_carrier="gas",
                groupby=["bus0"] + grouper,
                aggregate_across_components=True,
            )
            .reindex(eu27_idx, level="bus0")
            .groupby(by=grouper)
            .sum()
        )

        df_eu = n.statistics.supply(
            comps=supply_comps,
            bus_carrier="gas",
            groupby=grouper,
            aggregate_across_components=True,
        ).loc[lambda s: (~s.index.get_level_values("carrier").isin(df_countries.index))]

        df = pd.concat([df_countries, df_eu])
    elif table == "hydrogen_supply":
        grouper = ["carrier"]
        exclusions = [
            "H2 pipeline",
            "H2 pipeline OH",
            "H2 cavern-storage discharger",
            "H2 tank-storage discharger",
        ]
        df = n.statistics.supply(
            comps=supply_comps,
            bus_carrier="H2",
            groupby=["bus"] + grouper,
            aggregate_across_components=True,
        ).loc[lambda df: ~df.index.get_level_values("carrier").isin(exclusions)]
    elif table == "biomass_supply":
        grouper = ["carrier"]
        df_fed_btl = n.statistics.withdrawal(
            comps=demand_comps,
            bus_carrier="solid biomass",
            groupby=grouper,
            aggregate_across_components=True,
        )

        # Biogas not upgraded to biomethane is part of the FED in Open-TYNDP
        biogas_not_upgraded = opt["biogas_not_upgraded"][planning_horizons] * 1e6
        df_fed_btl.loc["biomass final energy demand"] -= biogas_not_upgraded
        df_fed_btl.loc["for biomethane"] = biogas_not_upgraded

        eff = float(opt["biomass_to_methane_efficiency"][planning_horizons])

        df_biogas = n.statistics.energy_balance(
            comps="Generator",
            bus_carrier="biogas",
            groupby=grouper,
            aggregate_across_components=True,
        ).div(eff)

        df = pd.concat([df_fed_btl, df_biogas])
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
                    "H2 cavern-storage discharger",
                    "H2 tank-storage discharger",
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
                    groupby_time=False,
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
    elif table in ["electricity_price", "hydrogen_price"]:
        carrier = "AC" if "electricity" in table else "H2"

        # Collect nodes with energy balance; nodes without balance will have the load shedding price
        idx_w_balance = n.statistics.energy_balance(
            bus_carrier=carrier, groupby="bus", aggregate_across_components=True
        ).index

        df = (
            n.statistics.prices(
                bus_carrier=carrier,
                weighting="time",
            )
            .to_frame("value")
            .rename_axis("bus")
            .assign(carrier=carrier)
            .set_index("carrier", append=True)
        )

        mask = ~df.index.get_level_values("bus").isin(idx_w_balance)
        df.loc[mask, :] = np.nan
    elif table in ["electricity_price_excl_shed", "hydrogen_price_excl_shed"]:
        carrier = "AC" if "electricity" in table else "H2"

        df = (
            n.statistics.prices(
                bus_carrier=carrier,
                weighting="time",
                groupby_time=False,
            )
            .pipe(lambda x: x.where(x < load_shedding.get(carrier, np.inf)))
            .mean(axis=1)
            .to_frame("value")
            .rename_axis("bus")
            .assign(carrier=carrier)
            .set_index("carrier", append=True)
        )
    elif table in ["crossborder_electricity", "crossborder_hydrogen"]:
        carrier = "AC" if "elec" in table else "H2"
        bus_carrier = [carrier, f"{carrier}_OH"]  # noqa F814
        if carrier == "H2":
            bus_carrier.extend(["import H2"])
        connector = " -> " if carrier == "H2" else "-"

        idx_b = n.buses.query("carrier.isin(@bus_carrier)").index  # noqa F841
        idx_l = n.links.query("bus0.isin(@idx_b) and bus1.isin(@idx_b)").index

        df = sws @ n.links_t.p0[idx_l]
        if carrier == "H2":
            df = df.rename(
                lambda x: re.sub(r"\b([A-Z]+)00\b", r"\1", x).replace("UK", "GB")
            )
        df = normalize_direction(
            df, buses_from_index=True, connector=connector, format_index=True
        )

        df = (
            df.to_frame("value")
            .assign(carrier=carrier)
            .set_index("carrier", append=True)
        )
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

    # Add EU27 (load-weighted average for prices)
    if "bus" in [c for c in ["bus", "carrier", "snapshot"] if c in df.columns]:
        if "price" in table:
            df_eu = df
            bus_name = "Pan-EU"
            weights = (
                n.statistics.withdrawal(
                    groupby="bus",
                    bus_carrier=carrier,
                    nice_names=False,
                )
                .groupby("bus")
                .sum()
                .T
            )
            normalizer = weights.sum()
        else:
            df_eu = df.loc[lambda x: x["bus"].isin(eu27_idx)]
            bus_name = "EU27"
            weights = pd.Series(1.0, index=df_eu.bus.unique())
            normalizer = 1.0

        df_eu = (
            df_eu.assign(value=lambda x: x.bus.map(weights).fillna(0) * x.value)
            .groupby(by=[c for c in ["carrier", "snapshot"] if c in df.columns])
            .value.sum()
            .div(normalizer)
            .reset_index()
            .assign(bus=bus_name)
        )
        df = pd.concat([df, df_eu])
    elif "border" not in df.columns:
        df = df.assign(bus="EU27")

    op = "sum" if "price" not in table else "mean"
    df = (
        df.groupby(
            by=[c for c in ["bus", "carrier", "snapshot", "border"] if c in df.columns]
        )
        .agg(op)
        .reset_index()
        .assign(
            table=table,
            unit=lambda x: "MWh"
            if opt["unit"] in ENERGY_UNITS
            else "MW"
            if opt["unit"] in POWER_UNITS
            else "EUR/MWh"
            if opt["unit"] in PRICE_UNITS
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
    load_shedding = snakemake.params["load_shedding"]
    low_voltage = snakemake.params["low_voltage"]
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    # Mapping from Open-TYNDP carrier names to benchmarking carrier names
    add_benchmarking_mappings(snakemake.input.carrier_mapping, options["tables"])

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
        planning_horizons=planning_horizons,
        load_shedding=load_shedding,
        low_voltage=low_voltage,
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
