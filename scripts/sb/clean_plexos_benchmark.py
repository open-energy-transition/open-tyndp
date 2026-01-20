# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the PLEXOS output data for benchmarking.

Reads PLEXOS MMStandardOutputFile xlsx files from both:
- "Yearly Outputs" sheet for power and electricity data
- "Yearly H2 Outputs" sheet for hydrogen-specific data

Note: Currently, only NT scenario processing is supported.
"""

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import country_converter as coco
import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


# Mapping from PLEXOS technology names to benchmark carrier names
PLEXOS_CARRIER_MAPPING = {
    # Power capacity and generation mappings
    "Nuclear": "nuclear",
    "Lignite old 1": "coal + other fossil",
    "Lignite old 2": "coal + other fossil",
    "Lignite new": "coal + other fossil",
    "Lignite CCS": "coal + other fossil",
    "Lignite biofuel": "biofuels",
    "Hard coal old 1": "coal + other fossil",
    "Hard coal old 2": "coal + other fossil",
    "Hard coal new": "coal + other fossil",
    "Hard coal CCS": "coal + other fossil",
    "Hard Coal biofuel": "biofuels",
    "Gas conventional old 1": "methane",
    "Gas conventional old 2": "methane",
    "Gas CCGT old 1": "methane",
    "Gas CCGT old 2": "methane",
    "Gas CCGT new": "methane",
    "Gas CCGT CCS": "methane",
    "Gas CCGT present 1": "methane",
    "Gas CCGT present 2": "methane",
    "Gas OCGT old": "methane",
    "Gas OCGT new": "methane",
    "Gas biofuel": "biofuels",
    "Light oil": "oil",
    "Heavy oil old 1": "oil",
    "Heavy oil old 2": "oil",
    "Light oil biofuel": "biofuels",
    "Heavy oil biofuel": "biofuels",
    "Oil shale old": "oil",
    "Oil shale new": "oil",
    "Oil shale biofuel": "biofuels",
    "Wind Onshore": "wind onshore",
    "Wind Offshore": "wind offshore",
    "Solar (Photovoltaic)": "solar",
    "Solar (Thermal)": "solar thermal",
    "Run-of-River": "hydro (exc. pump storage)",
    "Reservoir": "hydro (exc. pump storage)",
    "Pondage": "hydro (exc. pump storage)",
    "Pump Storage - Open Loop (turbine)": "hydro and pumped storage",
    "Pump Storage - Closed Loop (turbine)": "hydro and pumped storage",
    "Pump Storage - Open Loop (pump)": "hydro and pumped storage",
    "Pump Storage - Closed Loop (pump)": "hydro and pumped storage",
    "Others renewable": "chp and small thermal",  # >TODO check if that could be small scale res
    "Others non-renewable": "other non-res",
    "Battery Storage discharge (gen.)": "battery",
    "Battery Storage charge (load)": "battery",
    "Hydrogen CCGT": "hydrogen",
    "Hydrogen Fuel Cell": "hydrogen",
    "Demand Side Response Explicit": "demand shedding",
    "Demand Side Response Implicit": "demand shedding",
    # Hydrogen_demand
    "Native Demand (excl. H2 storage charge) [GWhH2]": "exogenous demand",
    # "power generation" - could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
    # "e-fuels" - not available in PLEXOS NT scenario, included in native demand
    # Hydrogen supply from H2 sheet
    "Electrolyser (gen.)": "p2g",
    "Steam methane reformer": "smr (grey)",
    # Note: PLEXOS doesn't distinguish between grey and blue SMR in the output files
    # "Exchanges with non-modeled nodes" → "imports (renewable & low carbon)" (handled separately)
    # NOT available in PLEXOS (will not appear in output):
    # - "ammonia imports"
    # - "bi product"
    # - "smr with ccs (blue)" - PLEXOS only has generic SMR
    # - "undefined for generation" . could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
}


# look up dictionary {name of plot: [sheet_name, output_type]}
LOOKUP_TABLES: dict[str, list[str]] = {
    "power_capacity": ["Yearly Outputs", "Installed Capacities [MW]"],
    "power_generation": ["Yearly Outputs", "Annual generation [GWh]"],
    "hydrogen_demand": [
        "Yearly H2 Outputs",
        "Native Demand (excl. H2 storage charge) [GWhH2]",
    ],
    "hydrogen_supply": ["Yearly H2 Outputs", "Annual generation [GWhH2]"],
}


def load_plexos_sheet(
    table_name: str,
    filepath: str | Path,
    eu27: list,
    skiprows: int = 5,
) -> pd.DataFrame:
    """
    Read PLEXOS sheet from xlsx file.

    Annual sheets have the same structure:
    - Rows 1-4: Metadata (Scenario, Simulator, Date, Status)
    - Row 5: Blank
    - Row 6: Headers (Output type, Output type, then country codes)
    - From row 7: Data rows with [output_type, technology, values...]

    Parameters
    ----------
    filepath : str or Path
        Path to the PLEXOS xlsx file.
    table_name : str
        Name of the table from LOOKUP_TABLES (e.g., "power_capacity").
    eu27 : list
        List of EU27 country codes.
    skiprows : int, default 5
        Number of metadata rows to skip.

    Returns
    -------
    pd.DataFrame
        DataFrame with multi-index rows [output_type, technology] and country columns.
    """
    sheet_name, output_type = LOOKUP_TABLES[table_name]

    df = pd.read_excel(
        filepath,
        sheet_name=sheet_name,
        skiprows=skiprows,
        header=0,
    )

    # Set multi-index
    level0 = df.iloc[:, 0].ffill()
    level1 = df.iloc[:, 1].fillna(level0)

    # Set as multiindex and keep remaining columns
    df.drop(df.columns[:2], inplace=True, axis=1)
    df.index = pd.MultiIndex.from_arrays([level0, level1])
    df.index.names = ["output_type", "carrier"]

    # Rename column names to country
    df.columns = df.columns.str[:2]

    # Only consider EU27 and sum
    df = df.T.groupby(df.columns).sum().reindex(eu27).sum()

    # Rename and group
    df = df.rename(index=PLEXOS_CARRIER_MAPPING, level=1).groupby(level=[0, 1]).sum()

    final = df.loc[output_type].rename("value").reset_index()
    final["table"] = table_name
    return final


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_plexos_benchmark",
            planning_horizons="2030",
            scenario="NT",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    scenario = snakemake.params["scenario"]
    planning_horizon = int(snakemake.wildcards.planning_horizons)

    # currently only implemented for NT
    if scenario != "NT":
        logger.warning(
            "Processing of TYNDP output files currently only implemented for NT scenario. Exporting empty Data Frame"
        )
        pd.DataFrame().to_csv(snakemake.output.benchmarks)

    # EU27 country codes
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()

    # PLEXOS input file
    plexos_file = snakemake.input.plexos_file

    # Plots for which Plexos output files provide data
    tables_to_process = [
        t for t in LOOKUP_TABLES.keys() if t in options["tables"].keys()
    ]

    logger.info(f"Processing tables: {", ".join(tables_to_process)}")

    tqdm_kwargs = {
        "ascii": False,
        "unit": " table",
        "total": len(tables_to_process),
        "desc": "Loading PLEXOS benchmark data",
    }

    func = partial(load_plexos_sheet, filepath=plexos_file, eu27=eu27, skiprows=5)

    # Process in parallel
    with mp.Pool(processes=snakemake.threads) as pool:
        benchmarks = list(tqdm(pool.imap(func, tables_to_process), **tqdm_kwargs))

    # Convert to dictionary
    plexos_data = pd.concat(benchmarks, ignore_index=True)

    plexos_data["scenario"] = f"TYNDP {scenario}"
    plexos_data["year"] = planning_horizon
    plexos_data["source"] = "Plexos output"

    # Save data
    plexos_data.to_csv(snakemake.output.benchmarks, index=False)
    logger.info(f"\nSaved to: {snakemake.output.benchmarks}")
