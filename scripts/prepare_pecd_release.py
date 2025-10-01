# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and filters the available raw PECD data for the subset of required climate years as specified in the configuration.

Outputs
-------
PECD prebuilt directory with filtered csv files including only relevant climate year data.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def process_pecd_files(
    pecd_file: str,
    dir_pecd: Path,
    output_dir: Path,
    cyears: pd.Series,
) -> pd.DataFrame:
    fn = Path(dir_pecd, pecd_file)

    # Malta CSP data file has an extra header row that must be skipped
    if pecd_file == "PECD_CSP_noStorage_2040_MT00_edition 2023.2.csv":
        skiprows = 11
    else:
        skiprows = 10

    if "xls" in pecd_file or "xlsx" in pecd_file:
        df = pd.read_excel(
            fn,
            skiprows=skiprows,  # first rows contain only file metadata
            usecols=lambda name: name == "Date"
            or name == "Hour"
            or name in cyears.values
            or name in cyears.astype(str).values
            or name in cyears.astype(float).astype(str).values,
            engine="openpyxl",
        ).rename(columns={str(float(cyear)): str(cyear) for cyear in cyears})
    else:
        df = pd.read_csv(
            fn,
            skiprows=skiprows,  # first rows contain only file metadata
            usecols=lambda name: name == "Date"
            or name == "Hour"
            or name in cyears.values
            or name in cyears.astype(str).values
            or name in cyears.astype(float).astype(str).values,
        ).rename(columns={str(float(cyear)): str(cyear) for cyear in cyears})

    output_file = Path(
        output_dir, pecd_file.replace(".xlsx", ".csv").replace(".xls", ".csv")
    )

    df.to_csv(output_file, index=False)

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_pecd_release",
            clusters="all",
            PECD_PREBUILT_VERSION="3.1+pre-built.0.1",
        )
    configure_logging(snakemake)  # pylint: disable=possibly-used-before-assignment
    set_scenario_config(snakemake)

    # Parameters
    # Climate year from snapshots
    cyears = pd.Series(snakemake.params.cyears).astype(int)
    available_cyears = np.arange(1982, 2020, 1)
    if set(cyears).difference(available_cyears):
        logger.warning(
            "Climate year doesn't match available TYNDP data. Only returning subset of available climate years."
        )
        cyears = pd.Series(list(set(cyears).intersection(available_cyears)))
    # Planning years for which PECD data is available for in the specified PECD version
    available_pyears = snakemake.params.available_pyears
    # Input and output directories and prebuilt version
    dir_pecd = snakemake.input.pecd_raw
    prebuilt_dir = snakemake.output.pecd_prebuilt
    prebuilt_version = snakemake.wildcards.pecd_prebuilt_version

    # Iterate over available planning years
    for year in available_pyears:
        dir_pecd_year = Path(dir_pecd, str(year))
        pecd_files = [
            f
            for f in os.listdir(dir_pecd_year)
            if os.path.isfile(Path(dir_pecd_year, f))
        ]
        output_dir = Path(prebuilt_dir, str(year))

        # create output directory to save new files in
        os.makedirs(output_dir, exist_ok=True)

        # Load and prep pecd data
        tqdm_kwargs = {
            "ascii": False,
            "unit": " nodes",
            "total": len(pecd_files),
            "desc": f"Building PECD prebuilt version {prebuilt_version} for year {year}",
        }

        func = partial(
            process_pecd_files,
            dir_pecd=dir_pecd_year,
            output_dir=output_dir,
            cyears=cyears,
        )

        with mp.Pool(processes=snakemake.threads) as pool:
            list(tqdm(pool.imap(func, pecd_files), **tqdm_kwargs))
