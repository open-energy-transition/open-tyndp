# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Loads and filters the available raw PECD data for the subset of required weather scenarios as specified in the configuration.

Outputs
-------
PECD prebuilt directory with filtered csv files including only relevant weather scenario data.
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

    skiprows = 10

    def _usecols(name):

        return (
            name in ("Date", "Hour")
            or name in cyears.values
            or name in cyears.astype(str).values
        )

    if "xls" in pecd_file or "xlsx" in pecd_file:
        df = pd.read_excel(
            fn,
            skiprows=skiprows,  # first rows contain only file metadata
            usecols=lambda name: _usecols(name),
            engine="openpyxl",
        )
    else:
        df = pd.read_csv(
            fn,
            skiprows=skiprows,  # first rows contain only file metadata
            usecols=lambda name: _usecols(name),
        )
    output_file = Path(output_dir, pecd_file).with_suffix(".csv")

    df.to_csv(output_file, index=False)

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_pecd_release",
            clusters="all",
        )
    configure_logging(snakemake)  # pylint: disable=possibly-used-before-assignment
    set_scenario_config(snakemake)

    # Parameters
    ############

    # Weather scenarios from snakemake params
    cyears = pd.Series(snakemake.params.cyears)
    available_cyears = [f"WS{x:03d}" for x in np.arange(1, 120, 1)]
    if set(cyears).difference(available_cyears):
        logger.warning(
            "Weather scenarios doesn't match available TYNDP data. Only returning subset of available weather scenarios."
        )
        cyears = pd.Series(list(set(cyears).intersection(available_cyears)))
    # Planning years for which PECD data is available for in the specified PECD version
    available_pyears = snakemake.params.available_pyears
    # Input and output directories and prebuilt version
    dir_pecd = snakemake.input.pecd_raw
    prebuilt_dir = snakemake.output.pecd_prebuilt

    # Iterate over available planning years
    #######################################
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
            "desc": f"Building PECD prebuilt for year {year}",
        }

        func = partial(
            process_pecd_files,
            dir_pecd=dir_pecd_year,
            output_dir=output_dir,
            cyears=cyears,
        )

        with mp.Pool(processes=snakemake.threads) as pool:
            list(tqdm(pool.imap(func, pecd_files), **tqdm_kwargs))
