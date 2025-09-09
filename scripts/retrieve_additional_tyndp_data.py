# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Retrieves additional TYNDP data not included in the Zenodo TYNDP data bundle.
Downloads the zip file from Google Drive and extracts it in the ``data/tyndp_2024_bundle``
subdirectory, such that all files of the TYNDP bundle are stored in it.
The original data is published by ENTSO-E and ENTSOG under Creative Commons Attribution 4.0 International License (CC-BY 4.0)
and can be found under https://2024.entsos-tyndp-scenarios.eu/download/.

Currently, this is used for two additional datasets:
* TYNDP PECD data: The TYNDP PECD v3.1 data contains input data for the 2024 TYNDP scenario building process.
* TYNDP hydro inflows: The TYNDP hydro inflow data from PEMMDB v2.4 contains hydro inflow data for different hydro technologies:
    * Run of River,
    * Pondage,
    * Reservoir,
    * PS Open,
    * PS Closed

**Outputs**

- ``data/tyndp_2024_bundle/<additional-dataset>``: Additional input dataset for TYNDP 2024 scenario building

"""

import logging
import os
import shutil
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# TODO: retrieve_additional_tyndp_data needs to be deprecated once all TYNDP data is added to the TYNDP data bundle


def retrieve_bundle(url: str, source: str, to_dir: str, disable_progress: bool = False):
    to_fn = Path(to_dir, Path(url).name)

    # download data
    logger.info(f"Downloading TYNDP {source} data from '{url}'.")
    progress_retrieve(url, to_fn, disable=disable_progress)

    # extract if needed
    if (
        zipfile.is_zipfile(to_fn) and to_fn.suffix.lower() != ".xlsx"
    ):  # xlsx are valid zip files
        logger.info(f"Extracting TYNDP {source} data.")
        with zipfile.ZipFile(to_fn, "r") as zip_ref:
            zip_ref.extractall(to_dir)
        os.remove(to_fn)
        shutil.rmtree(Path(to_dir, "__MACOSX"), ignore_errors=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_hydro_inflow")

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    url = snakemake.params.url
    source = snakemake.params.source
    to_dir = snakemake.output.dir

    retrieve_bundle(url, source, to_dir, disable_progress)

    logger.info(f"TYNDP {source} data available in '{to_dir}'.")
