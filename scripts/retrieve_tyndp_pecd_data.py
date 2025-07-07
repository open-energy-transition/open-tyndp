# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
The TYNDP PECD data contains input data for the 2024 TYNDP scenario building process.

This rule downloads the TYNDP PECD v3.1 data from Google Drive and extracts it in the ``data/tyndp_2024_bundle``
subdirectory, such that all files of the TYNDP bundle are stored in it.

**Outputs**

- ``data/tyndp_2024_bundle/PECD``: PECD input data for TYNDP 2024 scenario building

"""

import logging
import os
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
# TODO: retrieve_tyndp_pecd_data needs to be deprecated once PECD data is added to the TYNDP data bundle
url = "https://storage.googleapis.com/open-tyndp-data-store/PECD.zip"

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_pecd_data")
        rootpath = ".."
    else:
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = snakemake.output.dir
    tyndp_bundle_fn = Path(rootpath, snakemake.params["tyndp_bundle"])
    to_fn_zp = to_fn + ".zip"

    # download .zip file
    logger.info(f"Downloading TYNDP PECD data from '{url}'.")
    progress_retrieve(url, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting TYNDP PECD data.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(tyndp_bundle_fn)

    # remove .zip file
    os.remove(to_fn_zp)

    logger.info(f"TYNDP PECD data available in '{to_fn}'.")
