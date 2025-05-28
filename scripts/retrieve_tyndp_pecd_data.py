# SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
The TYNDP PECD data contains input data for the 2024 TYNDP scenario building process.

This rule downloads the TYNDP PECD data from `a gdrive and extracts it in the ``data/tyndp_2024_bundle``
subdirectory, such that all files of the TYNDP bundle are stored in it.

**Outputs**

- ``data/tyndp_2024_bundle/PECD``: PECD input data for TYNDP 2024 scenario building

"""

import logging
import os
import zipfile

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
url = "https://drive.usercontent.google.com/download?id=15YF5_DIhIKrkJvhwbzbj3FKTyTzrtfda&export=download&authuser=0&confirm=t"

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_bundle")
        rootpath = ".."
    else:
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = snakemake.output.dir
    tyndp_bundle_fn = snakemake.params["tyndp_bundle"]
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
