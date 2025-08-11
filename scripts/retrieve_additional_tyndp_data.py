# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Retrieves additional TYNDP data not included in the Zenodo TYNDP data bundle .
Downloads the zip file and extracts it in the ``data/tyndp_2024_bundle``
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
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# TODO: retrieve_additional_tyndp_data needs to be deprecated once all TYNDP data is added to the TYNDP data bundle

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

    # define url
    url = snakemake.params.url

    # download .zip file
    logger.info(f"Downloading additional TYNDP data from '{url}'.")
    progress_retrieve(url, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting additional TYNDP data.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(tyndp_bundle_fn)

    # remove .zip file
    os.remove(to_fn_zp)

    logger.info(f"Additional TYNDP data available in '{to_fn}'.")
