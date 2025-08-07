# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script retrieves the TYNDP 2024 Scenarios Report Data Figures package for benchmarking purposes.

This rule downloads the package from the `ENTSOs website
<https://2024.entsos-tyndp-scenarios.eu/download/>` and extracts it in the
``data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package`` subdirectory.

**Output**

- ``data/tyndp_2024_bundle/TYNDP-2024-Scenarios-Package/TYNDP_2024-Scenario-Report-Data-Figures_240522.xslx``:
  Excel sheet with the data and calculation used to produce the figures in the Scenarios Report

"""

import logging
import os
import shutil
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
# TODO: retrieve_tyndp_benchmark needs to be deprecated once the package is added to the TYNDP data bundle
url_package = "https://2024-data.entsos-tyndp-scenarios.eu/files/reports/TYNDP-2024-Scenarios-Package-20250128.zip"


def retrieve_bundle(url: str, to_fn: str, disable_progress: bool = False):
    to_fn_zp = Path(to_fn, Path(url).name)

    # download .zip file
    logger.info(f"Downloading TYNDP 2024 Scenarios Report Data Figures from '{url}'.")
    progress_retrieve(url, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting TYNDP 2024 Scenarios Report Data Figures package.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    # remove .zip file and __MACOSX
    os.remove(to_fn_zp)
    shutil.rmtree(Path(to_fn, "__MACOSX"), ignore_errors=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_benchmark")

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = Path(snakemake.output.scenarios_figures).parents[1]

    retrieve_bundle(url_package, to_fn, disable_progress)

    logger.info(f"TYNDP 2024 Scenarios Report Data Figures available in '{to_fn}'.")
