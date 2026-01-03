# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.14230568.svg
  :target: https://doi.org/10.5281/zenodo.14230568

The data bundle contains input data for the 2024 TYNDP scenario building process.

This rule downloads the TYNDP data bundle from `zenodo
<https://doi.org/10.5281/zenodo.14230568>` and extracts it in the ``data``
subdirectory, such that all files of the bundle are stored in the
``data/tyndp_2024_bundle`` subdirectory.

**Outputs**

- ``data/tyndp_2024_bundle``: input data for TYNDP 2024 scenario building

"""

import logging
import zipfile

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_bundle")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # extract
    logger.info("Extracting TYNDP data bundle.")
    with zipfile.ZipFile(snakemake.input.zip_file, "r") as zip_ref:
        zip_ref.extractall(snakemake.output.dir)

    logger.info(f"TYNDP data bundle available in '{snakemake.output.dir}'.")
