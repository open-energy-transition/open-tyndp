# SPDX-FileCopyrightText: Contributors to NGV-IEM project
#
# SPDX-License-Identifier: MIT
"""
Extract results from uncertainty scenarios and consolidates them.

Consolidation yields single values for all interconnections, to be used for restricting flows exogenously.
"""

import importlib
import logging
import os
import re
import sys
from functools import partial
from typing import Any

import linopy
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import yaml
from linopy.remote.oetc import OetcCredentials, OetcHandler, OetcSettings
from pypsa.descriptors import get_activity_mask
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "extract_uncertainty_results",
            clusters=all,
            sector_opts="",
            planning_horizons=2030,
            configfiles="config/test/config.tyndp.yaml",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    dispatches = []
    for n_fp in snakemake.input.networks:
        # n = pypsa.Network(n_fp)
        n = pypsa.Network(
            "/home/user/Documents/GitHub/NGV-IEM/results/tyndp/NT/low-demand_high-renewables/networks/base_s_all___2030.nc"
        )

        # Extract all relevant links (DC links from or to GB00)
        links_s = n.components.links.static
        relevant_links = links_s.loc[
            (links_s["carrier"] == "DC")
            & ((links_s["bus0"] == "GB00") | (links_s["bus1"] == "GB00"))
        ].index

        # Get dispatches for relevant links
        links_d = n.components.links.dynamic
        links_d = links_d["p0"][relevant_links]

        dispatches.append(links_d)

    # Consolidate the results of the different scenarios by averaging
    combined = pd.concat(dispatches, axis=1)
    consolidated = combined.T.groupby(level=0).mean().T

    # Save to CSV
    consolidated.to_csv(snakemake.output["line_limits"])
