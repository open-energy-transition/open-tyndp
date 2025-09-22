# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Collects and bundles the available PEMMDB v2.4 capacities and profiles for different PEMMDB technologies from TYNDP data bundle for a given planning horizon.

Outputs
-------
Cleaned csv file with all NT capacities (p_nom) in long format and netcdf file with must run obligations (p_min_pu) and availability (p_max_pu) for all different PEMMDB technologies.

- ``resources/pemdb_capacties_{planning_horizon}.csv`` in long format
- ``resources/pemmdb_profiles_{planning_horizon}.nc`` with the following structure:

    ===================  ====================  =========================================================
    Field                Dimensions            Description
    ===================  ====================  =========================================================
    profile              time, bus, carrier,   the per unit hourly availability and must-run obligations
                         type                  for each bus and PEMMDB technology
    ===================  ====================  =========================================================
"""

import logging

import pandas as pd
import xarray as xr
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_pemmdb_data",
            clusters="all",
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameter
    pemmdb_techs = snakemake.params.pemmdb_techs

    # Capacities
    capacities = []
    logger.info("Load PEMMDB capacities for ...")
    for tech in tqdm(pemmdb_techs):
        capacity = pd.read_csv(
            snakemake.input[f"pemmdb_capacities_{tech}"], index_col=0
        )
        capacities.append(capacity)

    capacities_df = pd.concat(capacities, axis=0).reset_index(drop=True)
    capacities_df.to_csv(snakemake.output.pemmdb_capacities, index=False)

    # Profiles
    profiles = []
    logger.info(
        f"Load PEMMDB must-run and availability profiles for {', '.join(pemmdb_techs)}..."
    )
    for tech in tqdm(pemmdb_techs):
        profile = xr.open_dataset(snakemake.input[f"pemmdb_profiles_{tech}"])
        profiles.append(profile)

    ds = xr.merge(profiles)
    ds.to_netcdf(snakemake.output.pemmdb_profiles)
