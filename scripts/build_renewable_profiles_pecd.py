# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Create renewable profiles for each region from PECD and TYNDP datasets, building on PECD climate data. The available
generation time series are read for each node across all renewable technologies, including onshore wind, AC-connected offshore wind,
DC-connected offshore wind, and solar PV generators.

.. note:: Hydroelectric profiles will be built in script :mod:`build_hydro_profiles_PECD`. Not yet implemented.

Outputs
-------

- ``resources/profile_pecd_{clusters}_{technology}.nc`` with the following structure

    ===================  ====================  =========================================================
    Field                Dimensions            Description
    ===================  ====================  =========================================================
    profile              year, bus, bin, time  the per unit hourly availability factors for each bus
    ===================  ====================  =========================================================
"""

import logging

import pandas as pd
import xarray as xr

from scripts._helpers import (
    configure_logging,
    safe_pyear,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_renewable_profiles_pecd",
            clusters="all",
            technology="Wind_Offshore",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    technology = snakemake.wildcards.technology
    pyears = snakemake.params.planning_horizons

    profiles = []

    for year in pyears:
        logger.info(
            f"Extract PECD capacity factor time series for year {year} for technology {technology}..."
        )
        year_i = year
        # falling back to latest available pyear if not in list of available years
        year = safe_pyear(
            int(year), available_years=snakemake.params.available_years, source="PECD"
        )
        # TODO: remove once PECD data is updated
        if year == 2050:
            logger.warning(
                "PECD input data for 2050 is incomplete. Falling back to 2040 PECD data."
            )
            year = 2040

        profile = (
            pd.read_csv(
                snakemake.input[f"pecd_data_{year}"], parse_dates=True, index_col=0
            )
            .rename_axis("time")
            .reset_index()
            .melt(id_vars=["time"], var_name="bus", value_name="profile")
            .assign(bin=0, year=year_i)
            .set_index(["time", "bus", "bin", "year"])
            .to_xarray()
        )

        profiles.append(profile)

    ds = xr.merge(profiles)
    ds.to_netcdf(snakemake.output.profile)
