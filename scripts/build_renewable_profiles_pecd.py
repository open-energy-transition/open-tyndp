# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Creates renewable profiles for each region from PECD and TYNDP input data containing the available
generation time series (based on PECD weather data) from the node for onshore wind, AC-connected offshore wind,
DC-connected offshore wind and solar PV generators.

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

import numpy as np
import pandas as pd
import xarray as xr

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_renewable_profiles_pecd",
            clusters="all",
            technology="offwind-ac",
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
        if int(year) not in [2030, 2040, 2050]:
            year = np.clip(10 * (year // 10), 2030, 2050)
            logger.warning(
                "Planning horizon doesn't match available TYNDP PECD data. "
                f"Falling back to previous available year {year}."
            )
        if year == 2050:
            logger.warning(
                "PECD input data for 2050 is incomplete. Falling back to 2040 PECD data instead."
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

    # TODO: Later on the max capacities for renewable technologies and regions can be added here
    # profiles = xr.merge(profiles)
    # p_nom_max = place_holder
    # average_distance = place_holder
    # ds = xr.merge(
    #     [
    #         profiles,
    #         p_nom_max.rename("p_nom_max"),
    #         average_distance.rename("average_distance"),
    #     ]
    # )

    ds.to_netcdf(snakemake.output.profile)
