# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds hydroelectric inflow time-series for each country based on TYNDP hydro inflow data.

Outputs
-------

- ``resources/profile_hydro_tyndp.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               bus, time,         Inflow to the state of charge (in MW),
                         year, hydro_tech   e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================
"""

import logging

import pandas as pd
import xarray as xr
from tqdm.contrib.itertools import product

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    safe_pyear,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_hydro_profile")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    years_in_time = pd.DatetimeIndex(time).year.unique()

    technologies = snakemake.params.technologies
    pyears = snakemake.params.planning_horizons

    inflows = []

    for year, technology in product(pyears, technologies):
        logger.info(
            f"Extracting hydro inflows for {technology} in {year}"
        )
        year_i = year
        # falling back to latest available pyear if not in list of available years
        year = safe_pyear(
            int(year),
            available_years=snakemake.params.available_years,
            source="PEMMDB hydro inflow",
        )

        inflow = (
            pd.read_csv(
                snakemake.input[f"hydro_inflow_tyndp_{technology}_{year}"],
                parse_dates=True,
                index_col=0,
            )
            .rename_axis("time")
            .reset_index()
            .melt(id_vars=["time"], var_name="bus", value_name="profile")
            .assign(year=year_i, hydro_tech=technology)
            .set_index(["bus", "time", "year", "hydro_tech"])
            .to_xarray()
        )

        inflows.append(inflow)

    ds = xr.merge(inflows)
    ds = ds.sel(time=time)

    ds.to_netcdf(snakemake.output.profile)
