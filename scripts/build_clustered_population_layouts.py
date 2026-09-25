# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build population layouts for all clustered model regions as total as well as
split by urban and rural population.
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, load_cutout, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_clustered_population_layouts", clusters=48)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    cutout = load_cutout(snakemake.input.cutout)

    clustered_regions = (
        gpd.read_file(snakemake.input.regions_onshore).set_index("name").buffer(0)
    )

    I = cutout.indicatormatrix(clustered_regions)  # noqa: E741

    pop = {}
    for item in ["total", "urban", "rural"]:
        pop_layout = xr.open_dataarray(snakemake.input[f"pop_layout_{item}"])
        pop[item] = I.dot(pop_layout.stack(spatial=("y", "x")))

    pop = pd.DataFrame(pop, index=clustered_regions.index)

    if snakemake.input.buses_tyndp:
        buses_tyndp = pd.read_csv(snakemake.input.buses_tyndp).set_index("bus_id")
        onshore_buses = buses_tyndp.index[buses_tyndp["category"] == "onshore"]
        pop = pop.reindex(onshore_buses, fill_value=0).sort_index()

    pop["ct"] = pop.index.str[:2]
    country_population = pop.total.groupby(pop.ct).sum()
    pop["fraction"] = (pop.total / pop.ct.map(country_population)).fillna(0)

    pop.to_csv(snakemake.output.clustered_pop_layout)
