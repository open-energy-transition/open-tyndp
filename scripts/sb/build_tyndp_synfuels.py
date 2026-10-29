# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds the TYNDP 2026 domestic synthetic fuel production capacities and the
synthetic fuel generator properties for one planning horizon.

Inputs
------

- `data/tyndp_2026/.../Hydrogen/SYNTHETIC FUEL LINES 2026.xlsx`: maximum H2 to
  synthetic fuel capacities per H2 Zone 2 node (sheet ``H2 LIMITS``).
- `data/tyndp_2026/.../Hydrogen/SYNTHETIC FUEL GENERATORS PROPERTIES.xlsx`:
  import offer prices of the synthetic fuels and annual limit of the biogenic
  CO2 available for domestic synthetic fuel production (sheet ``TEMPLATE``).

Outputs
-------

- `resources/synfuel_links_tyndp_{planning_horizons}.csv`: production link
  capacities with columns ``bus0`` (H2 Zone 2 node), ``fuel`` (``e-liquids`` or
  ``sng``) and ``p_nom`` (MW_H2).
- `resources/synfuel_generators_tyndp_{planning_horizons}.csv`: generator
  properties indexed by generator, one per synthetic fuel import (e.g.
  ``e-liquids``, ``sng``) and ``CO2-Synfuels`` for the biogenic CO2, with
  columns ``marginal_cost`` (offer price, EUR/MWh_H2) and ``e_sum_max``
  (annual maximum, MWh_H2 for the fuels and t_CO2 for ``CO2-Synfuels``).
"""

import logging

import pandas as pd

from scripts._helpers import (
    configure_logging,
    format_bz_names,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_synfuel_links(
    fn: str, planning_horizon: int, fuels: list[str]
) -> pd.DataFrame:
    """
    Load the TYNDP H2 to synthetic fuel capacities of one planning horizon.

    Parameters
    ----------
    fn : str
        Path to the TYNDP synthetic fuel lines Excel file.
    planning_horizon : int
        Planning horizon to read the capacities for.
    fuels : list[str]
        Synthetic fuel names, e.g. ``["e-liquids", "sng"]``.

    Returns
    -------
    pd.DataFrame
        Capacities with columns ``bus0``, ``fuel`` and ``p_nom`` (MW_H2),
        indexed by link name.
    """
    column_dict = {
        "YEAR": "year",
        "NODE FROM": "bus0",
        "SYNTHETIC FUEL (E-FUEL / SNG)": "fuel",
        "H2 - SYNTHETIC FUEL MAX CAPACITY [MW_h2]": "p_nom",
    }

    links = (
        pd.read_excel(fn, sheet_name="H2 LIMITS", engine="calamine")
        .rename(columns=column_dict)
        .query("year == @planning_horizon and p_nom > 0")
        .assign(bus0=lambda df: format_bz_names(df.bus0))
    )

    if not (unknown := links.fuel[~links.fuel.isin(fuels)]).empty:
        raise ValueError(f"Unknown TYNDP synthetic fuels: {unknown.unique()}")

    links.index = links.bus0 + " H2 to " + links.fuel
    links.index.name = "name"

    return links[["bus0", "fuel", "p_nom"]]


def load_synfuel_generators(
    fn: str, planning_horizon: int, fuels: list[str]
) -> pd.DataFrame:
    """
    Load the TYNDP synthetic fuel generator properties of one planning horizon.

    Parameters
    ----------
    fn : str
        Path to the TYNDP synthetic fuel generators properties Excel file.
    planning_horizon : int
        Planning horizon to read the properties for.
    fuels : list[str]
        Synthetic fuels, e.g. ``["e-liquids", "sng"]``.

    Returns
    -------
    pd.DataFrame
        Offer price ``marginal_cost`` (EUR/MWh_H2) and annual maximum
        ``e_sum_max`` (MWh_H2, t_CO2 for ``CO2-Synfuels``), indexed by generator.
    """
    column_dict = {
        "YEAR": "year",
        "GENERATOR (CO2/E-FUELI IMPORTS/SNG IMPORTS)": "generator",
        "OFFER PRICE [€/MWh_h2]": "marginal_cost",
        "MAX ENERGY YEAR [GWh] (CO2: ktonnes)": "e_sum_max",
    }
    names = [*fuels, "CO2-Synfuels"]

    generators = (
        pd.read_excel(fn, sheet_name="TEMPLATE", engine="calamine")
        .rename(columns=column_dict)
        .query("year == @planning_horizon and generator in @names")
        .set_index("generator")[["marginal_cost", "e_sum_max"]]
        .assign(e_sum_max=lambda df: df.e_sum_max * 1e3)  # GWh to MWh, kt to t
    )

    return generators.loc[names]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_synfuels",
            planning_horizons=2030,
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.planning_horizons)
    fuels = list(snakemake.params.synfuels)

    links = load_synfuel_links(snakemake.input.synfuel_lines, planning_horizon, fuels)
    generators = load_synfuel_generators(
        snakemake.input.synfuel_generators, planning_horizon, fuels
    )

    links.to_csv(snakemake.output.synfuel_links)
    generators.to_csv(snakemake.output.synfuel_generators)
