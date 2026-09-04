# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script extracts the TYNDP 2026 electricity reference grid's NTC (link
``p_nom``, plus ``bus0``/``bus1``/``length`` for links not yet in the network)
for a given wildcard planning horizon.

``build_tyndp_network.py`` builds the static electricity base network topology
once, from a single, fixed reference year (``electricity: tyndp_reference_year``)
- the border topology is identical across all available planning horizons in
the reference grid, only the NTC capacities differ. This script re-reads the
same reference grid's ``Year_{planning_horizons}`` sheet to extract the
per-horizon NTC, later overlaid onto the network in `prepare_sector_network`
(see `apply_tyndp_electricity_ntc`) whenever `electricity: base_network` is
`tyndp`.
"""

import logging

import geopandas as gpd
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_tyndp_network import build_links

logger = logging.getLogger(__name__)

NTC_COLUMNS = ["bus0", "bus1", "p_nom", "length"]

# Matches the constants `build_tyndp_network.add_links_missing_attributes`
# hardcodes for every TYNDP electricity link.
VOLTAGE = 380  # TODO Improve assumption, see build_tyndp_network.py
UNDERGROUND = True
UNDER_CONSTRUCTION = False


def apply_tyndp_electricity_ntc(
    n: pypsa.Network,
    ntc_fn: str,
    investment_year: int,
) -> None:
    """
    Overlay TYNDP electricity reference-grid NTC onto the network for a given
    planning horizon.

    Existing links (already part of the network's fixed-reference-year
    topology, see `build_tyndp_network.build_links`) get their `p_nom`
    updated to this horizon's NTC. Links that only appear in this horizon's
    reference grid (e.g. interconnectors commissioned after the reference
    year) are added as new links.

    Parameters
    ----------
    n : pypsa.Network
        Network to overlay the NTC onto, modified in place.
    ntc_fn : str
        Path to the per-horizon NTC CSV produced by this script's `__main__`.
    investment_year : int
        Planning horizon the NTC data applies to, used for logging only.
    """
    ntc = pd.read_csv(ntc_fn, index_col=0)

    known_links = n.links.index.intersection(ntc.index)
    n.links.loc[known_links, "p_nom"] = ntc.loc[known_links, "p_nom"]

    new_links = ntc.index.difference(n.links.index)
    if new_links.empty:
        return

    logger.info(
        f"Adding new TYNDP electricity links for planning horizon {investment_year}: "
        f"{', '.join(sorted(new_links))}"
    )
    p_max_pu = n.links["p_max_pu"].iloc[0] if not n.links.empty else 1.0
    p_min_pu = n.links["p_min_pu"].iloc[0] if not n.links.empty else -p_max_pu
    bus0 = ntc.loc[new_links, "bus0"]
    bus1 = ntc.loc[new_links, "bus1"]
    n.add(
        "Link",
        new_links,
        bus0=bus0,
        bus1=bus1,
        voltage=VOLTAGE,
        p_nom=ntc.loc[new_links, "p_nom"],
        length=ntc.loc[new_links, "length"] / 1e3,  # m -> km
        underground=UNDERGROUND,
        under_construction=UNDER_CONSTRUCTION,
        tags=bus0 + " -> " + bus1,
        carrier="DC",
        dc=True,
        p_max_pu=p_max_pu,
        p_min_pu=p_min_pu,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_electricity_ntc",
            planning_horizons=2030,
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pyear = int(snakemake.wildcards.planning_horizons)

    buses = gpd.read_file(snakemake.input.buses).set_index("bus_id")
    links = build_links(
        snakemake.input.elec_reference_grid, buses, reference_year=pyear
    )

    links[NTC_COLUMNS].to_csv(snakemake.output.ntc)
