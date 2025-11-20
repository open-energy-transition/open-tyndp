# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Simplify SB (Scenario Building) network for CBA analysis.

Extracts a planning horizon from the optimized network and applies
simplifications needed for CBA reference network.
"""

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("simplify_sb_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the solved base network
    # The wildcard has been expanded to {clusters}_{opts}_{sector_opts}_{planning_horizons}
    n = pypsa.Network(snakemake.input.network)

    # TODO: in the case of a perfect foresight network we need to extract a single planning horizon here

    # TODO: Apply CBA-specific simplifications
    # For now, pass through the network as-is

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
