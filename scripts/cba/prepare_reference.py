# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare CBA reference network with all TOOT projects included.

Placeholder: currently passes the simplified network through unchanged.
Once published TYNDP 2026 CBA results are available to build against, the
logic for adding TOOT/PINT projects and reference-grid corrections to the
SB base network belongs here.
"""

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_reference",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    # TODO: apply CBA project-addition / reference-grid corrections here
    # once published TYNDP 2026 CBA results are available to build against.

    n.export_to_netcdf(snakemake.output.network)
