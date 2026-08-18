# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Placeholder for capacity corrections aligning the SB reference grid with the
CBA reference grid of a given planning horizon.

SB and CBA are prepared at different stages, and the TYNDP 2026 reference
grid already gives the correct grid directly for every planning horizon (see
`build_tyndp_network.build_links`), so the SB base network fed into
`prepare_reference` is currently assumed to already match the target
horizon's reference grid: no correction is applied here for now, and this
script emits an empty corrections table.

The rule and script are kept as a placeholder rather than removed, since once
the CBA methodology is more mature it may need to patch the SB input network
independently (e.g. project-specific capacity overrides not reflected in the
raw reference grid). At that point, this script is where that correction
logic should live again.
"""

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

CORRECTIONS_COLUMNS = ["bus0", "bus1", "p_nom"]

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "fix_reference_sb_to_cba",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
            planning_horizons="2040",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pd.DataFrame(columns=CORRECTIONS_COLUMNS).to_csv(
        snakemake.output.corrections, quotechar="'"
    )
