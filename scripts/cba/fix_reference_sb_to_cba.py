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

The rule and script are kept as a placeholder rather than removed: once
published TYNDP 2026 CBA results are available to build against, patching the
SB input network independently (e.g. project-specific capacity overrides not
reflected in the raw reference grid) may again be needed, and this script is
where that correction logic should live.
"""

import pandas as pd

<<<<<<< HEAD
from scripts._helpers import configure_logging, set_scenario_config
=======
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)
from scripts.build_tyndp_network import add_links_missing_attributes
from scripts.sb._helpers import extract_grid_data_tyndp
from scripts.sb.build_tyndp_transmission_projects import read_invest_file


def read_guidelines(path: Path) -> pd.DataFrame:
    query = "In_ref_grid_2030=='no' and In_ref_grid_2040=='yes'"
    df = (
        pd.read_csv(path)
        .query(query)
        .rename(
            columns={
                "A_B_MW": "Summary Direction 1",
                "B_A_MW": "Summary Direction 2",
            }
        )
        .replace(
            {
                "ITS1-GR00": "GR00-ITS1",
                "MT00-ITSI": "ITSI-MT00",
                "UKNI-UK00": "UK00-UKNI",
            }
        )
        .replace({"UK": "GB"}, regex=True)
        .set_index("Border")[["Summary Direction 1", "Summary Direction 2"]]
        .dropna(how="all")
        .astype({"Summary Direction 1": float, "Summary Direction 2": float})
        .groupby("Border")
        .sum()
    )
    return df[~df.index.str.contains("internal|Offshore|OBZ")]


def build_links(projects_fix: pd.DataFrame, buses_fn: str) -> pd.DataFrame:
    links_fix = extract_grid_data_tyndp(
        projects_fix.reset_index(),
        expand_from_index=True,
        idx_connector="-",
        idx_suffix="-DC",
        idx_separator="",
    )
    buses = gpd.read_file(buses_fn).set_index("bus_id")
    links_fix = add_links_missing_attributes(links_fix, buses)

    links_fix["carrier"] = "DC"
    links_fix["length"] /= 1e3
    links_fix[["underground", "under_construction"]] = (
        links_fix[["underground", "under_construction"]] == "t"
    )

    return links_fix
>>>>>>> origin/master

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

# TODO: Adjust or deprecate reference fixes once once published
# TYNDP 2026 CBA results are available to build against.
pd.DataFrame(columns=CORRECTIONS_COLUMNS).to_csv(
    snakemake.output.corrections, quotechar="'"
)
