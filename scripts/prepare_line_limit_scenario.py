# SPDX-FileCopyrightText: Contributors to NGV-IEM project
#
# SPDX-License-Identifier: MIT
"""
Uses an already solved network and prepares it for uncertainty analysis by:
* Copying optimised capacities (`p_nom_opt`, `e_nom_opt`, ...) to nominal capacities (`p_nom`, `e_nom`, ...)
* Setting capacity extendable flags to False
* Modifying the network according to the uncertainty scenario, e.g., changing the demand or availability of renewables
"""

import logging

import pandas as pd
import pypsa
import numpy as np

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

def add_electrolysis_constraints(n):
    electrolysis_i = n.links[n.links.carrier=="H2 Electrolysis"].index
    n.links_t.p_set.loc[:,electrolysis_i] = abs(n.links_t.p0.loc[:,electrolysis_i])
    return n 

def remove_components_added_in_solve_network_py(n: pypsa.Network) -> pypsa.Network:
    """Removes components that were added in solve_network.py; we're planing on running this network through the same step again and want to avoid adding the components again."""

    logger.info("Removing components added in solve_network.py")

    # These components are not always part of the network, so
    # we check for their existence first
    if "co2_sequestration_limit" in n.global_constraints.index:
        n.remove(
            class_name="GlobalConstraint",
            name="co2_sequestration_limit",
        )

    if "load" in n.carriers.index:
        n.remove(
            class_name="Carrier",
            name="load",
        )
        gens_i = n.generators.query("`name`.str.endswith(' load')").index
        n.remove(
            class_name="Generator",
            name=gens_i,
        )

    if "curtailment" in n.carriers.index:
        n.remove(
            class_name="Carrier",
            name="curtailment",
        )
        gens_i = n.generators.query("`name`.str.endswith(' curtailment')").index
        n.remove(
            class_name="Generator",
            name=gens_i,
        )

    return n

def restrict_elec_flows(n: pypsa.Network, line_limits_fp: str) -> pypsa.Network:
    """
    Restrict electricity flows based on pre-calculated hourly line limits from and to GB.

    Restrictions are put in place by limiting `p_min_pu` and `p_max_pu` of each line connected to GB.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network instance
    line_limits_fp : str
        File path to CSV containing line limits

    Returns
    -------
    pypsa.Network
        PyPSA network instance with restricted line flows
    """
    logger.info(
        "Restricting electricity flows based on line limits from uncertainty scenarios."
    )
    line_limits = pd.read_csv(line_limits_fp, index_col=0, parse_dates=True)
    line_p_max_pu = 1.05 * n.components.links.dynamic["p_max_pu"]
    line_p_min_pu = 0.95 * n.components.links.dynamic["p_max_pu"]

    # Ensure that all lines for which line limits are provided exist in the network
    # (If not, then we are using the wrong input either for the network or the line limits)
    missing_lines = line_limits.columns.difference(n.components.links.static.index)
    if not missing_lines.empty:
        raise ValueError(
            f"The following lines from the line limits file are missing in the network: {missing_lines.tolist()}"
        )

    # Remove existing restrictions that are also part of the `line_limits` if there are any
    # This is not problematic, as the new restrictions are build upon the old restrictions,
    # i.e. the most restrictive limits will apply
    existing_restricted_links = line_limits.columns.intersection(line_p_max_pu.columns)
    if any(existing_restricted_links):
        logger.info(
            f"Removing existing link flow restrictions for GB-connected lines: {existing_restricted_links.tolist()}"
        )
        line_p_max_pu = line_p_max_pu.drop(columns=existing_restricted_links)
        line_p_min_pu = line_p_min_pu.drop(columns=existing_restricted_links)

    # Add new restrictions
    n.components.links.dynamic["p_max_pu"] = pd.concat(
        [line_p_max_pu, line_limits], axis="columns"
    )
    n.components.links.dynamic["p_min_pu"] = pd.concat(
        [line_p_min_pu, line_limits], axis="columns"
    )

    return n

#%%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_line_limit_scenarios",
            opts="",
            clusters="all",
            configfiles="config/config.tyndp.yaml",
            sector_opts="",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    
    n.optimize.fix_optimal_capacities()               
    n = remove_components_added_in_solve_network_py(n)
    n = add_electrolysis_constraints(n)

    n = restrict_elec_flows(n, snakemake.input.line_limits)

    n.export_to_netcdf(snakemake.output.network)
    