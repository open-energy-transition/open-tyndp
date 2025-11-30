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

import pypsa

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def turn_optimisation_to_dispatch(n: pypsa.Network) -> pypsa.Network:
    """Disable all capacity expansion by copying optimized capacities to nominal capacities"""

    logger.info("Turning optimisation problem into dispatch problem.")

    for c in n.components[["Generator", "Link", "Store", "StorageUnit", "Line"]]:
        ext_flag = [col for col in c.static.columns if "extendable" in col]
        if len(ext_flag) != 1:
            logger.error("Expected exactly one extendable flag column")
        ext_flag = ext_flag[0]
        capacity = ext_flag.replace("_extendable", "")
        capacity_extendable = f"{capacity}_extendable"

        # Set all components to not be extendable
        c.static[capacity_extendable] = False

        # Copy the optimized capacities to the nominal capacity columns
        c.static[capacity] = c.static[f"{capacity}_opt"]

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


def add_scenario_uncertainty(n: pypsa.Network, scenario_name: str) -> pypsa.Network:
    """Creates a new network that is modified according to the uncertainty scenario."""

    # Work on a copy of the network
    n = n.copy()

    # Change name
    n.name = f"{n.name} (sensitivity: {scenario_name})"

    # Modify demand and renewable availability according to the scenario
    all_factors = {
        "low-demand_low-renewables": (0.995, 1.0),
        "low-demand_high-renewables": (0.995, 1.005),
        "high-demand_low-renewables": (1.005, 0.995),
        "high-demand_high-renewables": (1.005, 1.005),
    }

    scenario_factors = all_factors[scenario_name]

    # Get all RES generators for the scenario and modify their availability (p_max_pu is upper bound for available dispatch)
    res_factor = scenario_factors[1]
    mask = n.generators.loc[
        n.generators["carrier"].str.startswith(("ror", "onwind", "solar-pv", "offwind"))
    ].index
    n.generators_t["p_max_pu"][mask] *= res_factor

    # Get all electricity loads for the scenario and modify their demand
    demand_factor = scenario_factors[0]
    mask = n.loads.loc[n.loads["carrier"] == "electricity"].index
    n.loads_t["p_set"][mask] *= demand_factor

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_uncertainty_scenarios",
            opts="",
            clusters="5",
            configfiles="config/test/config.overnight.yaml",
            sector_opts="",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)

    n = turn_optimisation_to_dispatch(n)
    n = remove_components_added_in_solve_network_py(n)

    for uncertainty_scenario in snakemake.params.uncertainty_scenarios:
        n_uncertain = add_scenario_uncertainty(n, uncertainty_scenario)

        output_path = [
            p for p in snakemake.output.networks if f"{uncertainty_scenario}" in p
        ]
        assert len(output_path) == 1, (
            "Expected exactly one output path for each uncertainty scenario"
        )
        output_path = output_path[0]

        logger.info(f"Saving uncertainty scenario network to {output_path}")
        n_uncertain.export_to_netcdf(output_path)
