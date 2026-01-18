# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Simplify SB (Scenario Building) network for CBA analysis.

Extracts a planning horizon from the optimized network and applies
simplifications needed for CBA reference network.
"""

import logging

import pandas as pd
import pypsa
from numpy import inf, isfinite

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def summarize_counts(s: pd.Series):
    ret = ""
    for key, count in s.value_counts().items():
        ret += f"- {key} ({count})\n"
    return ret


def disable_volume_limits(n: pypsa.Network):
    """
    Disable volume limits (e_sum_min) for generators and links.

    Volume limits constrain minimum energy production over the optimization period.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    for c in n.components[{"Generator", "Link"}]:
        has_e_sum_min = isfinite(c.static.get("e_sum_min", []))
        if has_e_sum_min.any():
            stats = summarize_counts(c.static.loc[has_e_sum_min, "carrier"])
            logger.info(f"Disabling e_sum_min volume limits of:\n{stats}")
            c.static.loc[has_e_sum_min, "e_sum_min"] = -inf


def disable_global_constraints(n: pypsa.Network):
    """
    Remove global constraints that are not needed for CBA analysis.

    Removes constraints on biomass sustainability and CO2 sequestration limits
    to allow unconstrained operation in the reference network.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    if "unsustainable biomass limit" in n.global_constraints.index:
        logger.info('Removing GlobalConstraint "unustainable biomass limit"')
        n.remove("GlobalConstraint", "unsustainable biomass limit")
    if "co2_sequestration_limit" in n.global_constraints.index:
        logger.info('Removing GlobalConstraint "co2_sequestration_limit"')
        n.remove("GlobalConstraint", "co2_sequestration_limit")


def disable_store_cyclicity(n: pypsa.Network, cyclic_carriers: list[str] | None = None):
    """
    Disable cyclic state of charge constraints for stores and storage units.

    Cyclic constraints force storage to end at the same state of charge as it started.
    For CBA, long-term storage should be free to have different start/end states,
    while short-term storage (e.g., batteries) should remain cyclic.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    cyclic_carriers : list[str], optional
        List of carrier names that should remain cyclic (e.g., ["battery", "home battery"]).
        If None, all storage cyclicity is disabled.
    """
    if cyclic_carriers is None:
        cyclic_carriers = []

    # Disable cyclicity for stores (except cyclic_carriers)
    has_e_cyclic = n.stores["e_cyclic"]
    is_cyclic_carrier = n.stores["carrier"].isin(cyclic_carriers)
    to_disable = has_e_cyclic & ~is_cyclic_carrier

    if to_disable.any():
        stats = summarize_counts(n.stores.loc[to_disable, "carrier"])
        logger.info(f"Disabling e_cyclic for stores:\n{stats}")
        n.stores.loc[to_disable, "e_cyclic"] = False

    if is_cyclic_carrier.any() and has_e_cyclic.any():
        kept_cyclic = has_e_cyclic & is_cyclic_carrier
        if kept_cyclic.any():
            stats = summarize_counts(n.stores.loc[kept_cyclic, "carrier"])
            logger.info(f"Keeping e_cyclic=True for short-term storage:\n{stats}")

    has_e_cyclic_per_period = n.stores["e_cyclic_per_period"]
    to_disable_per_period = has_e_cyclic_per_period & ~is_cyclic_carrier

    if to_disable_per_period.any():
        stats = summarize_counts(n.stores.loc[to_disable_per_period, "carrier"])
        logger.info(f"Disabling e_cyclic_per_period for stores:\n{stats}")
        n.stores.loc[to_disable_per_period, "e_cyclic_per_period"] = False

    # Disable cyclicity for storage units (except cyclic_carriers)
    has_cyclic_soc = n.storage_units["cyclic_state_of_charge"]
    is_cyclic_carrier_su = n.storage_units["carrier"].isin(cyclic_carriers)
    to_disable_su = has_cyclic_soc & ~is_cyclic_carrier_su

    if to_disable_su.any():
        stats = summarize_counts(n.storage_units.loc[to_disable_su, "carrier"])
        logger.info(f"Disabling cyclic_state_of_charge for storage units:\n{stats}")
        n.storage_units.loc[to_disable_su, "cyclic_state_of_charge"] = False

    if is_cyclic_carrier_su.any() and has_cyclic_soc.any():
        kept_cyclic_su = has_cyclic_soc & is_cyclic_carrier_su
        if kept_cyclic_su.any():
            stats = summarize_counts(n.storage_units.loc[kept_cyclic_su, "carrier"])
            logger.info(
                f"Keeping cyclic_state_of_charge=True for short-term storage:\n{stats}"
            )


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

    n.optimize.fix_optimal_capacities()

    disable_volume_limits(n)

    # Get cyclic carriers from config (short-term storage that should remain cyclic)
    cyclic_carriers = snakemake.params.get("cyclic_carriers", [])
    disable_store_cyclicity(n, cyclic_carriers=cyclic_carriers)

    disable_global_constraints(n)

    # Hurdle costs
    n.links.loc[n.links.carrier == "DC", "marginal_cost"] = (
        snakemake.params.hurdle_costs
    )

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
