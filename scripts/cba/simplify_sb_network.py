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


def extend_primary_fuel_sources(n: pypsa.Network, tyndp_conventional_carriers: list):
    """
    Remove capacity constraints on primary fuel source generators for CBA rolling horizon.

    When using capacities fixed from Scenario Building in rolling horizon optimization,
    peak fuel production can be artificially limited. Since primary fuel sources incur no
    capital costs, unlimited capacity ensures sufficient fuel supply without changing the
    objective function.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    tyndp_conventional_carriers : list
        List of conventional carrier names from TYNDP data, which may include
        fuel sub-types (e.g., 'oil-light', 'oil-heavy'). These are grouped by
        their base fuel type (e.g., 'oil').

    Returns
    -------
    pypsa.Network
        Modified network with infinite capacity for primary fuel sources
    """
    # split for 'oil-light', 'oil-heavy', 'oil-shale' -> 'oil'
    primary_fuel_carriers = pd.Series(
        [carrier.split("-")[0] for carrier in tyndp_conventional_carriers]
    ).unique()
    gen_i = n.generators[n.generators.carrier.isin(primary_fuel_carriers)]
    n.generators.loc[gen_i, "p_nom"] = inf
    return n


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


def disable_store_cyclicity(n: pypsa.Network):
    """
    Disable cyclic state of charge constraints for stores and storage units.

    Cyclic constraints force storage to end at the same state of charge as it started.
    For CBA storage should be free to have different start/end states.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    # Disable cyclicity for stores
    has_e_cyclic = n.stores["e_cyclic"]
    if has_e_cyclic.any():
        stats = summarize_counts(n.stores.loc[has_e_cyclic, "carrier"])
        logger.info(f"Disabling e_cyclic for stores:\n{stats}")
        n.stores["e_cyclic"] = False

    has_e_cyclic_per_period = n.stores["e_cyclic_per_period"]
    if has_e_cyclic_per_period.any():
        stats = summarize_counts(n.stores.loc[has_e_cyclic_per_period, "carrier"])
        logger.info(f"Disabling e_cyclic_per_period for stores:\n{stats}")
        n.stores["e_cyclic_per_period"] = False

    # Disable cyclicity for storage units
    has_cyclic_soc = n.storage_units["cyclic_state_of_charge"]
    if has_cyclic_soc.any():
        stats = summarize_counts(n.storage_units.loc[has_cyclic_soc, "carrier"])
        logger.info(f"Disabling cyclic_state_of_charge for storage units:\n{stats}")
        n.storage_units["cyclic_state_of_charge"] = False


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

    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers
    n = extend_primary_fuel_sources(n, tyndp_conventional_carriers)

    disable_volume_limits(n)

    disable_store_cyclicity(n)

    disable_global_constraints(n)

    # Hurdle costs
    n.links.loc[n.links.carrier == "DC", "marginal_cost"] = (
        snakemake.params.hurdle_costs
    )

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
