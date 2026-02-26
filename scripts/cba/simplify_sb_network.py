# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Simplify scenario building network for CBA analysis.

Creates the base CBA network with:
- Fixed optimal capacities from scenario building
- Hurdle costs on DC links
- Extended primary fuel source capacities

**Inputs**

- Solved network from scenario building workflow

**Outputs**

- ``resources/cba/networks/simple_{planning_horizons}.nc``: Simplified network for CBA
"""

import logging

import pandas as pd
import pypsa
from numpy import inf

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def extend_primary_fuel_sources(n: pypsa.Network, tyndp_conventional_carriers: list):
    """
    Set infinite capacity for primary fuel source generators.

    Rolling horizon with fixed capacities can artificially limit peak fuel production.
    Primary fuel sources have no capital costs, so unlimited capacity ensures
    sufficient supply without affecting the objective function.

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
    None
        Network is modified in place.
    """
    # split for 'oil-light', 'oil-heavy', 'oil-shale' -> 'oil'
    primary_fuel_carriers = pd.Series(
        [carrier.split("-")[0] for carrier in tyndp_conventional_carriers]
    ).unique()
    mask = n.generators.carrier.str.contains("|".join(primary_fuel_carriers))
    gen_i = n.generators[mask].index
    n.generators.loc[gen_i, "p_nom"] = inf
    return n


def remove_noisy_marginal_costs(
    n: pypsa.Network,
    threshold: float = 0.02,
    components: list[str] | None = None,
):
    """
    Remove small noisy marginal costs for selected components (Storage, StorageUnit).

    This removes the random noise injected in solve_network.prepare_network
    (approx 0.01 +/- 0.001) by zeroing marginal_cost values below the threshold.
    """
    if components is None:
        components = ["Store", "StorageUnit"]
    for comp_name in components:
        comp = n.components.get(comp_name)
        if comp is None:
            continue
        if "marginal_cost" not in comp.df:
            continue
        mc = comp.df["marginal_cost"]
        cleaned = mc.where(mc.abs() >= threshold, 0.0)
        changed = (cleaned != mc).sum()
        if changed:
            logger.info(
                "Removed noisy marginal_costs for %s: %d entries set to 0 (threshold=%s)",
                comp_name,
                changed,
                threshold,
            )
        comp.df["marginal_cost"] = cleaned
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

        snakemake = mock_snakemake(
            "simplify_sb_network",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the solved network from scenario building
    n = pypsa.Network(snakemake.input.network)

    # Extend primary fuel sources capacity
    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers
    extend_primary_fuel_sources(n, tyndp_conventional_carriers)

    # TODO: in the case of a perfect foresight network we need to extract a single planning horizon here

    # Fix optimal capacities from scenario building
    n.optimize.fix_optimal_capacities()

    # Add hurdle costs to DC links
    # Hurdle costs: 0.01 €/MWh (p.20, 104 TYNDP 2024 CBA implementation guidelines)
    hurdle_costs = snakemake.params.hurdle_costs
    n.links.loc[n.links.carrier == "DC", "marginal_cost"] = hurdle_costs
    logger.info(f"Applied hurdle costs of {hurdle_costs} EUR/MWh to DC links")
    # TODO: for DE/GA add merging of the two H2 zones
    # TODO: for DE/GA add EV electricity consumption from SB as fixed demand

    # Save base network
    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers
    n = extend_primary_fuel_sources(n, tyndp_conventional_carriers)

    if snakemake.params.get("remove_noisy_costs", False):
        n = remove_noisy_marginal_costs(
            n,
            threshold=snakemake.params.get("noisy_costs_threshold", 0.02),
            components=["Store", "StorageUnit"],
        )

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
    logger.info(f"Saved CBA base network to {snakemake.output.network}")
