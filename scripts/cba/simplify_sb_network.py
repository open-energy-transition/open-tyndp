# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Simplify CBA base network for rolling horizon dispatch.

Takes the CBA base network (with fixed capacities and hurdle costs) and applies
simplifications needed for CBA rolling horizon optimization:
- Extend primary fuel source capacities
- Disable volume limits
- Disable cyclicity for seasonal stores (they will receive MSV instead)
- Keep cyclicity for short-term storage (batteries)
- Disable global constraints
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
    mask = n.generators.carrier.str.contains("|".join(primary_fuel_carriers))
    gen_i = n.generators[mask].index
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

    Removes constraints on biomass sustainability to allow unconstrained
    operation in the reference network.

    Note: co2_sequestration_limit is already removed in prepare_cba_base.py

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    if "unsustainable biomass limit" in n.global_constraints.index:
        logger.info('Removing GlobalConstraint "unustainable biomass limit"')
        n.remove("GlobalConstraint", "unsustainable biomass limit")


def disable_store_cyclicity(
    n: pypsa.Network,
    cyclic_carriers: list[str] | None = None,
    seasonal_carriers: list[str] | None = None,
):
    """
    Disable cyclic state of charge constraints for stores and storage units.

    Cyclic constraints force storage to end at the same state of charge as it started.
    For CBA rolling horizon, seasonal storage should be non-cyclic (they receive MSV instead),
    while short-term storage (e.g., batteries) should remain cyclic.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    cyclic_carriers : list[str], optional
        List of carrier names that should remain cyclic (e.g., ["battery", "home battery"]).
        If None, all storage cyclicity is disabled.
    seasonal_carriers : list[str], optional
        List of carrier names that are seasonal (e.g., ["H2 Store", "gas"]).
        These will have cyclicity disabled and will receive MSV during solving.
        Used for logging purposes.
    """
    if cyclic_carriers is None:
        cyclic_carriers = []
    if seasonal_carriers is None:
        seasonal_carriers = []

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

    # Log seasonal carriers that will receive MSV
    if seasonal_carriers:
        is_seasonal = n.stores["carrier"].isin(seasonal_carriers)
        if is_seasonal.any():
            stats = summarize_counts(n.stores.loc[is_seasonal, "carrier"])
            logger.info(f"Seasonal stores (will receive MSV):\n{stats}")

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

    # Load the CBA base network (already has fixed capacities and hurdle costs)
    n = pypsa.Network(snakemake.input.network)

    # Extend primary fuel sources capacity
    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers
    n = extend_primary_fuel_sources(n, tyndp_conventional_carriers)

    disable_volume_limits(n)

    # Get storage carrier settings from config
    cyclic_carriers = snakemake.params.get("cyclic_carriers", [])
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    disable_store_cyclicity(
        n, cyclic_carriers=cyclic_carriers, seasonal_carriers=seasonal_carriers
    )

    disable_global_constraints(n)

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
