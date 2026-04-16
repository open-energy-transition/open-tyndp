# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Prepare the scenario-building network input used for the CBA workflow.

The CBA workflow reuses:
- the solved SB network from the base ``cba.sb_scenario`` as the source of
  optimized capacities and static topology, and
- a climate-year-specific pre-solve network as the source of snapshots and
  weather-dependent input data.

For the first planning horizon, the climate-year source is the output of
``add_existing_baseyear``.

For later planning horizons, the climate-year source is the output of
``prepare_sector_network``. Those later-horizon source networks use the
pre-brownfield naming convention, so this script aligns new assets with the
solved myopic naming/build-year convention before copying climate-dependent
inputs onto the reused solved base network.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


RESULT_KEYS = {
    "e",
    "marginal_price",
    "mu_energy_balance",
    "n_mod",
    "p",
    "p0",
    "p1",
    "p2",
    "p3",
    "q",
    "spill",
    "state_of_charge",
    "status",
}


def add_build_year_to_new_assets(
    n: pypsa.Network,
    baseyear: int,
    components: tuple[str, ...] = ("links", "generators", "stores"),
) -> None:
    """
    Add build-year suffixes to new assets in a pre-solved source network.

    This is an adaptation of ``add_existing_baseyear.add_build_year_to_new_assets.``:
    - assets with finite lifetime and build_year == 0 get the target year
    - asset names get suffixed with -{baseyear}
    - dynamic series columns get renamed to match

    ``prepare_sector_network`` outputs later-horizon assets without the build-year
    suffix that appears in solved myopic networks. To align the climate-year
    source network with the solved base network, we reproduce the relevant part
    of the SB naming logic here for the affected component types.
    """
    for component_name in components:
        static = getattr(n, component_name)
        dynamic = getattr(n, f"{component_name}_t")
        if static.empty:
            continue

        assets = static.index[(static.lifetime != np.inf) & (static.build_year == 0)]
        if assets.empty:
            continue

        static.loc[assets, "build_year"] = baseyear
        rename = {asset: f"{asset}-{baseyear}" for asset in assets}
        static.rename(index=rename, inplace=True)

        for key in list(dynamic.keys()):
            if not dynamic[key].empty:
                dynamic[key] = dynamic[key].rename(columns=rename)

        logger.info(
            "Added build-year suffix %s to %s %s assets",
            baseyear,
            len(assets),
            component_name,
        )


def overwrite_static_attr(
    base: pypsa.Network, source: pypsa.Network, component: str, attr: str
) -> None:
    """
    Overwrite a static component attribute with explicit intersection checks.

    This only overwrites attributes that are intended to be climate-year-dependent.
    The helper logs any non-shared rows so the source/base mismatch stays
    visible, but it only touches the shared components explicitly.
    """
    base_static = getattr(base, component)
    source_static = getattr(source, component)

    shared = base_static.index.intersection(source_static.index)
    if shared.empty:
        raise ValueError(
            f"Cannot overwrite {component}.{attr}: no shared component rows found."
        )

    base_static.loc[shared, attr] = source_static.loc[shared, attr].values
    missing = base_static.index.difference(source_static.index)
    extra = source_static.index.difference(base_static.index)
    logger.info(
        "Overwrote static %s.%s for %s shared rows (missing=%s, extra=%s)",
        component,
        attr,
        len(shared),
        len(missing),
        len(extra),
    )


def overwrite_dynamic_attr(
    base: pypsa.Network, source: pypsa.Network, component: str, attr: str
) -> None:
    """
    Overwrite a dynamic component attribute with explicit snapshot/column checks.

    This borrows the alignment idea from PyPSA's internal dynamic-data helpers, but slightly adapted for our use case:
    - snapshots must match the already-prepared base network snapshots
    - overwrite only shared columns
    - keep base-only columns untouched
    - no silent zero-filling
    """
    base_dynamic = getattr(base, f"{component}_t")
    source_dynamic = getattr(source, f"{component}_t")

    df = source_dynamic[attr].copy()

    if not df.index.equals(base.snapshots):
        missing = base.snapshots.difference(df.index)
        extra = df.index.difference(base.snapshots)
        raise ValueError(
            f"Cannot overwrite {component}_t.{attr}: snapshot mismatch "
            f"(missing={len(missing)}, extra={len(extra)})."
        )

    expected_columns = getattr(base, component).index
    shared = expected_columns.intersection(df.columns)
    if shared.empty:
        raise ValueError(
            f"Cannot overwrite {component}_t.{attr}: no shared columns found."
        )

    # Start from the base network's existing data for attributes/columns that do
    # not exist in the climate-year source (for example, base-only load shedding
    # assets). Then overwrite the shared columns explicitly from the source.
    result = base_dynamic[attr].copy()
    result.loc[:, shared] = df.loc[:, shared].copy()
    base_dynamic[attr] = result

    missing = expected_columns.difference(df.columns)
    extra = df.columns.difference(expected_columns)
    logger.info(
        "Overwrote dynamic %s_t.%s for %s shared columns (missing=%s, extra=%s)",
        component,
        attr,
        len(shared),
        len(missing),
        len(extra),
    )


def clear_result_tables(base: pypsa.Network) -> None:
    """
    Clear stale solved dispatch and dual outputs from the reused solved network.

    After we replace snapshots and exogenous climate-year-dependent inputs, the
    old solved time series from the base SB run are no longer meaningful. The
    static optimized capacities (``*_opt``) remain untouched and are still used
    later by ``fix_optimal_capacities()`` in ``simplify_sb_network``.
    """
    for component_name in [
        "generators_t",
        "links_t",
        "loads_t",
        "stores_t",
        "storage_units_t",
    ]:
        component = getattr(base, component_name)
        for key in list(component.keys()):
            if key.startswith("mu_") or key in RESULT_KEYS:
                component[key] = pd.DataFrame(
                    index=base.snapshots,
                    columns=component[key].columns,
                    dtype=float,
                )
                logger.info("Cleared %s.%s", component_name, key)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_cba_sb_input_network",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the solved SB network from the base scenario. This network provides
    # the static topology and optimized capacities that we want to reuse for CBA.
    base = pypsa.Network(snakemake.input.solved_network).copy()
    logger.info("Loaded solved SB base network from %s", snakemake.input.solved_network)

    # Load the climate-year-specific pre-solve network. For the first planning
    # horizon this is the output of add_existing_baseyear. For later horizons it
    # is the output of prepare_sector_network.
    climate_source = pypsa.Network(snakemake.input.climate_source_network)
    logger.info(
        "Loaded climate-year source network from %s",
        snakemake.input.climate_source_network,
    )

    # For later horizons, the source network comes from prepare_sector_network.
    # That stage still uses unsuffixed names for new assets, so align them with
    # the solved myopic naming convention before copying climate-year-dependent
    # data onto the reused solved base network.
    if not snakemake.params.is_first_horizon:
        add_build_year_to_new_assets(
            climate_source,
            int(snakemake.wildcards.planning_horizons),
        )

    # Replace the solved base network snapshots with the target climate-year
    # snapshots. This resets the time-dependent table indices to the target year.
    base.set_snapshots(climate_source.snapshots)
    base.snapshot_weightings = climate_source.snapshot_weightings.copy()
    logger.info(
        "Copied %s snapshots and snapshot_weightings from climate-year source",
        len(climate_source.snapshots),
    )

    # Overwrite static climate-year-dependent inputs. At the first horizon, the
    # only required static field is loads.p_set.
    overwrite_static_attr(base, climate_source, "loads", "p_set")

    # Overwrite dynamic climate-year-dependent inputs. These are the fields we
    # currently know must follow the target climate year:
    # - renewable availability through generators_t.p_max_pu
    # - demand through loads_t.p_set
    overwrite_dynamic_attr(base, climate_source, "generators", "p_max_pu")
    overwrite_dynamic_attr(base, climate_source, "loads", "p_set")

    # Remove solved dispatch/dual time series inherited from the base solved SB
    # run. They no longer correspond to the modified target-year network.
    clear_result_tables(base)

    base.export_to_netcdf(snakemake.output.network)
    logger.info("Saved prepared CBA SB input network to %s", snakemake.output.network)
