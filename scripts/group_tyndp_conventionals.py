# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Groups TYNDP conventional carriers according to a pre-determined mapping.

This script takes PEMMDB capacities and profiles and groups the TYNDP conventional carriers according
to a pre-determined mapping. For capacity and energy values, the sum of the grouped components is calculated.
For the aggregated efficiency value, the mean is taken. This script will only be executed if activated in the
configuration and will be executed once for each planning horizon.

Inputs
------

- ``pemmdb_capacities_{planning_horizon}.csv``: Processed PEMMDB capacities for the given planning_horizon.
- ``pemmdb_profiles_{planning_horizon}.nc``: Processed PEMMDB must-run and availability profiles for the given
    planning_horizon.
- ``data/tyndp_technology_map.csv``: TYNDP technology mapping used for the grouping.

Outputs
-------

- ``pemmdb_capacities_{planning_horizon}_grouped.csv``: Grouped PEMMDB capacities for the given planning_horizon.
- ``pemmdb_profiles_{planning_horizon}_grouped.nc``: Grouped PEMMDB must-run and availability profiles for the given
    planning_horizon.
"""

import logging

import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

# Aggregation strategies for specific columns
AGG_STRATEGIES = {
    "p_nom": "sum",
    "e_nom": "sum",
    "p_min": "sum",
    "p_max": "sum",
    "efficiency": "mean",
}


def _group_conv(
    df: pd.DataFrame, groupby: list[str], agg_strategies: dict[str, str]
) -> pd.DataFrame:
    """
    Group conventional carriers based on specified columns.
    """

    # Aggregation dict: use specific strategy if defined, else "first"
    agg_dict = {
        col: agg_strategies.get(col, "first")
        for col in df.columns
        if col not in groupby
    }

    return (
        df.groupby(groupby, dropna=False)
        .agg(agg_dict)
        .reset_index()
        .drop(columns="index_carrier")
        .rename(
            columns={
                "open_tyndp_type": "index_carrier",
            },
            errors="ignore",
        )
    )


def _group_capacities(
    caps: pd.DataFrame,
    conventional_carriers: list[str],
) -> pd.DataFrame:
    """
    Group conventional capacities by bus and open_tyndp_type.
    """
    # Split into conventional and other
    caps_conv = caps.query("carrier in @conventional_carriers")
    caps_other = caps.query("carrier not in @conventional_carriers")

    # Group conventional capacities together
    caps_conv_grouped = _group_conv(
        caps_conv, groupby=["bus", "open_tyndp_type"], agg_strategies=AGG_STRATEGIES
    )

    # Recombine
    return pd.concat([caps_conv_grouped, caps_other], ignore_index=True)


def _group_profiles(
    profiles: pd.DataFrame,
    caps_conv: pd.DataFrame,
    caps_grouped: pd.DataFrame,
    conventional_carriers: list[str],
) -> pd.DataFrame:
    """
    Group conventional profiles by time, bus, and open_tyndp_type.
    """
    # Split into conventional and other
    profiles_conv = profiles.query("carrier in @conventional_carriers")
    profiles_other = profiles.query("carrier not in @conventional_carriers")

    # Calculate absolute values for p_min_t and p_max_t
    profiles_conv_abs = profiles_conv.merge(
        caps_conv[["bus", "index_carrier", "p_nom"]],
        on=["bus", "index_carrier"],
        how="left",
    ).assign(
        p_min_t=lambda df: df.p_min_pu * df.p_nom,
        p_max_t=lambda df: df.p_max_pu * df.p_nom,
    )

    # Group absolute values together
    profiles_conv_grouped = _group_conv(
        profiles_conv_abs,
        groupby=["time", "bus", "open_tyndp_type"],
        agg_strategies=AGG_STRATEGIES,
    )

    # Convert must-run and availability values back to per-unit values using grouped capacities.
    # Some sub-types within a grouped carrier lack profiles (e.g., "gas-ccgt-old1" within "gas-ccgt" in BE00).
    # Solution: (1) merge total grouped capacity from before, (2) add missing capacity to p_max
    # assuming default p_max_pu=1.0, (3) recalculate per-unit values.
    profiles_conv_pu = (
        profiles_conv_grouped.rename(columns={"p_nom": "p_nom_profiles"})
        .merge(
            caps_grouped[["bus", "index_carrier", "p_nom"]],
            on=["bus", "index_carrier"],
            how="left",
        )
        .assign(
            p_max_t=lambda df: df.p_max_t
            + (df.p_nom - df.p_nom_profiles),  # add missing capacities to p_max
            p_min_pu=lambda df: df.p_min_t / df.p_nom,
            p_max_pu=lambda df: df.p_max_t / df.p_nom,
            open_tyndp_type=lambda df: df.index_carrier,
        )
        .drop(
            columns=["p_min_t", "p_max_t", "p_nom", "p_nom_profiles"], errors="ignore"
        )
    )

    # Recombine
    return pd.concat([profiles_conv_pu, profiles_other], ignore_index=True)


def group_tyndp_conventionals(
    pemmdb_capacities: pd.DataFrame,
    pemmdb_profiles: pd.DataFrame,
    tyndp_conventional_carriers: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Group conventional thermal generation technologies by TYNDP type.

    Groups capacities and profiles for conventional carriers according to the
    open_tyndp_type mapping, combining multiple carriers into single entries.

    Parameters
    ----------
    pemmdb_capacities : pd.DataFrame
        All PEMMDB capacities.
    pemmdb_profiles : pd.DataFrame
        All PEMMDB must-run and availability profiles.
    tyndp_conventional_carriers : list[str]
        List of TYNDP conventional carriers to group

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Grouped capacities and profiles.
    """
    logger.info("Grouping TYNDP conventional thermal generation technologies.")

    # Store original conventional capacities before grouping
    caps_conv_original = pemmdb_capacities.query(
        "carrier in @tyndp_conventional_carriers"
    )

    # Group capacities
    capacities_grouped = _group_capacities(
        pemmdb_capacities, tyndp_conventional_carriers
    )

    # Group profiles
    profiles_grouped = _group_profiles(
        pemmdb_profiles,
        caps_conv_original,
        capacities_grouped,
        tyndp_conventional_carriers,
    )

    return capacities_grouped, profiles_grouped


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "group_tyndp_conventionals",
            planning_horizon="2030",
        )

    configure_logging(snakemake)  # pylint: disable=E0606

    # Parameters
    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers

    # Read in PEMMDB data and trajectories
    pemmdb_capacities = pd.read_csv(snakemake.input.pemmdb_capacities)
    pemmdb_profiles = xr.open_dataset(snakemake.input.pemmdb_profiles).to_dataframe()

    # Read in TYNDP carrier mapping and conventional carriers
    conventional_carrier_mapping = (
        pd.read_csv(snakemake.input.carrier_mapping)
        .set_index("open_tyndp_index")[
            ["open_tyndp_carrier", "open_tyndp_type", "pypsa_eur_carrier"]
        ]
        .query("open_tyndp_carrier in @tyndp_conventional_carriers")
        .replace(
            {"open_tyndp_carrier": ["oil-light", "oil-heavy", "oil-shale"]}, "oil"
        )  # TODO To remove once the three carriers have been implemented
    )

    # Group TYNDP conventionals
    pemmdb_capacities, pemmdb_profiles = group_tyndp_conventionals(
        pemmdb_capacities=pemmdb_capacities,
        pemmdb_profiles=pemmdb_profiles,
        tyndp_conventional_carriers=tyndp_conventional_carriers,
    )

    # Save grouped capacities and profiles
    pemmdb_capacities.to_csv(snakemake.output.pemmdb_capacities_grouped, index=False)
    ds_grouped = xr.Dataset(
        {
            "p_min_pu": (["sample"], pemmdb_profiles["p_min_pu"]),
            "p_max_pu": (["sample"], pemmdb_profiles["p_max_pu"]),
        },
        coords={
            column: (["sample"], pemmdb_profiles[column])
            for column in pemmdb_profiles.columns
            if column not in ["p_min_pu", "p_max_pu"]
        },
    )
    ds_grouped.to_netcdf(snakemake.output.pemmdb_profiles_grouped)
