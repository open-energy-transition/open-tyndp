# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP 2026 H2 import corridor potentials and profiles for a given
planning horizon.

Inputs
------

- `H2 IMPORTS GENERATORS PROPERTIES.xlsx`: one row per corridor/band/year with
  ``bus0``, ``bus1``, static or time-series-flagged capacity, and offer price.
- `H2 IMPORT PROFILES.xlsx`: hourly (8760h, no weather-year dimension) import
  capacity per corridor/band/year, for the subset of corridors flagged
  ``"TIME-SERIES-DATA"`` in the properties file.

Outputs
-------

- `resources/h2_import_potentials_{planning_horizons}.csv`: one row per
  corridor, with ``p_nom`` set to the static capacity, or to the profile's
  maximum for time-series corridors.
- `resources/h2_import_profiles_{planning_horizons}.csv`: hourly ``p_max_pu``
  (normalized to each corridor's own maximum) for the time-series corridors
  only.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config
from scripts.sb.build_tyndp_demand import check_snapshot_year

logger = logging.getLogger(__name__)

TIME_SERIES_FLAG = "TIME-SERIES-DATA"


def load_import_potentials(fn: str, pyear: int) -> pd.DataFrame:
    """
    Load and clean the TYNDP 2026 H2 import corridor properties for one planning horizon.

    Parameters
    ----------
    fn : str
        Path to the TYNDP 2026 H2 import generator properties Excel file
        ("H2 IMPORTS GENERATORS PROPERTIES.xlsx").
    pyear : int
        Planning horizon to filter for.

    Returns
    -------
    pd.DataFrame
        Cleaned TYNDP H2 import corridor properties, indexed by ``Corridor``,
        with a boolean ``has_profile`` column marking corridors whose
        capacity is given by an hourly profile instead of a static value.
    """
    column_dict = {
        "CORRIDOR": "Corridor",
        "NODE FROM": "bus0",
        "NODE TO": "bus1",
        "MAX CAPACITY [MW]": "p_nom",
        "OFFER PRICE [€/MWh]": "marginal_cost",
    }

    imports = (
        pd.read_excel(fn, engine="calamine")
        .query("YEAR == @pyear")
        .rename(columns=column_dict)
        .replace({"Type": {"Lh2": "LH2"}})
        .set_index("Corridor")
    )
    imports["marginal_cost"] = pd.to_numeric(
        imports["marginal_cost"], errors="coerce"
    ).fillna(0.0)
    imports["has_profile"] = imports["p_nom"] == TIME_SERIES_FLAG
    imports["p_nom"] = pd.to_numeric(
        imports["p_nom"].where(~imports["has_profile"]), errors="coerce"
    )

    return imports[
        ["bus0", "bus1", "Type", "Fuel", "p_nom", "marginal_cost", "has_profile"]
    ]


def load_import_profiles(
    fn: str, pyear: int, corridors: pd.Index, year: int
) -> pd.DataFrame:
    """
    Load hourly TYNDP 2026 H2 import profiles for a set of corridors.

    Parameters
    ----------
    fn : str
        Path to the TYNDP 2026 H2 import profiles Excel file
        ("H2 IMPORT PROFILES.xlsx").
    pyear : int
        Planning horizon; profile columns are named ``"{Corridor}-{pyear}"``.
    corridors : pd.Index
        Corridors to load profiles for.
    year : int
        Year to assign to the resulting DatetimeIndex.

    Returns
    -------
    pd.DataFrame
        Hourly import capacity per corridor, columns named by ``Corridor``.
    """
    profiles = pd.read_excel(fn, engine="calamine", index_col="Hour")

    columns = {f"{corridor}-{pyear}": corridor for corridor in corridors}
    missing = set(columns) - set(profiles.columns)
    if missing:
        raise ValueError(
            f"H2 import profiles missing expected columns for {pyear}: {sorted(missing)}"
        )

    profiles = profiles[list(columns)].rename(columns=columns)
    profiles.index = pd.date_range(
        start=f"{year}-01-01", periods=len(profiles), freq="h"
    )
    profiles.index.name = "datetime"

    return profiles


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_imports",
            planning_horizons=2030,
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    year = snapshots[0].year
    check_snapshot_year(year, snakemake.params.drop_leap_day)

    # Load corridor properties
    import_potentials = load_import_potentials(
        snakemake.input.import_potentials_raw, pyear
    )

    profile_corridors = import_potentials.index[import_potentials.has_profile]
    if not profile_corridors.empty:
        import_profiles = load_import_profiles(
            snakemake.input.import_profiles_raw, pyear, profile_corridors, year
        )
        maxima = import_profiles.max()
        import_potentials.loc[profile_corridors, "p_nom"] = maxima
        import_profiles = import_profiles.div(maxima).fillna(0.0)
    else:
        import_profiles = pd.DataFrame(index=pd.DatetimeIndex([], name="datetime"))

    # Save corridor potentials and profiles
    import_potentials.drop(columns="has_profile").to_csv(
        snakemake.output.import_potentials
    )
    import_profiles.to_csv(snakemake.output.import_profiles)
