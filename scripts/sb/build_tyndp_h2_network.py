# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script loads and cleans the TYNDP H2 reference grid and interzonal connections for a given wildcard planning horizon
and TYNDP scenario as defined in the config file. The reference grid contains data for the TYNDP planning year 2030,
while depending on the scenario, different planning years (`pyear`) are available for the interzonal connections.
DE and GA are defined for 2030, 2035, 2040, 2045 and 2050. For the NT scenario no interzonal capacities are defined.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    extract_grid_data_tyndp,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def normalize_starting_grid_h2_nodes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize node IDs from the newer H2 starting grid workbook to the country-level
    H2 node IDs currently used in the TYNDP H2 topology.

    This is needed because the newer H2 starting grid workbook contains node IDs with trailing zeros (e.g. ``DE00``)
    and some exceptions (e.g. ``UK00`` instead of ``GB00``) that need to be normalized to match the country-level node IDs
    used in the TYNDP H2 topology (e.g. ``DE``, ``GB``).

    Examples:
    - ``AT00`` -> ``AT``
    - ``IBIT00`` -> ``IBIT``
    - ``UK00`` -> ``GB``
    """

    df = df.copy()
    for col in ["bus0", "bus1"]:
        df[col] = (
            df[col]
            .astype(str)
            .str.replace(r"00$", "", regex=True)
            .str.replace(r"^UK$", "GB", regex=True)
        )
    return df


def load_h2_interzonal_connections(fn, scenario="GA", pyear=2030):
    """
    Load and clean H2 interzonal connections.
    Returns the cleaned interzonal connections as dataframe.

    Parameters
    ----------
    fn : str
        Path to Excel file containing H2 interzonal data.
    scenario : str
        TYNDP scenario to use for interzonal connection data.
        Possible options are:
        - 'GA'
        - 'DE'
        - 'NT'
    pyear : int
        TYNDP planning horizon to use for interzonal connection data.
        Possible options are:
        - 2030
        - 2035
        - 2040
        - 2045
        - 2050

    Returns
    -------
    pd.DataFrame
        The function returns cleaned TYNDP H2 interzonal connections.
    """

    if scenario in ["DE", "GA"]:
        interzonal_raw = pd.read_excel(fn, sheet_name="Hydrogen_Interzonal")
        interzonal_raw.columns = interzonal_raw.columns.str.title()

        if int(pyear) not in [2030, 2035, 2040, 2045, 2050]:
            logger.warning(
                "Planning horizon doesn't match available TYNDP data. "
                "Falling back to closest available year between 2030 and 2050."
            )
            pyear = np.clip(5 * (pyear // 5), 2030, 2050)
        scenario_dict = {v: k for k, v in SCENARIO_DICT.items()}
        scenario = scenario_dict[scenario]
        interzonal_filtered = interzonal_raw.query(
            "Scenario == @scenario and Year == @pyear "
        )

        interzonal = extract_grid_data_tyndp(
            interzonal_filtered, idx_prefix="H2 pipeline", idx_connector="->"
        )
        # convert from GW to PyPSA base unit MW as raw H2 reference grid data is given in GW
        interzonal["p_nom"] = interzonal.p_nom.mul(1e3)

    elif scenario == "NT":
        logger.info(
            "No interzonal capacities for 'NT' scenario. Saving empty file for interzonal capacities."
        )
        interzonal = pd.DataFrame()
    else:
        raise ValueError(
            "Unknown scenario requested. Please, choose from 'GA', 'DE' or 'NT'."
        )

    return interzonal


def load_h2_grid_entsoe(fn_grid: str, fn_projects: str | None) -> pd.DataFrame:
    """
    Load and clean the ENTSO-E H2 reference grid and format data.

    Adds H2 investment candidate projects to the reference grid and returns the
    cleaned reference grid as dataframe.

    Parameters
    ----------
    fn_grid : str
        Path to Excel file containing the ENTSO-E H2 reference grid data.
    fn_projects : str
        Path to CSV file containing H2 projects data.

    Returns
    -------
    pd.DataFrame
        The function returns the cleaned TYNDP H2 reference grid.
    """

    h2_grid_raw = pd.read_excel(fn_grid)
    h2_grid = extract_grid_data_tyndp(
        h2_grid_raw, idx_prefix="H2 pipeline", idx_connector="->"
    )
    # convert from GW to PyPSA base unit MW as raw H2 reference grid data is given in GW
    h2_grid["p_nom"] = h2_grid.p_nom.mul(1e3)

    # add projects to the grid
    if fn_projects:
        h2_projects = pd.read_csv(fn_projects, index_col=0)
        h2_grid = (
            pd.concat([h2_grid, h2_projects])
            .groupby(level=0)
            .agg({"bus0": "first", "bus1": "first", "p_nom": sum})
        )

    return h2_grid


def load_h2_grid_entsos(fn_grid: str, pyear: int) -> pd.DataFrame:
    """
    Load and clean the ENTSO-E/ENTSOG joint scenarios portal H2 starting grid.

    The workbook contains planning-horizon-specific sheets (e.g. ``H_2030``,
    ``H_2040``) already resolved to electricity-node IDs such as ``DE00``.
    """

    available_years = [2030, 2040, 2050]
    if pyear not in available_years:
        fallback = min(available_years, key=lambda y: abs(y - pyear))
        logger.warning(
            "Planning horizon %s is not available in StartingGrid2030.xlsx. "
            "Falling back to H_%s.",
            pyear,
            fallback,
        )
        pyear = fallback

    sheet_name = f"H_{pyear}"
    h2_grid_raw = pd.read_excel(fn_grid, sheet_name=sheet_name)
    required_cols = {"Border", "Summary Direction 1", "Summary Direction 2"}
    missing_cols = required_cols.difference(h2_grid_raw.columns)
    if missing_cols:
        raise KeyError(
            f"Sheet '{sheet_name}' in {fn_grid} is missing required columns: {sorted(missing_cols)}"
        )

    h2_grid = extract_grid_data_tyndp(
        h2_grid_raw, idx_prefix="H2 pipeline", idx_connector="->"
    )
    h2_grid = normalize_starting_grid_h2_nodes(h2_grid)
    h2_grid["p_nom"] = h2_grid.p_nom.mul(1e3)

    return h2_grid


def load_h2_grid(
    source: str,
    fn_grid_entsoe: str,
    fn_grid_entsos: str,
    fn_projects: str | None,
    pyear: int,
) -> pd.DataFrame:
    """
    Load the corresponding H2 grid based on the source.
    """

    if source == "entsoe":
        return load_h2_grid_entsoe(fn_grid_entsoe, fn_projects)
    if source == "entsos":
        return load_h2_grid_entsos(fn_grid_entsos, pyear)

    raise ValueError(
        f"Unknown H2 reference grid source '{source}'. Expected 'entsoe' or 'entsos'."
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_network",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params.scenario
    source = snakemake.params.h2_reference_grid_source
    pyear = int(snakemake.wildcards.planning_horizons)
    cyear = get_snapshots(snakemake.params.snapshots)[0].year

    # Load and prep H2 reference grid and interzonal pipeline capacities
    h2_grid = load_h2_grid(
        source=source,
        fn_grid_entsoe=snakemake.input.h2_reference_grid_entsoe,
        fn_grid_entsos=snakemake.input.h2_reference_grid_entsos,
        fn_projects=snakemake.input.h2_projects,
        pyear=pyear,
    )
    interzonal = load_h2_interzonal_connections(
        fn=snakemake.input.h2_reference_grid_entsoe, scenario=scenario, pyear=pyear
    )

    # Save prepped H2 grid and interzonal
    h2_grid.to_csv(snakemake.output.h2_grid_prepped)
    interzonal.to_csv(snakemake.output.interzonal_prepped)
