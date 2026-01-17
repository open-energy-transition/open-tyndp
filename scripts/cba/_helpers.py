# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import re
from pathlib import Path

import pandas as pd


def load_method_assignment(method_assignment_path: str | Path) -> pd.DataFrame:
    """
    Load CBA method assignments from CSV file.

    Parameters
    ----------
    method_assignment_path : str or Path
        Path to method_assignment.csv

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: project_id, scenario, planning_horizon, cba_method
    """
    return pd.read_csv(method_assignment_path)


def filter_projects_by_method(
    projects: pd.DataFrame,
    method: str,
    planning_horizon: int,
    scenario: str,
    method_assignment: pd.DataFrame,
) -> pd.DataFrame:
    """
    Filter projects based on CBA method, planning horizon, and scenario.

    Parameters
    ----------
    projects : pd.DataFrame
        DataFrame with project data
    method : str
        CBA method ('toot' or 'pint')
    planning_horizon : int
        Planning horizon year (e.g., 2030, 2040)
    scenario : str
        Scenario name (e.g., 'NT', 'DE')
    method_assignment : pd.DataFrame
        DataFrame with method assignments (from method_assignment.csv)

    Returns
    -------
    pd.DataFrame
        Filtered projects matching the given method, horizon, and scenario
    """
    method = method.lower()
    if method not in ("toot", "pint"):
        raise ValueError(f"Unknown CBA method: {method}. Valid: toot, pint")

    # Filter method assignments for this scenario/horizon/method
    mask = (
        (method_assignment["scenario"] == scenario)
        & (method_assignment["planning_horizon"] == planning_horizon)
        & (method_assignment["cba_method"].str.lower() == method)
    )
    valid_project_ids = method_assignment.loc[mask, "project_id"].unique()

    # Return projects matching the valid IDs
    return projects[projects["project_id"].isin(valid_project_ids)]


def get_toot_projects(
    projects: pd.DataFrame,
    planning_horizon: int,
    scenario: str,
    method_assignment: pd.DataFrame,
) -> pd.DataFrame:
    """
    Get all TOOT projects for a given scenario and planning horizon.

    Used for building the unified CBA reference network.

    Parameters
    ----------
    projects : pd.DataFrame
        DataFrame with project data
    planning_horizon : int
        Planning horizon year (e.g., 2030, 2040)
    scenario : str
        Scenario name (e.g., 'NT', 'DE')
    method_assignment : pd.DataFrame
        DataFrame with method assignments

    Returns
    -------
    pd.DataFrame
        TOOT projects for this scenario/horizon
    """
    return filter_projects_by_method(
        projects, "toot", planning_horizon, scenario, method_assignment
    )


def get_pint_projects(
    projects: pd.DataFrame,
    planning_horizon: int,
    scenario: str,
    method_assignment: pd.DataFrame,
) -> pd.DataFrame:
    """
    Get all PINT projects for a given scenario and planning horizon.

    Parameters
    ----------
    projects : pd.DataFrame
        DataFrame with project data
    planning_horizon : int
        Planning horizon year (e.g., 2030, 2040)
    scenario : str
        Scenario name (e.g., 'NT', 'DE')
    method_assignment : pd.DataFrame
        DataFrame with method assignments

    Returns
    -------
    pd.DataFrame
        PINT projects for this scenario/horizon
    """
    return filter_projects_by_method(
        projects, "pint", planning_horizon, scenario, method_assignment
    )


def filter_projects_by_specs(
    project_list: list[str], spec_list: list[str] | str | None
) -> list[str]:
    """
    Filter projects based on specifications with inclusions and exclusions.

    Supports:
    - Single projects: 't1', 's4'
    - Ranges: 't20-t25', 's4-s6'
    - Exclusions: '-t22', '-s5'
    - Exclusion ranges: '-t22-t25'

    The function operates in two modes:
    - Inclusion mode (default): Start with empty set, add specified projects
    - Removal mode: Start with all projects, remove specified ones (when first spec starts with '-')

    Parameters
    ----------
    project_list : list[str]
        List of all available project names to filter from
    spec_list : list[str], str, or None
        List of specifications, a single specification string, or None to return all projects

    Returns
    -------
    list[str]
        Filtered list of projects, preserving order from project_list

    Examples
    --------
    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], ['t20-t22'])
    ['t20', 't21', 't22']

    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], 't20-t22')
    ['t20', 't21', 't22']

    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], ['t20-t23', '-t22'])
    ['t20', 't21', 't23']

    >>> filter_projects_by_specs(['t1', 't2', 't3', 't4'], ['-t2', '-t3'])
    ['t1', 't4']
    """

    if not spec_list:
        return project_list

    if isinstance(spec_list, str):
        spec_list = [spec_list]

    projects = set()
    range_pattern = re.compile(r"^([a-z])(\d+)-\1(\d+)$")

    removals = spec_list[0].startswith("-")

    for spec in spec_list:
        # Check if this is an exclusion
        if spec.startswith("-"):
            spec = spec[1:]
            op = projects.discard if not removals else projects.add
        else:
            op = projects.add if not removals else projects.discard

        # Try to match range pattern
        match = range_pattern.match(spec)
        if match:
            prefix = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))

            for i in range(start, end + 1):
                op(f"{prefix}{i}")
        else:
            # Single project
            op(spec)

    if not removals:
        return [p for p in project_list if p in projects]
    else:
        return [p for p in project_list if p not in projects]
