# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import re

import pandas as pd


def filter_projects_by_method(
    projects: pd.DataFrame, method: str, planning_horizon: int
) -> pd.DataFrame:
    """
    Filter projects based on CBA method and planning horizon.

    Parameters
    ----------
    projects : pd.DataFrame
        DataFrame with project data including in_reference columns
    method : str
        CBA method ('toot' or 'pint')
    planning_horizon : int
        Planning horizon year (e.g., 2030, 2040)

    Returns
    -------
    pd.DataFrame
        Filtered projects appropriate for the given method
    """
    method = method.lower()
    if method not in METHOD_PROJECT_FILTERS:
        raise ValueError(
            f"Unknown CBA method: {method}. Valid: {list(METHOD_PROJECT_FILTERS.keys())}"
        )

    return METHOD_PROJECT_FILTERS[method](projects, planning_horizon)


def _filter_toot(projects: pd.DataFrame, planning_horizon: int) -> pd.DataFrame:
    """TOOT: all projects with data for this planning horizon are candidates."""
    in_ref_col = f"in_reference{planning_horizon}"
    if in_ref_col not in projects.columns:
        raise ValueError(f"Column {in_ref_col} not found in projects DataFrame")
    return projects.loc[projects[in_ref_col].notna()]


def _filter_pint(projects: pd.DataFrame, planning_horizon: int) -> pd.DataFrame:
    """PINT: only projects with in_reference=False are candidates."""
    in_ref_col = f"in_reference{planning_horizon}"
    if in_ref_col not in projects.columns:
        raise ValueError(f"Column {in_ref_col} not found in projects DataFrame")
    return projects.loc[projects[in_ref_col] == False]


# Registry of method-specific project filters
METHOD_PROJECT_FILTERS = {
    "toot": _filter_toot,
    "pint": _filter_pint,
}


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
