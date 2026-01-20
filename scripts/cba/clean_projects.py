# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extracts and cleans CBA transmission and storage projects from Excel exports.

Reads transmission projects from the "Trans.Projects" sheet of the CBA projects Excel file.
For projects with multiple borders (newline-separated in the Excel), the script explodes
these into separate rows, creating one row per border. Bus codes are extracted from the
border strings (expected format: "BUS0-BUS1") and projects that don't match this format
are filtered out with a warning.

Also extracts CBA method assignments (TOOT/PINT) from scenario sheets (2030NT, 2040NT, 2040DE).

Storage project extraction is not yet implemented and returns an empty DataFrame.

**Inputs**

- ``data/tyndp_2024_bundle/cba_projects/20250312_export_transmission.xlsx``: Excel file containing CBA transmission projects
- ``data/tyndp_2024_bundle/cba_projects/20250312_export_storage.xlsx``: Excel file containing CBA storage projects (not yet processed)

**Outputs**

- ``resources/cba/transmission_projects.csv``: Cleaned CSV with columns:
  - ``project_id``: Integer project identifier
  - ``project_name``: Project name
  - ``in_reference2030``: Boolean indicating if project is in 2030 reference grid
  - ``in_reference2040``: Boolean indicating if project is in 2035 reference grid (used for 2040 horizon)
  - ``border``: Border string in format "BUS0-BUS1"
  - ``p_nom 0->1``: Transfer capacity increase from bus0 to bus1 (MW)
  - ``p_nom 1->0``: Transfer capacity increase from bus1 to bus0 (MW)
  - ``bus0``: Source bus code (4 alphanumeric characters)
  - ``bus1``: Destination bus code (4 alphanumeric characters)
  - ``length``: Total route length in km (from Trans.Investments sheet)

- ``resources/cba/storage_projects.csv``: Empty CSV with columns project_id and project_name (stub implementation)

- ``resources/cba/method_assignment.csv``: CBA method assignments with columns:
  - ``project_id``: Project identifier
  - ``scenario``: Scenario name (NT, DE)
  - ``planning_horizon``: Year (2030, 2040)
  - ``cba_method``: Assessment method (toot, pint)

"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

TRANSMISSION_PROJECTS_COLUMN_MAP = {
    "Project ID": "project_id",
    "Project Name": "project_name",
    "Is the project in the 2030 reference grid?": "in_reference2030",
    # 2035 reference grid is used for assessing CBA projects in planning horizon 2040:
    "Is the project in the 2035 reference grid?": "in_reference2040",
    "Is the project cross-border?": "is_crossborder",
    "Border": "border",
    "Transfer capacity increase A-B (MW)": "p_nom 0->1",
    "Transfer capacity increase B-A (MW)": "p_nom 1->0",
}

# Columns from Trans.Investments sheet for length data
INVESTMENT_LENGTH_COLUMNS = {
    "This investment belongs to project number…": "project_id",
    "Total route length (km)": "length",
}

# Mapping of Excel sheet names to (scenario, planning_horizon)
SCENARIO_SHEETS = {
    "2030NT": ("NT", 2030),
    "2040NT": ("NT", 2040),
    "2040DE": ("DE", 2040),
}


def extract_transmission_projects(
    excel_path: Path, existing_buses: pd.Index
) -> pd.DataFrame:
    projects = (
        pd.read_excel(
            excel_path,
            sheet_name="Trans.Projects",
            skiprows=1,
            usecols=list(TRANSMISSION_PROJECTS_COLUMN_MAP),
            na_values=["  "],
            true_values=["True", "WAHR"],
            false_values=["False", "FALSCH"],
        )
        .rename(columns=TRANSMISSION_PROJECTS_COLUMN_MAP)
        .dropna(how="all")  # contains many all nan rows
    )
    projects["project_id"] = projects["project_id"].astype(int)

    # The multiple lines columns in the "Expected transfer capacity increase" section
    # are exploded into separate rows
    multiline_columns = ["border", "p_nom 0->1", "p_nom 1->0"]
    projects = projects.assign(
        **{c: projects[c].str.split("\n") for c in multiline_columns}
    ).explode(multiline_columns)

    # Note: This dictionary is also defined in build_tyndp_network.py#build_links
    replace_dict = {"UK": "GB"}
    projects = projects.assign(
        **projects["border"]
        .replace(replace_dict, regex=True)
        .str.extract(r"(?P<bus0>[A-Za-z0-9]{4,}) ?- ?(?P<bus1>[A-Za-z0-9]{4,})$")
    )

    unclear_border = ~(
        projects["bus0"].isin(existing_buses) & projects["bus1"].isin(existing_buses)
    )
    logger.warning(
        "%d out of %d extensions do not follow the simple <bus0>-<bus1> format or are not defined in the base network, ignoring them:\n%s",
        unclear_border.sum(),
        len(unclear_border),
        projects.loc[
            unclear_border, ["project_id", "project_name", "border"]
        ].to_string(index=False, max_colwidth=40, line_width=100),
    )

    empty_capacity = projects["p_nom 0->1"].isna() & projects["p_nom 1->0"].isna()
    logger.warning(
        "%d out of %d extensions have no capacity, ignoring them:\n%s",
        empty_capacity.sum(),
        len(empty_capacity),
        projects.loc[
            empty_capacity, ["project_id", "project_name", "border"]
        ].to_string(index=False, max_colwidth=40, line_width=100),
    )

    projects = projects.loc[~(empty_capacity | unclear_border)]

    # Several projects have capacities with "Up to ..."
    up_to_projects = set()
    for col in ["p_nom 0->1", "p_nom 1->0"]:
        up_to = projects[col].str.startswith("Up to ")
        if up_to.any():
            projects.loc[up_to, col] = projects.loc[up_to, col].str[len("Up to ") :]
            up_to_projects.update(projects.loc[up_to, "project_name"])
    if up_to_projects:
        logger.info(
            f"Removed 'Up to ' capacity prefix from {len(up_to_projects)} projects:\n"
            + ", ".join(up_to_projects)
        )

    return projects


def extract_project_lengths(excel_path: Path) -> pd.DataFrame:
    """
    Extract total route length per project from Trans.Investments sheet.

    Aggregates length from all investments belonging to each project.

    Parameters
    ----------
    excel_path : Path
        Path to the transmission projects Excel file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: project_id, length (in km)
    """
    investments = pd.read_excel(
        excel_path,
        sheet_name="Trans.Investments",
        skiprows=1,
        usecols=list(INVESTMENT_LENGTH_COLUMNS.keys()),
    ).rename(columns=INVESTMENT_LENGTH_COLUMNS)

    investments = investments.dropna(subset=["project_id"])
    investments["project_id"] = investments["project_id"].astype(int)

    # Sum length per project (multiple investments per project)
    length_by_project = investments.groupby("project_id")["length"].sum().reset_index()

    logger.info(f"Extracted length data for {len(length_by_project)} projects")

    return length_by_project


def extract_method_assignment_from_excel(excel_path: Path) -> pd.DataFrame:
    """
    Extract CBA method assignments (TOOT/PINT) from Excel scenario sheets.

    Reads the "Chosen Approach" column from sheets 2030NT, 2040NT, 2040DE.
    Only includes market-relevant cross-border projects (those with a method assigned).

    Parameters
    ----------
    excel_path : Path
        Path to the transmission projects Excel file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: project_id, scenario, planning_horizon, cba_method
    """
    all_assignments = []

    for sheet_name, (scenario, planning_horizon) in SCENARIO_SHEETS.items():
        try:
            # Row 0 has column names, rows 1-2 are sub-headers, data starts row 3
            # Use usecols by index since skiprows shifts the header
            df = pd.read_excel(
                excel_path,
                sheet_name=sheet_name,
                skiprows=2,  # Skip the sub-header rows (keep row 0 as header)
                usecols=[0, 2],  # Project ID (col 0) and Chosen Approach (col 2)
            )
            df.columns = ["project_id", "cba_method"]

            # Filter out rows without method assignment (non-market-relevant projects)
            df = df.dropna(subset=["project_id", "cba_method"])
            df["project_id"] = df["project_id"].astype(int)

            # Normalize method names (TOOT, Light TOOT -> toot; PINT, Light PINT -> pint)
            df["cba_method"] = df["cba_method"].str.lower().str.replace("light ", "")

            # Add scenario and planning horizon
            df["scenario"] = scenario
            df["planning_horizon"] = planning_horizon

            all_assignments.append(
                df[["project_id", "scenario", "planning_horizon", "cba_method"]]
            )

            logger.info(
                f"Extracted {len(df)} method assignments from sheet {sheet_name}"
            )

        except Exception as e:
            logger.warning(f"Could not read sheet {sheet_name}: {e}")

    if not all_assignments:
        logger.warning("No method assignments found in any scenario sheet")
        return pd.DataFrame(
            columns=["project_id", "scenario", "planning_horizon", "cba_method"]
        )

    return pd.concat(all_assignments, ignore_index=True)


def extract_method_assignment_from_file(file_path: Path) -> pd.DataFrame:
    """
    Read existing method assignment from CSV file.

    Used when method_assignment mode is 'custom_labels'.

    Parameters
    ----------
    file_path : Path
        Path to the method_assignment.csv file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: project_id, scenario, planning_horizon, cba_method
    """
    return pd.read_csv(file_path)


def extract_storage_projects(
    excel_path: Path, existing_buses: pd.Index
) -> pd.DataFrame:
    """
    Stub method to extract storage projects.

    Returns an empty DataFrame with the expected column structure.
    TODO: Implement actual storage project extraction from Excel file.
    """
    logger.info(
        "Storage project extraction not yet implemented, returning empty DataFrame"
    )
    return pd.DataFrame(columns=["project_id", "project_name"])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_projects", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    network = pypsa.Network(snakemake.input.network)
    existing_buses = network.buses.index[network.buses.carrier == "AC"].unique()

    excel_path = Path(snakemake.input.dir) / "20250312_export_transmission.xlsx"

    # Extract transmission projects
    transmission_projects = extract_transmission_projects(excel_path, existing_buses)

    # Extract and merge project lengths from Trans.Investments sheet
    project_lengths = extract_project_lengths(excel_path)
    transmission_projects = transmission_projects.merge(
        project_lengths, on="project_id", how="left"
    )

    transmission_projects.to_csv(snakemake.output.transmission_projects, index=False)

    # Extract storage projects
    storage_projects = extract_storage_projects(
        Path(snakemake.input.dir) / "20250312_export_storage.xlsx",
        existing_buses,
    )
    storage_projects.to_csv(snakemake.output.storage_projects, index=False)

    # Extract method assignments based on config
    method_assignment_mode = snakemake.params.method_assignment

    if method_assignment_mode == "entsoe_labels":
        method_assignment = extract_method_assignment_from_excel(excel_path)
    elif method_assignment_mode == "custom_labels":
        # Read from existing file (user provides custom assignments)
        method_assignment = extract_method_assignment_from_file(
            snakemake.output.method_assignment
        )
    elif method_assignment_mode == "algorithmic":
        # TODO: Implement algorithmic method assignment
        logger.warning(
            "Algorithmic method assignment not yet implemented, "
            "falling back to entsoe_labels"
        )
        method_assignment = extract_method_assignment_from_excel(excel_path)
    else:
        raise ValueError(
            f"Unknown method_assignment mode: {method_assignment_mode}. "
            f"Valid options: entsoe_labels, custom_labels, algorithmic"
        )

    method_assignment.to_csv(snakemake.output.method_assignment, index=False)
    logger.info(
        f"Saved method assignments: {len(method_assignment)} entries "
        f"({len(method_assignment['project_id'].unique())} unique projects)"
    )
