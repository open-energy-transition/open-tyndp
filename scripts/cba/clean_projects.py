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

Storage project extraction is not yet implemented and returns an empty DataFrame.

**Inputs**

- ``data/tyndp_2024_bundle/cba_projects/20250312_export_transmission.xlsx``: Excel file containing CBA transmission projects
- ``data/tyndp_2024_bundle/cba_projects/20250312_export_storage.xlsx``: Excel file containing CBA storage projects (not yet processed)

**Outputs**

- ``resources/cba/transmission_projects.csv``: Cleaned CSV with columns:
  - ``project_id``: Integer project identifier
  - ``project_name``: Project name
  - ``in_reference2030``: Boolean indicating if project is in 2030 reference grid
  - ``in_reference2035``: Boolean indicating if project is in 2035 reference grid
  - ``border``: Border string in format "BUS0-BUS1"
  - ``p_nom 0->1``: Transfer capacity increase from bus0 to bus1 (MW)
  - ``p_nom 1->0``: Transfer capacity increase from bus1 to bus0 (MW)
  - ``bus0``: Source bus code (4 alphanumeric characters)
  - ``bus1``: Destination bus code (4 alphanumeric characters)

- ``resources/cba/storage_projects.csv``: Empty CSV with columns project_id and project_name (stub implementation)

"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

TRANSMISSION_PROJECTS_COLUMN_MAP = {
    "Project ID": "project_id",
    "Project Name": "project_name",
    "Is the project in the 2030 reference grid?": "in_reference2030",
    "Is the project in the 2035 reference grid?": "in_reference2035",
    "Border": "border",
    "Transfer capacity increase A-B (MW)": "p_nom 0->1",
    "Transfer capacity increase B-A (MW)": "p_nom 1->0",
}


def extract_transmission_projects(excel_path: Path) -> pd.DataFrame:
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

    projects = projects.assign(
        **projects["border"].str.extract(
            r"(?P<bus0>[A-Z0-9]{4})-(?P<bus1>[A-Z0-9]{4})$"
        )
    )

    unclear_border = projects["bus0"].isna() | projects["bus1"].isna()
    logger.warning(
        "%d out of %d projects do not follow the simple <bus0>-<bus1> format, ignoring them:\n%s",
        unclear_border.sum(),
        len(unclear_border),
        projects.loc[
            unclear_border, ["project_id", "project_name", "border"]
        ].to_string(index=False, max_colwidth=40, line_width=100),
    )
    projects = projects.loc[~unclear_border]
    return projects


def extract_storage_projects(excel_path: Path) -> pd.DataFrame:
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

    transmission_projects = extract_transmission_projects(
        Path(snakemake.input.dir) / "20250312_export_transmission.xlsx"
    )
    transmission_projects.to_csv(snakemake.output.transmission_projects, index=False)

    storage_projects = extract_storage_projects(
        Path(snakemake.input.dir) / "20250312_export_storage.xlsx"
    )
    storage_projects.to_csv(snakemake.output.storage_projects, index=False)
