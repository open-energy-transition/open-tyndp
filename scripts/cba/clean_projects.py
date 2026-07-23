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

- `data/tyndp_2024_bundle/cba_projects/20250312_export_transmission.xlsx`: Excel file containing CBA transmission projects
- `data/tyndp_2024_bundle/cba_projects/20250312_export_storage.xlsx`: Excel file containing CBA storage projects (not yet processed)

**Outputs**

- `resources/cba/transmission_projects.csv`: Cleaned CSV with columns:
  - `project_id`: Integer project identifier
  - `project_name`: Project name
  - `border`: Border string in format "BUS0-BUS1"
  - `p_nom 0->1`: Transfer capacity increase from bus0 to bus1 (MW)
  - `p_nom 1->0`: Transfer capacity increase from bus1 to bus0 (MW)
  - `bus0`: Source bus code (4 alphanumeric characters)
  - `bus1`: Destination bus code (4 alphanumeric characters)
  - `length_km`: Total route length in km (from Trans.Investments)
  - `capex_meur`: Total estimated CAPEX in MEUR (from Trans.Investments)
  - `underwater_fraction`: Fraction of route that is offshore cable

- `resources/cba/storage_projects.csv`: Empty CSV with columns project_id and project_name (stub implementation)

"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

TRANSMISSION_PROJECTS_COLUMN_MAP = {
    "Project ID": "project_id",
    "Project Name": "project_name",
    "Is the project cross-border?": "is_crossborder",
    "Border": "border",
    "Transfer capacity increase A-B (MW)": "p_nom 0->1",
    "Transfer capacity increase B-A (MW)": "p_nom 1->0",
}

OFFSHORE_ELEMENT_TYPES = {
    "OffshoreDCTransmissionCable",
    "OffshoreACTransmissionCable",
}


def read_tyndp_electricity_buses(buses_fn: str):
    """
    Read node list for electricity from tyndp data input.

    Parameters
    ----------
        - buses_fn (str): Path to "LIST OF NODES.xlsx" from tyndp bundle

    Returns
    -------
        - buses: Index of electricity buses as used in Open-TYNDP

    See Also
    --------
        build_tyndp_network.py : build_buses
    """
    buses = pd.Index(
        pd.read_excel(buses_fn)
        .replace("UK", "GB", regex=True)
        .rename({"NODE": "bus_id"}, axis=1)["bus_id"]
    )

    # Manually add Italian virtual nodes and Corsica
    buses = buses.union(["ITCO", "ITVI", "FR15"])

    return buses


def remove_unclear_border(
    projects: pd.DataFrame, existing_buses: pd.DataFrame
) -> pd.DataFrame:
    """
    Remove projects defined by unclear borders from the list of projects.

    Parameters
    ----------
    projects : pd.DataFrame
        List of projects to assess.
    existing_buses : pd.DataFrame
        List of existing buses.

    Returns
    -------
    pd.DataFrame
        Curated list of projects that only use existing buses.
    """
    if projects.empty:
        return projects

    unclear_border = ~(
        projects["bus0"].isin(existing_buses) & projects["bus1"].isin(existing_buses)
    )
    if unclear_border.sum() > 0:
        logger.warning(
            "%d out of %d extensions do not follow the simple <bus0>-<bus1> format or are not defined in the base network, ignoring them:\n%s",
            unclear_border.sum(),
            len(unclear_border),
            projects.loc[
                unclear_border, ["project_id", "project_name", "border"]
            ].to_string(index=False, max_colwidth=40, line_width=100),
        )

    return projects.loc[~unclear_border]


def remove_no_capacity(projects: pd.DataFrame) -> pd.DataFrame:
    """
    Remove projects with no capacity from the list of projects.

    Parameters
    ----------
    projects : pd.DataFrame
        List of projects to clean.

    Returns
    -------
    pd.DataFrame
        Curated list of projects with a defined capacity
    """
    if projects.empty:
        return projects

    empty_capacity = projects["p_nom 0->1"].isna() & projects["p_nom 1->0"].isna()
    if empty_capacity.sum() > 0:
        logger.warning(
            "%d out of %d extensions have no capacity, ignoring them:\n%s",
            empty_capacity.sum(),
            len(empty_capacity),
            projects.loc[
                empty_capacity, ["project_id", "project_name", "border"]
            ].to_string(index=False, max_colwidth=40, line_width=100),
        )

    return projects.loc[~empty_capacity]


def extract_transmission_projects(
    excel_path: Path, existing_buses: pd.Index
) -> pd.DataFrame:
    """
    Extract transmission projects.

    Parameters
    ----------
    transmission_path : Path
        File path to transmission projects.
    existing_buses : pd.Index
        List of existing buses.

    Returns
    -------
    pd.DataFrame
        Curated list of projects.
    """
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

    # For project t339, the border is given as ITCS-ITSI and ITSA-ITSI, when it should be connected to the virtual node ITVI instead: ITCS-ITVI and ITSA-ITVI
    # Manually fixing this here (changing the border and bus1 columns)
    t339_mask = projects.loc[projects["project_id"] == 339].index
    projects.loc[t339_mask, ["border", "bus1"]] = projects.loc[
        t339_mask, ["border", "bus1"]
    ].replace(
        {
            "border": {"ITCS-ITSI": "ITCS-ITVI", "ITSA-ITSI": "ITSA-ITVI"},
            "bus1": {"ITSI": "ITVI"},
        }
    )

    # Clean the project list by removing unclear border and projects with no capacity
    projects = remove_unclear_border(projects, existing_buses)
    projects = remove_no_capacity(projects)

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


def extract_custom_transmission_projects(
    custom_transmission_path: Path, existing_buses: pd.Index
) -> pd.DataFrame:
    """
    Extract custom transmission projects.

    Parameters
    ----------
    custom_transmission_path : Path
        File path to custom transmission projects.
    existing_buses : pd.Index
        List of existing buses.

    Returns
    -------
    pd.DataFrame
        Curated list of custom projects.
    """
    custom_transmission_projects = (
        pd.read_csv(custom_transmission_path)
        .assign(border=lambda df: df.bus0 + "-" + df.bus1)
        .drop(["source", "further description"], axis=1)
    )

    custom_transmission_projects = remove_unclear_border(
        custom_transmission_projects, existing_buses
    )
    custom_transmission_projects = remove_no_capacity(custom_transmission_projects)

    return custom_transmission_projects


def extract_investment_attributes(excel_path: Path) -> pd.DataFrame:
    """
    Extract length, CAPEX, and underwater fraction from Trans.Investments sheet.

    Aggregates investment-level data to the project level by summing route
    lengths and CAPEX, and computing the underwater fraction from offshore
    cable lengths.
    """
    inv = pd.read_excel(
        excel_path,
        sheet_name="Trans.Investments",
        skiprows=1,
        usecols=[
            "This investment belongs to project number…",
            "Total route length (km)",
            "Estimated CAPEX (MEUR)",
            "Type of Element",
        ],
    ).rename(
        columns={
            "This investment belongs to project number…": "project_id",
            "Total route length (km)": "length_km",
            "Estimated CAPEX (MEUR)": "capex_meur",
            "Type of Element": "element_type",
        }
    )

    is_offshore = inv["element_type"].isin(OFFSHORE_ELEMENT_TYPES)

    agg = inv.groupby("project_id").agg(
        length_km=("length_km", "sum"),
        capex_meur=("capex_meur", "sum"),
    )
    offshore_km = inv.loc[is_offshore].groupby("project_id")["length_km"].sum()
    agg["underwater_fraction"] = (offshore_km / agg["length_km"]).fillna(0).round(3)

    return agg


def overwrite_projects(
    projects: pd.DataFrame, custom_projects: pd.DataFrame
) -> pd.DataFrame:
    """
    Apply custom project modifications to the list of projects.

    Parameters
    ----------
    projects : pd.DataFrame
        Base list of projects.
    custom_projects : pd.DataFrame
        Custom project modifications.

    Returns
    -------
    pd.DataFrame
        Updated list of projects with custom project modifications applied (if applicable).
    """
    if custom_projects.empty:
        return projects

    idx = ["project_id", "bus0", "bus1"]
    custom_isnull = custom_projects[idx].isnull().any(axis=1)
    if custom_isnull.any():
        malformed = custom_projects.loc[custom_isnull, idx].to_string(index=False)
        raise ValueError(
            f"Custom projects must define project_id, bus0 and bus1, but the "
            f"following rows have missing values:\n{malformed}"
        )

    custom_projects = custom_projects.set_index(idx, verify_integrity=True).sort_index()
    projects = projects.set_index(idx).sort_index()

    # Identify existing and new projects
    new_projects = custom_projects.index.difference(projects.index)
    existing_projects = custom_projects.index.intersection(projects.index)

    # Fill missing values with existing projects
    custom_projects = custom_projects.reindex(columns=projects.columns).fillna(projects)
    custom_projects["is_crossborder"] = (
        custom_projects["is_crossborder"].fillna(True).astype(bool)
    )
    custom_projects = custom_projects.fillna(0).infer_objects(copy=False)

    # Overwrite unique pairs of (project_id, bus0, bus1)
    projects.loc[existing_projects] = custom_projects.loc[existing_projects]

    # Add projects that don't already exist
    if len(new_projects) > 0:
        projects = pd.concat([projects, custom_projects.loc[new_projects]])

    return projects.reset_index()


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


def normalize_yes_no(value: str) -> str:
    return str(value).strip().lower()


def compute_method(flag: str) -> str:
    return "TOOT" if flag == "yes" else "PINT"


def build_method_assignments(
    guidelines: pd.DataFrame, projects: pd.DataFrame
) -> pd.DataFrame:
    guidelines = guidelines.rename(
        columns={
            "ID": "project_id",
            "Project_name": "project_name",
            "In_ref_grid_2030": "in_ref_2030",
            "In_ref_grid_2040": "in_ref_2040",
        }
    )

    for col in ["in_ref_2030", "in_ref_2040"]:
        if col in guidelines.columns:
            guidelines[col] = guidelines[col].map(normalize_yes_no)

    base = guidelines[["project_id", "project_name", "in_ref_2030", "in_ref_2040"]]
    base = base.dropna(subset=["project_id"])

    agg = base.groupby("project_id", as_index=False).agg(
        project_name=("project_name", "first"),
        in_ref_2030=("in_ref_2030", lambda s: "yes" if (s == "yes").any() else "no"),
        in_ref_2040=("in_ref_2040", lambda s: "yes" if (s == "yes").any() else "no"),
    )

    all_project_ids = projects["project_id"].unique()
    assigned = []
    for horizon, col in [(2030, "in_ref_2030"), (2040, "in_ref_2040")]:
        rows = agg[["project_id", "in_ref_2030", "in_ref_2040"]].copy()
        rows["planning_horizon"] = horizon
        rows["method"] = rows[col].map(compute_method)
        rows = rows.rename(
            columns={
                "in_ref_2030": "in_ref_grid_2030",
                "in_ref_2040": "in_ref_grid_2040",
            }
        )

        missing_ids = set(all_project_ids) - set(rows["project_id"])
        if missing_ids:
            missing_rows = pd.DataFrame(
                {
                    "project_id": list(missing_ids),
                    "in_ref_grid_2030": "no",
                    "in_ref_grid_2040": "no",
                    "planning_horizon": horizon,
                    "method": "PINT",
                }
            )
            rows = pd.concat([rows, missing_rows], ignore_index=True)
        assigned.append(rows)

    assigned = pd.concat(assigned, ignore_index=True)
    return projects.merge(assigned, on="project_id", how="left")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_projects", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    transmission_path = Path(snakemake.input.dir) / "20250312_export_transmission.xlsx"
    storage_path = Path(snakemake.input.dir) / "20250312_export_storage.xlsx"
    custom_transmission_path = Path(snakemake.input.custom_transmission)

    existing_buses = read_tyndp_electricity_buses(snakemake.input.buses)

    # Transmission projects
    transmission_projects = extract_transmission_projects(
        transmission_path, existing_buses
    )
    custom_transmission_projects = extract_custom_transmission_projects(
        custom_transmission_path, existing_buses
    )

    investment_attrs = extract_investment_attributes(transmission_path)
    transmission_projects = transmission_projects.merge(
        investment_attrs, on="project_id", how="left"
    )

    transmission_projects = overwrite_projects(
        transmission_projects, custom_transmission_projects
    )

    transmission_projects.to_csv(snakemake.output.transmission_projects, index=False)

    # Storage projects
    storage_projects = extract_storage_projects(storage_path, existing_buses)
    # TODO Add overwrite_projects for storages
    storage_projects.to_csv(snakemake.output.storage_projects, index=False)

    # Method definition (PINT / TOOT)
    guidelines = pd.read_csv(snakemake.input.guidelines)
    methods = build_method_assignments(guidelines, transmission_projects)
    methods.to_csv(snakemake.output.methods, index=False)
