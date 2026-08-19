# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extracts and cleans CBA transmission and storage projects from Excel exports.

Reads transmission projects from the "Trans.Projects" sheet of the CBA projects Excel file.
For projects with multiple borders (newline-separated in the Excel), the script explodes
these into separate rows, creating one row per border. Bus codes are extracted from the
border strings (expected format: "BUS0-BUS1") and rows are filtered out with a warning
when the format does not match, when either bus is absent from the TYNDP node list, or
when no capacity is reported in either direction.

Reads storage projects from the "Stor.Projects" sheet. Each project's country is mapped
to a representative electricity bus. Storage projects are assigned TOOT or PINT based on
whether they are already part of the reference grid, same as transmission projects; TOOT
storage projects are not yet supported downstream (see `prepare_project.py`).

Custom PINT transmission projects can be configured using `data/custom_cba_transmission_projects.csv`. With it,
the user can modify existing projects and add new ones. Transmission capacities are in MW.

- Using an existing PINT combination (`project_id`, `bus0`, `bus1`), the user can overwrite any
  field of an existing project with a custom value. Not all fields need to be specified; leaving
  a field empty keeps its existing value.

- New projects must use a new `project_id` and are added as PINT links.

**Inputs**

- `data/tyndp_2024_bundle/cba_projects/20250312_export_transmission.xlsx`: Excel file containing CBA transmission projects
- `data/tyndp_2024_bundle/cba_projects/20250312_export_storage.xlsx`: Excel file containing CBA storage projects
- `rules.retrieve_tyndp.output.nodes`: TYNDP electricity node list used to validate borders
- `rules.retrieve_cba_guidelines_reference_projects.output.file`: Table of projects as defined in the Implementation Guidelines Appendix B.1
- `data/cba_guidelines_reference_projects/.../table_B1_CBA_Implementations_Guidelines_TYNDP2024.csv`: CBA guidelines reference table, used to assign the TOOT/PINT method per project and planning horizon
- `data/cba/cba_project_corrections.csv`: Manually curated bus0/bus1/p_nom corrections for CBA projects, applied in place of the corresponding raw Excel entries
- `data/custom_cba_transmission_projects.csv`: File used to configure custom transmission projects. With it, the user can modify existing projects and add new ones.

**Outputs**

- `resources/cba/transmission_projects.csv`: Cleaned CSV with columns:
  - `project_id`: Integer project identifier
  - `project_name`: Project name
  - `is_crossborder`: Whether the project is reported as cross-border
  - `border`: Border string in format "BUS0-BUS1"
  - `p_nom 0->1`: Transfer capacity increase from bus0 to bus1 (MW)
  - `p_nom 1->0`: Transfer capacity increase from bus1 to bus0 (MW)
  - `bus0`: Source bus code (4 alphanumeric characters)
  - `bus1`: Destination bus code (4 alphanumeric characters)
  - `length_km`: Total route length in km (from Trans.Investments)
  - `capex_meur`: Total estimated CAPEX in MEUR (from Trans.Investments)
  - `underwater_fraction`: Fraction of route that is offshore cable

- `resources/cba/storage_projects.csv`: Cleaned CSV with columns:
  - `project_id`: Integer project identifier
  - `project_name`: Project name
  - `capex_meur`: Estimated total CAPEX in MEUR
  - `opex_meur_per_year`: Estimated annual OPEX in MEUR/year
  - `p_nom_discharge`: Total turbining/discharge capacity (MW)
  - `e_nom_gwh`: Storage capacity (GWh)
  - `p_nom_charge`: Total pumping/charge capacity (MW)
  - `roundtrip_efficiency`: Roundtrip efficiency (fraction)
  - `lifetime_years`: Operational lifetime (years); missing or zero values are replaced with `cba.storage.default_lifetime`
  - `carrier`: PyPSA carrier name for the storage technology
  - `bus`: Representative electricity bus for the project's country

- `resources/cba/cba_project_methods.csv`: Table defining the assignment method of each project.

"""

import logging
from pathlib import Path

import country_converter as coco
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_tyndp_network import AC_VIRTUAL_NODES_IT

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

STORAGE_PROJECTS_COLUMN_MAP = {
    "Project ID": "project_id",
    "Project Name": "project_name",
    "Country": "country",
    "Storage technology": "technology",
    "Estimated CAPEX (MEUR)": "capex_meur",
    "Estimated OPEX (MEUR/year)": "opex_meur_per_year",
    "Total turbining capacity (MW)": "p_nom_discharge",
    "Total Pumping Capacity of Pump Storage (MW)": "p_nom_charge",
    "Storage Capacity (GWh)": "e_nom_gwh",
    "Roundtrip efficiency (%)": "roundtrip_efficiency",
    "Operational lifetime (Years)": "lifetime_years",
    "Is the project in the reference grid 2030?": "in_ref_grid_2030",
    "Is the project in the reference grid 2035?": "in_ref_grid_2035",
}


def load_storage_technology_carrier_map(carrier_mapping_fn: str) -> pd.Series:
    """
    Load the storage technology -> carrier mapping from the technology map CSV.

    Parameters
    ----------
    carrier_mapping_fn : str
        Path to `data/tyndp_technology_map.csv`.

    Returns
    -------
    pd.Series
        Mapping from `cba_storage_technology` (as reported in the CBA storage
        projects Excel export) to `open_tyndp_carrier`.
    """
    tech_map = pd.read_csv(carrier_mapping_fn)
    tech_map = tech_map.dropna(subset=["cba_storage_technology"])
    return tech_map.set_index("cba_storage_technology")["open_tyndp_carrier"]


def read_tyndp_electricity_buses(
    buses_fn: str, col_name: str, virtual_buses: list[str] | None = None
) -> pd.Index:
    """
    Read node list for electricity from tyndp data input.

    Parameters
    ----------
    buses_fn : str
        Path to a TYNDP node list Excel file ("LIST OF NODES.xlsx"
        or offshore hubs "NODE.xlsx").
    col_name : str
        Column which is selected from the dataframe.
    virtual_buses : list
        List of virtual buses to add, not present in the raw node list
        (e.g. added later in build_tyndp_network.py).

    Returns
    -------
    buses : pandas.Index
        Index of electricity buses as used in Open-TYNDP.

    See Also
    --------
    build_tyndp_network.py : build_buses
    build_tyndp_offshore_hubs.py : load_offshore_hubs
    """
    virtual_buses = virtual_buses if virtual_buses is not None else []

    buses = pd.read_excel(buses_fn).replace("UK", "GB", regex=True).set_index(col_name)

    if "OFFSHORE_NODE_TYPE" in buses.columns:
        # drop radial offshore nodes, which are not built as hub buses
        buses = buses[buses.OFFSHORE_NODE_TYPE != "Radial"]

    return buses.index.union(virtual_buses)


def get_existing_buses(buses_fn: str, offshore_buses_fn: str | list | bool) -> pd.Index:
    """
    Return the electricity buses used to validate CBA project borders.
    Combines onshore buses with offshore hub buses, if provided, into a
    single index of existing bus names.

    Parameters
    ----------
    buses_fn : str
        Path to the file defining TYNDP electricity buses.
    offshore_buses_fn : str | list | bool
        Path(s) to the file(s) defining offshore hub buses. If falsy, only
        onshore buses are returned.

    Returns
    -------
    pd.Index
        Existing bus names, combining onshore and offshore buses.
    """
    virtual_buses = [x for pair in AC_VIRTUAL_NODES_IT.items() for x in pair]
    existing_buses = read_tyndp_electricity_buses(
        buses_fn, col_name="NODE", virtual_buses=virtual_buses
    )

    if offshore_buses_fn:
        existing_oh_buses = read_tyndp_electricity_buses(
            offshore_buses_fn, col_name="OFFSHORE_NODE"
        )
        existing_buses = existing_buses.union(existing_oh_buses)

    return existing_buses.rename("bus_ID")


def apply_cba_project_corrections(
    corrections_path: Path, projects: pd.DataFrame
) -> pd.DataFrame:
    """
    Replace bus0/bus1/p_nom of select projects with manually curated corrections.

    Parameters
    ----------
    corrections_path : Path
        Path to the file containing manual corrections.
    projects : pd.DataFrame
        List of transmission projects with their detailed characteristics.

    Returns
    -------
    pd.DataFrame
        List of transmission projects with corrected bus0, bus1 and p_nom values
        for the projects covered by the corrections.
    """

    # Read in CBA project corrections
    corrections = pd.read_csv(corrections_path)
    corrected_ids = corrections["project_id"].unique()

    logger.info(
        "Applying CBA project corrections for %d projects with project ID:\n%s",
        len(corrected_ids),
        ", ".join(corrected_ids.astype(str)),
    )

    # Add border column to corrections
    corrections["border"] = corrections.bus0 + "-" + corrections.bus1

    # Drop projects which are in corrections
    projects_reduced = projects[~projects.project_id.isin(corrected_ids)]

    # Add corrections
    return pd.concat(
        [projects_reduced, corrections.drop("notes", axis=1)], ignore_index=True
    ).sort_values(by="project_id")


def remove_unclear_border(
    projects: pd.DataFrame, existing_buses: pd.Index
) -> pd.DataFrame:
    """
    Remove projects defined by unclear borders from the list of projects.

    Parameters
    ----------
    projects : pd.DataFrame
        List of projects to assess.
    existing_buses : pd.Index
        List of existing buses.

    Returns
    -------
    pd.DataFrame
        Curated list of projects that only use existing buses.
    """
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
        Curated list of projects with a defined capacity.
    """
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
    transmission_path: Path, corrections_path: Path, existing_buses: pd.Index
) -> pd.DataFrame:
    """
    Read and clean the transmission projects from the "Trans.Projects" sheet.

    Projects reporting several expected capacity increases are exploded into one row per
    border. Rows are dropped when the border cannot be parsed or when no capacity is
    reported in either direction.

    Parameters
    ----------
    transmission_path : Path
        Path to the Excel export defining the transmission projects.
    corrections_path : Path
        Path to the file containing manual corrections.
    existing_buses : pd.Index
        Electricity buses as used in Open-TYNDP.

    Returns
    -------
    pd.DataFrame
        List of projects with their detailed characteristics. One row per project and border.
    """
    projects = (
        pd.read_excel(
            transmission_path,
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

    # Apply manual CBA project corrections
    projects = apply_cba_project_corrections(corrections_path, projects)

    # Clean the project list by removing projects with unclear borders or no capacity
    projects = remove_unclear_border(projects, existing_buses)
    projects = remove_no_capacity(projects)

    # Several projects have capacities with "Up to ..."
    cols = ["p_nom 0->1", "p_nom 1->0"]
    up_to_projects = set()
    for col in cols:
        up_to = projects[col].str.startswith("Up to ", na=False)
        if up_to.any():
            projects.loc[up_to, col] = projects.loc[up_to, col].str[len("Up to ") :]
            up_to_projects.update(projects.loc[up_to, "project_name"])
    if up_to_projects:
        logger.info(
            f"Removed 'Up to ' capacity prefix from {len(up_to_projects)} projects:\n"
            + ", ".join(up_to_projects)
        )

    # convert to numeric
    projects[cols] = projects[cols].apply(pd.to_numeric, errors="coerce")

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
        pd.read_csv(
            custom_transmission_path,
        )
        .assign(border=lambda df: df.bus0 + "-" + df.bus1)
        .drop(["source", "further description"], axis=1, errors="ignore")
    )

    # Remove unclear borders
    custom_transmission_projects = remove_unclear_border(
        custom_transmission_projects, existing_buses
    )

    # Remove projects without capacity
    custom_transmission_projects = remove_no_capacity(custom_transmission_projects)

    # Remove projects without project_id
    mask_null = custom_transmission_projects.project_id.notnull()
    null_projects = custom_transmission_projects[~mask_null]
    if not null_projects.empty:
        logger.warning(
            f"Some custom projects have no project_id (mandatory field), ignoring them:\n{null_projects.to_string(index=False)}"
        )
    custom_transmission_projects = custom_transmission_projects[mask_null]
    custom_transmission_projects["project_id"] = custom_transmission_projects[
        "project_id"
    ].astype(int)

    # Remove projects with two identical buses
    mask_dup_buses = (
        custom_transmission_projects.bus0 != custom_transmission_projects.bus1
    )
    dup_projects = custom_transmission_projects[~mask_dup_buses]
    if not dup_projects.empty:
        logger.warning(
            f"Some custom projects have identical bus0 and bus1, ignoring them:\n{dup_projects.to_string(index=False)}"
        )
    custom_transmission_projects = custom_transmission_projects[mask_dup_buses]

    return custom_transmission_projects


def extract_investment_attributes(transmission_path: Path) -> pd.DataFrame:
    """
    Extract length, CAPEX, and underwater fraction from Trans.Investments sheet.

    Aggregates investment-level data to the project level by summing route
    lengths and CAPEX, and computing the underwater fraction from offshore
    cable lengths.

    Parameters
    ----------
    transmission_path : Path
       Path to the Excel export defining the transmission projects and their investment attributes.

    Returns
    -------
    pd.DataFrame
       Route length, CAPEX and underwater fraction per project, indexed by ``project_id``.

    """
    inv = pd.read_excel(
        transmission_path,
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


def assign_country_bus(countries: pd.Series, existing_buses: pd.Index) -> pd.Series:
    """
    Assign each project's country name to a representative electrical bus.

    Country names (e.g. "Austria") are converted to ISO2 codes and matched
    against the bus prefix (e.g. "AT00"). If a country has several buses,
    the alphabetically first one is used and a warning is logged.

    Parameters
    ----------
    countries : pd.Series
        Country names as given in the storage projects Excel sheet.
    existing_buses : pd.Index
        Electricity buses as used in Open-TYNDP.

    Returns
    -------
    pd.Series
        Bus code per project, aligned with ``countries``.
    """
    iso2 = pd.Index(coco.convert(countries.unique().tolist(), to="ISO2"))
    country_to_iso2 = dict(zip(countries.unique(), iso2))

    buses_by_country = existing_buses.to_series().groupby(existing_buses.str[:2])
    representative_bus = buses_by_country.first()

    multi_bus_countries = buses_by_country.nunique().loc[lambda s: s > 1]
    if not multi_bus_countries.empty:
        logger.warning(
            "%d countries have several electricity buses; assigning storage "
            "projects to the first bus alphabetically:\n%s",
            len(multi_bus_countries),
            representative_bus.loc[multi_bus_countries.index].to_string(),
        )

    bus = countries.map(country_to_iso2).map(representative_bus)
    unmatched = bus.isna()
    if unmatched.any():
        logger.warning(
            "%d storage projects could not be matched to a bus for countries:\n%s",
            unmatched.sum(),
            countries.loc[unmatched].unique(),
        )
    return bus


def extract_storage_projects(
    excel_path: Path,
    existing_buses: pd.Index,
    default_lifetime: float,
    technology_carrier_map: pd.Series,
) -> pd.DataFrame:
    """
    Read and clean the storage projects from the "Stor.Projects" sheet.

    Each project's country is mapped to a representative electricity bus.
    Rows with an unknown technology or that cannot be matched to a bus are
    dropped with a warning. Projects with a missing or zero operational
    lifetime are overwritten with ``default_lifetime``.

    Parameters
    ----------
    excel_path : Path
        Path to the Excel export defining the storage projects.
    existing_buses : pd.Index
        Electricity buses as used in Open-TYNDP.
    default_lifetime : float
        Lifetime (years) used to replace missing or zero operational
        lifetimes reported in the Excel export.
    technology_carrier_map : pd.Series
        Mapping from `cba_storage_technology` to `open_tyndp_carrier`, as
        returned by `load_storage_technology_carrier_map`.

    Returns
    -------
    pd.DataFrame
        List of storage projects with their detailed characteristics, one row per project.
    """
    projects = (
        pd.read_excel(
            excel_path,
            sheet_name="Stor.Projects",
            skiprows=1,
            usecols=list(STORAGE_PROJECTS_COLUMN_MAP),
        )
        .rename(columns=STORAGE_PROJECTS_COLUMN_MAP)
        .dropna(subset=["project_id"])
    )
    projects["project_id"] = projects["project_id"].astype(int)
    projects["roundtrip_efficiency"] = projects["roundtrip_efficiency"] / 100.0

    unknown_technology = ~projects["technology"].isin(technology_carrier_map.index)
    if unknown_technology.any():
        logger.warning(
            "%d out of %d storage projects have an unknown storage technology, ignoring them:\n%s",
            unknown_technology.sum(),
            len(projects),
            projects.loc[
                unknown_technology, ["project_id", "project_name", "technology"]
            ].to_string(index=False, max_colwidth=40),
        )
    projects = projects.loc[~unknown_technology]
    projects["carrier"] = projects["technology"].map(technology_carrier_map)

    projects["bus"] = assign_country_bus(projects["country"], existing_buses)
    projects = projects.dropna(subset=["bus"])

    # Data quirk: project 1064 reports a negative pumping capacity
    projects["p_nom_charge"] = projects["p_nom_charge"].abs()

    # Data quirk: some projects report a missing or zero operational lifetime
    # (e.g. project 1076); fall back to default_lifetime in that case
    missing_lifetime = projects["lifetime_years"].isna() | (
        projects["lifetime_years"] <= 0
    )
    projects.loc[missing_lifetime, "lifetime_years"] = default_lifetime

    return projects.drop(columns=["technology", "country"])


# The Stor.Projects sheet only reports reference-grid flags for 2030 and 2035;
# 2035 is used as the nearest available proxy for the 2040 planning horizon.
STORAGE_REF_GRID_HORIZON_COLUMN = {2030: "in_ref_grid_2030", 2040: "in_ref_grid_2035"}


def build_storage_method_assignments(
    storage_projects: pd.DataFrame, planning_horizons: list[int]
) -> pd.DataFrame:
    """
    Define the assignment method of storage projects.

    A storage project is assigned TOOT if it is reported as already part of
    the reference grid for the given planning horizon, and PINT otherwise
    (i.e. it represents a new capacity addition). The reference-grid flag is
    only reported for 2030 and 2035 in the Excel; 2035 is used as a proxy for
    the 2040 planning horizon.

    Parameters
    ----------
    storage_projects : pd.DataFrame
        List of storage projects with their detailed characteristics,
        including the in_ref_grid_2030/in_ref_grid_2035 boolean columns.
    planning_horizons : list[int]
        Planning horizons for which to assign a method.

    Returns
    -------
    pd.DataFrame
        Table defining the assignment method of each storage project.
    """
    rows = []
    for horizon in planning_horizons:
        rows.append(
            storage_projects[["project_id", "project_name"]].assign(
                planning_horizon=horizon,
                method=storage_projects[STORAGE_REF_GRID_HORIZON_COLUMN[horizon]].map(
                    {True: "toot", False: "pint"}
                ),
            )
        )
    methods = pd.concat(rows, ignore_index=True)
    methods["project_type"] = "storage"
    return methods


def split_investment_attributes_per_line(
    investment_attrs: pd.DataFrame, transmission_projects: pd.DataFrame
) -> pd.DataFrame:
    """
    Split investment costs and length evenly across transmission lines.

    Investment costs and length are given per project and not per
    transmission line, therefore these attributes need to be split before
    merging.

    Parameters
    ----------
    investment_attrs : pd.DataFrame
        Investment attributes indexed by project_id.
    transmission_projects : pd.DataFrame
        Transmission projects with a project_id column.

    Returns
    -------
    pd.DataFrame
        investment_attrs with length_km and capex_meur divided by the number
        of lines per project.
    """
    link_counts = transmission_projects.groupby("project_id").size()
    return investment_attrs.assign(
        length_km=lambda d: d.length_km / d.index.map(link_counts).fillna(1),
        capex_meur=lambda d: d.capex_meur / d.index.map(link_counts).fillna(1),
    )


def normalize_yes_no(value: str) -> str:
    return str(value).strip().lower()


def compute_method(flag: str) -> str:
    return "toot" if flag == "yes" else "pint"


def build_method_assignments(
    guidelines_fn: str,
    projects: pd.DataFrame,
    custom_transmission_projects: pd.DataFrame,
) -> pd.DataFrame:
    """
    Determine the CBA assessment method for each project. The method is PINT (default) or TOOT and
    depends on the planning horizon (2030 or 2040).

    Parameters
    ----------
    guidelines_fn : str
        Path to the table of projects as defined in the Implementation Guidelines Appendix B.1.
    projects: pd.DataFrame
        List of projects with their detailed characteristics.
    custom_transmission_projects: pd.DataFrame
        List of custom projects with their detailed characteristics.

    Returns
    -------
    pd.DataFrame
        Table defining the assignment method of each project.
    """
    guidelines = pd.read_csv(guidelines_fn).rename(
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

    all_project_ids = set(projects["project_id"]).union(
        set(custom_transmission_projects["project_id"])
    )
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

        missing_ids = all_project_ids - set(rows["project_id"])
        if missing_ids:
            missing_rows = pd.DataFrame(
                {
                    "project_id": list(missing_ids),
                    "in_ref_grid_2030": "no",
                    "in_ref_grid_2040": "no",
                    "planning_horizon": horizon,
                    "method": "pint",
                }
            )
            rows = pd.concat([rows, missing_rows], ignore_index=True)
        assigned.append(rows)

    assigned = pd.concat(assigned, ignore_index=True).query(
        "project_id in @projects.project_id or project_id in @custom_transmission_projects.project_id"
    )
    assigned["project_type"] = "transmission"
    return assigned


def apply_custom_projects(
    projects: pd.DataFrame, custom_projects: pd.DataFrame, methods: pd.DataFrame
) -> pd.DataFrame:
    """
    Apply custom PINT project modifications to a list of projects.

    Parameters
    ----------
    projects : pd.DataFrame
        Base list of projects.
    custom_projects : pd.DataFrame
        Custom project modifications.
    methods : pd.DataFrame
        DataFrame of projects with the corresponding method to apply.

    Returns
    -------
    pd.DataFrame
        Updated list of projects with custom project modifications applied (if applicable).
    """

    # Validate inputs
    if custom_projects.empty:
        return projects

    idx = ["project_id", "bus0", "bus1"]

    mask_toot = custom_projects["project_id"].isin(
        methods.query("method=='toot'").project_id
    )
    custom_toot = custom_projects[mask_toot]
    custom_projects = custom_projects[~mask_toot]
    if not custom_toot.empty:
        logger.warning(
            f"Custom projects must refer to PINT projects. The following rows are ignored "
            f"because they refer to TOOT projects in at least one planning horizon:\n"
            f"{custom_toot[idx].to_string(index=False)}"
        )

    def set_valid_index(df: pd.DataFrame, is_custom: bool) -> pd.DataFrame:
        try:
            return df.set_index(idx, verify_integrity=True).sort_index()
        except ValueError as e:
            malformed = df.loc[df[idx].duplicated(), idx].to_string(index=False)
            label = "Custom projects" if is_custom else "Projects"
            raise ValueError(
                f"{label} must define a unique set of project_id, bus0 and bus1, but the "
                f"following rows have duplicated keys:\n{malformed}"
            ) from e

    custom_projects = set_valid_index(custom_projects, is_custom=True)
    projects = set_valid_index(projects, is_custom=False)

    # Identify existing and new projects
    new_projects = custom_projects.index.difference(projects.index)
    existing_projects = custom_projects.index.intersection(projects.index)

    # Fill missing values using existing projects
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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_projects", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # File paths
    transmission_path = Path(snakemake.input.dir, "20250312_export_transmission.xlsx")
    storage_path = Path(snakemake.input.dir, "20250312_export_storage.xlsx")
    custom_transmission_path = Path(snakemake.input.custom_transmission)
    corrections_path = snakemake.input.cba_project_corrections

    # Get existing buses
    existing_buses = get_existing_buses(
        snakemake.input.buses, snakemake.input.offshore_buses
    )

    # Transmission projects
    transmission_projects = extract_transmission_projects(
        transmission_path, corrections_path, existing_buses
    )

    # Custom transmission projects
    custom_transmission_projects = extract_custom_transmission_projects(
        custom_transmission_path, existing_buses
    )

    # Investment costs and length transmission
    investment_attrs = extract_investment_attributes(transmission_path)
    investment_attrs_per_line = split_investment_attributes_per_line(
        investment_attrs, transmission_projects
    )

    transmission_projects = transmission_projects.merge(
        investment_attrs_per_line, on="project_id", how="left"
    )

    # Method definition (PINT / TOOT) for transmission projects
    transmission_methods = build_method_assignments(
        snakemake.input.guidelines, transmission_projects, custom_transmission_projects
    )

    # Apply custom projects
    transmission_projects = apply_custom_projects(
        transmission_projects, custom_transmission_projects, transmission_methods
    )

    transmission_projects.to_csv(snakemake.output.transmission_projects, index=False)

    # Storage projects
    # TODO Add overwrite_projects for storage projects
    technology_carrier_map = load_storage_technology_carrier_map(
        snakemake.input.carrier_mapping
    )
    storage_projects = extract_storage_projects(
        storage_path,
        existing_buses,
        snakemake.params.storage_default_lifetime,
        technology_carrier_map,
    )
    storage_projects.to_csv(snakemake.output.storage_projects, index=False)

    # Method definition (PINT / TOOT) for storage projects
    storage_methods = build_storage_method_assignments(
        storage_projects, snakemake.params.planning_horizons
    )

    methods = pd.concat([transmission_methods, storage_methods], ignore_index=True)
    methods.to_csv(snakemake.output.methods, index=False)
