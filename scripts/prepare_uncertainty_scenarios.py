# SPDX-FileCopyrightText: Contributors to NGV-IEM project
#
# SPDX-License-Identifier: MIT
"""
Uses an already solved network and prepares it for uncertainty analysis by:
* Copying optimised capacities (`p_nom_opt`, `e_nom_opt`, ...) to nominal capacities (`p_nom`, `e_nom`, ...)
* Setting capacity extendable flags to False
* Modifying the network according to the uncertainty scenario, e.g., changing the demand or availability of renewables
"""

import logging

import pandas as pd
import pypsa
import numpy as np

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

def add_electrolysis_constraints(n):
    "Enforce the electrolysis dispatch to the optimal dispatch found in the solved network."
    electrolysis_i = n.links[n.links.carrier=="H2 Electrolysis"].index
    n.links_t.p_set.loc[:,electrolysis_i] = n.links_t.p0.loc[:,electrolysis_i]
    return n

def remove_components_added_in_solve_network_py(n: pypsa.Network) -> pypsa.Network:
    """Removes components that were added in solve_network.py; we're planing on running this network through the same step again and want to avoid adding the components again."""

    logger.info("Removing components added in solve_network.py")

    # These components are not always part of the network, so
    # we check for their existence first
    if "co2_sequestration_limit" in n.global_constraints.index:
        n.remove(
            class_name="GlobalConstraint",
            name="co2_sequestration_limit",
        )

    if "load" in n.carriers.index:
        n.remove(
            class_name="Carrier",
            name="load",
        )
        gens_i = n.generators.query("`name`.str.endswith(' load')").index
        n.remove(
            class_name="Generator",
            name=gens_i,
        )

    if "curtailment" in n.carriers.index:
        n.remove(
            class_name="Carrier",
            name="curtailment",
        )
        gens_i = n.generators.query("`name`.str.endswith(' curtailment')").index
        n.remove(
            class_name="Generator",
            name=gens_i,
        )

    return n

def extend_primary_fuel_sources(n):
    primary_fuel_sources = ['EU lignite', 'EU coal', 'EU oil primary', 'EU uranium', 'EU gas']
    n.generators.loc[primary_fuel_sources,'p_nom_extendable'] = True
    return n

def add_scenario_uncertainty(
    n: pypsa.Network, scenario_name: str, error_fp: str = None
) -> pypsa.Network:
    """
    Creates a new network that is modified according to the uncertainty scenario.

    Parameters
    ----------
    n : pypsa.Network
        The input network used as a base for the uncertainty scenario.
        This network will not be modified, a modified copy is returned.
    scenario_name : str
        The name of the uncertainty scenario to apply. Currently supported:
        - "trader-forecast": Applies forecast errors to renewable generation availability factors from file
    error_fp : str, optional
        File path to the parquet file containing the relative errors for the "trader-forecast" scenario.
        Required if scenario_name is "trader-forecast".

    Returns
    -------
    pypsa.Network
        A new network modified according to the specified uncertainty scenario.
    """

    # Work on a copy of the network
    n = n.copy()

    # Change name
    n.name = f"{n.name} (sensitivity: {scenario_name})"

    # Use errors from file
    if scenario_name == "trader-forecast":
        relative_errors = pd.read_parquet(error_fp)

        # Datetimeindex contains one entry too many (8761 instead of 8760), remove the last
        relative_errors = relative_errors.iloc[:8760]

        # Realign the datetime index to be for 2009 (hourly)
        relative_errors.index = pd.date_range(
            start="2009-01-01", periods=8760, freq="h"
        )

        # Rename index name to 'snapshot' for consistency with PyPSA
        relative_errors.index.name = "snapshot"

        # Manually curated mapping between bidding zones used in the error data
        # and the node names used in the TYNDP model
        bz_mapping = {
            # "AL": None, # missing
            # "BA": None, # missing
            # "BG": None, # missing
            # "EE": None, # missing
            # "GR": None, # missing
            # "HU": None, # missing
            # "ITVI": None, # missing
            # "LT": None, # missing
            # "LU": None, # missing
            # "LV": None, # missing
            # "ME": None, # missing
            # "MK": None, # missing
            # "MT": None, # missing
            # "NO2": None, # not mapped
            # "RO": None, # missing
            # "RS": None, # missing
            # "SI": None, # missing
            # "SK": None, # missing
            "AT": ["AT00"],
            "BE": ["BEIOH01", "BEOH001", "BE00"],
            "CH": ["CH00"],
            "CZ": ["CZ00"],
            "DE": ["DEOH002", "DEOR001", "DE00", "DEOH001"],
            "DK1": ["DKE1", "DKE0R01"],
            "DK2": ["DKW1", "DKW0R01"],
            "ES": ["ES00", "ESOH001", "ESOH003", "ESOH002"],
            "FI": ["FI00", "FIOH001"],
            "FR": ["FR00", "FR15", "FROH001", "FROH002", "FROH003"],
            "GB": [
                "GB00",
                "GBNI",
                "GBNIOR1",
                "GBOH001",
                "GBOH002",
                "GBOH003",
                "GBOH004",
                "GBOH005",
                "GBOH006",
                "GBOR001",
            ],
            "HR": ["HR00", "HROH001"],
            "IE": ["IEOH001", "IE00"],
            "ITCA": ["ITCA"],
            "ITCN": ["ITCN", "ITCNOR1"],
            "ITCS": ["ITCS", "ITCSOR1"],
            "ITN1": ["ITN1", "ITN1OR1"],
            "ITS1": ["ITS1", "ITS1OR1"],
            "ITSA": ["ITSA", "ITSAOR1"],
            "ITSI": ["ITSI", "ITSIOR1"],
            "NL": ["NL00", "NLOR001", "NLOH001"],
            "NO1": ["NOS0", "NOSOH01", "NOSOH02", "NOSOR01"],
            "NO3": ["NOM1", "NOMOH01"],
            "NO4": ["NON1", "NONOH01"],
            "PL": ["PL00", "PLOH001"],
            "PT": ["PT00", "PTOH001"],
            "SE1": ["SE01", "SEOR001"],
            "SE2": ["SE02", "SEOH001"],
            "SE3": ["SE03", "SEOH002"],
            "SE4": ["SE04", "SEOH003"],
        }

        ## Restructure the error dataframe to match the nodes names and component names from the TYNDP model
        # Duplicate the error data for each mapped node
        expanded_errors_l: list = []
        for bz, nodes in bz_mapping.items():
            bz_cols = [col for col in relative_errors.columns if col.startswith(bz)]

            for node in nodes:
                node_cols = [
                    col.replace(bz, node) for col in bz_cols if node is not None
                ]
                expanded_errors_l.append(
                    relative_errors[bz_cols].rename(
                        columns=dict(zip(bz_cols, node_cols))
                    )
                )

        expanded_errors: pd.DataFrame = pd.concat(expanded_errors_l, axis=1)

        # Expand errors per technology, duplicating the dataframes
        tech_mapping = {
            "_Wind": ["|onwind", "|offwind"],
            "_Solar": ["|solar-pv-utility", "|solar-pv-rooftop"],
            "_Load": ["|load"],
        }
        expanded_errors_l: list = []
        for suffix_old, suffixes_new in tech_mapping.items():
            tech_cols = [col for col in expanded_errors.columns if suffix_old in col]
            for suffix_new in suffixes_new:
                new_cols = {
                    col: col.replace(suffix_old, suffix_new) for col in tech_cols
                }
                expanded_errors_l.append(
                    expanded_errors[tech_cols].rename(columns=new_cols)
                )

        expanded_errors: pd.DataFrame = pd.concat(expanded_errors_l, axis="columns")
        # Take the column names, split them on " " and turn the split into a multiindex
        multiindex_tuples = [col.split("|") for col in expanded_errors.columns]
        expanded_errors.columns = pd.MultiIndex.from_tuples(
            multiindex_tuples, names=["bus", "technology"]
        )

        ## Combining the errors with the time-series data from the TYNDP model
        # Generators
        for bus, tech in expanded_errors.columns:
            if "load" in tech.casefold():
                comp = n.components["loads"].dynamic["p_set"]
                cols = [col for col in comp.columns if col == bus]
            else:
                comp = n.components["generators"].dynamic["p_max_pu"]
                cols = [
                    col for col in comp.columns if col.startswith(bus) and tech in col
                ]

            if not cols:
                continue
            logger.info(f"Applying errors to {bus} {tech} for columns: {cols}")

            # Apply the errors onto all columns from generators[col]
            new_p_max_pu = comp[cols].multiply(
                1 + expanded_errors.loc[:, (bus, tech)], axis="index"
            )

            # Errors may cause values below which is unrealistic, so clip accordingly
            # We could also clip > 1, but then we need to differentiate between
            # loads (absolute timeseries) and generators (pu timeseries)
            new_p_max_pu = new_p_max_pu.clip(lower=0) #, upper=max_value)

            # Assign the new values back to the generators dataframe
            # (this propagates to the network object n because it is a reference, not a copy)
            # Make sure to align snapshots first
            new_p_max_pu = new_p_max_pu.loc[comp.index]
            comp[cols] = new_p_max_pu
    return n

#%%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_uncertainty_scenarios",
            opts="",
            clusters="all",
            configfiles="config/config.tyndp.yaml",
            sector_opts="",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    
    n.optimize.fix_optimal_capacities()               
    n = remove_components_added_in_solve_network_py(n)
    n = add_electrolysis_constraints(n)
    n = extend_primary_fuel_sources(n)

    for uncertainty_scenario in snakemake.params.uncertainty_scenarios:
        n_uncertain = add_scenario_uncertainty(
            n, scenario_name=uncertainty_scenario, error_fp=snakemake.input.errors
        )

        output_path = [
            p for p in snakemake.output.networks if f"{uncertainty_scenario}" in p
        ]
        assert len(output_path) == 1, (
            "Expected exactly one output path for each uncertainty scenario"
        )
        output_path = output_path[0]

        # Better now than later
        n.consistency_check()
        n_uncertain.consistency_check()

        logger.info(f"Saving uncertainty scenario network to {output_path}")
        n_uncertain.optimize.create_model()
        n_uncertain.model.to_file(output_path[:-3] + '.lp', explicit_coordinate_names=True)
        n_uncertain.export_to_netcdf(output_path)
