# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Calculate CBA indicators by comparing reference and project scenarios.

This script computes the B1 indicator (Total System Cost difference) and other
CBA metrics by analyzing the solved networks for reference and project cases.

PINT (Put In at a Time):
    - Reference: Network WITHOUT any projects
    - Project: Network WITH the specific project added
    - B1 = Cost(reference) - Cost(with project)

TOOT (Take Out One at a Time):
    - Reference: Network WITH all projects (current plan)
    - Project: Network WITHOUT the specific project (removed)
    - B1 = Cost(without project) - Cost(reference)


References:
- CBA guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/news/2024/entso-e_4th_CBA_Guideline_240409.pdf
    - section 3.2.2: TOOT and PINT, page 23-24
- CBA implementation guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/tyndp-documents/TYNDP2024/foropinion/CBA_Implementation_Guidelines.pdf
    - section 5.1: B1 - SEW, page 58-59

"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)

INDICATOR_UNITS = {
    "B1_total_system_cost_change": "EUR/year",
    "cost_reference": "EUR/year",
    "capex_reference": "EUR/year",
    "opex_reference": "EUR/year",
    "cost_project": "EUR/year",
    "capex_project": "EUR/year",
    "opex_project": "EUR/year",
    "capex_change": "EUR/year",
    "opex_change": "EUR/year",
    "co2_variation": "t/year",
    "co2_ets_price": "EUR/t",
    "co2_societal_cost": "EUR/t",
    "B2_societal_cost_variation": "EUR/year",
    "B3_res_capacity_change_mw": "MW",
    "B3_res_generation_change_mwh": "MWh/year",
    "B3_annual_avoided_curtailment_mwh": "MWh/year",
    "B4a_nox": "kg/year",
    "B4b_nh3": "kg/year",
    "B4c_sox": "kg/year",
    "B4d_pm25": "kg/year",
    "B4e_pm10": "kg/year",
    "B4f_nmvoc": "kg/year",
}

CARRIER_TO_EMISSION_FACTORS = {
    "gas": ("Gas", "Gas biofuel"),
    "coal": ("Hard coal", "Hard Coal biofuel"),
    "lignite": ("Lignite", "Lignite biofuel"),
    "uranium": ("Nuclear", None),
    "oil": ("Light Oil", "Light Oil biofuel"),  # all oil fuels have same factors
}


def calculate_total_system_cost(n):
    """
    Calculate total annualized system cost using PyPSA built-in statistics.

    This implementation handles:
    - Annualized capital costs
    - Time-aggregated operational costs
    - All component types (generators, links, storage, etc.)

    Args:
        n: PyPSA Network (must be solved)

    Returns:
        float: Total system cost in currency units (Euros)
    """
    if not n.is_solved:
        raise ValueError("Network must be solved before calculating costs")

    # Use PyPSA's built-in statistics methods
    capex = n.statistics.capex().sum()
    opex = n.statistics.opex(aggregate_time="sum").sum()
    total = capex + opex
    return {
        "total": total,
        "capex": capex,
        "opex": opex,
    }


def check_method(method: str) -> str:
    """
    Normalize and validate the CBA method name.

    If the method is not recognized as either "pint" or "toot", a ValueError is raised.
    """
    method = method.lower()
    if method not in ["pint", "toot"]:
        raise ValueError(f"Method must be 'pint' or 'toot', got: {method}")
    return method


def calculate_co2_emissions_per_carrier(n: pypsa.Network) -> float:
    """
    Calculate net CO2 emissions using the final snapshot of the CO2 store.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object for which to calculate CO2 emissions.

    Returns
    -------
    float
        Net CO2 emissions at the final snapshot.
    """
    stores_by_carrier = n.stores_t.e.T.groupby(n.stores.carrier).sum().T
    net_co2 = stores_by_carrier["co2"].iloc[-1]  # get final snapshot value
    return float(net_co2)


def load_non_co2_emission_factors(path: str) -> pd.DataFrame:
    """Load non-CO2 emission factors and compute mean values."""
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df[df["Fuel"].notna()]
    pollutant_cols = {
        "NOX emission factor": "B4a_nox",
        "NH3 emission factor": "B4b_nh3",
        "SO emission factor": "B4c_sox",
        "PM2.5 and smaller emission factor": "B4d_pm25",
        "PM10 emission factor": "B4e_pm10",
        "NMVOC emission factor": "B4f_nmvoc",
    }
    factors = pd.DataFrame(
        {
            key: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            for col, key in pollutant_cols.items()
        }
    )
    factors["Fuel"] = df["Fuel"]
    df_mean = factors.groupby("Fuel")[list(pollutant_cols.values())].mean()
    return df_mean


def calculate_res_capacity_per_carrier(
    n: pypsa.Network, res_carriers: list[str]
) -> pd.Series:
    """
    Calculate total RES capacity per carrier (MW) using optimal capacity.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.

    Returns
    -------
    pandas.Series
        A Series containing RES capacity per carrier (MW).
    """
    capacity = n.statistics.optimal_capacity(nice_names=False).groupby("carrier").sum()
    res_capacity = capacity.reindex(res_carriers).dropna()
    return res_capacity


def calculate_res_generation_per_carrier(
    n: pypsa.Network, res_carriers: list[str]
) -> pd.Series:
    """
    Calculate total RES generation per carrier (MWh/year) using energy balance.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.

    Returns
    -------
    pandas.Series
        A Series containing RES generation per carrier (MWh/year).
    """
    energy_balance = (
        n.statistics.energy_balance(aggregate_time="sum", nice_names=False)
        .groupby("carrier")
        .sum()
    )
    res_generation = energy_balance.reindex(res_carriers).dropna()
    return res_generation


def calculate_res_dump_per_carrier(
    n: pypsa.Network, res_carriers: list[str]
) -> pd.Series:
    """
    Calculate total RES dumped energy per carrier (MWh/year) using curtailment.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.

    Returns
    -------
    pandas.Series
        A Series containing RES dumped energy per carrier (MWh/year).
    """
    curtailed = n.statistics.curtailment(nice_names=False).groupby("carrier").sum()
    res_dump = curtailed.reindex(res_carriers).dropna()
    return res_dump


def get_co2_ets_price(config, planning_horizon) -> float:
    """
    Retrieve the CO2 ETS price for a given planning horizon from the configuration.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing emission prices under the "costs" key.
    planning_horizon : int or str
        The year or period for which the CO2 ETS price is requested.

    Returns
    -------
    float
        The CO2 ETS price for the specified planning horizon.
    """
    emission_prices = config.get("costs", {}).get("emission_prices", {})
    if not emission_prices.get("enable", False):
        raise KeyError("Emission prices are not enabled in the config")

    co2_prices = emission_prices.get("co2", {})
    price = co2_prices.get(planning_horizon, co2_prices.get(str(planning_horizon)))
    if price is None:
        raise KeyError(f"Missing CO2 ETS price for {planning_horizon}")
    return float(price)


def calculate_b1_indicator(n_reference, n_project, method="pint"):
    """
    Calculate B1 indicator: change in total system cost.

    The interpretation depends on the method:
    - PINT: positive B1 means beneficial (project reduces costs)
    - TOOT: positive B1 means beneficial (removing project increases costs)

    Args:
        n_reference: Reference network
        n_project: Project network
        method: Either "pint" or "toot" (case-insensitive)

    Returns:
        dict: Dictionary with B1 and component costs
    """
    # Calculate costs for both scenarios
    cost_reference = calculate_total_system_cost(n_reference)
    cost_project = calculate_total_system_cost(n_project)

    if method == "pint":
        # PINT: positive B1 means beneficial (project reduces costs)
        # Reference is without project
        # Project is with project
        b1 = cost_reference["total"] - cost_project["total"]
    else:  # toot
        # TOOT: positive B1 means beneficial (removing project increases costs)
        # Reference is with all projects
        # Project is without project
        b1 = cost_project["total"] - cost_reference["total"]

    is_beneficial = b1 > 0

    if method == "pint":
        # Speak in terms of sew.
        if is_beneficial:
            interpretation = "The project reduces costs compared to the reference scenario without the project."
        else:
            interpretation = "The project increases costs compared to the reference scenario without the project."
    else:  # toot
        if is_beneficial:
            interpretation = "The project is beneficial as removing it increases costs compared to the reference scenario with all projects."
        else:
            interpretation = "The project increases costs as removing it decreases costs compared to the reference scenario with all projects."

    results = {}
    results["B1_total_system_cost_change"] = b1  # in Euros. Positive is beneficial.
    results["is_beneficial"] = is_beneficial
    results["interpretation"] = interpretation
    results["cost_reference"] = cost_reference["total"]
    results["capex_reference"] = cost_reference["capex"]
    results["opex_reference"] = cost_reference["opex"]
    results["cost_project"] = cost_project["total"]
    results["capex_project"] = cost_project["capex"]
    results["opex_project"] = cost_project["opex"]

    if method == "pint":
        results["capex_change"] = cost_reference["capex"] - cost_project["capex"]
        results["opex_change"] = cost_reference["opex"] - cost_project["opex"]
    else:  # toot
        results["capex_change"] = cost_project["capex"] - cost_reference["capex"]
        results["opex_change"] = cost_project["opex"] - cost_reference["opex"]

    return results


def calculate_b2_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    co2_societal_costs: dict,
    co2_ets_price: float,
) -> dict:
    """
    Calculate B2 indicator: change in CO2 emissions and societal cost.

    Returns totals for CO2 (t) and societal cost (EUR/year) for low/central/high
    societal cost assumptions.
    """

    co2_reference = calculate_co2_emissions_per_carrier(n_reference)
    co2_project = calculate_co2_emissions_per_carrier(n_project)

    if method == "pint":
        # Reference is without project, project is with project
        co2_diff = co2_reference - co2_project
    else:  # toot
        # Reference is with all projects, project is without project
        co2_diff = co2_project - co2_reference

    results = {
        "co2_variation": co2_diff,
        "co2_ets_price": co2_ets_price,
        "co2_societal_cost_low": co2_societal_costs["low"],
        "co2_societal_cost_central": co2_societal_costs["central"],
        "co2_societal_cost_high": co2_societal_costs["high"],
    }

    for level in ["low", "central", "high"]:
        b2_val = co2_diff * (co2_societal_costs[level] - co2_ets_price)
        results[f"B2_societal_cost_variation_{level}"] = b2_val

    return results


def calculate_b3_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    res_carriers: list[str],
) -> dict:
    """
    Calculate B3 indicator: change in RES capacity (MW) and generation (MWh/year).

    Parameters
    ----------
    n_reference : pypsa.Network
        Reference network.
    n_project : pypsa.Network
        Project network.
    method : str
        Either "pint" or "toot" (case-insensitive).

    Returns
    -------
    dict
        Dictionary with B3 indicators:
        - B3_res_capacity_change_mw
        - B3_res_generation_change_mwh
        - B3_res_dump_change_mwh
    """
    n_with, n_without = (
        (n_project, n_reference) if method == "pint" else (n_reference, n_project)
    )

    capacity_with = calculate_res_capacity_per_carrier(n_with, res_carriers)
    capacity_without = calculate_res_capacity_per_carrier(n_without, res_carriers)
    generation_with = calculate_res_generation_per_carrier(n_with, res_carriers)
    generation_without = calculate_res_generation_per_carrier(n_without, res_carriers)
    dump_with = calculate_res_dump_per_carrier(n_with, res_carriers)
    dump_without = calculate_res_dump_per_carrier(n_without, res_carriers)

    capacity_diff = capacity_with.sum() - capacity_without.sum()
    generation_diff = generation_with.sum() - generation_without.sum()
    dump_diff = dump_with.sum() - dump_without.sum()

    return {
        "B3_res_capacity_change_mw": capacity_diff,
        "B3_res_generation_change_mwh": generation_diff,
        "B3_annual_avoided_curtailment_mwh": dump_diff,
    }


def calculate_b4_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    emission_factors: pd.DataFrame,
    conventional_carriers: list[str],
) -> dict:
    """
    Calculate B4 indicator: non-CO2 emissions (tons/year).

    Parameters
    ----------
    n_reference : pypsa.Network
        Reference network.
    n_project : pypsa.Network
        Project network.
    method : str
        Either "pint" or "toot" (case-insensitive).
    emission_factors : pd.DataFrame
        DataFrame with non-CO2 emission factors (kg/MWh) indexed by carrier and
        with columns for different pollutants and statistics (min, mean, max).

    Returns
    -------
    dict
        Dictionary with B4 indicators for each pollutant:
        - B4{sub}_{pollutant}
    """

    pollutant_keys = list(emission_factors.columns)

    normalized = []
    for c in conventional_carriers:
        normalized.append("oil" if c.startswith("oil-") else c)
    conventional_carriers = sorted(set(normalized))

    # initialize emissions dictionary
    ref_emissions = {key: 0.0 for key in pollutant_keys}
    proj_emissions = {key: 0.0 for key in pollutant_keys}

    # function to check if link is biofuel/biomass/biogas link
    def is_biofuel_link(link_carrier: str, network: pypsa.Network) -> bool:
        links = network.links[network.links.carrier == link_carrier]
        if links.empty:
            return False
        buses = links[["bus0", "bus1"]].astype(str).agg(" ".join, axis=1).str.lower()
        return buses.str.contains("biofuel|biomass|biogas", regex=True).any()

    for carrier in conventional_carriers:
        if carrier not in CARRIER_TO_EMISSION_FACTORS:
            raise ValueError(
                f"Carrier '{carrier}' missing in CARRIER_TO_EMISSION_FACTORS"
            )
        regular_fuel, biofuel = CARRIER_TO_EMISSION_FACTORS[carrier]

        ref_balance = n_reference.statistics.energy_balance(
            nice_names=False, bus_carrier=carrier
        )
        proj_balance = n_project.statistics.energy_balance(
            nice_names=False, bus_carrier=carrier
        )

        # use only positive values
        ref_balance = ref_balance[ref_balance > 0]
        proj_balance = proj_balance[proj_balance > 0]

        # if regular fuel, use regular emission factors
        # if biofuel, use biofuel emission factors
        fuel_key = str(regular_fuel).strip()
        if fuel_key not in emission_factors.index:
            raise ValueError(f"No emission factors found for fuel '{regular_fuel}'")
        regular_factors = emission_factors.loc[fuel_key]
        biofuel_factors = None
        if biofuel:
            biofuel_key = str(biofuel).strip()
            if biofuel_key not in emission_factors.index:
                raise ValueError(f"No emission factors found for fuel '{biofuel}'")
            biofuel_factors = emission_factors.loc[biofuel_key]

        def accumulate(balance, emissions, network):
            for (component, item, _), value in balance.items():
                if component == "Generator":
                    factors = regular_factors
                elif component == "Link" and biofuel_factors is not None:
                    factors = (
                        biofuel_factors
                        if is_biofuel_link(item, network)
                        else regular_factors
                    )
                else:
                    factors = regular_factors
                for pollutant_key, kg_value in factors.items():
                    emissions[pollutant_key] += float(value) * kg_value

        accumulate(ref_balance, ref_emissions, n_reference)
        accumulate(proj_balance, proj_emissions, n_project)

    results = {}
    for pollutant_key in pollutant_keys:
        ref_val = ref_emissions.get(pollutant_key, 0.0)
        proj_val = proj_emissions.get(pollutant_key, 0.0)
        if method == "pint":
            diff = ref_val - proj_val
        else:
            diff = proj_val - ref_val
        results[pollutant_key] = diff

    return results


def apply_indicator_units(df: pd.DataFrame) -> pd.DataFrame:
    indicators = set(df["indicator"])
    missing = indicators - set(INDICATOR_UNITS)
    if missing:
        raise ValueError(
            f"Missing units for indicators: {sorted(missing)}. "
            "Update INDICATOR_UNITS to continue."
        )
    extra = set(INDICATOR_UNITS) - indicators
    if extra:
        logger.warning("Unused indicator units defined: %s", sorted(extra))
    df = df.copy()
    df["units"] = df["indicator"].map(INDICATOR_UNITS)
    return df


def build_long_indicators(indicators: dict) -> pd.DataFrame:
    meta = {
        "project_id": indicators.get("project_id"),
        "method": indicators.get("cba_method"),
        "is_beneficial": indicators.get("is_beneficial"),
        "interpretation": indicators.get("interpretation"),
    }

    rows = []
    skip_keys = set(meta.keys()) | {"cba_method"}
    for key, value in indicators.items():
        if key in skip_keys:
            continue

        indicator = key
        subindex = ""
        for suffix in ["_min", "_mean", "_max", "_low", "_central", "_high"]:
            if indicator.endswith(suffix):
                indicator = indicator[: -len(suffix)]
                subindex = suffix.lstrip("_")
                break

        rows.append(
            {
                **meta,
                "indicator": indicator,
                "subindex": subindex,
                "value": value,
            }
        )

    return pd.DataFrame(rows)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("make_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load both networks
    n_reference = pypsa.Network(snakemake.input.reference)
    n_project = pypsa.Network(snakemake.input.project)

    # Validate networks are solved
    if not n_reference.is_solved:
        raise ValueError("Reference network is not solved")
    if not n_project.is_solved:
        raise ValueError("Project network is not solved")

    # Detect method from wildcards (toot or pint)
    method = check_method(snakemake.wildcards.cba_method)
    planning_horizon = int(snakemake.wildcards.planning_horizons)

    # Calculate indicators
    indicators = calculate_b1_indicator(n_reference, n_project, method=method)

    co2_societal_costs_map = snakemake.config["cba"]["co2_societal_cost"]
    co2_societal_costs = get(co2_societal_costs_map, planning_horizon)

    co2_ets_price = get_co2_ets_price(snakemake.config, planning_horizon)
    b2_indicators = calculate_b2_indicator(
        n_reference,
        n_project,
        method=method,
        co2_societal_costs=co2_societal_costs,
        co2_ets_price=co2_ets_price,
    )
    indicators.update(b2_indicators)

    res_carriers = snakemake.config.get("electricity", {}).get(
        "tyndp_renewable_carriers"
    )
    b3_indicators = calculate_b3_indicator(
        n_reference,
        n_project,
        method=method,
        res_carriers=res_carriers,
    )
    indicators.update(b3_indicators)

    emission_factors = load_non_co2_emission_factors(snakemake.input.non_co2_emissions)
    conventional_carriers = snakemake.config.get("electricity", {}).get(
        "tyndp_conventional_carriers", []
    )
    b4_indicators = calculate_b4_indicator(
        n_reference,
        n_project,
        method=method,
        emission_factors=emission_factors,
        conventional_carriers=conventional_carriers,
    )
    indicators.update(b4_indicators)

    # Add project metadata
    cba_project = snakemake.wildcards.cba_project
    indicators["project_id"] = int(cba_project[1:])  # assuming format 't123'
    indicators["cba_method"] = method.upper()

    logger.info(
        f"Project {indicators['project_id']} is {'beneficial' if indicators['is_beneficial'] else 'not beneficial'} for {indicators['cba_method']}. B1 indicator: {indicators['B1_total_system_cost_change']} Euros"
    )

    # Convert to DataFrame and save
    df = build_long_indicators(indicators)
    df = apply_indicator_units(df)
    df.to_csv(snakemake.output.indicators, index=False)
