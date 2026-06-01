# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Calculate CBA indicators by comparing reference and project scenarios.

This script computes the B1 indicator (operational cost difference)
and other CBA metrics by analyzing the solved networks for reference and project
cases.

PINT (Put In at a Time):
    - Reference: Network WITHOUT any projects
    - Project: Network WITH the specific project added
    - B1 = OPEX(reference) - OPEX(with project)

TOOT (Take Out One at a Time):
    - Reference: Network WITH all projects (current plan)
    - Project: Network WITHOUT the specific project (removed)
    - B1 = OPEX(without project) - OPEX(reference)


References:
- CBA guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/news/2024/entso-e_4th_CBA_Guideline_240409.pdf
    - section 3.2.2: TOOT and PINT, page 23-24
- CBA implementation guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/tyndp-documents/TYNDP2024/foropinion/CBA_Implementation_Guidelines.pdf
    - section 5.1: B1 - SEW, page 58-59

"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)

INDICATOR_UNITS = {
    "B1_total_system_cost_change": "Meuro/year",
    "cost_reference": "Meuro/year",
    "capex_reference": "Meuro/year",
    "opex_reference": "Meuro/year",
    "cost_project": "Meuro/year",
    "capex_project": "Meuro/year",
    "opex_project": "Meuro/year",
    "capex_change": "Meuro/year",
    "opex_change": "Meuro/year",
    "B2a_co2_variation": "t/year",
    "co2_ets_price": "EUR/t",
    "co2_societal_cost": "EUR/t",
    "B2a_societal_cost_variation": "Meuro/year",
    "B3a_res_capacity_change": "MW",
    "B3_res_generation_change": "MWh/year",
    "B3_annual_avoided_curtailment": "MWh/year",
    "B4a_nox": "ton/year",
    "B4b_nh3": "ton/year",
    "B4c_sox": "ton/year",
    "B4d_pm25": "ton/year",
    "B4e_pm10": "ton/year",
    "B4f_nmvoc": "ton/year",
}

CARRIER_TO_EMISSION_FACTORS = {
    "gas": ("Gas", "Gas biofuel"),
    "coal": ("Hard coal", "Hard Coal biofuel"),
    "lignite": ("Lignite", "Lignite biofuel"),
    "uranium": ("Nuclear", None),
    "oil": ("Light Oil", "Light Oil biofuel"),  # all oil fuels have same factors
}


def _apply_original_costs(n, remove_noisy_costs: bool) -> None:
    if not remove_noisy_costs:
        return
    for t in n.iterate_components():
        if "marginal_cost_original" in t.static:
            mask = t.static["marginal_cost_original"].notna()
            t.static.loc[mask, "marginal_cost"] = t.static.loc[
                mask, "marginal_cost_original"
            ].astype(t.static["marginal_cost"].dtype)

    for t in n.iterate_components(["Line", "Link"]):
        if "capital_cost_original" in t.static:
            mask = t.static["capital_cost_original"].notna()
            t.static.loc[mask, "capital_cost"] = t.static.loc[
                mask, "capital_cost_original"
            ].astype(t.static["capital_cost"].dtype)


def calculate_total_system_cost(n, remove_noisy_costs: bool = False):
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

    _apply_original_costs(n, remove_noisy_costs)

    # Use PyPSA's built-in statistics methods
    capex = n.statistics.capex().sum() / 1e6
    opex = n.statistics.opex(aggregate_time="sum").sum() / 1e6
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


def difference_by_method(reference: float, project: float, method: str) -> float:
    """
    Calculate differences for indicator values by method (TOOT / PINT).

    This maintains consistency in sign conventions across indicators (B2 and B4), for example.

    - PINT: reference is without project, project is with project, so difference = reference - project
    - TOOT: reference is with project, project is without project, so difference = project - reference

    Using this method, for the B2 and B4 indicators, a positive difference will always indicate a beneficial
    impact of the project (e.g. emissions reduction, even if the reference and project are swapped
    between TOOT and PINT).

    Parameters
    ----------
    reference : float
        The value from the reference scenario.
    project : float
        The value from the project scenario.
    method : str
        The CBA method being used, either "pint" or "toot".

    Returns
    -------
    float
        The calculated difference according to the method's sign convention.
    """
    if method == "pint":
        return reference - project
    return project - reference


def _normalize_conventional_carrier(carrier: str) -> str:
    """
    Map carrier aliases to the conventional fuel names used in CBA.

    This maps all the oil carriers (such as "oil-light", "oil-heavy", and "oil-shale") to "oil" to use the same emission factors, while keeping other carriers unchanged.

    Parameters
    ----------
    carrier : str
        The original carrier name from the network.

    Returns
    -------
    str
        The normalized carrier name for matching with emission factors.
    """
    return "oil" if carrier.startswith("oil-") else carrier


def _get_snapshot_weightings(n: pypsa.Network) -> pd.Series:
    return n.snapshot_weightings.generators.reindex(n.snapshots).fillna(1.0)


def get_ac_electricity_producing_assets(n: pypsa.Network) -> pd.Series:
    """
    Get list of assets that produce electricity.

    Check energy balance for components with bus_carrier "AC" and positive net balance.
    This represents generators and links that produce electricity on AC buses.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.

    Returns
    -------
    pandas.Series
        A Series containing the electricity-producing assets on AC buses.
    """
    balance = n.statistics.energy_balance(
        groupby_time="sum",
        nice_names=False,
        bus_carrier="AC",
        groupby=["name", "carrier", "bus_carrier"],
    )
    balance = balance[balance > 0]

    assets = balance[
        balance.index.get_level_values("component").isin(["Generator", "Link"])
    ]
    return assets


def get_ac_energy_balance(
    n: pypsa.Network,
    assets: pd.Series,
    *,
    bus_carrier: str | None = None,
) -> pd.DataFrame:
    """
    Return energy balance for provided assets.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.
    assets : pandas.Series
        Assets for which to calculate energy balance.
    bus_carrier : str, optional
        If set, filter energy balance to this bus carrier (e.g. "co2").
    """
    balance = n.statistics.energy_balance(
        groupby_time=False,
        nice_names=False,
        bus_carrier=bus_carrier,
        groupby=["name", "carrier"],
    )

    return balance.reindex(assets.index.droplevel("bus_carrier")).dropna(how="all")


def calculate_power_sector_co2_emissions(
    n: pypsa.Network, ac_assets: pd.Series | None = None
) -> float:
    """
    Calculate annual power-sector CO2 emissions for assets producing on AC buses.

    Generators use carrier-specific ``co2_emissions`` intensities. Links use the
    explicit CO2 port flows of the electricity-producing asset.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object.
    ac_assets : pandas.Series, optional
        Pre-filtered Series of electricity-producing assets on AC buses.
        If not provided, it will be computed within the function.
        However, calculating it before calling this function can improve performance.

    Returns
    -------
    float
        Total annual CO2 emissions.
    """
    if ac_assets is None:
        ac_assets = get_ac_electricity_producing_assets(n)
    weights = _get_snapshot_weightings(n)
    co2_balance = get_ac_energy_balance(n, ac_assets, bus_carrier="co2")
    co2_emissions_per_h = co2_balance.sum()
    total_emissions = co2_emissions_per_h.mul(weights).sum()

    return float(total_emissions)


def load_non_co2_emission_factors(path: str) -> pd.DataFrame:
    """
    Load non-CO2 emission factors from ENTSO-E data and compute mean values.

    Parameters
    ----------
    path : str
        Path to the CSV file containing non-CO2 emission factors.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the mean non-CO2 emission factors for each fuel type.
    """
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


def calculate_b1_indicator(
    n_reference, n_project, method="pint", remove_noisy_costs: bool = False
):
    """
    Calculate B1 indicator.

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
    # Calculate full cost breakdowns for reporting
    cost_reference = calculate_total_system_cost(n_reference, remove_noisy_costs)
    cost_project = calculate_total_system_cost(n_project, remove_noisy_costs)

    # For B1 calculation, use only OPEX
    if method == "pint":
        # PINT: positive B1 means beneficial (project reduces costs)
        # Reference is without project
        # Project is with project
        b1 = cost_reference["opex"] - cost_project["opex"]
    else:  # toot
        # TOOT: positive B1 means beneficial (removing project increases costs)
        # Reference is with all projects
        # Project is without project
        b1 = cost_project["opex"] - cost_reference["opex"]

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

    results = {
        "B1_total_system_cost_change": b1,
        "is_beneficial": is_beneficial,
        "interpretation": interpretation,
        "cost_reference": cost_reference["total"],
        "capex_reference": cost_reference["capex"],
        "opex_reference": cost_reference["opex"],
        "cost_project": cost_project["total"],
        "capex_project": cost_project["capex"],
        "opex_project": cost_project["opex"],
    }
    units = {
        "B1_total_system_cost_change": "Meuro/year",
        "cost_reference": "Meuro/year",
        "capex_reference": "Meuro/year",
        "opex_reference": "Meuro/year",
        "cost_project": "Meuro/year",
        "capex_project": "Meuro/year",
        "opex_project": "Meuro/year",
    }

    if method == "pint":
        results["capex_change"] = cost_reference["capex"] - cost_project["capex"]
        results["opex_change"] = cost_reference["opex"] - cost_project["opex"]
    else:  # toot
        results["capex_change"] = cost_project["capex"] - cost_reference["capex"]
        results["opex_change"] = cost_project["opex"] - cost_reference["opex"]

    units["capex_change"] = "Meuro/year"
    units["opex_change"] = "Meuro/year"

    return results, units


def calculate_b2_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    co2_societal_costs: dict,
    co2_ets_price: float,
    ac_assets_reference: pd.Series | None = None,
    ac_assets_project: pd.Series | None = None,
) -> tuple[dict, dict]:
    """
    Calculate B2 indicator: change in CO2 emissions and societal cost.

    Returns totals for CO2 (t) and societal cost (EUR/year) for low/central/high
    societal cost assumptions.

    Parameters
    ----------
    n_reference : pypsa.Network
        Reference network.
    n_project : pypsa.Network
        Project network.
    method : str
        Either "pint" or "toot" (case-insensitive).
    co2_societal_costs : dict
        Dictionary with keys "low", "central", "high" for societal cost of CO2 in EUR/t.
    co2_ets_price : float
        The CO2 ETS price in EUR/t for the relevant planning horizon.
    ac_assets_reference : pandas.Series, optional
        Pre-filtered Series of electricity-producing assets on AC buses for the reference network.
        If not provided, it will be computed within the function.
        Providing these can improve performance by avoiding redundant calculations.
    ac_assets_project : pandas.Series, optional
        Pre-filtered Series of electricity-producing assets on AC buses for the project network.
        If not provided, it will be computed within the function.
        Providing these can improve performance by avoiding redundant calculations.

    Returns
    -------
    tuple[dict, dict]
        A tuple containing:
        - A dictionary with B2 indicators:
            - B2a_co2_variation (ktonnes/year)
            - co2_ets_price (EUR/t)
            - co2_societal_cost_low (EUR/t)
            - co2_societal_cost_central (EUR/t)
            - co2_societal_cost_high (EUR/t)
            - B2a_societal_cost_variation_{level} (Meuro/year), where the level is low, central, or high
    """

    co2_reference = calculate_power_sector_co2_emissions(
        n_reference, ac_assets=ac_assets_reference
    )
    co2_project = calculate_power_sector_co2_emissions(
        n_project, ac_assets=ac_assets_project
    )

    # co2_diff is avoided emissions: (positive = beneficial) for monetisation
    # co2_diff is calculating the difference of emissions without project minus with project, regardless of method
    co2_diff = difference_by_method(co2_reference, co2_project, method)

    co2_diff_ktonnes = co2_diff / 1000.0
    results = {
        "B2a_co2_variation": -co2_diff_ktonnes,  # Report B2a with guideline sign convention for CO2 variation
        "co2_ets_price": co2_ets_price,
        "co2_societal_cost_low": co2_societal_costs["low"],
        "co2_societal_cost_central": co2_societal_costs["central"],
        "co2_societal_cost_high": co2_societal_costs["high"],
    }
    units = {
        "B2a_co2_variation": "ktonnes/year",
        "co2_ets_price": "EUR/t",
        "co2_societal_cost_low": "EUR/t",
        "co2_societal_cost_central": "EUR/t",
        "co2_societal_cost_high": "EUR/t",
    }

    for level in ["low", "central", "high"]:
        b2_val = co2_diff * (co2_societal_costs[level] - co2_ets_price)
        results[f"B2a_societal_cost_variation_{level}"] = b2_val / 1e6
        units[f"B2a_societal_cost_variation_{level}"] = "Meuro/year"

    return results, units


def calculate_b3_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    res_carriers: list[str],
) -> tuple[dict, dict]:
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
    tuple[dict, dict]
        Dictionary with B3 indicators:
        - B3a_res_capacity_change
        - B3_res_generation_change
        - B3_annual_avoided_curtailment
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
    generation_diff_mwh = generation_with.sum() - generation_without.sum()
    dump_diff_mwh = dump_with.sum() - dump_without.sum()
    avoided_curtailment_gwh = dump_diff_mwh / 1000.0

    results = {
        "B3a_res_capacity_change": capacity_diff,
        "B3_res_generation_change": generation_diff_mwh / 1000.0,
        "B3_annual_avoided_curtailment": avoided_curtailment_gwh,
    }
    units = {
        "B3a_res_capacity_change": "MW",
        "B3_res_generation_change": "GWh/year",
        "B3_annual_avoided_curtailment": "GWh/year",
    }
    return results, units


def calculate_b4_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    emission_factors: pd.DataFrame,
    conventional_carriers: list[str],
    ac_assets_reference: pd.Series | None = None,
    ac_assets_project: pd.Series | None = None,
) -> tuple[dict, dict]:
    """
    Calculate B4 indicator: non-CO2 emissions (ton/year).

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
    conventional_carriers = sorted(
        {_normalize_conventional_carrier(carrier) for carrier in conventional_carriers}
    )

    ref_emissions = {key: 0.0 for key in pollutant_keys}
    proj_emissions = {key: 0.0 for key in pollutant_keys}

    def factors_for_fuel(carrier: str) -> tuple[pd.Series, pd.Series | None]:
        if carrier not in CARRIER_TO_EMISSION_FACTORS:
            raise ValueError(
                f"Carrier '{carrier}' missing in CARRIER_TO_EMISSION_FACTORS"
            )
        regular_fuel, biofuel = CARRIER_TO_EMISSION_FACTORS[carrier]

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

        return regular_factors, biofuel_factors

    factors_by_fuel = {
        carrier: factors_for_fuel(carrier) for carrier in conventional_carriers
    }
    for network, assets, emissions in [
        (n_reference, ac_assets_reference, ref_emissions),
        (n_project, ac_assets_project, proj_emissions),
    ]:
        if assets is None:
            assets = get_ac_electricity_producing_assets(network)
        weights = _get_snapshot_weightings(network)

        for fuel in conventional_carriers:
            bal = get_ac_energy_balance(network, assets, bus_carrier=fuel)

            # function to check if component is biofuel/biomass/biogas
            is_bio = bal.index.get_level_values("name").astype(str).str.contains(
                "bio", case=False, na=False
            ) | bal.index.get_level_values("carrier").astype(str).str.contains(
                "bio", case=False, na=False
            )

            fuel_use = -bal.clip(
                upper=0.0
            )  # convert negative consumption values to positive
            regular_use = fuel_use.loc[
                ~is_bio
            ].sum()  # sum non-biomass fuel consumption
            bio_use = fuel_use.loc[is_bio].sum()  # sum biomass fuel consumption

            regular_total = regular_use.mul(weights).sum()
            bio_total = bio_use.mul(weights).sum()
            regular_factors, bio_factors = factors_by_fuel[fuel]

            for pollutant_key, kg_value in regular_factors.items():
                emissions[pollutant_key] += float(regular_total) * kg_value
            if bio_factors is not None:
                for pollutant_key, kg_value in bio_factors.items():
                    emissions[pollutant_key] += float(bio_total) * kg_value

    results = {}
    units = {}
    for pollutant_key in pollutant_keys:
        ref_val = ref_emissions.get(pollutant_key, 0.0)
        proj_val = proj_emissions.get(pollutant_key, 0.0)
        diff = difference_by_method(ref_val, proj_val, method)
        # Convert kg to tons
        results[pollutant_key] = (
            -diff / 1000.0
        )  # report indicators to align with guideline sign convention (negative means reduction in emissions)
        units[pollutant_key] = "ton/year"

    return results, units


def apply_indicator_units(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply units to indicators based on INDICATOR_UNITS mapping.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing an 'indicator' column and optionally a 'units' column.

    Returns
    -------
    pandas.DataFrame
        DataFrame with 'units' column filled in based on INDICATOR_UNITS mapping for
        indicators that were missing units.
    """
    df = df.copy()
    missing_units = df["units"].isna() | (df["units"] == "")
    df.loc[missing_units, "units"] = df.loc[missing_units, "indicator"].map(
        INDICATOR_UNITS
    )
    unresolved = (
        df.loc[df["units"].isna() | (df["units"] == ""), "indicator"]
        .dropna()
        .unique()
        .tolist()
    )
    if unresolved:
        logger.warning("Skipping units for unmapped indicators: %s", sorted(unresolved))
    extra = set(INDICATOR_UNITS) - set(df["indicator"])
    if extra:
        logger.warning("Unused indicator units defined: %s", sorted(extra))
    return df


def build_long_indicators(indicators: dict, units: dict) -> pd.DataFrame:
    """
    Create long-format DataFrame of indicators

    Parameters
    ----------
    indicators : dict
        Dictionary of indicator values, where keys are indicator names (sometimes with suffixes like _min, _mean, _max)
        and values are the corresponding indicator values.
    units : dict
        Dictionary mapping indicator names (without suffixes) to their units. This is used to assign
        units to indicators, even if the indicator keys in the indicators dict have suffixes.

    Returns
    -------
    pandas.DataFrame
        A long-format DataFrame with columns for meta information (planning horizon, project id, etc.) and
        columns for indicator, subindex (e.g. min/mean/max), units, and value.
    """
    meta = {
        "planning_horizon": indicators.get("planning_horizon"),
        "project_id": indicators.get("project_id"),
        "project_code": indicators.get("project_code"),
        "project_type": indicators.get("project_type"),
        "method": indicators.get("cba_method"),
        "is_beneficial": indicators.get("is_beneficial"),
        "interpretation": indicators.get("interpretation"),
        "source": indicators.get("source", "Open-TYNDP"),
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

        unit = units.get(key, "")
        if not unit:
            unit = units.get(indicator, "")

        rows.append(
            {
                **meta,
                "indicator": indicator,
                "subindex": subindex,
                "units": unit,
                "value": value,
            }
        )

    return pd.DataFrame(rows)


def load_benchmark_rows(
    benchmark_path: str | Path,
    project_id: int,
    planning_horizon: int,
    scenario: str | None,
    project_type: str | None,
) -> pd.DataFrame:
    """
    Load benchmark rows for a given project and planning horizon from the ENTSOE TYNDP 2024 results.

    This function looks through the TYNDP 2024 results CSV file for rows matching the specified project ID,
    planning horizon, scenario, and project type.

    The function applies some transformations to align the benchmark data with the calculated indicators,
    such as adjusting the planning horizon to 2040 if benchmarks are only available for 2030 and 2040,
    and mapping any "indicator_mapped" values to the "indicator" column.

    Parameters
    ----------
    benchmark_path : str or Path
        Path to the ENTSOE TYNDP 2024 results CSV file.
    project_id : int
        ID of the project for which to load benchmarks.
    planning_horizon : int
        Planning horizon for which to load benchmarks.
    scenario : str or None
        Scenario for which to load benchmarks.
    project_type : str or None
        Type of the project for which to load benchmarks.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the benchmark rows for the specified project and planning horizon.
    """
    benchmark_horizon = planning_horizon
    if benchmark_horizon not in [2030, 2040]:
        logger.warning(
            "Benchmark data only available for 2030/2040. Using 2040 benchmarks for planning horizon %s.",
            planning_horizon,
        )
        benchmark_horizon = 2040
    path = Path(benchmark_path)
    if not path.exists():
        logger.warning("Benchmark file %s not found, skipping", path)
        return pd.DataFrame()

    benchmark = pd.read_csv(path)
    if benchmark.empty:
        return pd.DataFrame()

    if scenario and not str(scenario)[:4].isdigit():
        scenario = f"{benchmark_horizon}{scenario}"

    benchmark = benchmark.loc[
        (benchmark["project_id"] == project_id)
        & (benchmark["planning_horizon"] == benchmark_horizon)
        & (benchmark["project_type"] == project_type)
        & (benchmark["scenario"] == scenario)
    ].copy()

    if benchmark.empty:
        logger.warning(
            "No benchmark rows found for project %s (horizon %s, scenario %s, type %s)",
            project_id,
            planning_horizon,
            scenario,
            project_type,
        )
        return pd.DataFrame()

    benchmark["planning_horizon"] = planning_horizon
    if "indicator_mapped" in benchmark.columns:
        benchmark["indicator"] = benchmark["indicator_mapped"].fillna(
            benchmark["indicator"]
        )

    benchmark["method"] = (
        benchmark.get("method", None).astype(str).str.upper()
        if "method" in benchmark.columns
        else None
    )
    benchmark["source"] = "TYNDP 2024"
    benchmark["is_beneficial"] = None
    benchmark["interpretation"] = None

    return benchmark[
        [
            "planning_horizon",
            "project_id",
            "project_code",
            "project_type",
            "method",
            "is_beneficial",
            "interpretation",
            "indicator",
            "subindex",
            "units",
            "value",
            "source",
        ]
    ]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_indicators",
            run="NT",
            cba_project="t4",
            planning_horizons="2030",
            configfiles=["config/config.tyndp.yaml"],
        )
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

    planning_horizon = int(snakemake.wildcards.planning_horizons)
    # Detect method from assignments (toot or pint)
    cba_project = snakemake.wildcards.cba_project
    project_id = int(cba_project[1:])
    methods = pd.read_csv(snakemake.input.methods)
    method_row = methods[
        (methods["project_id"] == project_id)
        & (methods["planning_horizon"] == planning_horizon)
    ]
    if method_row.empty:
        raise ValueError(
            f"Missing CBA method for project {project_id} and horizon {planning_horizon}"
        )
    method = check_method(method_row["method"].iloc[0])

    # Calculate indicators
    indicators = {}
    units = {}

    noisy_costs_option = snakemake.config["cba"].get("remove_noisy_costs", False)
    b1_indicators, b1_units = calculate_b1_indicator(
        n_reference,
        n_project,
        method=method,
        remove_noisy_costs=noisy_costs_option,
    )
    indicators.update(b1_indicators)
    units.update(b1_units)

    co2_societal_costs_map = snakemake.config["cba"]["co2_societal_cost"]
    co2_cost_horizon = planning_horizon
    if co2_cost_horizon not in [2030, 2040]:
        logger.warning(
            "CO2 societal costs only available for 2030/2040. Using 2040 for planning horizon %s.",
            planning_horizon,
        )
        co2_cost_horizon = 2040
    co2_societal_costs = get(co2_societal_costs_map, co2_cost_horizon)

    co2_ets_price = get_co2_ets_price(snakemake.config, planning_horizon)
    ac_assets_reference = get_ac_electricity_producing_assets(n_reference)
    ac_assets_project = get_ac_electricity_producing_assets(n_project)

    b2_indicators, b2_units = calculate_b2_indicator(
        n_reference,
        n_project,
        method=method,
        co2_societal_costs=co2_societal_costs,
        co2_ets_price=co2_ets_price,
        ac_assets_reference=ac_assets_reference,
        ac_assets_project=ac_assets_project,
    )
    indicators.update(b2_indicators)
    units.update(b2_units)

    res_carriers = snakemake.config.get("electricity", {}).get(
        "tyndp_renewable_carriers"
    )
    b3_indicators, b3_units = calculate_b3_indicator(
        n_reference,
        n_project,
        method=method,
        res_carriers=res_carriers,
    )
    indicators.update(b3_indicators)
    units.update(b3_units)

    emission_factors = load_non_co2_emission_factors(snakemake.input.non_co2_emissions)
    conventional_carriers = snakemake.config.get("electricity", {}).get(
        "tyndp_conventional_carriers", []
    )
    b4_indicators, b4_units = calculate_b4_indicator(
        n_reference,
        n_project,
        method=method,
        emission_factors=emission_factors,
        conventional_carriers=conventional_carriers,
        ac_assets_reference=ac_assets_reference,
        ac_assets_project=ac_assets_project,
    )
    indicators.update(b4_indicators)
    units.update(b4_units)

    # Add project metadata
    project_type = "storage" if cba_project.startswith("s") else "transmission"
    indicators["planning_horizon"] = planning_horizon
    indicators["project_id"] = project_id  # numeric id
    indicators["project_code"] = cba_project
    indicators["project_type"] = project_type
    indicators["cba_method"] = method.upper()
    logger.info(
        f"Project {indicators['project_id']} is {'beneficial' if indicators['is_beneficial'] else 'not beneficial'} for {indicators['cba_method']}. B1 indicator: {indicators['B1_total_system_cost_change']:.2f} Meuro/year"
    )

    # Convert to DataFrame and save
    df_model = build_long_indicators(indicators, units)

    benchmark_scenario = snakemake.config.get("cba", {}).get(
        "sb_scenario"
    ) or snakemake.config.get("tyndp_scenario")
    benchmark_rows = load_benchmark_rows(
        snakemake.input.benchmark,
        indicators["project_id"],
        planning_horizon,
        benchmark_scenario,
        project_type,
    )

    if not benchmark_rows.empty:
        benchmark_rows = benchmark_rows.reindex(
            columns=df_model.columns, fill_value=None
        )
        df = pd.concat([df_model, benchmark_rows], ignore_index=True)
    else:
        df = df_model

    df = apply_indicator_units(df)

    reference_cyears = pd.DatetimeIndex(n_reference.snapshots).year.unique()
    project_cyears = pd.DatetimeIndex(n_project.snapshots).year.unique()

    # check that both reference and project networks contain exactly one climate year
    if len(reference_cyears) != 1 or len(project_cyears) != 1:
        raise ValueError(
            "More than one climate year found in reference or project snapshots."
        )

    # check that both reference and project networks use the same climate year
    if reference_cyears[0] != project_cyears[0]:
        raise ValueError(
            f"Reference and project networks use different climate years: "
            f"{reference_cyears[0]} != {project_cyears[0]}"
        )

    df["cyear"] = int(reference_cyears[0])
    df.to_csv(snakemake.output.indicators, index=False)
