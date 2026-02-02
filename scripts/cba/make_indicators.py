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
from pathlib import Path

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

CARRIER_BUS_COMPONENT = {
    "OCGT": ("gas", "Link"),
    "gas-ccgt": ("gas", "Link"),
    "gas-ccgt-ccs": ("gas", "Link"),
    "gas-conv": ("gas", "Link"),
    "gas-ocgt": ("gas", "Link"),
    "coal": ("coal", "Generator"),
    "lignite": ("lignite", "Generator"),
    "uranium": ("uranium", "Link"),
    "oil-heavy": ("oil", "Link"),
    "oil-light": ("oil", "Link"),
    "oil-shale": ("oil", "Link"),
}

CARRIER_TO_EMISSION_FACTORS = {
    "OCGT": ("Gas", ["OCGT old", "OCGT new"]),
    "gas-ccgt": (
        "Gas",
        ["CCGT old 1", "CCGT old 2", "CCGT present 1", "CCGT present 2", "CCGT new"],
    ),
    "gas-ccgt-ccs": (
        "Gas",
        ["CCGT CCS"],
    ),
    "gas-conv": ("Gas", ["conventional old 1", "conventional old 2", "OCGT old"]),
    "gas-ocgt": ("Gas", ["OCGT old", "OCGT new"]),
    "coal": ("Hard coal", ["old 1", "old 2", "new"]),
    "lignite": ("Lignite", ["old 1", "old 2", "new"]),
    "uranium": ("Nuclear", []),
    "oil-heavy": ("Heavy Oil", ["old 1", "old 2"]),
    "oil-light": ("Light Oil", []),
    "oil-shale": ("Oil shale", ["old", "new"]),
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
    """
    Load non-CO2 emission factors and convert to kg/MWh.

    The input factors are in units of kg/GJ (thermal) and are converted using the
    standard efficiency (NCV) and 1 MWh = 3.6 GJ.
    """
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df[df["Fuel"].notna()].copy()
    df["Fuel_norm"] = df["Fuel"].astype(str).str.strip().str.lower()
    df["Type_norm"] = df["Type"].fillna("").astype(str).str.strip().str.lower()
    efficiency = (
        df["Standard efficiency in NCV terms"]
        .astype(str)
        .str.replace("%", "", regex=False)
        .pipe(pd.to_numeric, errors="coerce")
        / 100.0
    )

    pollutant_cols = {
        "NOX emission factor": "B4a_nox",
        "NH3 emission factor": "B4b_nh3",
        "SO emission factor": "B4c_sox",
        "PM2.5 and smaller emission factor": "B4d_pm25",
        "PM10 emission factor": "B4e_pm10",
        "NMVOC emission factor": "B4f_nmvoc",
    }

    factors = {}
    for col, key in pollutant_cols.items():
        # kg/GJ_el = kg/GJ_th / efficiency
        # kg/MWh_el = kg/GJ_el * 3.6
        values = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        factors[key] = values / efficiency * 3.6

    factor_df = pd.DataFrame(factors)
    factor_df["Fuel_norm"] = df["Fuel_norm"]
    factor_df["Type_norm"] = df["Type_norm"]
    factor_df = factor_df[efficiency.notna() & (efficiency > 0)]
    return factor_df


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

    results = {
        "B1_total_system_cost_change": b1 / 1e6,  # convert to Meuro/year
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
        "cost_reference": "EUR/year",
        "capex_reference": "EUR/year",
        "opex_reference": "EUR/year",
        "cost_project": "EUR/year",
        "capex_project": "EUR/year",
        "opex_project": "EUR/year",
    }

    if method == "pint":
        results["capex_change"] = cost_reference["capex"] - cost_project["capex"]
        results["opex_change"] = cost_reference["opex"] - cost_project["opex"]
    else:  # toot
        results["capex_change"] = cost_project["capex"] - cost_reference["capex"]
        results["opex_change"] = cost_project["opex"] - cost_reference["opex"]

    units["capex_change"] = "EUR/year"
    units["opex_change"] = "EUR/year"

    return results, units


def calculate_b2_indicator(
    n_reference: pypsa.Network,
    n_project: pypsa.Network,
    method: str,
    co2_societal_costs: dict,
    co2_ets_price: float,
) -> tuple[dict, dict]:
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

    co2_diff_ktonnes = co2_diff / 1000.0
    results = {
        "B2a_co2_variation": co2_diff_ktonnes,
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
) -> tuple[dict, dict]:
    """
    Calculate B4 indicator: non-CO2 emissions (kg/year).

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
        Dictionary with B4 indicators for each pollutant and statistic (stat in [min, mean, max]):
        - {pollutant}_{stat}
    """

    pollutant_keys = [
        col for col in emission_factors.columns if col not in ("Fuel_norm", "Type_norm")
    ]

    # initialize emissions dictionary
    ref_emissions = {key: 0.0 for key in pollutant_keys}
    proj_emissions = {key: 0.0 for key in pollutant_keys}

    for carrier in conventional_carriers:
        # if carrier from conventional_carriers is 'gas', expand to gas sub-carriers:
        # OCGT, gas-ccgt, gas-ccgt-ccs, gas-conv, gas-ocgt
        carriers = (
            ["OCGT", "gas-ccgt", "gas-ccgt-ccs", "gas-conv", "gas-ocgt"]
            if carrier == "gas"
            else [carrier]
        )
        for subcarrier in carriers:
            if subcarrier not in CARRIER_BUS_COMPONENT:
                raise ValueError(
                    f"Carrier '{subcarrier}' missing in CARRIER_BUS_COMPONENT"
                )
            if subcarrier not in CARRIER_TO_EMISSION_FACTORS:
                raise ValueError(
                    f"Carrier '{subcarrier}' missing in CARRIER_TO_EMISSION_FACTORS"
                )
            bus_carrier, component = CARRIER_BUS_COMPONENT[subcarrier]
            fuel, types = CARRIER_TO_EMISSION_FACTORS[subcarrier]

            # get the energy balance for carrier in reference and project networks
            # get the absolute value and take the sum
            ref_mwh = (
                n_reference.statistics.energy_balance(
                    nice_names=False,
                    carrier=subcarrier,
                    bus_carrier=bus_carrier,
                    components=[component],
                )
                .abs()
                .sum()
            )
            proj_mwh = (
                n_project.statistics.energy_balance(
                    nice_names=False,
                    carrier=subcarrier,
                    bus_carrier=bus_carrier,
                    components=[component],
                )
                .abs()
                .sum()
            )

            fuel_norm = str(fuel).strip().lower()
            type_norms = [t.strip().lower() for t in types]
            factors = emission_factors[emission_factors["Fuel_norm"] == fuel_norm]
            # if there are specific types for a given fuel, filter for those
            # otherwise, filter for empty type
            if type_norms:
                factors = factors[factors["Type_norm"].isin(type_norms)]
            else:
                factors = factors[factors["Type_norm"] == ""]
            if factors.empty:
                raise ValueError(
                    f"No emission factors found for carrier '{subcarrier}'"
                )

            # if multiple factors exist for a carrier, take the mean
            kg_per_mwh = factors[pollutant_keys].agg("mean")
            for pollutant_key, kg_value in kg_per_mwh.items():
                ref_emissions[pollutant_key] += float(ref_mwh) * kg_value
                proj_emissions[pollutant_key] += float(proj_mwh) * kg_value

    results = {}
    units = {}
    for pollutant_key in pollutant_keys:
        ref_val = ref_emissions.get(pollutant_key, 0.0)
        proj_val = proj_emissions.get(pollutant_key, 0.0)
        if method == "pint":
            diff = ref_val - proj_val
        else:
            diff = proj_val - ref_val
        results[pollutant_key] = diff
        units[pollutant_key] = "kg/year"

    return results, units


def apply_indicator_units(df: pd.DataFrame) -> pd.DataFrame:
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
        logger.warning(
            "Skipping units for unmapped indicators: %s", sorted(unresolved)
        )
    extra = set(INDICATOR_UNITS) - set(df["indicator"])
    if extra:
        logger.warning("Unused indicator units defined: %s", sorted(extra))
    return df


def build_long_indicators(indicators: dict, units: dict) -> pd.DataFrame:
    meta = {
        "project_id": indicators.get("project_id"),
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
    path = Path(benchmark_path)
    if not path.exists():
        logger.warning("Benchmark file %s not found, skipping", path)
        return pd.DataFrame()

    benchmark = pd.read_csv(path)
    if benchmark.empty:
        return pd.DataFrame()

    if scenario and not str(scenario)[:4].isdigit():
        scenario = f"{planning_horizon}{scenario}"

    benchmark = benchmark.loc[
        (benchmark["project_id"] == project_id)
        & (benchmark["planning_horizon"] == planning_horizon)
        & (benchmark["type"] == project_type)
        & (benchmark["scenario"] == scenario)
    ].copy()

    if benchmark.empty:
        raise ValueError("No benchmark rows found")

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
            "project_id",
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
    indicators = {}
    units = {}

    b1_indicators, b1_units = calculate_b1_indicator(
        n_reference, n_project, method=method
    )
    indicators.update(b1_indicators)
    units.update(b1_units)

    co2_societal_costs_map = snakemake.config["cba"]["co2_societal_cost"]
    co2_societal_costs = get(co2_societal_costs_map, planning_horizon)

    co2_ets_price = get_co2_ets_price(snakemake.config, planning_horizon)
    b2_indicators, b2_units = calculate_b2_indicator(
        n_reference,
        n_project,
        method=method,
        co2_societal_costs=co2_societal_costs,
        co2_ets_price=co2_ets_price,
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
    )
    indicators.update(b4_indicators)
    units.update(b4_units)

    # Add project metadata
    cba_project = snakemake.wildcards.cba_project
    indicators["project_id"] = int(cba_project[1:])  # assuming format 't123'
    indicators["cba_method"] = method.upper()
    indicators.update(units)
    logger.info(
        f"Project {indicators['project_id']} is {'beneficial' if indicators['is_beneficial'] else 'not beneficial'} for {indicators['cba_method']}. B1 indicator: {indicators['B1_total_system_cost_change']} Euros"
    )

    # Convert to DataFrame and save
    df_model = build_long_indicators(indicators, units)

    scenario = snakemake.config.get("tyndp_scenario") or snakemake.config.get(
        "run", {}
    ).get("name")
    project_type = "storage" if cba_project.startswith("s") else "transmission"
    benchmark_rows = load_benchmark_rows(
        snakemake.input.benchmark,
        indicators["project_id"],
        planning_horizon,
        scenario,
        project_type,
    )

    if not benchmark_rows.empty:
        benchmark_rows = benchmark_rows.reindex(
            columns=df_model.columns, fill_value=None
        )
        df = pd.concat([df_model, benchmark_rows], ignore_index=True)
    else:
        df = df_model

    print(df)
    df = apply_indicator_units(df)
    df.to_csv(snakemake.output.indicators, index=False)
