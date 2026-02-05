# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Cost-Benefit Analysis (CBA) configuration.

See docs in https://open-tyndp.readthedocs.io/en/latest/configuration.html#cba
"""

from typing import Any, Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _CbaStorageConfig(ConfigModel):
    """Configuration for `cba.storage` settings."""

    cyclic_carriers: list[str] = Field(
        default_factory=lambda: ["battery", "home battery"],
        description="Carriers that should remain cyclic (short-term storage). These will maintain e_cyclic=True in prepare_rolling_horizon.",
    )
    seasonal_carriers: list[str] = Field(
        default_factory=lambda: ["H2 Store", "gas", "co2 sequestered"],
        description="Seasonal carriers that will receive MSV (long-term storage). These will have cyclicity disabled and MSV applied as marginal_cost.",
    )


class _CbaMsvExtractionConfig(ConfigModel):
    """Configuration for `cba.msv_extraction` settings."""

    resolution: bool | str = Field(
        default=False,
        description="Temporal resolution for MSV extraction solve. False uses native resolution, or a string like '24H', '48H' for faster solve.",
    )
    resample_method: Literal["ffill", "interpolate"] = Field(
        default="ffill",
        description="Method for resampling MSV to target network resolution.",
    )


class _CbaSolvingConfig(ConfigModel):
    """Configuration for `cba.solving` settings."""

    options: dict[str, Any] = Field(
        default_factory=lambda: {
            "horizon": 168,
            "overlap": 0,
            "load_shedding": True,
            "io_api": "direct",
        },
        description="CBA-specific solving options.",
    )
    solver: dict[str, str] = Field(
        default_factory=lambda: {"name": "highs", "options": "highs-simplex"},
        description="CBA solver configuration.",
    )
    solver_options: dict[str, dict[str, Any]] = Field(
        default_factory=dict,
        description="CBA solver options.",
    )


class CbaConfig(BaseModel):
    """Configuration for top level `cba` (Cost-Benefit Analysis) settings."""

    hurdle_costs: float = Field(
        0.01,
        description="Marginal cost for transmission lines in Cost-Benefit Analysis networks (EUR/MWh). Used to represent congestion costs and market inefficiencies.",
    )
    co2_societal_cost: dict[int, dict[str, float]] = Field(
        default_factory=dict,
        description="Dictionary mapping planning horizons to societal costs of CO2 emissions (EUR/t) for 'low', 'central', and 'high' scenarios. Used in the CBA workflow to assess the societal cost of CO2 emissions (B2 indicator).",
    )
    planning_horizons: list[int] = Field(
        default_factory=list,
        description="List of planning horizons for which to run CBA.",
    )
    methods: list[Literal["toot", "pint"]] = Field(
        default_factory=lambda: ["toot"],
        description="CBA methodologies to apply. Options include 'toot' (Take One Out at a Time) and 'pint' (Put In at a Time).",
    )
    projects: list[str] = Field(
        default_factory=list,
        description="List of TYNDP project identifiers to evaluate in the CBA workflow. Format follows TYNDP project naming conventions (e.g., 't1-t35').",
    )
    storage: _CbaStorageConfig = Field(
        default_factory=_CbaStorageConfig,
        description="Storage configuration for CBA workflow.",
    )
    msv_extraction: _CbaMsvExtractionConfig = Field(
        default_factory=_CbaMsvExtractionConfig,
        description="MSV extraction settings for seasonal storage dispatch.",
    )
    solving: _CbaSolvingConfig = Field(
        default_factory=_CbaSolvingConfig,
        description="Configuration options specific to CBA network optimization using rolling horizon approach. Uses the same structure as the top-level `solving` section with CBA-specific defaults for `options: horizon` (168 hours), `options: io_api` (direct), and `solver: name` (highs).",
    )
