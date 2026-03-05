# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Cost-benefit analysis configuration.

See docs in https://open-tyndp.readthedocs.io/en/latest/configuration.html#cba
"""

from typing import Any, Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _CbaStorageConfig(ConfigModel):
    """Configuration for `cba.storage` settings."""

    cyclic_carriers: list[str] = Field(
        default_factory=lambda: ["battery", "home battery"],
        description="Carriers that should remain cyclic (short-term storage). All other store and storage unit carriers automatically receive marginal storage value and have cyclicity disabled.",
    )
    soc_boundary_carriers: list[str] = Field(
        default_factory=lambda: ["hydro-reservoir"],
        description="Storage unit carriers for which the state of charge is pinned at the boundaries between rolling horizon windows, using values pre-computed from the perfect foresight (full-year) optimisation.",
    )


class _CbaMsvExtractionConfig(ConfigModel):
    """Configuration for `cba.msv_extraction` settings."""

    resolution: bool | str = Field(
        default=False,
        description="Temporal resolution for extraction solve. False uses native resolution, or a string like '24H', '48H' for faster solve.",
    )
    resample_method: Literal["ffill", "interpolate"] = Field(
        default="ffill",
        description="Method for resampling marginal storage value to target network resolution.",
    )


class _CbaSolvingConfig(ConfigModel):
    """Configuration for `cba.solving` settings."""

    options: dict[str, Any] = Field(
        default_factory=lambda: {
            "horizon": 168,
            "overlap": 1,
            "load_shedding": {"enable": True},
            "remove_noisy_costs": True,
            "io_api": "direct",
        },
        description="Solving options for rolling horizon dispatch.",
    )
    solver: dict[str, str] = Field(
        default_factory=lambda: {"name": "highs", "options": "highs-simplex"},
        description="Solver configuration.",
    )
    solver_options: dict[str, dict[str, Any]] = Field(
        default_factory=dict,
        description="Solver-specific options.",
    )


class _CbaSbToCbaConfig(ConfigModel):
    """Configuration for using pre-solved SB networks in the CBA workflow."""

    use_presolved: bool = Field(
        False,
        description="If true, use a pre-solved SB network from an external archive instead of running the SB workflow.",
    )
    sb_version: str = Field(
        "latest",
        description="Version tag for the pre-solved SB network archive (used when cba_scenario_input.use_presolved is true).",
    )


class CbaConfig(BaseModel):
    """Configuration for top level `cba` (cost-benefit analysis) settings."""

    hurdle_costs: float = Field(
        0.01,
        description="Marginal cost for transmission lines in cost-benefit analysis networks (EUR/MWh).",
    )
    co2_societal_cost: dict[int, dict[str, float]] = Field(
        default_factory=dict,
        description="Dictionary mapping planning horizons to societal costs of CO2 emissions (EUR/t) for 'low', 'central', and 'high' scenarios.",
    )
    planning_horizons: list[int] = Field(
        default_factory=list,
        description="List of planning horizons for which to run cost-benefit analysis.",
    )
    cba_scenario_input: _CbaSbToCbaConfig = Field(
        default_factory=_CbaSbToCbaConfig,
        description="Settings for using pre-solved SB networks as inputs to the CBA workflow.",
    )
    methods: list[Literal["toot", "pint"]] = Field(
        default_factory=lambda: ["toot"],
        description="Methodologies to apply: 'toot' (take one out at a time) and 'pint' (put in one at a time).",
    )
    projects: list[str] = Field(
        default_factory=list,
        description="List of project identifiers to evaluate (e.g., 't1-t35').",
    )
    area: Literal["tyndp", "entso-e", "eu27"] = Field(
        default="tyndp",
        description="Geographical area for cost-benefit analysis. Options include 'tyndp', 'entso-e', and 'eu27'.",
    )
    storage: _CbaStorageConfig = Field(
        default_factory=_CbaStorageConfig,
        description="Storage configuration for the cost-benefit analysis workflow.",
    )
    msv_extraction: _CbaMsvExtractionConfig = Field(
        default_factory=_CbaMsvExtractionConfig,
        description="Marginal storage value extraction settings for seasonal storage dispatch.",
    )
    solving: _CbaSolvingConfig = Field(
        default_factory=_CbaSolvingConfig,
        description="Configuration for rolling horizon network optimization. Uses the same structure as the top-level `solving` section.",
    )
