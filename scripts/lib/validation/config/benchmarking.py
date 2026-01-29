# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Benchmarking configuration.

See docs in https://open-tyndp.readthedocs.io/en/latest/configuration.html#benchmarking
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _BenchmarkingTableNamesConfig(ConfigModel):
    """Configuration for `benchmarking.tables.{table_name}.names` settings."""

    index: list[str] = Field(
        default_factory=list,
        description="Names for index levels corresponding to `index_col`.",
    )
    column: list[str] = Field(
        default_factory=list,
        description="Names for column levels corresponding to `header`.",
    )


class _BenchmarkingTableConfig(ConfigModel):
    """Configuration for `benchmarking.tables.{table_name}` settings."""

    table_type: Literal["scenario_comparison", "time_series"] = Field(
        ...,
        description="Table type configuration to use.",
    )
    sheet_name: str = Field(
        ...,
        description="Excel sheet name containing reference data.",
    )
    unit: str = Field(
        ...,
        description="Data unit (e.g. TWh or GW) for conversion to the base unit using the unit conversion factors defined in `benchmarking:unit_conversion`.",
    )
    index_col: int | list[int] = Field(
        ...,
        description="Column index or list of column indices to use as row index when reading tables. This is an Excel-like index, starting at 1.",
    )
    header: int | list[int] = Field(
        ...,
        description="Row number or list of row numbers to use as column headers when reading tables. This is an Excel-like index, starting at 1.",
    )
    nrows: int | None = Field(
        None,
        description="Number of rows to read from reference data (optional).",
    )
    ncolumns: int | None = Field(
        None,
        description="Number of columns to read from reference data (optional).",
    )
    names: _BenchmarkingTableNamesConfig = Field(
        default_factory=_BenchmarkingTableNamesConfig,
        description="Names configuration for index and column levels.",
    )
    mapping: dict[str, str] = Field(
        default_factory=dict,
        description="Mapping from model carrier names to reference data names.",
    )
    biomass_to_methane_efficiency: dict[int, float] | None = Field(
        None,
        description="Efficiency to convert biomass to biomethane. Dictionary with planning horizons as keys. Only used for biomass supply table benchmarking.",
    )


class BenchmarkingConfig(BaseModel):
    """Configuration for top level `benchmarking` settings."""

    enable: bool = Field(
        False,
        description="Switch to enable or disable benchmarking functionality for comparing model results with reference data.",
    )
    remove_last_day: bool = Field(
        False,
        description="Switch to optionally remove the last day of the year to ensure the benchmarked values have exactly 52 weeks.",
    )
    tables: dict[str, _BenchmarkingTableConfig] = Field(
        default_factory=dict,
        description="Tables configuration for benchmarking. Each key is a table name, and the value contains the table configuration.",
    )
