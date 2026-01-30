# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Electricity configuration.

See docs in https://open-tyndp.readthedocs.io/en/latest/configuration.html#electricity
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from scripts.lib.validation.config._base import ConfigModel


class _OperationalReserveConfig(ConfigModel):
    """Configuration for `electricity.operational_reserve` settings."""

    activate: bool = Field(
        False,
        description="Whether to take operational reserve requirements into account during optimisation.",
    )
    epsilon_load: float = Field(
        0.02,
        description="share of total load.",
    )
    epsilon_vres: float = Field(
        0.02,
        description="share of total renewable supply.",
    )
    contingency: float = Field(
        4000,
        description="Fixed reserve capacity (MW).",
    )


class _MaxHoursConfig(BaseModel):
    """Configuration for `electricity.max_hours` settings."""

    battery: float = Field(
        6,
        description="Maximum state of charge capacity of the battery in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    H2: float = Field(
        168,
        description="Maximum state of charge capacity of the hydrogen storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )


class _ExtendableCarriersConfig(BaseModel):
    """Configuration for `electricity.extendable_carriers` settings."""

    Generator: list[str] = Field(
        default_factory=lambda: [
            "solar",
            "solar-hsat",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "offwind-float",
            "OCGT",
            "CCGT",
        ],
        description="Defines existing or non-existing conventional and renewable power plants to be extendable during the optimization. Conventional generators can only be built/expanded where already existent today. If a listed conventional carrier is not included in the `conventional_carriers` list, the lower limit of the capacity expansion is set to 0.",
    )
    StorageUnit: list[str] = Field(
        default_factory=list,
        description="Adds extendable storage units (battery and/or hydrogen) at every node/bus after clustering without capacity limits and with zero initial capacity.",
    )
    Store: list[str] = Field(
        default_factory=lambda: ["battery", "H2"],
        description="Adds extendable storage units (battery and/or hydrogen) at every node/bus after clustering without capacity limits and with zero initial capacity.",
    )
    Link: list[str] = Field(
        default_factory=list,
        description="Adds extendable links (H2 pipelines only) at every connection where there are lines or HVDC links without capacity limits and with zero initial capacity. Hydrogen pipelines require hydrogen storage to be modelled as `Store`.",
    )


class _TechnologyMappingConfig(BaseModel):
    """Configuration for `electricity.estimate_renewable_capacities.technology_mapping` settings."""

    Offshore: str = Field(
        "offwind-ac",
        description="PyPSA-Eur carrier that is considered for existing offshore wind technology (IRENA, GEM).",
    )
    Onshore: str = Field(
        "onwind",
        description="PyPSA-Eur carrier that is considered for existing onshore wind capacities (IRENA, GEM).",
    )
    PV: str = Field(
        "solar",
        description="PyPSA-Eur carrier that is considered for existing solar PV capacities (IRENA, GEM).",
    )


class _EstimateRenewableCapacitiesConfig(BaseModel):
    """Configuration for `electricity.estimate_renewable_capacities` settings."""

    enable: bool = Field(
        True,
        description="Activate routine to estimate renewable capacities in rule `add_electricity`. This option should not be used in combination with pathway planning `foresight: myopic` or `foresight: perfect` as renewable capacities are added differently in `add_existing_baseyear`.",
    )
    from_gem: bool = Field(
        True,
        description="Add renewable capacities from `Global Energy Monitor's Global Solar Power Tracker <https://globalenergymonitor.org/projects/global-solar-power-tracker/>`_ and `Global Energy Monitor's Global Wind Power Tracker <https://globalenergymonitor.org/projects/global-wind-power-tracker/>`_.",
    )
    year: int = Field(
        2020,
        description="Renewable capacities are based on existing capacities reported by IRENA (IRENASTAT) for the specified year.",
    )
    expansion_limit: float | bool = Field(
        False,
        description="Artificially limit maximum IRENA capacities to a factor. For example, an `expansion_limit: 1.1` means 110% of capacities. If false are chosen, the estimated renewable potentials determine by the workflow are used.",
    )
    technology_mapping: _TechnologyMappingConfig = Field(
        default_factory=_TechnologyMappingConfig,
        description="Mapping between PyPSA-Eur and powerplantmatching technology names.",
    )


class _AutarkyConfig(BaseModel):
    """Configuration for `electricity.autarky` settings."""

    enable: bool = Field(
        False,
        description="Require each node to be autarkic by removing all lines and links.",
    )
    by_country: bool = Field(
        False,
        description="Require each country to be autarkic by removing all cross-border lines and links. `electricity: autarky` must be enabled.",
    )


class _PecdPreBuiltConfig(BaseModel):
    """Configuration for `electricity.pecd_renewable_profiles.pre_built` settings."""

    cyears: list[int] = Field(
        default_factory=lambda: [1995, 2008, 2009],
        description="List of climate years to filter for when creating the PECD pre-built. Can be any year between 1982 and 2019.",
    )

    @field_validator("cyears")
    @classmethod
    def validate_cyears(cls, v: list[int]) -> list[int]:
        """Validate that climate years are between 1982 and 2019."""
        for year in v:
            if not 1982 <= year <= 2019:
                raise ValueError(f"Climate year {year} must be between 1982 and 2019")
        return v


class _PecdRenewableProfilesConfig(BaseModel):
    """Configuration for `electricity.pecd_renewable_profiles` settings."""

    enable: bool = Field(
        False,
        description="Activate PECD renewable profiles from TYNDP 2024 instead of default renewable profiles for specified renewable technologies below.",
    )
    pre_built: _PecdPreBuiltConfig = Field(
        default_factory=_PecdPreBuiltConfig,
        description="A new pre-built version is prepared when the configured PECD version (`tyndp_pecd`) is not already set to a pre-built version.",
    )
    fill_gaps_method: str = Field(
        "zero",
        description="The chosen method for filling gaps in the PECD data to the modelled nodes. Can be either `zero` to fill with zero values or any other aggregation method such as `mean`, `median`, `max` or similar.",
    )
    available_years: list[int] = Field(
        default_factory=lambda: [2030, 2040, 2050],
        description="List of years for which PECD data is available.",
    )
    technologies: list[
        Literal[
            "Wind_Offshore",
            "LFSolarPVRooftop",
            "LFSolarPVUtility",
            "CSP_noStorage",
            "CSP_withStorage",
            "Wind_Onshore",
        ]
    ] = Field(
        default_factory=lambda: [
            "Wind_Offshore",
            "LFSolarPVRooftop",
            "LFSolarPVUtility",
            "Wind_Onshore",
        ],
        description="The PECD technologies whose PECD profile is used. These PECD technologies and their profiles are mapped to `tyndp_renewable_carriers`. Please make sure that the mapping is included in the `tyndp_technology_map.csv` file.",
    )


class _PemmdbHydroProfilesConfig(BaseModel):
    """Configuration for `electricity.pemmdb_hydro_profiles` settings."""

    enable: bool = Field(
        False,
        description="Activate PEMMDB hydro inflow profiles from 2024 TYNDP instead of default hydro profiles.",
    )
    available_years: list[int] = Field(
        default_factory=lambda: [2030, 2040, 2050],
        description="List of years for which PEMMDB hydro inflows data is available.",
    )
    technologies: list[
        Literal["Run_of_River", "Pondage", "Reservoir", "PS_Open", "PS_Closed"]
    ] = Field(
        default_factory=lambda: [
            "Run_of_River",
            "Pondage",
            "Reservoir",
            "PS_Open",
            "PS_Closed",
        ],
        description="The hydro technologies for which PEMMDB hydro inflows data is used.",
    )


class _PemmdbCapacitiesConfig(BaseModel):
    """Configuration for `electricity.pemmdb_capacities` settings."""

    enable: bool = Field(
        False,
        description="Activate PEMMDB capacities, must-runs and availabilities from TYNDP 2024.",
    )
    nprocesses: int = Field(
        4,
        description="Number of parallel processes when reading pemmdb data.",
    )
    available_years: list[int] = Field(
        default_factory=lambda: [2030, 2040, 2050],
        description="List of years for which PEMMDB data is available.",
    )
    technologies: list[
        Literal[
            "Solar",
            "Wind",
            "Hydro",
            "Other_RES",
            "Gas",
            "Nuclear",
            "Hard_coal",
            "Lignite",
            "Light_oil",
            "Heavy_oil",
            "Oil_shale",
            "Other_Non-RES",
            "Hydrogen",
            "Electrolyser",
            "Battery",
        ]
    ] = Field(
        default_factory=lambda: [
            "Solar",
            "Wind",
            "Hydro",
            "Other_RES",
            "Gas",
            "Nuclear",
            "Hard_coal",
            "Lignite",
            "Light_oil",
            "Heavy_oil",
            "Oil_shale",
            "Other_Non-RES",
            "Hydrogen",
            "Electrolyser",
            "Battery",
        ],
        description="The PEMMDB technologies which are available and which data shall be used.",
    )


class ElectricityConfig(BaseModel):
    """Configuration for `electricity` settings."""

    voltages: list[float] = Field(
        default_factory=lambda: [220.0, 300.0, 330.0, 380.0, 400.0, 500.0, 750.0],
        description="Voltage levels to consider.",
    )
    base_network: Literal["entsoegridkit", "osm", "tyndp"] = Field(
        "osm",
        description="Specify the underlying base network, i.e. GridKit (based on ENTSO-E web map extract), OpenStreetMap (OSM), or TYNDP.",
    )
    gaslimit_enable: bool = Field(
        False,
        description="Add an overall absolute gas limit configured in `electricity: gaslimit`.",
    )
    gaslimit: float | bool = Field(
        False,
        description="Global gas usage limit.",
    )
    co2limit_enable: bool = Field(
        False,
        description="Add an overall absolute carbon-dioxide emissions limit configured in `electricity: co2limit` in `prepare_network`. **Warning:** This option should currently only be used with electricity-only networks, not for sector-coupled networks.",
    )
    co2limit: float = Field(
        7.75e7,
        description="Cap on total annual system carbon dioxide emissions.",
    )
    co2base: float = Field(
        1.487e9,
        description="Reference value of total annual system carbon dioxide emissions if relative emission reduction target is specified in `{opts}` wildcard.",
    )
    operational_reserve: _OperationalReserveConfig = Field(
        default_factory=_OperationalReserveConfig,
        description="Settings for reserve requirements following `GenX <https://genxproject.github.io/GenX/dev/core/#Reserves>`_.",
    )
    max_hours: _MaxHoursConfig = Field(
        default_factory=_MaxHoursConfig,
        description="Maximum storage hours configuration.",
    )
    extendable_carriers: _ExtendableCarriersConfig = Field(
        default_factory=_ExtendableCarriersConfig,
        description="Defines which carriers are extendable during optimization.",
    )
    powerplants_filter: str | bool = Field(
        "(DateOut >= 2024 or DateOut != DateOut) and not (Country == 'Germany' and Fueltype == 'Nuclear')",
        description="Filter query for the default powerplant database.",
    )
    custom_powerplants: str | bool = Field(
        False,
        description="Filter query for the custom powerplant database.",
    )
    everywhere_powerplants: list[str] = Field(
        default_factory=list,
        description="List of conventional power plants to add to every node in the model with zero initial capacity. To be used in combination with `extendable_carriers` to allow for building conventional powerplants irrespective of existing locations.",
    )
    conventional_carriers: list[str] = Field(
        default_factory=lambda: [
            "nuclear",
            "oil",
            "OCGT",
            "CCGT",
            "coal",
            "lignite",
            "geothermal",
            "biomass",
        ],
        description="List of conventional power plants to include in the model from `resources/powerplants_s_{clusters}.csv`. If an included carrier is also listed in `extendable_carriers`, the capacity is taken as a lower bound.",
    )
    tyndp_conventional_carriers: list[
        Literal[
            "gas",
            "coal",
            "lignite",
            "uranium",
            "oil-light",
            "oil-heavy",
            "oil-shale",
        ]
    ] = Field(
        default_factory=list,
        description="List of TYNDP conventional power plants to include in the model.",
    )
    group_tyndp_conventionals: bool = Field(
        True,
        description="Whether to group together different TYNDP conventional technologies according to a predetermined mapping specified in `data/tyndp_technology_map.csv`.",
    )
    renewable_carriers: list[str] = Field(
        default_factory=lambda: [
            "solar",
            "solar-hsat",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "offwind-float",
            "hydro",
        ],
        description="List of renewable generators to include in the model.",
    )
    tyndp_renewable_carriers: list[
        Literal[
            "solar-pv",
            "solar-pv-utility",
            "solar-pv-rooftop",
            "onwind",
            "offwind-ac-fb-r",
            "offwind-ac-fl-r",
            "offwind-dc-fb-r",
            "offwind-dc-fl-r",
            "offwind-dc-fb-oh",
            "offwind-dc-fl-oh",
            "offwind-h2-fb-oh",
            "offwind-h2-fl-oh",
            "other-res",
            "hydro-ror",
            "hydro-reservoir",
            "hydro-pondage",
            "hydro-phs",
            "hydro-phs-pure",
        ]
    ] = Field(
        default_factory=list,
        description="List of TYNDP renewable generators to include in the model. Technologies covered by specified TYNDP renewable carriers need to be removed from `estimate_renewable_carriers` technology list.",
    )
    tyndp_stores: list[Literal["battery", "h2_cavern", "h2_tank"]] = Field(
        default_factory=list,
        description="List of TYNDP storage units (battery and/or hydrogen).",
    )
    pecd_renewable_profiles: _PecdRenewableProfilesConfig = Field(
        default_factory=_PecdRenewableProfilesConfig,
        description="PECD (Pan-European Climate Database) renewable profiles configuration.",
    )
    pemmdb_hydro_profiles: _PemmdbHydroProfilesConfig = Field(
        default_factory=_PemmdbHydroProfilesConfig,
        description="PEMMDB hydro profiles configuration.",
    )
    pemmdb_capacities: _PemmdbCapacitiesConfig = Field(
        default_factory=_PemmdbCapacitiesConfig,
        description="PEMMDB capacities configuration.",
    )
    estimate_renewable_capacities: _EstimateRenewableCapacitiesConfig = Field(
        default_factory=_EstimateRenewableCapacitiesConfig,
        description="Configuration for estimating renewable capacities.",
    )
    autarky: _AutarkyConfig = Field(
        default_factory=_AutarkyConfig,
        description="Autarky configuration.",
    )
    transmission_limit: str = Field(
        "vopt",
        description="Limit on transmission expansion. The first part can be `v` (for setting a limit on line volume) or `c` (for setting a limit on line cost). The second part can be `opt` or a float bigger than one (e.g. 1.25). If `opt` is chosen line expansion is optimised according to its capital cost (where the choice `v` only considers overhead costs for HVDC transmission lines, while `c` uses more accurate costs distinguishing between overhead and underwater sections and including inverter pairs). The setting `v1.25` will limit the total volume of line expansion to 25% of currently installed capacities weighted by individual line lengths. The setting `c1.25` will allow to build a transmission network that costs no more than 25 % more than the current system.",
    )

    model_config = ConfigDict(populate_by_name=True)
