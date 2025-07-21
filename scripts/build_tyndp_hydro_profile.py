# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Build hydroelectric inflow time-series for each country based on TYNDP hydro inflow data.

Outputs
-------

- ``resources/profile_hydro_tyndp.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               countries, time,   Inflow to the state of charge (in MW),
                         year, hydro_tech   e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================
"""
