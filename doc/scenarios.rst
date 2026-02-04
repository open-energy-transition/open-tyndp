.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _scenarios:

####################
Open TYNDP Scenarios
####################

The modelling for the reference, National Trends (NT), 
and the deviation, Distributed Energy (DE) and Global Ambition (GA), scenarios
differ in terms of assumptions, storyline and method (capacity expansion vs dispatch modelling).

Here we explain these modelling differences and these differences have been implemented into open-tyndp.
We discuss the relevant configuration settings and the implications of these implementation decisions.

Background information can be found in the report from ENTSO-E and ENTSO-G "TYNDP 2024 Scenarios Methdoology Report".

National Trends+
~~~~~~~~~~~~~~~~

- Uses predefined final energy demands and generation capacity collected from transmission system operators.
- All energy carriers included in final energy demand, not just electricity and gas (as with TYDNP < 2024)
- Run for 2030 and 2040

Deviation Scenarios
~~~~~~~~~~~~~~~~~~~

- Begins with National Trends+ 2030 scenario results for grid topology, generation capacity and final energy demands
- Uses predefined final energy demands for 2040 and 2050 collected from transmission system operators and public consultation
- Capacity expansion model is used to identify where and when investment in generation capacity is required
- Exogenous constraints on generation investment are imposed to force:
    - no new nuclear power capacity
    - all existing fossil gas plants are decommissioned
- upper and lower bounds are also imposed to force trajectories for
    - the expansion of solar, wind, prosumer batteries and large scale batteries
    - H2 import potentials