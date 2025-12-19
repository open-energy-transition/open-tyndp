..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

.. _innovation_roadmap:

##########################################
Innovation Roadmap
##########################################

The [TYNDP Innovation Roadmap]()

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - TYNDP Innovation Roadmap
     - Open TYNDP Innovation Roadmap
   * - TYNDP 2024 brings together a suite of tools into a toolchain to develop the scenarios. These include:
       
       - Energy Transition Model
       - Supply Tool
       - DFT
       - PLEXOS
       - Visualisation Platform
       - Data files

     - Open TYNDP integrates all these steps into a single integrated snakemake workflow. This ensures full reproducibility of the scenario development process and allows for easy adaptation and extension of the workflow for future TYNDPs.
   

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Energy Transition Model
     - Open TYNDP Innovation Roadmap   
   * - Dashboard: All graphs for the TYNDP 2026 report will be integrated into the dashboard
     - Open TYNDP will feature a dashboard that allows users to explore the results of the scenarios interactively.
   * - Improvement of Reference Values: Update residential space heating and hot water technology shares from national statistics rather than 2019 EUROSTAT energy statistics
     - 
   * - Addition of climate year functionality: Integrate weather years into energy demand scenarios including heating demand (electricity demand from heat pumps and boilers) and all final energy demands. Produce a set of sectoral energy demand profiles for each scenario.
     - 
   * - Include missing non-EU countries such as Norway, Switzerland, and Serbia. Split the UK into separate datasets for Great Britain and Northern Ireland.
     - 
   * - Stable ETM server for 2026 cycle - ensure that a stable version of the ETM is available with consistent data and features for the duration of the 2026 cycle
     - label
   * - Integrate supply tool features in ETM 
     - Open TYNDP already integrates supply and demand features into a full integrated workflow 
   * - Demand profile modelling in ETM
         - Model hourly methane and hydrogen profiles instead of using internal ENTSOG tool ensuring full consistency with scenario assumptions
         - Model hourly electricity demand profiles, creating a strong link to scenario parameters and demand profiles while leaving flexibility for TSOs to choose adoption.
     - Open TYNDP already integrates hourly methane, hydrogen and electricity demand profiles, building on rigorous open-source methodologies

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   