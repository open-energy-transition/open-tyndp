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
   * - Addition of climate year functionality
       
       - Integrate weather years into energy demand scenarios including heating demand (electricity demand from heat pumps and boilers) and all final energy demands. 
       - Produce a set of sectoral energy demand profiles for each scenario.
     - 
   * - Include missing non-EU countries 
         - Include missing non-EU countries such as Norway, Switzerland, and Serbia. 
         - Split the UK into separate datasets for Great Britain and Northern Ireland.
     - 
   * - Stable ETM server for 2026 cycle
         - ensure that a stable version of the ETM is available with consistent data and features for the duration of the 2026 cycle
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

   * - Pan-European Market Modelling Database App
     - Open TYNDP Innovation Roadmap
   * - The PEMMDB app will provide an API to allow efficient data transfer into the PLEXOS model 
     - The Open TYNDP workflow already automates data integration though REST API calls to various sources, as well as automatically downloading, extracting, transforming and filtering data for all of its inputs.

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Quality Control
     - Open TYNDP Innovation Roadmap
   * - Perform sanity checks on each scenario, for example preventing simultaneous dispatch of electrolysers and H2/gas-fired plants, ensuring reasonable levels of curtailment and Energy Not Served, comparing generator margins, corss-sector assets and storage investment costs
     - <Validation?>

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - System Modelling Innovations
     - Open TYNDP Innovation Roadmap
   * - Hydrogen storage
         - differentiate between short-term and seasonal storage needs
         - better representation of operational constraints of storage facilities such as salt caverns and aquifers
         - align national and modelling studies on storage capacities
         - incorporate techno-economic constraints of hydrogen supply
         - improve pipelines modelling to reflect transport flexibility
     - <TODO>
   * - Integration of Hybrid heat pumps
         - Ensure hybrid systems are correctly sized for applications, considering peak demand scenarios
         - Ensure assumptions incorporate both economics and behavioural considerations
     - <TODO>
   * - Grid topology
         - Incorporate more detailed representation of H2 network topology that approximates physical H2 flow to ensure key H2 corridors and infrastructure are well represented
     - <TODO>
   * - Gas turbine usage and peaking unit utilisation 
         - Explore dynamic operational needs of gas turbines given increasing reliance of variable renewable resources
         - Update assumptions around CH4 to H2 retrofitting projects, which have struggled to compete in current markets
     - <TODO>
   * - EV Modelling
         - Ensure that electricity flows follow charge/discharge cycles of EV batteries
     - <TODO>
   * - Economic Assessment
         - incorporate economic assessment of key technologies such as Steam Methane Reformers (SMR), nuclear plants, ammonia regasification terminals which are currently represented without economic attributes
     - <TODO>
   * - Methane pricing structure and formation
     -
   * - Ammonia Import costs
     -
   * - Distinction in Hydrogen usage
     -
   * - Flexibility of heat pumps
     -
   * - Modelling of E-fuels 
     -
   * - Higher granularity topology
     -
   * - Improved modelling of prosumer demand
     -
   * - Consider peaking units as expansion candidates
     -
   * - Check on remaining CO2 emissions in 2050
     -
   * - Implementation of hybrid electrolyser plants
     -
   * - Hydrogen imports and pipeline assessment
     -
   * - Geographical correlation in hydrogen production
     -

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Stakeholder Reference Group
     - Open TYNDP Innovation Roadmap
   * - Synthetic fuels
     -
   * - Climate Variability
     -
   * - Industrial applications
     -
   * - Sector-specific modelling
     -
   * - EV modelling techniques
     -
   * - District heating
     -
   * - Liquified hydrogen
     -
   * - H2 Import Quotas
     -
   * - Electric heat pumps
     -
   * - Optimisation across energy vectors
     -
   * - Transmission system Losses
     -
   * - Additional hydrogen production pathways
     -
   * - Flexibility in modelling
     -
   * - Price setting for hydrogen
     -
   * - Sensitivity to commodity prices
     -
   * - Inclusion of emerging technologies
     -
   * - Out of scope Innovations
     -
