..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

.. _innovation_roadmap:

##########################################
Innovation Roadmap
##########################################

The tables below compare the current and upcoming features of the Open TYNDP workflow (and those already implemented in PyPSA-Eur) to the
`TYNDP Innovation Roadmap <https://tyndp.entsoe.eu/resources/tyndp-scenarios-innovation-roadmap>`_ which lists desirable features for the 2026 TYNDP cycle.

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

Innovations on the Energy Transition Model (ETM)
===============================================

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Innovation Roadmap Details
     - Comparison with PyPSA-Eur / Open-TYNDP   
   * - Dashboard: All graphs for the TYNDP 2026 report will be integrated into a reactive dashboard
     - Open TYNDP will feature a dashboard that allows users to explore the results of the scenarios interactively. PyPSA-Eur has introduced automated plotting of energy balance maps and heatmap time series. Recent updates include interactive bus-balance plots and heat-source maps.
   * - Improvement of Reference Values
           - Update residential space heating and hot water technology shares from national statistics rather than 2019 EUROSTAT energy statistics
     - PyPSA-Eur recently updated its energy balances to JRC-IDEES-2021, switching the reference year to 2019. Open-TYNDP has refined electricity demand and biomass potentials to achieve an exact match with reference values.
   * - Addition of climate year functionality
         - Integrate weather years into energy demand scenarios including heating demand (electricity demand from heat pumps and boilers) and all final energy demands. 
         - Produce a set of sectoral energy demand profiles for each scenario.
     - PyPSA-Eur now supports spanning multiple consecutive or meteorological weather years in a single optimization. It provides pre-built cutouts for a wide range of years (e.g., 1996, 2010–2023).
   * - Include missing non-EU countries 
         - Include missing non-EU countries such as Norway, Switzerland, and Serbia. 
         - Split the UK into separate datasets for Great Britain and Northern Ireland.
     - PyPSA-Eur covers the full ENTSO-E area and has recently integrated Ukraine, Moldova, and Kosovo. It supports NUTS-level clustering across these regions.
   * - Stable ETM server for 2026 cycle
         - ensure that a stable version of the ETM is available with consistent data and features for the duration of the 2026 cycle
     - Open-source frameworks achieve stability through version control. To manage computational environment, conda-lock files are used for dependency management and Snakemake version requirements 
       (e.g., minimum version 9.0) to ensure cross-platform reproducibility.
   * - Integrate supply tool features in ETM 
            - Merging supply modeling features (CCS, imports, biomass) directly into the ETM to reduce interfaces and provide a coherent energy system representation
     - Open TYNDP already integrates supply and demand features into a full integrated workflow. PyPSA-Eur and Open-TYNDP are inherently integrated cross-sectoral models. 
       They already include rules for CO2 sequestration potentials (CO2Stop), biomass transport costs, and automated import price configurations.
   * - Demand profile modelling in ETM
         - Model hourly methane and hydrogen profiles instead of using internal ENTSOG tool ensuring full consistency with scenario assumptions
         - Model hourly electricity demand profiles, creating a strong link to scenario parameters and demand profiles while leaving flexibility for TSOs to choose adoption.
     - Open-TYNDP has implemented the attachment of exogenous TYNDP gas and hydrogen demands to the network. PyPSA-Eur generates hourly profiles for all energy carriers as a core function of its sector-coupled workflow.

Pan-European Market Modelling Database App
==========================================

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Innovation Roadmap Details
     - Comparison with Open-TYNDP
   * - The PEMMDB app will provide an API to allow efficient data transfer into the PLEXOS model 
     - The Open TYNDP workflow already automates data integration though REST API calls to various sources, as well as automatically downloading, extracting, transforming and filtering data for all of its inputs.

Quality Control
===============

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Innovation Roadmap Details
     - Comparison with Open-TYNDP 
   * - Perform sanity checks on each scenario, for example preventing simultaneous dispatch of electrolysers and H2/gas-fired plants, ensuring reasonable levels of curtailment and Energy Not Served, comparing generator margins, corss-sector assets and storage investment costs
     - <Validation?>

System Modelling Innovations
============================

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Innovation Roadmap Details
     - Comparison with Open-TYNDP 
   * - Hydrogen storage
         - differentiate between short-term and seasonal storage needs
         - better representation of operational constraints of storage facilities such as salt caverns and aquifers
         - align national and modelling studies on storage capacities
         - incorporate techno-economic constraints of hydrogen supply
         - improve pipelines modelling to reflect transport flexibility
     - PyPSA-Eur recently added aquifer thermal energy storage (ATES) and fixed bugs related to underground H2 cavern creation. Open-TYNDP includes regionalized H2 salt cavern potentials.
   * - Integration of Hybrid heat pumps
         - Ensure hybrid systems are correctly sized for applications, considering peak demand scenarios
         - Ensure assumptions incorporate both economics and behavioural considerations
     - PyPSA-Eur refined heat pump CAPEX allocations and COP approximations. It also introduced river- and sea-water sourced heat pumps.
   * - Grid topology
         - Incorporate more detailed representation of H2 network topology that approximates physical H2 flow to ensure key H2 corridors and infrastructure are well represented
     - PyPSA-Eur supports administrative clustering (NUTS0 to NUTS3), allowing the network to be resolved at highly granular levels. Open-TYNDP uses the TYNDP H2 topology, including Z1 and Z2 setup.
   * - Gas turbine usage and peaking unit utilisation 
         - Explore dynamic operational needs of gas turbines given increasing reliance of variable renewable resources
         - Update assumptions around CH4 to H2 retrofitting projects, which have struggled to compete in current markets
     - <TODO>
   * - EV Modelling
         - Ensure that electricity flows follow charge/discharge cycles of EV batteries
     - PyPSA-Eur refined Vehicle-to-Grid (V2G) dispatch capacity and temperature-dependent energy demand correction factors for EVs.
   * - Economic Assessment
         - incorporate economic assessment of key technologies such as Steam Methane Reformers (SMR), nuclear plants, ammonia regasification terminals which are currently represented without economic attributes
     - <TODO>
   * - Methane pricing structure and formation
     -
   * - Ammonia Import costs
     -
   * - Distinction in Hydrogen usage
           - Separating hydrogen used directly as a gas from hydrogen used as a feedstock for producing synthetic fuels.
     - PyPSA-Eur added a suite of technologies for methanol-to-power, reforming, and kerosene, and updated locations/capacities for ammonia plants to accurately distribute demand.
   * - Flexibility of heat pumps
     -
   * - Modelling of E-fuels 
     -
   * - Higher granularity topology
           - Adding more nodes per country and differentiating between prosumer and non-prosumer households for accurate grid interaction modeling
     - PyPSA-Eur supports multiple resource classes for wind and solar per region to improve accuracy at low spatial resolutions. It also allows for behind-the-meter rooftop PV modeling.
   * - Improved modelling of prosumer demand
     -
   * - Consider peaking units as expansion candidates
     -
   * - Check on remaining CO2 emissions in 2050
     -
   * - Implementation of hybrid electrolyser plants
           - Implementing plants connected to both dedicated renewables and the grid to optimize production and market coupling
     - Open-TYNDP implements offshore wind hubs where wind farms can connect to both the network and P2G units for H2 production. PyPSA-Eur added minimum unit dispatch settings for electrolysis.
   * - Hydrogen imports and pipeline assessment
           - Modeling practical volumes and prices for pipeline imports to evaluate energy security and dependence.
     - PyPSA-Eur implemented renewable energy imports for H2, ammonia, methanol, and oil with configurable prices and volume limits.
   * - Geographical correlation in hydrogen production
     -

Stakeholder Reference Group (SRG) Proposals
===========================================

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Stakeholder Reference Group
     - Open TYNDP Innovation Roadmap
   * - Synthetic Fuels: Incorporating methanol for the maritime sector to align with decarbonization goals and identify infrastructure needs.
     - PyPSA-Eur introduced methanol-based technologies (e.g., biomass-to-methanol) in its 2024.09 release. Open-TYNDP defaults maritime demand to methanol.
   * - Climatic Variability: Suggesting models run with three different climatic years to assess impact on energy security.
     - PyPSA-Eur is designed for this; it integrates with atlite to process multi-year datasets and supports spanning these in a single model.
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
   * - Optimisation Across Energy Vectors
           - Expanding PLEXOS to integrate electricity, gas, and hydrogen systems holistically.
     - While the roadmap lists this as a goal, PyPSA-Eur and Open-TYNDP use the Linopy backend to optimize these vectors simultaneously by default.
   * - Transmission System Losses
           - Reassessing losses to reflect actual power flow dynamics more accurately.
     - PyPSA-Eur allows for piecewise linear approximation of transmission losses and provides the option to disable efficiency losses for specific carriers.
   * - Additional hydrogen production pathways
           - Integrating methane pyrolysis and waste-to-hydrogen processes.
     - PyPSA-Eur has already integrated biomass-to-hydrogen (with or without carbon capture) and supports custom technology adjustments via configuration.
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