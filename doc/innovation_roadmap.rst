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

     **Current Implementation Scope:**
While PyPSA-Eur theoretically has the capability to cover all components of the TYNDP toolchain (including calculating capacity factors, heat demand time series, and total annual demands), the current Open-TYNDP implementation for Scenario Building focuses on:

- Replicating the functionalities of PLEXOS as the core market simulation tool
- Providing automated workflow for input data processing via Snakemake
- Implementing visualization of the results 
We use existing outputs from Supply Tool and DFT as inputs, rather than replacing these tools in the Open-TYNDP.

**Potential for Full Integration:**
Supply Tool and DFT could be replaced within the Open-TYNDP framework, but this would require:

- Code adaptations to be integrated back into PyPSA-Eur core functionality
- Explicit assumptions about technology specifications (e.g., types of onshore wind farms, PV panel characteristics) to generate accurate capacity factors
- Either adopting PyPSA-Eur's standard assumptions or documenting a new set of transparent assumptions (noting that existing TYNDP assumptions behind Supply Tool/DFT are not publicly available)

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
       They already include rules for CO2 sequestration potentials (CO2Stop), CCS technologies, carbon networks, regional biomass potentials and transport costs, and automated import price configurations.
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
     - PyPSA-Eur recently fixed bugs related to underground H2 cavern creation. Open-TYNDP includes regionalized H2 salt cavern potentials.
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
     - Open TYNDP allows for endogenous retrofitting of CH4 pipeline capacity to H2.
   * - EV Modelling
         - Ensure that electricity flows follow charge/discharge cycles of EV batteries
     - PyPSA-Eur refined Vehicle-to-Grid (V2G) dispatch capacity and temperature-dependent energy demand correction factors for EVs.
   * - Economic Assessment
         - incorporate economic assessment of key technologies such as Steam Methane Reformers (SMR), nuclear plants, ammonia regasification terminals which are currently represented without economic attributes
     - Open TYNDP incorporates default capital investment and fixed and variable operation cost assumptions for all technologies
   * - Methane pricing structure and formation
            - Resolve pricing inconsistences between synthetic natural gas, biomethane and hydrogen
     - TODO
   * - Ammonia Import costs
           - Ensure ammonia import costs reflect the entire supply chain
     - Open TYNDP includes location and capacities of European ammonia plants. Ammonia import prices and volumes can be configured.
   * - Distinction in Hydrogen usage
           - Separating hydrogen used directly as a gas from hydrogen used as a feedstock for producing synthetic fuels.
     - PyPSA-Eur added a suite of technologies for methanol-to-power, reforming, and kerosene, and updated locations/capacities for ammonia plants to accurately distribute demand.
   * - Flexibility of heat pumps
           - Incorporate heat-pump modelling into PLEXOS to better represent flexibility, thermal inertia and heat storage
     - Open TYNDP uses an energy-to-power ratio constraiint for thermal energy storage. Proportional sizing is used for chargers and dischargers.
       Includes an option to calculate dynamic storage capacities for thermal energy storage. Also includes aquifer thermal energy storage, and supplemental heating such as booster heat pumps.
   * - Modelling of E-fuels 
           - Allow transportation of e-methanol, e-methane and e-kerosene via pipelines, ships or tanker trucks
           - Incorporate full load hours into spatial optimisation of e-fuel refineries and supply infrastructure
     - <TODO>
   * - Higher granularity topology
           - Adding more nodes per country and differentiating between prosumer and non-prosumer households for accurate grid interaction modeling
     - PyPSA-Eur supports multiple resource classes for wind and solar per region to improve accuracy at low spatial resolutions. It also allows for behind-the-meter rooftop PV modeling.
   * - Improved modelling of prosumer demand
     - Open TYNDP allows connection of microgeneration e.g. residential solar PV to be connected to low voltage buses.
       Residential and utility scale PV are treated separately, with separate rules in the workflow to build and cluster rooftop potentials.
       Open TYNDP also distinguishes between stationary battery storage and EV batteries.
   * - Consider peaking units as expansion candidates
     - PyPSA-Eur is a capacity expansion model by nature and can be set to build new peaking units whenever they are the cost-optimal way to ensure reliability.
   * - Check on remaining CO2 emissions in 2050
     - TODO
   * - Implementation of hybrid electrolyser plants
           - Implementing plants connected to both dedicated renewables and the grid to optimize production and market coupling
     - Open-TYNDP implements offshore wind hubs where wind farms can connect to both the network and P2G units for H2 production. PyPSA-Eur added minimum unit dispatch settings for electrolysis.
   * - Hydrogen imports and pipeline assessment
           - Modeling practical volumes and prices for pipeline imports to evaluate energy security and dependence.
     - PyPSA-Eur implemented renewable energy imports for H2, ammonia, methanol, and oil with configurable prices and volume limits.
   * - Geographical correlation in hydrogen production
     - Open-TYNDP ensures geographical correlation by attaching planning-year dependent renewable profiles from the PECD to specific generators within interconnected zones

Stakeholder Reference Group (SRG) Proposals
===========================================

.. list-table::
   :widths: 50 50
   :header-rows: 1    

   * - Stakeholder Reference Group
     - Open TYNDP Innovation Roadmap
   * - Synthetic Fuels
           - Incorporating methanol for the maritime sector to align with decarbonization goals and identify infrastructure needs.
     - PyPSA-Eur introduced methanol-based technologies (e.g., biomass-to-methanol) in its 2024.09 release. Open-TYNDP defaults maritime demand to methanol.
   * - Climatic Variability
           - Suggesting models run with three different climatic years to assess impact on energy security.
     - PyPSA-Eur is designed for this; it integrates with atlite to process multi-year datasets and supports spanning these in a single model.
   * - Industrial applications
           - Verify technical and commercial viability of converting industrial gas offtakes to H2 or other carriers before grid expansion
     - Open TYNDP uses a script to interpolate industry sector transition pathways, gradually switching processes from status quo to best-in-class energy consumption per ton of material output
   * - Sector-specific modelling
           - Discuss the Z1 Z2 concept to streamline management across gas, electricity, and hydrogen vectors
     - Open-TYNDP has already introduced the TYNDP H2 topology, which specifically includes the H2 Z1 and Z2 setup, production, and storage technologies
   * - EV modelling techniques
           - refine assumptions on EV charging behavior and their impact on potential grid bottlenecks
     - PyPSA-Eur now limits Vehicle-to-Grid (V2G) dispatch capacity based on the fraction of vehicles participating in demand-side management. It also refines temperature-dependent correction factors for EV energy demand.
   * - District heating
           - Create a dedicated tool distinct from the ETM to simulate production from biomass, geothermal, and other sources
     - PyPSA-Eur already features a highly detailed district heating module. Recent additions include geothermal district heating, aquifer thermal energy storage (ATES), and booster heat pumps for supplemental heating
   * - Liquified hydrogen
           - Explore LH2 import methods to understand logistical, storage, and cost constraints
     - PyPSA-Eur implements a "H2 liquid" bus at each location to specifically handle hydrogen liquefaction costs for shipping demand
   * - H2 Import Quotas
           - Align hydrogen import quotas with RepowerEU targets to avoid overestimating domestic production
     - 
   * - Electric heat pumps
           - Move heat pump modeling to PLEXOS to better capture thermal inertia and load management
     - PyPSA-Eur manages this via an energy-to-power ratio constraint for thermal storage, linking storage capacity to charger capacity to simulate load shifting
   * - Optimisation Across Energy Vectors
           - Expanding PLEXOS to integrate electricity, gas, and hydrogen systems holistically.
     - PyPSA-Eur and Open-TYNDP optimize these vectors simultaneously by default.
   * - Transmission System Losses
           - Reassessing losses to reflect actual power flow dynamics more accurately.
     - PyPSA-Eur allows for piecewise linear approximation of transmission losses and provides the option to disable efficiency losses for specific carriers.
   * - Additional hydrogen production pathways
           - Integrating methane pyrolysis and waste-to-hydrogen processes.
     - PyPSA-Eur has already integrated biomass-to-hydrogen (with or without carbon capture) and supports custom technology adjustments via configuration.
   * - Flexibility in modelling
           - Focus on load displacement (shifting demand) rather than just load reduction to maximize renewable use
     - Open TYNDP allows for demand side management functionalities to shift loads in time as well as in magnitude.
   * - Price setting for hydrogen
           - Revise methodology to reflect real-world contracts and costs like dehydrogenation
     - As an integrated sector-coupled model, endogenous pricing of hydrogen includes all represented upstream processes
   * - Sensitivity to commodity prices
           - Conduct sensitivity analyses on price fluctuations (gas, oil, H2) to understand investment risks
     - Workflow management tool snakemake enables the simultaneous execution of multiple scenarios with single calls and configuration overrides
   * - Inclusion of emerging technologies
     - TODO
   * - Out of scope Innovations
     - TODO
