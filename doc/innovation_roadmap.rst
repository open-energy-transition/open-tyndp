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
   :widths: 33 33 34
   :header-rows: 1

   * - TYNDP 2024
     - Current Open TYNDP Implementation  
     - Comparison with PyPSA-Eur
   * - TYNDP 2024 brings together a suite of tools into a toolchain to develop the scenarios. These include:
           - Energy Transition Model
           - Supply Tool
           - DFT
           - PLEXOS
           - Visualisation Platform
           - Data files
           - Antares*
           - Promed (internal)*
           - APG Tool (internal)*
       Asterix (*) marks tools used only for Cost Benefit Analysis
     - The current Open-TYNDP implementation for Scenario Building focuses on:
           - Replicating the functionalities of PLEXOS as the core market simulation tool
           - Providing automated workflow for input data processing via Snakemake
           - Implementing visualization of the results

       We use existing outputs from Supply Tool and DFT as inputs, rather than replacing these tools in the Open-TYNDP.

     - PyPSA-Eur has the capability to cover all components of the TYNDP toolchain (including calculating capacity factors, heat demand time series, and total annual demands).
       
       Supply Tool and DFT could be replaced within the Open-TYNDP framework, but this would require:

       - Code adaptations to be integrated back into PyPSA-Eur core functionality
       - Explicit assumptions about technology specifications (e.g., types of onshore wind farms, PV panel characteristics) to generate accurate capacity factors
       - Either adopting PyPSA-Eur's standard assumptions or documenting a new set of transparent assumptions (noting that existing TYNDP assumptions behind Supply Tool/DFT are not publicly available)

Innovations on the Energy Transition Model (ETM)
================================================

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Innovation Roadmap Details
     - Current Open TYNDP Implementation  
     - Features available in PyPSA-Eur
   * - 5.1 Dashboard
           - All graphs for the TYNDP 2026 report will be integrated into a reactive dashboard
     - Open TYNDP will feature a dashboard that allows users to explore the results of the scenarios interactively
     - PyPSA-Eur offers automated plotting of energy balance maps and heatmap time series. Recent updates include interactive bus-balance plots and heat-source maps
   * - 5.2 Improvement of Reference Values
           - Update residential space heating and hot water technology shares using more granular national data rather than 2019 EUROSTAT energy statistics
     - Open-TYNDP has refined electricity demand and biomass potentials to achieve an exact match with reference values
     - This challenge is shared by Open TYNDP. For the existing heating technologies in terms of capacities per country the `DG ENER mapping analysis <https://energy.ec.europa.eu/publications/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment-fossilrenewables-1_en>`_ is used, an update is needed.
       For district heating shares `Fraunhofer ISI data <https://fordatis.fraunhofer.de/handle/fordatis/341.5>`_.
       This data source is also used for geothermal heat potentials.
       PyPSA-Eur recently updated its energy balances to JRC-IDEES-2021, switching the reference year to 2019.
   * - 5.3 Addition of climate year functionality
         - Integrate weather years into energy demand scenarios including heating demand (electricity demand from heat pumps and boilers) and all final energy demands 
         - Produce a set of sectoral energy demand profiles for each scenario
     - The final energy demand is based on historical data in PyPSA-Eur from Eurostat
     - PyPSA-Eur now supports spanning multiple consecutive or meteorological weather years in a single optimization. It provides pre-built cutouts for a wide range of years (e.g., 1996, 2010–2023). Total electricity demand is endogenous in the optimisation, as it is calculated from energy service demands and the endogenous choice of demand technologies (such as heat-pumps or gas boilers). The performance of the demand technologies as well as the demand itself is therefore a function of weather parameters calculated from climate scenarios. This could be integrated into the automated Open TYNDP workflow instead of using exogenous demands for e.g. natural gas.
   * - 5.4 Include missing non-EU countries 
           - Include missing non-EU countries such as Norway, Switzerland, and Serbia
           - Split the UK into separate datasets for Great Britain and Northern Ireland
     - Open TYNDP covers the full ENTSO-E area and has recently integrated Ukraine, Moldova, and Kosovo
     - PyPSA-Eur supports NUTS-level clustering across these regions. This makes it relatively easy to create a spatial representation including Norway, Switzerland and Serbia (or any other countries).
   * - 5.5 Stable ETM server for 2026 cycle
           - ensure that a stable version of the ETM is available with consistent data and features for the duration of the 2026 cycle
     - Open-source frameworks achieve stability through version control. To manage computational environment, conda-lock files are used for dependency management and Snakemake version requirements 
       (e.g., minimum version 9.0) to ensure cross-platform reproducibility.
     - ditto
   * - 5.6 Integrate supply tool features in ETM 
            - Merging supply modeling features (CCS, imports, biomass) directly into the ETM to reduce interfaces and provide a coherent energy system representation
     - Open TYNDP currently uses exogenous demand assumptions in line with TYNDP 2024
     - PyPSA-Eur already integrates supply and demand features into a full integrated workflow.
       PyPSA-Eur and Open-TYNDP are inherently integrated cross-sectoral models. 
       They already include CO_2 sequestration potentials (``CO2Stop``), and optimising CO2 storage, CCS technologies, carbon networks, regional biomass potentials and transport costs, 
       and automated import price configurations.
   * - 5.7 Demand profile modelling in ETM
           - Model hourly methane and hydrogen profiles instead of using internal ENTSOG tool ensuring full consistency with scenario assumptions
           - Model hourly electricity demand profiles, creating a strong link to scenario parameters and demand profiles while leaving flexibility for TSOs to choose adoption
     - Open-TYNDP has implemented the attachment of exogenous TYNDP gas and hydrogen demands to the network.
     - Unlike TYNDP's fixed hydrogen and methane demands, PyPSA-Eur models final energy service demands (industrial process heat, heating, etc.) and endogenously optimises which energy carrier (hydrogen, methane, electricity) meets each demand. This allows the model to capture technology competition and fuel switching based on relative costs and availability.

Pan-European Market Modelling Database App
==========================================

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Innovation Roadmap Details
     - Current Open TYNDP Implementation  
     - Features available in PyPSA-Eur
   * - The PEMMDB app will provide an API to allow efficient data transfer into the PLEXOS model 
     - The Open TYNDP workflow already automates data integration through REST API calls to various sources, as well as automatically downloading, extracting, transforming and filtering data for all of its inputs.
     - All features are already implemented in Open TYNDP

Quality Control
===============

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Innovation Roadmap Details
     - Current Open TYNDP Implementation  
     - Features available in PyPSA-Eur
   * - Perform sanity checks on each scenario, for example preventing simultaneous dispatch of electrolysers and H2/gas-fired plants, ensuring reasonable levels of curtailment and Energy Not Served, 
       comparing generator margins, corss-sector assets and storage investment costs
     - Open TYNDP already offers:
           - Automated plotting of energy balances for each carrier, enabling quick identification of potential issues such as simultaneous electrolyser operation and hydrogen consumption
           - Curtailment with technology-specific values
           - Comprehensive n.statistics module for rapid evaluation key metrics
           - Basic consistency checking via n.consistency_check() to validate network topology and parameter ranges
     - Further development work could allow automated sanity checks with warnings for scenario-specific thresholds (e.g., flagging curtailment levels above 15%)

System Modelling Innovations
============================

The bulk of innovations are listed under system modelling.

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Innovation Roadmap Details
     - Current Open TYNDP Implementation  
     - Features available in PyPSA-Eur
   * - 8.1 Hydrogen storage
           - differentiate between short-term and seasonal storage needs
           - better representation of operational constraints of storage facilities such as salt caverns and aquifers
           - align national and modelling studies on storage capacities
           - incorporate techno-economic constraints of hydrogen supply
           - improve pipelines modelling to reflect transport flexibility
     - Open-TYNDP follows the same methodology as TYNDP for hydrogen infrastructure.
     - One could adapt back to PyPSA-Eur assumptions with some small modifications. PyPSA-Eur already implements:
           - Long-term storage in salt caverns and aquifers, and short-term storage in medium or high-pressure steel tanks
           - Technical potential constraints - Underground storage capacities in salt caverns are limited based on technical potential estimations from Caglayan et al. (2020)
           - Environmental considerations - Default configuration includes only nearshore underground storage potential (within 50 km of shore) due to environmental concerns regarding highly saline brine disposal. 
             However, the model offers flexibility to select from three storage types: nearshore (<50 km from sea), onshore (>50 km from sea), or offshore
           - Pipeline modeling - Endogenous optimization of methane-to-hydrogen pipeline retrofitting, hydrogen transport losses, and transport flow modeling in the pipeline network
       The following could be implemented with minor code extensions:
           - Minimum hydrogen storage fill levels to account for working gas requirements
           - Electrolyser ramp-up and ramp-down constraints reflecting technical flexibility limits
           - Minimum run requirements for electrolysers (must-run constraints)
           - Line packing in H2 pipelines: Could be modeled as additional storage capacity
       While this would be hard to implement:
           - Physical hydrogen flows with pressure modeling: Full representation of pressure-dependent flow dynamics would require fundamental changes to the transport model architecture and substantially increase computational complexity
   * - 8.2 Integration of Hybrid heat pumps
         - Ensure hybrid systems are correctly sized for applications, considering peak demand scenarios
         - Ensure assumptions incorporate both economics and behavioural considerations
     - The current implementation supports hybrid heating systems as an investment option alongside standalone heat pumps, resistive heaters, and gas boilers. Various heat pump types and the calculation of the corresponding COP are included. 
       Heating capacity sizing is endogenously optimized based on full-year hourly heating demand profiles, ensuring adequate capacity for peak demand periods while minimising total system costs.
       
     - In PyPSA-Eur, all investment decisions are purely economics-driven based on cost optimization. Behavioral considerations, such as consumer preferences are not currently incorporated.
       Behavioral constraints could be integrated if formulated as explicit technical or policy constraints, for example, 
       maximum deployment rates to limit annual heat pump installations to realistic adoption curves (e.g., maximum 5% annual increase in heat pump penetration per region)
   * - 8.3 Grid topology
         - Incorporate more detailed representation of H2 network topology that approximates physical H2 flow to ensure key H2 corridors and infrastructure are well represented
     - Detailed H2 network topology is already implemented. Open-TYNDP uses the TYNDP H2 topology, including Z1 and Z2 setup. H2 flow are represented as a linear transport model. 
     - Line packing could be added easily.
       Physical flows including pressure drops would be harder, and result in a non-linear optimisation.
       PyPSA-Eur supports administrative clustering (NUTS0 to NUTS3), allowing the network to be resolved at highly granular levels.
   * - 8.4 Gas turbine usage and peaking unit utilisation 
         - Explore dynamic operational needs of gas turbines given increasing reliance of variable renewable resources
         - Update assumptions around CH4 to H2 retrofitting projects, which have struggled to compete in current markets
     - Open TYNDP allows for endogenous retrofitting of CH4 plants to operate with H2, or retrofitting of gas boilers to run with H2
     - All relevant features of already implemented
   * - 8.5 EV Modelling
         - Ensure that electricity flows follow charge/discharge cycles of EV batteries
     - Simultaneous charging/discharging is avoided in PyPSA framework used by Open TYNDP by the use of small marginal costs.
     - PyPSA-Eur refines Vehicle-to-Grid (V2G) dispatch capacity and temperature-dependent energy demand correction factors for EVs. These will be incorporated in Open TYNDP during a future update.
   * - 8.6 Economic Assessment
         - incorporate economic assessment of key technologies such as Steam Methane Reformers (SMR), nuclear plants, ammonia regasification terminals which are currently represented without economic attributes
     - Open TYNDP incorporates default capital investment and fixed and variable operation cost assumptions for all technologies based on the open-source `technology-database <https://github.com/PyPSA/technology-data>`_
     - All relevant features of already implemented
   * - 8.7 Methane pricing structure and formation
            - Resolve pricing inconsistences between synthetic natural gas, biomethane and hydrogen
     - High consistency is assured through using open licensed databases which are checked by multiple people. The automated workflow guarantees same currency year and units
     - All relevant features of already implemented
   * - 8.8 Ammonia Import costs
           - Ensure ammonia import costs reflect the entire supply chain
     - Open TYNDP includes location and capacities of European ammonia plants. Ammonia import prices and volumes can be configured.
       Supply chain costs can be modelled directly for all green carriers. See `Neumann et al. (2025) <https://www.nature.com/articles/s41467-025-60652-1>`_
     - All relevant features of already implemented
   * - 8.9 Distinguish between hydrogen use as energy or feedstock 
           - Separating hydrogen used directly as a gas from hydrogen used as a feedstock for producing synthetic fuels
     - H2 use for energy and feedstock purposes is already clearly distinguished in Open TYNDP which includes a suite of technologies for methanol-to-power, reforming, and kerosene, and updated locations/capacities for ammonia plants to accurately distribute demand.
     - All relevant features of already implemented
   * - 8.10 Flexibility of heat pumps
           - Incorporate heat-pump modelling into PLEXOS to better represent flexibility, thermal inertia and heat storage
     - Open TYNDP implements the TYNDP methodology for modelling hybrid heating
     - PyPSA-Eur incorporates heat pump modelling for various types of heat pumps, including the calculation of the COP
           - differentiates between rural/urban, service/residential, district/individual heating
           - for district heating different types of thermal storage are modelling, including constraints on the energy-to-power ration, booster heat pumps to get to the necessary temperature level of the district heating network
           - for individual heating small thermal storage in form from water tanks is included
           - endogenous decision on building renovation can be modelled, with country-specific building stock data and renovation costs
           - thermal inertia of buildings is not modelled, could be presented by a free thermal storage
             Includes an option to calculate dynamic storage capacities for thermal energy storage. Also includes aquifer thermal energy storage, and supplemental heating such as booster heat pumps.
   * - 8.11 Modelling of E-fuels 
           - Allow transportation of e-methanol, e-methane and e-kerosene via pipelines, ships or tanker trucks
           - Incorporate full load hours into spatial optimisation of e-fuel refineries and supply infrastructure
     - <TODO>
     - <TODO>
   * - 8.12 Higher granularity topology
           - Adding more nodes per country and differentiating between prosumer and non-prosumer households for accurate grid interaction modeling
     - Open TYNDP builds on PyPSA Eur which has enables highly spatial and temporally resolved modelling
     - PyPSA-Eur was designed for highly resolved spatial and temporal modelling. It supports multiple resource classes for wind and solar per region to improve accuracy at low spatial resolutions. It also allows for behind-the-meter rooftop PV modeling.
           - higher spatial resolution is possible, also the data can be added on a higher spatial resolution with finer granularity
           - distribution grid with associated costs and losses can be modelled, but data for existing infrastructure of the distribution grid as well as regional investment costs is missing
           - on the distribution level connected are EV charging, rooftop PV, electricity demand for individual heating technologies
           - different types of PV panels can be modelled on the distribution grid level (see `Rahdan et al. (2025) <https://www.nature.com/articles/s41467-025-61492-9>`_, and `Rahdan et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0306261924001041>`_
   * - 8.13 Improved modelling of prosumer demand
     - Open TYNDP allows connection of microgeneration e.g. residential solar PV to be connected to low voltage buses.
       Residential and utility scale PV are treated separately, with separate rules in the workflow to build and cluster rooftop potentials.
       Open TYNDP also distinguishes between stationary battery storage and EV batteries.
     - All relevant features of already implemented
   * - 8.14 Consider peaking units as expansion candidates
     - Open TYNDP is a capacity expansion model by nature and can be set to build new peaking units whenever they are the cost-optimal way to ensure reliability.
     - All relevant features of already implemented
   * - 8.15 Check on remaining CO2 emissions in 2050
     - Open TYNDP models the full CO2 management (not only emissions, but also CCS, storage, CO2 networks). Remaining emissions can be checked in the final csvs and automated plots.
     - A simple extension is to add an extra constraint on top of the CO2 price, e.g. that CO2 emissions have to be net-zero in 2050
   * - 8.16 Implementation of hybrid electrolyser plants
           - Implementing plants connected to both dedicated renewables and the grid to optimize production and market coupling
     - Open-TYNDP implements offshore wind hubs where wind farms can connect to both the network and P2G units for H2 production
     - Existing studies such as `Zeyen et al. (2024) <https://iopscience.iop.org/article/10.1088/1748-9326/ad2239>`_ have investigated this with PyPSA
       showing that implementation is possible for Open TYNDP model.
   * - 8.17 Hydrogen imports and pipeline assessment
           - Modeling practical volumes and prices for pipeline imports to evaluate energy security and dependence.
     - <TODO>
     - PyPSA-Eur implemented renewable energy imports for H2, ammonia, methanol, and oil with configurable prices and volume limits.
   * - 8.18 Geographical correlation in hydrogen production
     - Open-TYNDP ensures geographical correlation by attaching planning-year dependent renewable profiles from the PECD to specific generators within interconnected zones
     - <TODO>

Stakeholder Reference Group (SRG) Proposals
===========================================

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Innovation Roadmap Details
     - Current Open TYNDP Implementation  
     - Features available in PyPSA-Eur
   * - 9.1 Synthetic Fuels
           - Incorporating methanol for the maritime sector to align with decarbonization goals and identify infrastructure needs
     - <TODO>
     - PyPSA-Eur introduced methanol-based technologies (e.g., biomass-to-methanol) in its 2024.09 release. PyPSA-Eur defaults maritime demand to methanol.
       Methanol can be used also in various other sectors (e.g. as back up power, in industry, as kerosene). See `Glaum et al. (2025) <https://arxiv.org/html/2505.09277v1>`_
   * - 9.2 Climatic Variability
           - Suggesting models run with three different climatic years to assess impact on energy security
     - <TODO>
     - PyPSA-Eur is designed for this; it integrates with atlite to process multi-year datasets and supports spanning these in a single model.
  
   * - 9.3 Industrial applications
           - Verify technical and commercial viability of converting industrial gas offtakes to H2 or other carriers before grid expansion
     - Open TYNDP uses a script to interpolate industry sector transition pathways, gradually switching processes from status quo to best-in-class energy consumption per ton of material output.
       One can also choose to model the supply of process heat for industry (split in low, medium, high) endogenously, so cost-optimal solution would be found for potential switch from methane to hydrogen/power/biomass
     - 
   * - 9.4 Sector-specific modelling
           - Discuss the Z1 Z2 concept to streamline management across gas, electricity, and hydrogen vectors
     - Open-TYNDP has already introduced the TYNDP H2 topology, which specifically includes the H2 Z1 and Z2 setup, production, and storage technologies
     -
   * - 9.5 EV modelling techniques
           - refine assumptions on EV charging behavior and their impact on potential grid bottlenecks
     - 
     - PyPSA-Eur now limits Vehicle-to-Grid (V2G) dispatch capacity based on the fraction of vehicles participating in demand-side management. It also refines temperature-dependent correction factors for EV energy demand.
       Since the distribution grid is modelled (even though without the corresponding topology), just as an capacity expansion with corresponding costs, one can investigate the relation between flexible EV charging and necessary distribution grid capacity
   * - 9.6 District heating
           - Create a dedicated tool distinct from the ETM to simulate production from biomass, geothermal, and other sources
     - 
     - PyPSA-Eur already features a highly detailed district heating module. Recent additions include geothermal district heating, aquifer thermal energy storage (ATES), and booster heat pumps for supplemental heating

   * - 9.7 Liquified hydrogen
           - Explore LH2 import methods to understand logistical, storage, and cost constraints
     -
     - PyPSA-Eur implements a "H2 liquid" bus at each location to specifically handle hydrogen liquefaction costs for shipping demand
   * - 9.8 H2 Import Quotas
           - Align hydrogen import quotas with RepowerEU targets to avoid overestimating domestic production
     - 
     -
   * - 9.9 Electric heat pumps
           - Move heat pump modeling to PLEXOS to better capture thermal inertia and load management
     - 
     - PyPSA-Eur manages this via an energy-to-power ratio constraint for thermal storage, linking storage capacity to charger capacity to simulate load shifting
   * - 9.10 Optimisation Across Energy Vectors
           - Expanding PLEXOS to integrate electricity, gas, and hydrogen systems holistically
     - Open-TYNDP optimizes these vectors simultaneously by default
     - PyPSA-Eur optimizes these vectors simultaneously by default
   * - 9.11 Transmission System Losses
           - Reassessing losses to reflect actual power flow dynamics more accurately
     -
     - PyPSA-Eur allows for piecewise linear approximation of transmission losses and provides the option to disable efficiency losses for specific carriers.
   * - 9.12 Additional hydrogen production pathways
           - Integrating methane pyrolysis and waste-to-hydrogen processes
     - 
     - PyPSA-Eur has already integrated biomass-to-hydrogen (with or without carbon capture) and supports custom technology adjustments via configuration.
   * - 9.13 Flexibility in modelling
           - Focus on load displacement (shifting demand) rather than just load reduction to maximize renewable use
     - Open TYNDP follows TYNDP methodology so exhibits the same limitation. 
     - PyPSA-Eur has higher flexibility since all sectors are modelled (instead of fixed exogenous demand for H2).
   * - 9.14 Price setting for hydrogen
           - Revise methodology to reflect real-world contracts and costs like dehydrogenation
     - As an integrated sector-coupled model, endogenous pricing of hydrogen includes all represented upstream processes
     - ditto
   * - 9.15 Sensitivity to commodity prices
           - Conduct sensitivity analyses on price fluctuations (gas, oil, H2) to understand investment risks
     - Workflow management tool snakemake enables the simultaneous execution of multiple scenarios with single calls and configuration overrides
     - ditto
   * - 9.16 Inclusion of emerging technologies
     - The collaborative approach offered by Open TYNDP provides a formal review process by which new technologies can be included in the analysis
     - 
   * - 9.17 Out of scope Innovations
           - Carbon Capture and utilisation
           - Innovative Grid Technologies
     - Open TYNDP includes a detailed representation of CCS and CDR technologies including carbon sequestration sites
     - 
