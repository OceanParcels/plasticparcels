---
title:  'plasticparcels: A python package for marine plastic dispersal simulations using Parcels'
tags:
  - Plastic dispersal simulations
  - Lagrangian oceanography
  - Fluid dynamics
authors:
  - name: Michael C. Denes
    orcid: 0000-0003-1868-0043
    affiliation: 1
  - name: Erik van Sebille
    orcid: 0000-0003-2041-0704
    affiliation: 1

affiliations:
  - name: Institute for Marine and Atmospheric Research, Utrecht University, the Netherlands
    index: 1
date: 24 April 2024
bibliography: paper.bib
link-bibliography: false
---

# Summary
`plasticparcels` is a python package for simulating the transport and dispersion of plastics in the ocean. The tool is based on `v3.0.2` of the `Parcels` computational Lagrangian ocean analysis framework [@Lange2017,@Delandmeter2019], providing a modular and customizable collection of methods, notebooks, and tutorials for advecting virtual plastic particles with a wide range of physical properties. The tool applies a collection of physical processes to the virtual particles, such as Stokes drift, wind-induced drift, biofouling, and turbulent mixing, via custom particle behaviour programmed in the form of `Kernels`. In addition to the fine-scale physics parameterisations, `plasticparcels` provides global particle initialisation maps that represent best estimates for plastic pollution emissions along coastlines [@Jambeck2015], from river sources [@Meijer2021], in the open-ocean from fishing-related activities [@Kroodsma2018], as well as a current best estimate of buoyant plastic concentrations globally [@Kaandorp2023]. We envisage plasticparcels as a tool for easy-to-run plastic dispersal simulations; as well as for rapid prototyping, development, and testing of new fine-scale physics parameterisations.

The current version supports nano- and microplastic behaviour, with support for macroplastics planned in the near-future. It has been designed for use with the Copernicus Marine Service platform [@CMEMS] **update .bib with full-details** <mark>EvS: I don't think they use the abbreviation CMEMS anymore</mark>, providing new plastic modelling capabilities as part of the NECCTON project. `plasticparcels` is easily adapted to run on local machines and high-performance computing (HPC) architecture with various hydrodynamic, biogeochemical, and other model fields as inputs. A future goal is to embed `plasticparcels` within a cloud platform to allow for even more rapid prototyping, development, and simulations.


# Statement of need
Marine plastic debris can be found almost everywhere in the ocean. A recent study estimates that there is approximately 3,200 kilotonnes of (initially) positively buoyant plastics in the global ocean in the year 2020 [@Kaandorp2023], where 59-62\% of these plastics are found at the ocean surface, 36-39\% within the deeper ocean, and 1.5-1.9\% along the coastline. They estimate that 500 kilotonnes of positively buoyant plastic enters the ocean each year, where 39-42\% originate from mismanaged waste along coastlines, 45-48\% originate from fishing-related activities (e.g. fishing lines, nets, traps, and crates), and 12-13\% from mismanaged waste entering the ocean via rivers. These estimates are for positively buoyant plastics, which make up only a fraction of the total production of virgin plastics each year. **(Citation for percentage breakdown of pos. buoyant vs. neg. buoyant?)**

Due to its durable, inert, and cheap-to-manufacture nature, plastic has become one of the most abundant manufactured synthetic materials on Earth. Between 1950 and 2017 an estimated 8,300 million tonnes [@Geyer2017] of virgin plastic was produced, with the rate of production only set to increase. Its durability is of primary concern to the marine environment, where, without intervention, plastics will remain for millennia to come **(find citation)**, and will likely degrade and fragment into smaller pieces that will disperse across ever larger distances. These plastics interact and interfere with marine wildlife, either entangling, or being inadvertently ingested, with documented cases affecting over 900 marine species so far [@Kuhn2020]. To better understand and predict the effects of plastic pollution on the marine environment, it is of paramount importance that we understand where and how plastic enters our ocean, and the pathways of transport, dispersal patterns, and ultimate fate of these plastics.

Lagrangian ocean analysis, where virtual particles are tracked in hydrodynamic flow fields, is widely used to uncover and investigate the pathways and timescales of the dispersion of plastic particulates in the ocean [@Lebreton2012,@Hardesty2017,@JalonRojas2019,@Chassignet2021,@Kaandorp2023]. However, two important questions arise when performing such Lagrangian simulations. Firstly, what physical processes drive the transport and dispersal of a plastic particle? The properties of plastic particles (e.g., size, shape, and density) determine what the dominant physical processes are at play, and due to the chaotic nature of the ocean, plastics of different properties will have unique dispersal patterns and transport behaviours. Current state-of-the-art ocean models are either too coarse in resolution to capture these processes, or disregard these processes entirely, and so parameterising these processes is important in modelling their effects. Secondly, what are the initial release locations and concentrations of marine plastic pollution? Forecasting spatial maps of near-future plastic concentrations is largely an initial value problem, relying on accurate initial conditions for a realistic simulation output.

The past decade has seen a growing number of community-developed software packages for performing Lagrangian simulations [@Paris2013,@Fredj2016,@Lange2017,@Doos2017,@Dagestad2018,@JalonRojas2019,@Delandmeter2019]. In many cases, these packages are specific to particular particle classes or hydrodynamic models, or are written and embedded in proprietary software languages, and can be inflexible or difficult to integrate into different applications. In the case of plastic dispersal simulations, where the physical processes are still researched and their implementation is development [@vanSebille2020], an open-source, flexible, and modular approach to performing Lagrangian simulations is necessary for rapid plastic dispersal simulations, as well as for prototyping, developing, and testing new physical process parameterisations. **improve/expand on this last sentence**

Here, we have developed `plasticparcels` to unify plastic dispersion modelling into one easy-to-use code. While `plasticparcels` has been designed for researchers who routinely perform plastic particle dispersion simulations, it is equally useful to novice users who are new to Lagrangian ocean analysis techniques.

# Description of the software
`plasticparcels` has been designed as a layer over the `Parcels` Lagrangian framework [@Lange2017,@Delandmeter2019]. The core functionality of `Parcels` are its `FieldSets`, `Kernels` and `ParticleSets`. These objects are designed to be as flexible and customisable as possible so that users can perform Lagrangian simulations of a wide variety of particulates, such as tuna, plastic, plankton, icebergs, turtles **(etc. etc. + CITATIONS)**. However, due to the flexible nature of the software, there is a steep learning curve for new users, who often find it difficult to setup their simulations in a rapid fashion due to the complexity of modern hydrodynamic model output. We have developed `plasticparcels` as user-friendly tool specifically designed for easy-to-generate plastic dispersal simulations. While `plasticparcels` is primarily designed for use in the cloud and in HPC environments (due to the ever increasing size of hydrodynamic datasets generated from ocean general circulation models), it can be easily installed and run on local machines.

![`plasticparcels` schematic.\label{fig:schematic}](schematic.png){width=80%}

The core features of plasticparcels are: 1) a user-friendly python notebook layer on top of `Parcels` that provides a streamlined workflow for performing plastic dispersal simulations, 2) custom `Parcels` kernels designed to simulate the fine-scale physical processes that influence the transport of nano- and microplastic particulates, and 3) global particle initialisation maps which represent the best estimate locations of plastic pollution emissions from coastal sources, river sources, open ocean fishing-related activity emission sources, and a current best estimate of buoyant plastic concentrations.

In addition, due to the flexibility of the package, users may use the functions and modular design of `plasticparcels` to enhance their existing `Parcels` simulations and workflow. For example, users can use the initialisation maps, associated `ParticleSet` creation methods, and/or the custom physics kernels with their own `Parcels` simulations. Post-processing and analysis of the generated trajectory datasets is purposefully left to the user, however some tutorials are provided in the [`plasticparcels` github repository](\url{https://github.com/OceanParcels/plasticparcels}), along with the tutorials on the [`Parcels` github repository](\url{https://github.com/OceanParcels/parcels}). Below we provide an example of how `plasticparcels` may be used, utilising a developed `ipynb` tutorial notebook.


![Particle initialisations maps based on a) current global surface concentrations, b) coastal plastic emissions, c) river plastic emissions, d) fishing activity plastic emissions.\label{fig:initialisation_maps}](initialisation_maps.png){width=80%}

# Usage Examples
**TODO fill this in once created**
The `plasticparcels` github repository provides several useful tutorials. Here, we briefly demonstrate how `plasticparcels` can be used for a microplastic dispersal simulation in the Mediterranean Sea. The tutorial can be found on the `plasticparcels` github repository, ...

[Describe a release location, advection time, and show 2 plots - a) dispersal pathways b) heatmap/concentrations]

# Acknowledgments
We would like to thank the OceanParcels team at Utrecht University for their helpful suggestions developing the tool. MCD is supported by the NECCTON project, which has received funding from Horizon Europe RIA under grant agreement No 101081273. EvS was supported by the project Tracing Marine Macroplastics by Unraveling the Ocean’s Multiscale Transport Processes with project number VI.C.222.025 of the research programme Talent Programme Vici 2022 which is financed by the Dutch Research Council (NWO).

# References