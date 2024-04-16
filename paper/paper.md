---
title:  'PlasticParcels: A python package for marine plastic dispersal simulations using Parcels'
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
date: 15 April 2024
bibliography: paper.bib
link-bibliography: false
---

# Summary
`PlasticParcels` is a python package for simulating the transport and dispersion of plastics in the ocean. The tool is based on `v3.0.2` of the `Parcels` computational Lagrangian ocean analysis framework [@Lange2017, @Delandmeter2019], providing a modular and customizable collection of methods, notebooks, and tutorials for advecting virtual plastic particles with a wide range of physical properties. The tool applies a collection of physical processes to the virtual particles, such as Stokes drift, wind-induced drift, biofouling, and turbulent mixing, via custom particle behaviour programmed in the form of `Kernels`. In addition to the fine-scale physics parameterisations, `PlasticParcels` provides global particle initialisation maps that represent the best estimates for plastic pollution emissions along coastlines [@Jambeck2015], from river sources [@Meijer2021], in the open-ocean from fishing-related activities [@Kroodsma2018], as well as a current best estimate of buoyant plastic concentrations globally [@Kaandorp2023]. We envisage PlasticParcels as a tool for rapid plastic dispersal simulations; as well as for rapid prototyping, development, and testing of new fine-scale physics parameterisations.

The current version supports nano- and microplastic behavior, with support for macroplastics planned in the near-future. It has been designed for use with the Copernicus Marine Service platform [@CMEMS], providing new plastic modelling capabilities as part of the NECCTON project. `PlasticParcels` is easily adapted to run on local machines and high-performance computing (HPC) architecture with various hydrodynamic, biogeochemical, and other model inputs. A future goal is to embed `PlasticParcels` within a cloud platform to allow for even more rapid prototyping, development, and simulations.


# Statement of need
**General statement on plastic in the ocean**
Marine plastic debris can be found almost anywhere in the ocean.
**Quantify and give location to the statement**
A recent study estimates that there is approximately 3,200 kilotonnes of (initially) positively buoyant plastics in the global ocean in the year 2020 [@Kaandorp2023], where 59-62\% of these plastics are found at the ocean surface, 36-39\% within the deeper ocean, and 1.5-1.9\% along the coastline. They estimate that 500 kilotonnes of positively buoyant plastic enters the ocean each year, where 39-42\% originate from mismanaged waste along coastlines, 45-48\% originate from fishing-related activities (e.g. fishing lines, nets, traps, and crates), and 12-13\% from mismanaged waste entering the ocean via rivers. These estimates are for positively buoyant plastics, which make up only a fraction of the total production of virgin plastics each year (Citation for percentage breakdown?).

**Description of the impacts to marine ecology**
Due to it's durable, inert, and cheap-to-manufacture nature, plastic has become one of the most abundant man-made materials on Earth. Between 1950 and 2017 an estimated 8,300 million tonnes [@Geyer2017] of virgin plastic was produced, with the rate of production only set to increase. It's durability is of primary concern to the marine environment, where, without intervention, plastics will remain for millenia to come, and will likely degrade and fragment into smaller pieces that will disperse across ever larger distances. These plastics interact and interfere with marine wildlife, either entangling, or being accidently ingested, with over 900 species documented so far [Kuhn2020].
**Statement on why understanding dispersal patterns, pathways, transport, fate, is important**
To better understand and predict the effects of plastic pollution on the marine environment, it is of paramount importance that we understand where and how plastic enters our oceans, and the pathways of transport, dispersal patterns, and ultimate fate of these plastics.



**Paragraph on approach to solve said problem**
Lagrangian ocean analysis, where virtual particles are tracked in hydrodynamic flow fields, is widely used to uncover and investigate the pathways and timescales of the dispersion of plastic particulates in the ocean [@Lebreton2012, @Hardesty2017, @JalonRojas2019, @Chassignet2021, @Kaandorp2023]. However, two important questions arise when performing such Lagrangian simulations. Firstly, what physical processes drive the transport and dispersal of a plastic particle? The properties of plastic particles (e.g., size, shape, and density) determine what the dominant physical processes are at play, and due to the chaotic nature of the ocean, plastics of different properties will have unique dispersal patterns and transport behaviours. Current state-of-the-art ocean models are either too coarse in resolution to capture these processes, or disregard these processes entirely, and so parameterising these processes is important in modelling their effects. Secondly, what are the initial release locations and concentrations of marine plastic pollution? Predicting spatial maps of future plastic concentrations is, in effect, an initial value problem, relying on accurate initial conditions for a realistic simulation output.

The past decade has seen a growing number of community-developed software packages for performing Lagrangian simulations [@Paris2013, @Fredj2016, @Lange2017, @Doos2017, @Dagestad2018, @JalonRojas2019, @Delandmeter2019]. In many cases, these packages are specific to particular particle classes, or hydrodynamic models, and can be inflexible or difficult to integrate into different applications. In the case of plastic dispersal simulations, where the physical processes are still under research and development [@vanSebille2020], a flexible and modular approach to performing Lagrangian simulations is necessary for rapid plastic disperal simulations, as well as for prototyping, developing, and testing new physical process parameterisations.

**Who is PP designed for (maybe for summary?)**
Here, we have developed `PlasticParcels` to unify plastic dispersion modelling into one easy-to-use (user friendly?) code. While `PlasticParcels` has been designed for researchers who perform plastic particle dispersion simulations, it is equally useful to ... (finish off paragraph).



# Description of the software
`PlasticParcels` has been designed as a layer over the `Parcels` Lagrangian framework [@Lange2017, @Delandmeter2019].
The core functionality of `Parcels` are it's `FieldSets`, `Kernels` and `ParticleSets`.
These objects are designed to be as flexible and customisable as possible so that users can perform Lagrangian simulations of a wide variety of particulates, such as tuna, plastic, plankton, etc. etc. + CITATIONS).
However, due to the flexible nature of the software, there is a steep learning curve for new users, who often find it difficult to setup their simulations in a rapid fashion. We have developed `PlasticParcels` as user-friendly tool specifically designed for easy to generate plastic dispersal simulations. While `PlasticParcels` is primarily designed for use in the cloud and in HPC environments (due to the ever increasing size of hydrodynamic datasets generated from ocean general circulation models), it can be easily installed and run on local machines.

The core features of PlasticParcels are: 1) a user-friendly python notebook layer on top of `Parcels` that provides a streamlined workflow for performing plastic dispersal simulations, 2) custom `Parcels` kernels designed to simulate the fine-scale physical forces that influence the transport of nano- and microplastic particulates, and 3) particle inisialisation maps which represent the best-guess locations of emissions from coastal sources, river sources, open ocean fishing activity related emission sources, and a current best estimate of buoyant plastic concentrations.

In addition, due to the flexibility of the package, users may use the the functions and modular design of `PlasticParcels` to enhance their existing `Parcels` simulations and workflow.
For example, users may use the particle initialisation maps and associated `ParticleSet` creation methods, and/or the custom physics kernels with their own `Parcels` simulations.
Post-processing and analysis of these trajectory datasets is purposefully left to the user, however some tutorials are provided in the `PlasticParcels` github repository, along with the tutorials on the `Parcels` github repository. Below we detail the specific kernels implemented, as well as describe how the particle initialisation maps are generated. We also provide an example of how `PlasticParcels` may be used, modelled from an `ipynb` tutorial notebook.

![`PlasticParcels` schematic](schematic.png "PlasticParcels` schematic")

## Physics kernels
The `Parcels` Lagrangian framework is a tool for advecting virtual particles. It works by integrating the velocity fields from a hydrodynamic model while including any additional \textit{behaviour} of the particle. Mathematically, particle trajectories are computed by solving the following equation:
$$
\begin{equation}
\mathbf{x}(t) = \mathbf{x}(0) + \int_{0}^{t} \mathbf{v}(\mathbf{x}(s), s) + \mathbf{B},
\end{equation}
$$
where $\mathbf{x}(t)$ describes the particle position at time $t$, $\mathbf{v}$ is the hydrodynamic model velocity field, and $\mathbf{B}$ describes any changes to the particle position due to any additional behaviour the particle exhibits. Examples of these additional behaviours are described below.

Numerically, we solve Eq (1) using a time-stepping approach, effectively solving

$$
\begin{equation}
\frac{\text{d}\mathbf{x}(t)}{\text{d}t} = \mathbf{v}(\mathbf{x}(t), t) + \frac{\text{d}\mathbf{B}}{\text{d}t}.
\end{equation}
$$
For simplicity, we use the default fourth-order Runge-Kutta scheme of parcels to solve the advection of the particle from the hydrodynamic model velocity field $\mathbf{v}$, and a Euler-forward scheme for all other additional behaviours $\mathbf{B}$. 


### Stokes Drift
We include a kernel to parameterise the effect of Stokes drift on a particle, based on the Phillips spectrum approximation developed in [@Breivik2016]. Specifically, we model the additional behaviour of a particle as $\mathbf{B}_{\text{Stokes}}$, where the change in the particle position is described by
$$
\begin{equation}
\frac{\text{d}\mathbf{B}_{\text{Stokes}}}{\text{d}t} = \mathbf{v}_{\text{Stokes}}(\mathbf{x}(t), t) =\mathbf{v}_{\text{Stokes}}(\mathbf{x}_{z=0}(t),t)\bigg(e^{2k_p z} - \beta\sqrt{-2\pi k_p z}\text{ erfc}(-2k_p z) \bigg).
\end{equation}
$$
Here, $z$ is the vertical component of the particle position, $\mathbf{v}_{\text{Stokes}}(\mathbf{x}_{z=0}(t),t)$ is the surface Stokes drift velocity, $\beta=1$ (as we assume a Phillips spectrum), and erfc is the complementary error function. $k_p$ is the peak wave number, computed as $k_p = \omega_{p}^2/9.81$, with $\omega_p = 2 \pi / T_p$, where $\omega_p$ is the peak wave freaquency computed from the peak wave period $T_p$. 

! Include where this kernel has been used before
[@Onink2021] uses Stokes drift at the surface only, and uses SummedFields (so RK4 approach), and not an explicit kernel.

### Wind-induced drift
$$\mathbf{B}_{\text{Wind}} = $$

! Include where this kernel has been used before

### Biofouling kernel
@Lobelle2021
@Fischer2022
@Kaandorp2023
based on @Kooi2017
$$\mathbf{B}_{\text{Biofouling}} = $$

This kernel has been used in multiple forms in ...

### Vertical mixing kernel
@Onink2022
$$\mathbf{B}_{\text{Vertical Mixing}} = $$



### Sea-ice capture

$$\mathbf{B}_{\text{Sea-ice Capture}} = $$

### In development:
-beaching, fragmentation, degradation, etc.

When performing a plastic dispersal simulation with `PlasticParcels`, users have the explicit option of choosing which additional behaviour to include. For example, when solving 



## Description of algorithms for particle release location maps
Datasets:
Coastal
@Jambeck2015
Description of creating coastal [@Jambeck2015] + image

Rivers
@Meijer2021
Description of creating rivers [@Meijer2021] + image

Fisheries
@Kroodsma2018
Description of creating fisheries [@Kroodsma2018] + image

Current concentrations
@Kaandorp2023
Description of creating current map + image

## Installation

1) Download software
2) Download <mark> Jambeck, Kroodsma, Meijer, Kaandorp datasets </mark>
3) Update `settings.txt` <mark> give more description </mark>
4) Run `create_masks.py` and `create_release_locations.py` <mark> wrap this into a single script? </mark>
5) Open and run `PlasticParcels_tutorial.ipynb`
6) Open, modify, and run `PlasticParcels_template.ipynb`


OR

The latest version of `PlasticParcels` can be installed directly from github, via:
```
git clone https://github.com/OceanParcels/PlasticParcels.git
cd PlasticParcels; pip install -r requirements.txt
python PlasticParcels/run_initialisation.py
export PYTHONPATH="$PYTHONPATH:$PWD"
```
Ensure that you have updated `settings.txt` with the required directories and filenames.

# Usage Examples
The `PlasticParcels` github repository provides several useful tutorials. Here, we briefly demonstrate how `PlasticParcels` can be used for a microplastic dispersal simulation in the Mediterranean Sea.

Dispersal of plastic in the Mediterranean? <mark> Include as an example, `PlasticParcels_Mediterranean_example.ipynb`? </mark>

2 plots - a) dispersal pathways b) heatmap/concentrations

# Acknowledgments

We would like to thank ...

# References