# PlasticParcels
`PlasticParcels` is a python package for simulating the transport and dispersion of plastics in the ocean. The tool is based on `v3.0.2` of the [`Parcels`](https://oceanparcels.org/) computational Lagrangian ocean analysis framework ([@Lange2017](http://dx.doi.org/10.5194/gmd-10-4175-2017) and [@Delandmeter2019](http://dx.doi.org/10.5194/gmd-12-3571-2019)), providing a modular and customisable collection of methods, notebooks, and tutorials for advecting virtual plastic particles with a wide range of physical properties.

# Table of contents
0. [Description of Software](#description)
1. [Installation](#installation)
2. [Physics Kernels](#physicskernels)
    1.  [Stokes Drift](#stokes)
    2.  [Wind-induced Drift / Leeway](#winddrift)
    3.  [Biofouling](#biofouling)
    4.  [Vertical Mixing](#verticalmixing)
    5.  [Sea-ice Capture](#seaice)
3. [Particle Initialisation Maps](#initialisationmaps)
    1. [Coastal mismanaged plastic waste emissions](#coastalrelease)
    2. [Riverine mismanaged plastic waste emissions](#riverrelease)
    3. [Open-sea fishing-related plastic emissions](#fishingrelease)
    4. [Current global ocean plastic concentrations](#staterelease)


## Description of software
An open-source article describing `PlasticParcels` can be found here **link to article**. The tool applies a collection of physical processes to the virtual particles, such as Stokes drift, wind-induced drift, biofouling, and turbulent mixing, via custom particle behaviour programmed in the form of `Kernels`. In addition to the fine-scale physics parameterisations, `PlasticParcels` provides global particle initialisation maps that represent best estimates for plastic pollution emissions along coastlines [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), from river sources [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803), in the open-ocean from fishing-related activities [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646), as well as a current best estimate of buoyant plastic concentrations globally [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0). We envisage PlasticParcels as a tool for easy-to-run plastic dispersal simulations; as well as for rapid prototyping, development, and testing of new fine-scale physics parameterisations.

The current version supports nano- and microplastic behaviour, with support for macroplastics planned in the near-future. It has been designed for use with the [Copernicus Marine Service platform](https://marine.copernicus.eu/), providing new plastic modelling capabilities as part of the [NECCTON](https://neccton.eu) project. `PlasticParcels` is easily adapted to run on local machines and high-performance computing (HPC) architecture with various hydrodynamic, biogeochemical, and other model fields as inputs. A future goal is to embed `PlasticParcels` within a cloud platform to allow for even more rapid prototyping, development, and simulations.

Below we provide instructions to install `PlasticParcels`, we detail the specific physics kernels implemented in the tool, and we describe how the particle initialisation maps are generated.

## Installation
**(Is this necessary? Or can be relegated to the github readme?)**

**EvS: Yes, I think this is necessary. And I think we should make a conda package, so that the installation is simply `conda install plasticparcels`.**

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
