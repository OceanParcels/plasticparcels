# PlasticParcels
`PlasticParcels` is a python package for simulating the transport and dispersion of plastics in the ocean. The tool is based on `v3.0.2` of the `Parcels` computational Lagrangian ocean analysis framework [@Lange2017,@Delandmeter2019], providing a modular and customizable collection of methods, notebooks, and tutorials for advecting virtual plastic particles with a wide range of physical properties.

# Table of contents
0. [Description of Software](#description)
1. [Installation](#installation)
2. [Physics Kernels](#physicskernels)
    a. [Stokes Drift](#stokes)
    b. [Wind-induced Drift / Leeway](#winddrift)
    c. [Biofouling](#biofouling)
    d. [Vertical Mixing](#verticalmixing)
    e. [Sea-ice Capture](#seaice)
3. [Particle Initialisation Maps](#initialisationmaps)
    a. [Coastal mismanaged plastic waste emissions](#coastalrelease)
    b. [Riverine mismanaged plastic waste emissions](#riverrelease)
    c. [Open-sea fishing-related plastic emissions](#fishingrelease)
    d. [Current global ocean plastic concentrations](#staterelease)


## Description of software
An open-source article describing `PlasticParcels` can be found here **link to article**. The tool applies a collection of physical processes to the virtual particles, such as Stokes drift, wind-induced drift, biofouling, and turbulent mixing, via custom particle behaviour programmed in the form of `Kernels`. In addition to the fine-scale physics parameterisations, `PlasticParcels` provides global particle initialisation maps that represent best estimates for plastic pollution emissions along coastlines [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), from river sources [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803), in the open-ocean from fishing-related activities [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646), as well as a current best estimate of buoyant plastic concentrations globally [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0). We envisage PlasticParcels as a tool for easy-to-run plastic dispersal simulations; as well as for rapid prototyping, development, and testing of new fine-scale physics parameterisations.

The current version supports nano- and microplastic behaviour, with support for macroplastics planned in the near-future. It has been designed for use with the Copernicus Marine Service platform [@CMEMS](https://marine.copernicus.eu/), providing new plastic modelling capabilities as part of the NECCTON project. `PlasticParcels` is easily adapted to run on local machines and high-performance computing (HPC) architecture with various hydrodynamic, biogeochemical, and other model fields as inputs. A future goal is to embed `PlasticParcels` within a cloud platform to allow for even more rapid prototyping, development, and simulations.

Below we detail the specific physics kernels implemented, as well as describe how the particle initialisation maps are generated.

## Installation
**(Is this necessary? Or can be relegated to the github readme?)**

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



## Physics kernels <a name="physicskernels"></a>

The `Parcels` Lagrangian framework is a tool for advecting virtual particles that are assumed to be spherical in shape. It works by numerically integrating the velocity fields from a hydrodynamic model while including any additional \textit{behaviour} of the particle. Mathematically, particle trajectories are computed by solving the following equation:

```math
\mathbf{x}(t) = \mathbf{x}(0) + \int_{0}^{t} \mathbf{v}(\mathbf{x}(s), s) + \mathbf{B}(\mathbf{x}(s),s) \text{d}s,
```

where $`\mathbf{x}(t)`$ describes the particle position at time $t$, $`\mathbf{v} = (u,v,w)`$ is the hydrodynamic model velocity field, and $`\mathbf{B}(\mathbf{x}(t),t)`$ describes any displacements to the particle position caused by additional behaviour the particle exhibits or experiences. When performing a plastic dispersal simulation with `PlasticParcels`, users have the explicit option of choosing which additional behaviour to include. Examples of these additional behaviours are described below.

Numerically, we solve the above equation using a time-stepping approach, where we compute the displacements in the particle position as

```math
\frac{\text{d}\mathbf{x}(t)}{\text{d}t} = \mathbf{v}(\mathbf{x}(t), t) + \mathbf{B}(\mathbf{x}(t), t),
```

and updating the particle position at each timestep. For simplicity, by default we use the fourth-order Runge-Kutta scheme of `Parcels` to solve the advection of the particle from the hydrodynamic model velocity field $`\mathbf{v}`$, and an Euler-forward scheme for all other additional behaviours realised in $`\mathbf{B}`$.


### Stokes Drift <a name="stokes"></a>

An important process that affects plastic particle dispersal in the upper ocean is the Stokes drift, whereby a particle subjected to a surface wave will experience a net displacement in the direction of wave propagation. We include a kernel to parameterise the effect of Stokes drift on a particle, based on the Phillips spectrum approximation developed in [@Breivik2016](http://dx.doi.org/10.1016/j.ocemod.2016.01.005). Specifically, we model this additional behaviour as $`\mathbf{B}_{\text{Stokes}}`$, where the change in the particle position is described by

```math
\mathbf{B}_{\text{Stokes}} := \mathbf{v}_{\text{Stokes}}(\mathbf{x}(t), t) =\mathbf{v}_{\text{Stokes}}(\mathbf{x}_{z=0}(t),t)\bigg(e^{2k_p z} - \beta \sqrt{-2\pi k_p z}\text{ erfc}(-2k_p z) \bigg).
```

Here, $z$ is the depth of the particle, $`\mathbf{v}_{\text{Stokes}}(\mathbf{x}_{z=0}(t),t)`$ is the surface Stokes drift velocity, $`\beta=1`$ (as we assume a Phillips spectrum), and erfc is the complementary error function. The peak wave number $`k_p`$ is computed as $`k_p = \omega_{p}^2/9.81`$, where $`\omega_p`$ is the peak wave frequency $`\omega_p = 2 \pi / T_p`$, using the peak wave period $`T_p = T_p(\mathbf{x}_{z=0}(t),t)`$.

Our particular implementation of the Stokes drift kernel requires a surface Stokes velocity field $`\mathbf{v}_{\text{Stokes}}(\mathbf{x}_{z=0}(t),t)`$, as well as a peak wave period field $`T_p(\mathbf{x}_{z=0}(t),t)`$. Earlier versions of this kernel have been used in the following research articles **(INCLUDE CITATIONS - but which?)**.

**-----[@Onink2021] uses Stokes drift at the surface only, and uses a `SummedFields` approach (= RK4 approach also), and not an explicit kernel.**

### Wind-induced drift / Leeway <a name="winddrift"></a>
Plastic particles at the ocean surface that are not completely submersed will experience a force from the relative wind due to a wind drag, leading to a wind-induced drift. This wind-induced drift of the particle is called leeway [@Allen1999], which can be decomposed into a downwind component (in the direction of the wind), and a crosswind component (which is typically non-zero for asymmetric objects). As we assume that each plastic particle is spherical, we can ignore the crosswind component of leeway, and only consider the downwind component of leeway. The downwind component follows an almost linear relationship with the relative 10m wind speed [@Allen2005], so we model the leeway as


```math
\mathbf{B}_{\text{Wind}} := c \cdot \big(\mathbf{v}_{\text{Wind}}(\mathbf{x}(t),t) - \mathbf{v}(\mathbf{x}(t),t)\big),
```


where $`\mathbf{v}_{\text{Wind}}`$ is the wind velocity 10m above sea level, and $`c`$ is the leeway rate (a percentage of wind speed, which we refer to in the code as the windage coefficient). Ignoring all additional behaviour of the particle, then $`\mathbf{v}_{\text{Wind}} - \mathbf{v}`$ is the relative wind acting on the particle.

**! Include where this kernel has been used before + references to literature on what percentages to use**

### Biofouling kernel <a name="biofouling"></a>
Plastic particles in the ocean can be a hotbed for the accumulation and growth of organisms, known as biofouling. The formation of a biofilm on the surface of a plastic particle can result in a density change, affecting the buoyancy of the particle. An initially buoyant particle may become negatively buoyant, and sink or settle, depending on the surrounding seawater density.

We model the biofouling of a plastic particle following the approach of [@Kooi2017](http://dx.doi.org/10.1021/acs.est.6b04702), where the settling velocity of a particle is computed from the relative density difference of the plastic particle and the surrounding seawater. Here, we assume the biofilm growth (and decay) is primarily microbial algae, and is distributed homogeneously over the particle surface. The density of the biofouled plastic particle depends on the radius and density of the particle, and the thickness and density of the algal biofilm. The primary component of the biofouling kernel is modelling the change in the number of attached algae (denoted by $A$) on the surface of the plastic particle. As in [@Kooi2017](http://dx.doi.org/10.1021/acs.est.6b04702), we model the attached algal growth as


```math
\frac{\text{d}A}{\text{d}t} := \underbrace{\frac{A_A \beta_A}{\theta_\text{Plastic}}}_{\text{Collisions}} + \overbrace{\mu_A A}^{\text{Algal growth}} - \underbrace{m_A A}_{\text{Mortality}} - \overbrace{Q_{10}^{(T-20)/10}R_{20}A}^{\text{Respiration}}.
```


The first term models growth of algae due to collisions of the particle with algae in the surrounding seawater, where $`A_A`$ is the ambient algal amount, $`\beta_A`$ is the encounter kernel rate, $`\theta_{\text{Plastic}}`$ is the surface area of the plastic particle. The second term models the growth of the biofilm, where the growth term $`\mu_A`$ is computed from the total productivity provided by model output. The third and fourth terms model the (grazing) mortality and respiration of the biofilm respectively. As in [@Kooi2017](http://dx.doi.org/10.1021/acs.est.6b04702), we use constant mortality $`\m_A`$ and respiration $`R_{20}`$ rates, with a temperature dependent term $`\big(Q_{10}^{(T-20)/10}\big)`$ included in the respiration component (see [@Kooi2017](http://dx.doi.org/10.1021/acs.est.6b04702) for more details).

As described above, the modelled attached algal growth drives a change in the settling velocity of the biofouled particle, $`\mathbf{v}_{\text{Biofouling}}`$. Hence, we model the additional behaviour of the particle due to biofouling as

```math
\mathbf{B}_{\text{Biofouling}} := \mathbf{v}_{\text{Biofouling}}.
```

**(Include where this has been used)** This kernel has been used in various forms in [@Lobelle2021](http://dx.doi.org/10.1029/2020JC017098)[@Fischer2022](http://dx.doi.org/10.5194/bg-19-2211-2022)[@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0).

### Vertical mixing kernel <a name="verticalmixing"></a>
An important process that is unresolved in even high-resolution ocean models is wind-driven turbulent mixing, which occurs at scales far smaller than a typical model ocean grid cell. In the vertical direction, this turbulent mixing can distribute even positively buoyant plastic particles throughout the mixed layer. To model this process, we take the approach of [@Onink2022](http://dx.doi.org/10.5194/gmd-15-1995-2022), by employing a Markov-0 styled stochastic parameterisation. 

Denote by $`K_z = K_z(\mathbf{x}(t))`$ the vertical diffusion coefficient profile based on a $`K`$-profile parameterisation (KPP) model [@Large1994](http://dx.doi.org/10.1029/94RG01872). Then the displacement of a particle (in the vertical direction) with a settling velocity $`w`$ can be modelled as an SDE [@Grawe2012](http://dx.doi.org/10.1007/s10236-012-0523-y),

```math
\mathbf{B}_{\text{Vertical Mixing}} := \bigg(w + \frac{\partial K_z}{\partial z}\bigg)\text{d}t + \sqrt{2 K_z}\text{d}W(t),
```

where $`\text{d}W(t)`$ is a Wiener noise increment with zero mean and a variance of $`\text{d}t`$. In our case, the displacement due to the settling velocity of a particle is already accounted for in the biofouling kernel, hence we only model the stochastic term (by setting $`w=0`$). To numerically solve this equation, we use the stochastic generalisation of the Euler-forward scheme, called the Euler-Maruyama scheme [@Maruyama1955](http://dx.doi.org/10.1007/BF02846028). 


### Sea-ice capture <a name="seaice"></a>
**TODO once implemented**

### In development:
**TODO -beaching, fragmentation, degradation, etc.???**

## Description of algorithms for particle initialisation maps
Included in the `PlasticParcels` package are four algorithms to create particle initialisation maps, which represent best estimates for plastic pollution emmisions along with the current state of plastic concentrations in our oceans globally. Below we describe each of these algorithms. Each initialisation map, however, requires that particles be placed in ocean grid cells, hence we provide algorithms to generate these intialisation maps, rather than the maps themselves. These maps are land mask dependent, we include scripts to generate a land mask file, as well as a coast mask file, if the model does not provide one.


### Coastal mismanaged plastic waste emissions <a name="coastalrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from coastal communities, we use a global mismanaged plastic waste dataset provided per country [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352). Specifically, for each country, we use the 'Mismanaged plastic waste [kg/person/day]' data to identify the amount of plastic entering the ocean along a coastline. The algorithm is as follows:


**Coastal emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**.
3. Load the Gridded Population of the World dataset [@NASA] **(add reference to .bib)**.
4. For each country in the country boundaries shapefile:
    a. Extract the coordinates of the vertices of the country border (border vertices).
    b. ~~Compute the distance between each border vertex, and every coastal model grid-cell center.~~
    c Create a list of coastal model grid-cells that are within $r$ km of a border vertex.
    d. For each identified coastal model grid-cell, identify the maximum population density from the GPW data within a specified distance $\phi$ (in degrees) north/south or east/west from the coastal model grid-cell center.
    e. Create an array with the coastal model grid-cell and it's associated area, the country name, continent name, region name, and subregion name from the shapefile, and the identified population density.
5. Combine all entries generated in Step 4) a) into one array.
6. Load the global mismanaged plastic waste data [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), and join it to the array generated in Step 5), by 'left joining' on country name$^*$. Create an additional column 'MPW_cell', mismanaged plastic waste across the grid cell by multiplying the mismanaged plastic waste per kilogram per day with the population density and the grid-cell area.
7. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.


$^*$We pre-process the country names in the [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352) data to account for small differences in the naming conventions of each country. Here, we use $r=50$ km, and $\phi$ is chosen as the model grid width in degrees. A sample plot of the initialisation map is shown in Figure X **add link**.

### Riverine mismanaged plastic waste emissions <a name="riverrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from river sources, we use a global riverine input dataset [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803). This dataset is provided in the form of a shapefile, providing a location (latitude and longitude) and amount of plastic emissions (in units of tonnes per year). The algorithm is as follows:

**Riverine emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**, and extract the coordinates of the vertices of every country border.
3. Load the riverine emissions shapefile.
4. For each location in the riverine emissions shapefile:
    a. Compute the distance from the emission source to the center of every coastal grid cell, and identify the closest coastal grid cell.
    b. Compute the distance from the emission source to every country border point, and identify the closest country border point.
    c. Create an array with the coastal model grid-cell, the country name, continent name, region name, and subregion name of the closest border point from the Natural Earth shapefile, and the associated emissions amount.
5. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.


### Open-sea fishing-related plastic emissions <a name="fishingrelease"></a>
To generate a particle initialisation map of plastic pollution emitted into the ocean from fishing-related activities, we use the global fishing watch dataset, first described in [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646). Assuming plastic emissions from fishing lines, trawlers, nets, (etc.?) are proportional to the amount of fishing hours in a given location, we generate a fishing-related plastic emissions initialisation map using the following algorithm:

**Fishing-related emissions initialisation map algorithm**
**(currently we store as monthly, to discuss what we should do, an average over the months? keep the seasonality? etc.)**
1. Load the global fishing watch dataset.
2. Aggregate the dataset by summing the fishing hours over the latitude, longitude, flag, geartype, and date columns, creating a daily dataset.
3. Create a new column containing the month and year of each row, and again aggregate over all columns, replacing the date with the month and year column.
4. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**, and extract a unique list of country name, continent name, region name, subregion name, country name, and country flag (using the ISO standard, or if this is missing/invalid, using the SU_A3 column).
5. Join the country information dataset from Step 4) onto the aggregated dataset from Step 3), by 'left joining' on the flag columns.
6. Load (or generate) the coast mask file from the selected ocean model.
7. Find the closest ocean grid cell for each entry in the aggregated dataset from Step 5) using a KD-Tree approach.
8. Aggregate the data by summing the fishing hours over the following columns: country name, continent name, flag, gear type, date (month and year), ocean grid cell.
9. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.



### Current global ocean plastic concentrations <a name="staterelease"></a>
**TODO once implemented [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0)**



![Particle initialisations maps based on a) current global surface concentrations, b) coastal plastic emissions, c) river plastic emissions, d) fishing activity plastic emissions.](paper/initialisation_maps.png){width=80%}



TODO:
1. change references to DOI Links
2. Update algorithms to be clearer
3. Include additional kernels explanations
4. Include release dataset explanations that are missing
5. Include installation instructions