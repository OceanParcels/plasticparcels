# Description of algorithms for particle initialisation maps
Included in the `PlasticParcels` package are four algorithms to create particle initialisation maps, which represent best estimates for plastic pollution emissions along with the current state of plastic concentrations in our oceans globally. Below we describe each of these algorithms. Each initialisation map, however, requires that particles be placed in ocean grid cells, hence we provide algorithms to generate these initialisation maps, rather than the maps themselves. These maps are land mask dependent, we include scripts to generate a land mask file, as well as a coast mask file, if the model does not provide one.


### Coastal mismanaged plastic waste emissions <a name="coastalrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from coastal communities, we use a global mismanaged plastic waste dataset provided per country [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352). Specifically, for each country, we use the 'Mismanaged plastic waste [kg/person/day]' data to identify the amount of plastic entering the ocean along a coastline. The algorithm is as follows:


**Coastal emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**.
3. Load the Gridded Population of the World dataset [@NASA] **(add reference to .bib)**.
4. For each country in the country boundaries shapefile:
    1. Extract the coordinates of the vertices of the country border (border vertices).
    2. ~~Compute the distance between each border vertex, and every coastal model grid-cell center.~~ <mark>Remove this?</mark>
    3. Create a list of coastal model grid-cells that are within $r$ km of a border vertex.
    4. For each identified coastal model grid-cell, identify the maximum population density from the GPW data within a specified distance $\phi$ (in degrees) north/south or east/west from the coastal model grid-cell center.
    5. Create an array with the coastal model grid-cell and its associated area, the country name, continent name, region name, and subregion name from the shapefile, and the identified population density.
5. Combine all entries generated in Step 4.i. into one array.
6. Load the global mismanaged plastic waste data [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), and join it to the array generated in Step 5, by 'left joining' on country name$`^*`$. Create an additional column 'MPW_cell', mismanaged plastic waste across the grid cell by multiplying the mismanaged plastic waste per kilogram per day with the population density and the grid-cell area.
7. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.


$`^*`$We pre-process the country names in the [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352) data to account for small differences in the naming conventions of each country. Here, we use $`r=50`$ km, and $`\phi`$ is chosen as the model grid width in degrees. A sample plot of the initialisation map is shown in Figure X **add link**.

### Riverine mismanaged plastic waste emissions <a name="riverrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from river sources, we use a global riverine input dataset [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803). This dataset is provided in the form of a shapefile, providing a location (latitude and longitude) and amount of plastic emissions (in units of tonnes per year). The algorithm is as follows:

**Riverine emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**, and extract the coordinates of the vertices of every country border.
3. Load the riverine emissions shapefile.
4. For each location in the riverine emissions shapefile:
    1. Compute the distance from the emission source to the center of every coastal grid cell, and identify the closest coastal grid cell.
    2. Compute the distance from the emission source to every country border point, and identify the closest country border point.
    3. Create an array with the coastal model grid-cell, the country name, continent name, region name, and subregion name of the closest border point from the Natural Earth shapefile, and the associated emissions amount.
5. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.


### Open-sea fishing-related plastic emissions <a name="fishingrelease"></a>
To generate a particle initialisation map of plastic pollution emitted into the ocean from fishing-related activities, we use the global fishing watch dataset, first described in [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646). Assuming plastic emissions from fishing activity is proportional to the amount of fishing hours in a given location, we generate a fishing-related plastic emissions initialisation map using the following algorithm:

**Fishing-related emissions initialisation map algorithm**
**(currently we store as monthly, to discuss what we should do, an average over the months? keep the seasonality? etc.)**
1. Load the global fishing watch dataset.
2. Aggregate the dataset by summing the fishing hours over the latitude, longitude, flag, geartype, and date columns, creating a daily dataset.
3. Create a new column containing the month and year of each row, and again aggregate over all columns, replacing the date with the month and year column.
4. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth] **(add reference to .bib)**, and extract a unique list of country name, continent name, region name, subregion name, country name, and country flag (using the ISO standard, or if this is missing/invalid, using the SU_A3 column).
5. Join the country information dataset from Step 4) onto the aggregated dataset from Step 3, by 'left joining' on the flag columns.
6. Load (or generate) the coast mask file from the selected ocean model.
7. Find the closest ocean grid cell for each entry in the aggregated dataset from Step 5 using a KD-Tree approach.
8. Aggregate the data by summing the fishing hours over the following columns: country name, continent name, flag, gear type, date (month and year), ocean grid cell.
9. Save the data into a `.csv` file, to be read and processed by `PlasticParcels`.



### Current global ocean plastic concentrations <a name="staterelease"></a>
**TODO once implemented [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0)**



![Particle initialisations maps based on a) current global surface concentrations, b) coastal plastic emissions, c) river plastic emissions, d) fishing activity plastic emissions.](paper/initialisation_maps.png){width=80%}



TODO:
1. ~~change references to DOI Links~~ -- add links to non-published papers
2. Update algorithms to be clearer
3. Include additional kernels explanations
4. Include release dataset explanations that are missing
5. Include installation instructions