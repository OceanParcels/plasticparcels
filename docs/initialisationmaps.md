# Description of algorithms for particle initialisation maps
Included in the `plasticparcels` package are four particle initialisation maps, along with the algorithms to create them. These maps represent best estimates for plastic pollution emissions along coastlines [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), from river sources [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803), in the open-ocean from fishing-related activities [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646), as well as a current best estimate of buoyant plastic concentrations globally [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0).
Each initialisation map, however, requires that particles be placed in ocean grid cells, so we also provide algorithms to generate these ocean masks too.

The code for these algorithims can be found in `plasticparcels/scripts/create_release_maps.py`. Below we describe each of the algorithms.


## Coastal mismanaged plastic waste emissions <a name="coastalrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from coastal communities, we use a global mismanaged plastic waste dataset provided per country [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352). Specifically, for each country, we use the 'Mismanaged plastic waste [kg/person/day]' (MPW) data to identify the amount of plastic entering the ocean along a coastline. We utilise a country border shapefile, a coastal ocean grid-cell mask, and the MPW dataset to create a `.csv` file. The algorithm is as follows:


**Coastal emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/).
3. Load the Gridded Population of the World dataset [@NASA](https://sedac.ciesin.columbia.edu/data/collection/gpw-v4).
4. For each country in the country boundaries shapefile:
    1. Extract the coordinates of the vertices of the country border (border vertices).
    2. Create a list of coastal model grid-cells that are within $r$ km of a border vertex.
    3. For each identified coastal model grid-cell, identify the maximum population density from the GPW data within a specified distance $\phi$ (in degrees) north/south or east/west from the coastal model grid-cell center.
    4. Create an array with the coastal model grid-cell and its associated area, the country name, continent name, region name, and subregion name from the shapefile, and the identified population density.
5. Combine all entries generated in Step 4.4. into one array.
6. Load the global mismanaged plastic waste data [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352), and join it to the array generated in Step 5, by 'left joining' on country name$^*$. Create an additional column 'MPW_cell', which represents the mismanaged plastic waste across the grid cell, by multiplying the mismanaged plastic waste per kilogram per day with the population density and the grid-cell area.
7. Save the data into a `.csv` file, to be read and processed by `plasticparcels`.


$^*$We pre-process the country names in the [@Jambeck2015](http://dx.doi.org/10.1126/science.1260352) data to account for small differences in the naming conventions of each country. We use $`r=50`$ km, and $`\phi`$ is chosen as the model grid width in degrees. A sample plot of the initialisation map is shown in Figure X **add link**.

## Riverine mismanaged plastic waste emissions <a name="riverrelease"></a>
To generate a particle initialisation map of plastic pollution that enters the ocean from river sources, we use a global riverine input dataset [@Meijer2021](http://dx.doi.org/10.1126/sciadv.aaz5803). This dataset is in the form of a shapefile, providing a location (latitude and longitude) and amount of plastic emissions (in units of tonnes per year). The algorithm is as follows:

**Riverine emissions initialisation map algorithm**
1. Load (or generate) the coast mask file from the selected ocean model.
2. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/), and extract the coordinates of the vertices of every country border.
3. Load the riverine emissions shapefile.
4. For each location in the riverine emissions shapefile:
    1. Compute the distance from the emission source to the center of every coastal grid cell, and identify the closest coastal grid cell.
    2. Compute the distance from the emission source to every country border point, and identify the closest country border point.
    3. Create an array with the coastal model grid-cell, the country name, continent name, region name, and subregion name of the closest border point from the Natural Earth shapefile, and the associated emissions amount.
5. Save the data into a `.csv` file, to be read and processed by `plasticparcels`.


## Open-sea fishing-related plastic emissions <a name="fishingrelease"></a>
To generate a particle initialisation map of plastic pollution emitted into the ocean from fishing-related activities, we use the global fishing watch dataset, first described in [@Kroodsma2018](http://dx.doi.org/10.1126/science.aao5646), which is provided as a set of daily csv files. Assuming plastic emissions from fishing activity is proportional to the amount of fishing hours in a given location, we generate a fishing-related plastic emissions initialisation map using the following algorithm:

**Fishing-related emissions initialisation map algorithm**
1. Load the global fishing watch dataset$^*$.
2. Aggregate the dataset by summing the fishing hours over the latitude, longitude, flag, geartype, and date columns, creating an aggregated daily dataset.
3. Create a new column containing the month and year of each row, and again aggregate over all columns, replacing the date with the month-year column.
4. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/), and extract a unique list of country name, continent name, region name, subregion name, country name, and country flag (using the ISO standard, ISO_A3 column, or if this is missing/invalid, using the SU_A3 column).
5. Join the country information dataset from Step 4. onto the aggregated dataset from Step 3., by 'left joining' matching the flag columns.
6. Load (or generate) the coast mask file from the selected ocean model.
7. Find the closest ocean grid cell for each entry in the aggregated dataset from Step 5. using a KD-Tree approach.
8. Aggregate the data by summing the fishing hours over the following columns: country name, continent name, flag, gear type, date (month and year), ocean grid cell.
9. Save the data into a `.csv` file, to be read and processed by `plasticparcels`.

$^*$We use the `fleet-daily-csvs-100-v2-2020` files, which are for the year 2020 only.


## Current global ocean plastic concentrations <a name="staterelease"></a>
To generate a particle initialisation map of the current best-estimate of global plastic concentrations in the ocean and along coastlines, we use an estimate produced in [@Kaandorp2023](http://dx.doi.org/10.1038/s41561-023-01216-0). This is a gridded dataset provided in a `netCDF` format. Specifically, we create a particle initialisation map based on the modelled global concentrations along coastlines in the year 2020 for all plastic types, as well as a map for the open ocean in the year 2020 for all plastic types in the near-surface ocean. To generate this initialisation map we use the following algorithm.

**Current global ocean plastic concentrations initialisation map algorithm**
1. Load the global concentrations dataset.
2. Extract two sets of data, the first for `concentration_beach_mass_log10`, containing only the data for `size_nominal='all'`, `time=2020`, the second for `concentration_mass_log10`, containing only the data for `size_nominal='all'`, `time=2020`, and `depth='0 - <5m'`.
3. Load (or generate) the coast mask file from the selected ocean model.
4. Load the Natural Earth country boundaries shapefile at 1:50m resolution [@NaturalEarth](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/), and extract the coordinates of the vertices of every country border.
5. For every coastal model cell, identify the closest `(lon_beach, lat_beach)` within a 50km radius from the `concentration_beach_mass_log10` dataset, and the closest country boundary vertex from the Natural Earth shapefile.
6. Create an array with the coastal model cell, the plastic concentration amount from the `concentration_beach_mass_log10` dataset (converting it into a mass instead of a `log10` mass), and the continent name, region name, subregion name, country name, and country flag from the Natural Earth shapefile.
7. Load (or generate) the land mask file from the selected ocean model.
8. Interpolate the `concentration_mass_log10` to the ocean-grid cells, using  an `RegularGridInterpolator` function from `scipy.interpolate`, with the grid and data being `(lon, lat)` and `concentration_mass_log10` from the `concentration_mass_log10` dataset.
9. For all valid concentrations identified in Step 8., identify the closest country boundary vertex from the Natural Earth shapefile.
10. Create an array with the ocean model cell, the interpolated plastic concentration amount (converting it into a mass insteaf of a `log10` mass), and the continent name, region name, subregion name, country name, and country flag from the Natural Earth shapefile.
11. Combine the arrays generated in Steps 6. and 10., and save the data as a `.csv` file, to be read and processed `plasticparcels`.