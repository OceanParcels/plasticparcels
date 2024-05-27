# Functions to create the particle release maps
# @author Michael Denes


# Library imports
import os
import xarray as xr
import numpy as np
import pandas as pd
import urllib.request as urlr
import cartopy.io.shapereader as shpreader
import geopandas as gpd
import glob
from scipy import spatial
from scipy.interpolate import RegularGridInterpolator

from utils import distance, get_coords_from_polygon


# Function definitions
def create_coastal_mpw_jambeck_release_map(mask_coast_filepath, coords_filepath, gpw_filepath,
                                           distance_threshhold=50., grid_range=0.083,
                                           gpw_column_name='Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes',
                                           gpw_raster_number=4,
                                           jambeck_filepath=None,
                                           jambeck_url='https://www.science.org/doi/suppl/10.1126/science.1260352/suppl_file/1260352_supportingfile_suppl._excel_seq1_v1.xlsx',
                                           jambeck_url_backup='http://jambeck.engr.uga.edu/wp-content/uploads/2015/01/JambeckSData.xlsx'
                                           ):
    """
    Description
    ----------

    A function to create a particle release map based on the Mismanaged Plastice Waste data from Jambeck et al. (2015).
    We first load the Natural Earth coastline vector dataset. For each vertex in the vector dataset, we identify any
    associated coastal cell within distance_threshhold km. We then assign a population density from the GPW dataset,
    computed as the maximum population density within grid_range degrees. We then left-join the Jambeck mismanaged
    plastic waste data to compute the total mismanaged plastic waste in kg/day per coastal grid cell.


    Parameters
    ----------

    mask_coast_filepath : str
        File path to the coastal mask generated in ...py TODO
    coords_filepath : str
        File path to the model grid coordinates.
    gpw_filepath : str
        File path to the 'Gridded Population of the World (GPW)' dataset.
    distance_threshhold : float, optional
        The threshhold distance used to find coastal grid cells associated to Natural Earth vector vertices, by default 50
    grid_range : float, optional
        The approximate grid width to search for the maximum population density surrounding a coastal grid cell, by default 0.083
    gpw_column_name : str, optional
        The column name of the GPW dataset that corresponds to the population density, by default 'Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'
    gpw_raster_number : int, optional
        The raster number of the dataset to use. When using the GPW v4.11 data, raster = 4 provides the 2020 population density.
    jambeck_filepath : str, optional
        The filepath to the Jambeck dataset if downloaded manually.
    jambeck_url : str, optional
        The URL of the Jambeck dataset from the Science article.
    jambeck_url_backup : str, optional
        The URL of the Jambeck dataset from Jenna Jambeck's website.

    Returns
    -------
    coastal_density_mpw_df
        A pandas dataframe containing the coastal grid cells, associated countries, and associated mismanaged plastic waste in kg/day.
    """

    # Open the Natural Earth coastline vector dataset
    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    # Open the GPW dataset, and set NaNs to zeros (i.e. zero population density in the ocean)
    gpw = xr.open_dataset(gpw_filepath)
    gpw = gpw.fillna(0.)

    # Load in coast mask and model coordinates
    data_mask_coast = xr.open_dataset(mask_coast_filepath)
    coords = xr.open_dataset(coords_filepath, decode_cf=False)

    lats_coast = data_mask_coast['lat'].data[np.where(data_mask_coast['mask_coast'])]
    lons_coast = data_mask_coast['lon'].data[np.where(data_mask_coast['mask_coast'])]

    # Compute the area of each grid cell
    cell_areas = coords['e1t'][0] * coords['e2t'][0]/10e6  # in km**2
    coastal_cell_areas = cell_areas.data[np.where(data_mask_coast['mask_coast'])]

    # Loop through all countries from Natural Earth dataset
    coastal_density_list = []
    for country in countries:
        # Country information
        continent = country.attributes['CONTINENT']
        region_un = country.attributes['REGION_UN']
        subregion = country.attributes['SUBREGION']
        country_name = country.attributes['NAME_LONG']

        # Extract longitude and latitude of the coastal point
        country_coords = get_coords_from_polygon(country.geometry)
        country_lons, country_lats = country_coords[:, 0], country_coords[:, 1]

        # Loop through country points to find coastal points within distance_threshhold km
        all_coastal_indices = []
        for i in range(len(country_lons)):
            distances = distance(np.repeat(country_lons[i], len(lons_coast)), np.repeat(country_lats[i], len(lats_coast)), lons_coast, lats_coast)
            coastal_indices = np.where(distances <= distance_threshhold)[0]  # Coastal indices are those that are within the thresshold distance
            all_coastal_indices.append(coastal_indices)

        # Concatenate into one list and identify the unique coastal cells
        all_coastal_indices = np.unique(np.hstack(all_coastal_indices))

        # For all coastal points assigned to the country, find the maximum population density around that point
        country_coastal_density_list = []
        for i in range(len(all_coastal_indices)):
            coastal_id = all_coastal_indices[i]
            coastal_point = [lons_coast[coastal_id], lats_coast[coastal_id]]
            coastal_cell_area = coastal_cell_areas[coastal_id]

            # Compute the population density as the maximum population density around the coastal cell
            # Because the grid is ordered longitude ascending, latitude descending, the slice ordering is swapped for lat
            population_density = gpw[gpw_column_name].sel(longitude=slice(coastal_point[0] - grid_range, coastal_point[0] + grid_range),
                                                          latitude=slice(coastal_point[1] + grid_range, coastal_point[1] - grid_range),
                                                          raster=gpw_raster_number).max()

            population_density = np.float32(population_density)

            country_coastal_density_list.append({'Continent': continent,
                                                 'Region': region_un,
                                                 'Subregion': subregion,
                                                 'Country': country_name,
                                                 'Longitude': coastal_point[0],
                                                 'Latitude': coastal_point[1],
                                                 'Area[km2]': coastal_cell_area,
                                                 'PopulationDensity': population_density})

        country_coastal_density_df = pd.DataFrame.from_records(country_coastal_density_list)
        coastal_density_list.append(country_coastal_density_df)

    # Concatenate all the information into one dataframe
    coastal_density_df = pd.concat(coastal_density_list)

    # Open the Jambeck dataset
    if jambeck_filepath is not None:
        MPW = pd.read_excel(jambeck_filepath, 'Jambeck et al. (2014)')
    else:
        try:
            socket = urlr.urlopen(jambeck_url)
        except:
            socket = urlr.urlopen(jambeck_url_backup)
        spreadsheet_from_url = pd.ExcelFile(socket.read())
        MPW = pd.read_excel(spreadsheet_from_url, 'Jambeck et al. (2014)')

    # Perform some data cleaning steps
    # 1. Rename the columns to make it easier to work with
    MPW = MPW.rename(columns={'Economic status1': 'Economic status',
                              'Mismanaged plastic waste [kg/person/day]7': 'Mismanaged plastic waste [kg/person/day]'})

    # 2. Set the bottom rows with black info to blank strings and only keep data with valid country names
    MPW = MPW.fillna('')
    MPW = MPW[MPW['Country'] != '']

    # 3. Remove all sub/superscripts and standardise the 'and' and ampersands.
    MPW['Country'] = MPW['Country'].replace(regex='[0-9]', value='').str.replace('&', 'and')

    # 4. Manually rename some countries to match the Natural Earth Dataset
    MPW['Country'] = MPW['Country'].replace('Brunei', 'Brunei Darussalam')
    MPW['Country'] = MPW['Country'].replace('Curacao', 'Curaçao')
    MPW['Country'] = MPW['Country'].replace('Congo, Dem rep. of', 'Democratic Republic of the Congo')
    MPW['Country'] = MPW['Country'].replace('Korea, North', 'Dem. Rep. Korea')
    MPW['Country'] = MPW['Country'].replace('Korea, South (Republic of Korea)', 'Republic of Korea')
    MPW['Country'] = MPW['Country'].replace('Burma/Myanmar', 'Myanmar')
    MPW['Country'] = MPW['Country'].replace('Micronesia', 'Federated States of Micronesia')
    MPW['Country'] = MPW['Country'].replace('Faroe Islands', 'Faeroe Islands')
    MPW['Country'] = MPW['Country'].replace('Falkland Islands', 'Falkland Islands / Malvinas')
    MPW['Country'] = MPW['Country'].replace('Cote d\'Ivoire', 'Côte d\'Ivoire')
    MPW['Country'] = MPW['Country'].replace('East Timor', 'Timor-Leste')
    MPW['Country'] = MPW['Country'].replace('Russia', 'Russian Federation')
    MPW['Country'] = MPW['Country'].replace('Saint Pierre', 'Saint Pierre and Miquelon')
    MPW['Country'] = MPW['Country'].replace('Congo Rep of', 'Republic of the Congo')
    MPW['Country'] = MPW['Country'].replace('Palestine (Gaza Strip is only part on the coast)', 'Palestine')
    MPW['Country'] = MPW['Country'].replace('Saint Maarten, DWI', 'Sint Maarten')
    MPW['Country'] = MPW['Country'].replace('USVI', 'United States Virgin Islands')
    MPW['Country'] = MPW['Country'].replace('Sao Tome and Principe', 'São Tomé and Principe')

    # TODO: Sort out these countries...
    # Below commented out because it obviously doubles up rows when left-joining.
    # MPW['Country'] = MPW['Country'].replace('Northern Cyprus', 'Cyprus') # Combines the two sets
    # MPW['Country'] = MPW['Country'].replace('Gibraltar', 'United Kingdom') # Add to UK
    # MPW['Country'] = MPW['Country'].replace('French Guiana', 'France')
    # MPW['Country'] = MPW['Country'].replace('Guadeloupe', 'France')
    # MPW['Country'] = MPW['Country'].replace('Martinique', 'France')
    # MPW['Country'] = MPW['Country'].replace('Christmas Island', 'Australia')
    # MPW['Country'] = MPW['Country'].replace('Reunion', 'France')
    # MPW['Country'] = MPW['Country'].replace('Netherlands Antilles', 'Netherlands')

    # Combine the coastal cell density data with the mismanaged plastic waste data, performing a left-join on the country column
    coastal_density_mpw_df = pd.merge(coastal_density_df, MPW[['Country', 'Economic status', 'Mismanaged plastic waste [kg/person/day]']], on="Country", how="left")

    # Compute the mismanaged plastic waste in units of kg/day
    coastal_density_mpw_df['MPW_Cell'] = coastal_density_mpw_df['Area[km2]']*coastal_density_mpw_df['PopulationDensity']*coastal_density_mpw_df['Mismanaged plastic waste [kg/person/day]']

    return coastal_density_mpw_df


def create_rivers_meijer_release_map(mask_coast_filepath, river_filepath):
    """
    Description
    ----------

    A function to create a particle release map based on the rivers mismanaged plastic waste data from Meijer et al. (2021).
    TODO: More description


    Parameters
    ----------
    mask_coast_filepath : str
        File path to the coastal mask generated in ...py TODO
    river_filepath : _type_
        File path to the Meijer et al. (2021) dataset. This dataset can be downloaded here: https://figshare.com/articles/dataset/Supplementary_data_for_More_than_1000_rivers_account_for_80_of_global_riverine_plsatic_emissions_into_the_ocean_/14515590

    Returns
    -------
    river_emissions_df
        A dataset containing the coastal cells associated with plastic emissions from rivers.
    """
    # Load Meijer Data
    data = gpd.read_file(river_filepath)

    lon_river = data.geometry.x
    lat_river = data.geometry.y
    output_river = data['dots_exten'].values  # in tonnes per year

    # Load in coast mask
    data_mask_coast = xr.open_dataset(mask_coast_filepath)
    lats_coast = data_mask_coast['lat'].data[np.where(data_mask_coast['mask_coast'])]
    lons_coast = data_mask_coast['lon'].data[np.where(data_mask_coast['mask_coast'])]

    # Load Natural Earth dataset for attaching country information to river source
    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    countries_list = []
    for country in countries:
        continent = country.attributes['CONTINENT']
        region_un = country.attributes['REGION_UN']
        subregion = country.attributes['SUBREGION']
        country_name = country.attributes['NAME_LONG']

        country_coords = get_coords_from_polygon(country.geometry)
        country_lons, country_lats = country_coords[:, 0], country_coords[:, 1]

        country_df = pd.DataFrame({'Continent': np.repeat(continent, len(country_lons)),
                                   'Region': np.repeat(region_un, len(country_lons)),
                                   'Subregion': np.repeat(subregion, len(country_lons)),
                                   'Country': np.repeat(country_name, len(country_lons)),
                                   'Longitude': country_lons,
                                   'Latitude': country_lats})
        countries_list.append(country_df)
    coastal_df = pd.concat(countries_list)

    # Create river emissions dataset
    river_emissions_list = []

    for i, (lon, lat, output) in enumerate(zip(lon_river, lat_river, output_river)):
        # Find closest coastal cell
        distances = distance(np.repeat(lon, len(lons_coast)), np.repeat(lat, len(lats_coast)), lons_coast, lats_coast)
        closest_coast_id = np.argmin(distances)

        # Find closest country point to river point to assign country information
        distances_country = distance(np.repeat(lon, len(coastal_df['Longitude'])),
                                     np.repeat(lat, len(coastal_df['Latitude'])),
                                     coastal_df['Longitude'],
                                     coastal_df['Latitude'])
        closest_country_id = np.argmin(distances_country)

        river_emissions_list.append({'Continent': coastal_df['Continent'].iloc[closest_country_id],
                                     'Region': coastal_df['Region'].iloc[closest_country_id],
                                     'Subregion': coastal_df['Subregion'].iloc[closest_country_id],
                                     'Country': coastal_df['Country'].iloc[closest_country_id],
                                     'Longitude': lons_coast[closest_coast_id],
                                     'Latitude': lats_coast[closest_coast_id],
                                     'Emissions': output})
    river_emissions_df = pd.DataFrame.from_records(river_emissions_list)

    return river_emissions_df


def create_fisheries_gfwv2_release_map(fisheries_filepath, mask_land_filepath):
    """
    Description
    ----------

    A function to create a particle release map based on the Global Fishing Watch v2 data from Kroodsma et al. (2018).
    TODO: More description


    Parameters
    ----------
    fisheries_filepath : str
        File path to the Global Fishing Watch v2 dataset.
    mask_land_filepath : str
        File path to the land mask generated in ...TODO

    Returns
    -------
    agg_data_fisheries_info
        A dataset containing fishing related plastic release locations.

    """
    # Sort all the fisheries files
    fisheries_files = sorted(glob.glob(fisheries_filepath + 'fleet*/*'))

    # Load Global Fishing Watch data and concatenate into one
    data_fisheries = []
    for file_ in fisheries_files:
        data_fisheries_day = pd.read_csv(file_)
        data_fisheries_day = data_fisheries_day[data_fisheries_day['fishing_hours'] > 0]
        data_fisheries.append(data_fisheries_day)

    data_fisheries = pd.concat(data_fisheries, axis=0, ignore_index=True)

    # Aggregate by location, flag, geartype, and date
    agg_data_fisheries = data_fisheries.groupby(['cell_ll_lat', 'cell_ll_lon', 'flag', 'geartype', 'date'])['fishing_hours'].agg('sum').reset_index()

    # Aggregate by location, flag, geartype, and month
    agg_data_fisheries['date'] = pd.to_datetime(agg_data_fisheries['date'])
    agg_data_fisheries['month'] = agg_data_fisheries['date'].values.astype('datetime64[M]')
    agg_data_fisheries = agg_data_fisheries.groupby(['cell_ll_lat', 'cell_ll_lon', 'flag', 'geartype', 'month'])['fishing_hours'].agg('sum').reset_index()

    # Load in the Natural Earth vector dataset
    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    countries_list = []
    for country in countries:
        continent = country.attributes['CONTINENT']
        region_un = country.attributes['REGION_UN']
        subregion = country.attributes['SUBREGION']
        country_name = country.attributes['NAME_LONG']
        country_flag = country.attributes['ISO_A3']

        # Where the Natural Earth dataset doesn't match the ISO standard (e.g. France), use the SU_A3 column
        if country_flag == '-99':
            country_flag = country.attributes['SU_A3']

        countries_list.append({'Continent': continent,
                               'Region': region_un,
                               'Subregion': subregion,
                               'Country': country_name,
                               'flag': country_flag})
    countries_df = pd.DataFrame.from_records(countries_list)

    # Match flag to country attributes from the Natural Earth dataset, and cleanup column names
    agg_data_fisheries_info = pd.merge(agg_data_fisheries, countries_df, on='flag', how='left')
    agg_data_fisheries_info = agg_data_fisheries_info.rename(columns={'cell_ll_lat': 'Latitude',
                                                                      'cell_ll_lon': 'Longitude',
                                                                      'flag': 'Flag',
                                                                      'geartype': 'Geartype',
                                                                      'month': "Month"})

    # Create a smaller dataset where the fishing hours are matched to the model grid
    model_agg_data_fisheries_info = agg_data_fisheries_info.copy(deep=True)
    data_mask_land = xr.open_dataset(mask_land_filepath)

    lats_ocean = data_mask_land['lat'].data[np.where(~data_mask_land['mask_land'])]
    lons_ocean = data_mask_land['lon'].data[np.where(~data_mask_land['mask_land'])]

    # A list of ocean points
    ocean_points = np.array([lons_ocean, lats_ocean]).T

    # Find closest ocean cell to the fisheries data
    fishing_points = np.array(model_agg_data_fisheries_info[['Longitude', 'Latitude']])
    distances_deg, indices = spatial.KDTree(ocean_points).query(fishing_points)

    mapped_ocean_points = ocean_points[indices]

    # Add these mapped points to the fisheries data as extra columns
    model_agg_data_fisheries_info.insert(0, "ModelLongitude", mapped_ocean_points[:, 0])
    model_agg_data_fisheries_info.insert(1, "ModelLatitude", mapped_ocean_points[:, 1])

    # Create a data set where longitude and latitude are from the model
    model_agg_data_fisheries_info = model_agg_data_fisheries_info.groupby(['ModelLongitude', 'ModelLatitude', 'Flag', 'Geartype', 'Month', 'Continent', 'Region', 'Subregion', 'Country'])['fishing_hours'].agg('sum').reset_index()
    model_agg_data_fisheries_info = model_agg_data_fisheries_info.rename(columns={'ModelLatitude': 'Latitude', 'ModelLongitude': 'Longitude'})

    # Return both raw and model datasets
    return agg_data_fisheries_info, model_agg_data_fisheries_info


def create_global_concentrations_kaandorp_release_map(mask_land_filepath, mask_coast_filepath, coords_filepath,
                                                      kaandorp_filepath, distance_thresshold=50.):
    """
    Description
    ----------

    A function to create a particle release map based on the global concentrations data produced by [@Kaandorp2023].
    We match this data to the model coastal and ocean cells for a better coverage release.
    We use only the 2020 data for all plastic class sizes, and for the ocean cells we only use the surface 0-5m depth.
    """

    # Load in the data and select year 2020, all plastic sizes, and for ocean concentrations take the surface layer
    ds = xr.open_dataset(kaandorp_filepath)
    beach = ds['concentration_beach_mass_log10'].sel({'size_nominal': 'all', 'time': 2020})
    ocean = ds['concentration_mass_log10'].sel({'size_nominal': 'all', 'time': 2020, 'depth': '0 - <5'})

    # Load in coast mask and model coordinates
    data_mask_coast = xr.open_dataset(mask_coast_filepath)

    lats_coast = data_mask_coast['lat'].data[np.where(data_mask_coast['mask_coast'])]
    lons_coast = data_mask_coast['lon'].data[np.where(data_mask_coast['mask_coast'])]

    # Tackle the beach concentrations first:
    # Create a list of lons, lats, concentration values
    lon_beach = beach.lon_beach.values
    lat_beach = beach.lat_beach.values
    conc_beach = np.power(10, beach.values)  # beach.values are in log10 form

    # Load Natural Earth dataset for attaching country information to beach source
    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    countries_list = []
    for country in countries:
        continent = country.attributes['CONTINENT']
        region_un = country.attributes['REGION_UN']
        subregion = country.attributes['SUBREGION']
        country_name = country.attributes['NAME_LONG']

        country_coords = get_coords_from_polygon(country.geometry)
        country_lons, country_lats = country_coords[:, 0], country_coords[:, 1]

        country_df = pd.DataFrame({'Continent': np.repeat(continent, len(country_lons)),
                                   'Region': np.repeat(region_un, len(country_lons)),
                                   'Subregion': np.repeat(subregion, len(country_lons)),
                                   'Country': np.repeat(country_name, len(country_lons)),
                                   'Longitude': country_lons,
                                   'Latitude': country_lats})
        countries_list.append(country_df)
    coastal_df = pd.concat(countries_list)

    # Create coastal concentrations dataset
    coast_concentration_list = []
    distance_threshhold = 50.

    for i, (lon, lat) in enumerate(zip(lons_coast, lats_coast)):
        # Find the closest beach concentration
        distances = distance(np.repeat(lon, len(lon_beach)), np.repeat(lat, len(lat_beach)), lon_beach, lat_beach)
        closest_beach_id = np.argmin(distances)
        if distances[closest_beach_id] > distance_threshhold:  # skip coastal grid cells not within a threshhold from a littered beach
            continue
        else:
            # Find the closest country point to the coastal cell to assign country information
            distances_country = distance(np.repeat(lon, len(coastal_df['Longitude'])),
                                         np.repeat(lat, len(coastal_df['Latitude'])),
                                         coastal_df['Longitude'],
                                         coastal_df['Latitude'])
            closest_country_id = np.argmin(distances_country)

            coast_concentration_list.append({'Continent': coastal_df['Continent'].iloc[closest_country_id],
                                             'Region': coastal_df['Region'].iloc[closest_country_id],
                                             'Subregion': coastal_df['Subregion'].iloc[closest_country_id],
                                             'Country': coastal_df['Country'].iloc[closest_country_id],
                                             'Longitude': lon,
                                             'Latitude': lat,
                                             'Concentration': conc_beach[closest_beach_id],
                                             'ConcentrationType': 'Beach'})

    coast_concentration_df = pd.DataFrame.from_records(coast_concentration_list)

    # Now tackle the surface ocean concentrations:
    conc_ocean = np.power(10, ocean.values)  # Values are in log10 space

    data_mask_land = xr.open_dataset(mask_land_filepath)
    lats_ocean = data_mask_land['lat'].data[np.where(~data_mask_land['mask_land'])]
    lons_ocean = data_mask_land['lon'].data[np.where(~data_mask_land['mask_land'])]

    # Function to interpolate the ocean concentrations
    f_interp_conc_ocean = RegularGridInterpolator((ocean.lon, ocean.lat), conc_ocean.T, method='nearest', bounds_error=False, fill_value=None)

    # Interpolate the ocean concentrations to the model grid cells
    interp_conc_ocean = f_interp_conc_ocean((lons_ocean, lats_ocean))

    # Create ocean concentration dataset where values are non-NaN
    non_nan_id = ~np.isnan(interp_conc_ocean)
    ocean_concentration_list = []
    for i in np.where(non_nan_id)[0]:
        ocean_concentration_list.append({'Continent': 'N/A',
                                         'Region': 'N/A',
                                         'Subregion': 'N/A',
                                         'Country': 'N/A',
                                         'Longitude': lons_ocean[i],
                                         'Latitude': lats_ocean[i],
                                         'Concentration': interp_conc_ocean[i],
                                         'ConcentrationType': 'Ocean'})
    ocean_concentration_df = pd.DataFrame.from_records(ocean_concentration_list)

    # Combine the two beach and ocean datasets
    concentration_df = pd.concat([ocean_concentration_df, coast_concentration_df])

    return concentration_df


output_data = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/release/generated_files/'
mask_land_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/output_data/masks/mask_land_NEMO0083.nc'
mask_coast_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/output_data/masks/mask_coast_NEMO0083.nc'
coords_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/input_data/MOi/domain_ORCA0083-N006/coordinates.nc'

# Create coastal release data
gpw_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/release/GPWv4/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc'
output_name = output_data + 'coastal_population_MPW_NEMO0083.csv'

if not os.path.isfile(output_name):
    coastal_dataset = create_coastal_mpw_jambeck_release_map(mask_coast_filepath=mask_coast_filepath,
                                                             coords_filepath=coords_filepath,
                                                             gpw_filepath=gpw_filepath)
    coastal_dataset.to_csv(output_name)
    print("Coastline mismanaged plastic waste file created:", output_name)
else:
    print("Coastline mismanaged plastic waste file already exists:", output_name)


# Create rivers release data
river_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/release/Meijer2021_midpoint_emissions/'
output_name = output_data + 'river_emissions_NEMO0083.csv'

if not os.path.isfile(output_name):
    river_dataset = create_rivers_meijer_release_map(mask_coast_filepath=mask_coast_filepath,
                                                     river_filepath=river_filepath)
    river_dataset.to_csv(output_name)
    print("River mismanaged plastic waste file created:", output_name)
else:
    print("River mismanaged plastic waste file already exists:", output_name)


# Create fisheries release data
fisheries_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/release/GlobalFishingWatch/'
output_name = output_data+'agg_data_fisheries_info_raw.csv'
output_on_model_name = output_data+'agg_data_fisheries_info_NEMO0083.csv'

if not os.path.isfile(output_name) and not os.path.isfile(output_on_model_name):
    fisheries_dataset, fisheries_on_model_grid = create_fisheries_gfwv2_release_map(fisheries_filepath=fisheries_filepath, mask_land_filepath=mask_land_filepath)
    fisheries_dataset.to_csv(output_name)
    fisheries_on_model_grid.to_csv(output_on_model_name)
    print("Fishing related plastic waste file created:", output_name, 'and', output_on_model_name)
else:
    print("Fishing related plastic waste file already exists:", output_name)

# Create current concentrations release data
kaandorp_filepath = '/Users/denes001/Research/Projects/PlasticParcels/PlasticParcels/data/release/Kaandorp2023/AtlantECO-MAPS_Global_plastic_mass_budget_Kaandorp_etal_2023.nc'
output_name = output_data+'global_concentrations_NEMO0083.csv'

if not os.path.isfile(output_name):
    concentration_dataset = create_global_concentrations_kaandorp_release_map(mask_land_filepath, mask_coast_filepath, coords_filepath, kaandorp_filepath, distance_thresshold=50.)
    concentration_dataset.to_csv(output_name)
    print("Concentration map file created:", output_name)
else:
    print("Concentration map file already exists:", output_name)
