import os
import numpy as np
import xarray as xr


import pandas as pd
from parcels import FieldSet, Field, ParticleSet, JITParticle, Variable, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.converters import Geographic, GeographicPolar
from datetime import datetime, timedelta

from PlasticParcels.kernels import PolyTEOS10_bsq, Stokes_drift, windage_drift, settling_velocity, biofouling, vertical_mixing, unbeaching, periodicBC, checkErrorThroughSurface, deleteParticle, checkThroughBathymetry
from PlasticParcels.utils import select_files, getclosest_ij


## Fieldset creation
def create_fieldset(model_settings, particle_settings):
    """ Helper function to create a Parcels.Fieldset

    Parameters
    ----------
    model_settings :
        A dictionary of model settings used to create the fieldset
    particle_settings :
        A dictionary of particle settings used to define some ....

    Returns
    -------
    fieldset
        A parcels.FieldSet object
    """     
    
    ## Directories and arrays of files
    dirread_model = os.path.join(model_settings['input_data_dir'], model_settings['ocean_dir'], model_settings['ocean_filename'])

    start_date = particle_settings['start_date']
    runtime = int(np.ceil(particle_settings['runtime'].total_seconds()/86400.)) # convert to days

    ## If running backward simulations, we need to change the start date and runtime for loading the fieldset
    #if np.sign(runtime) < 0:
    #    start_date = start_date + particle_settings['runtime']
    #    runtime = -runtime
    
    
    ## Mesh masks
    ocean_mesh = os.path.join(model_settings['input_data_dir'], model_settings['ocean_mesh'])      #mesh_mask

    ## Setup input for fieldset creation
    ufiles = select_files(dirread_model,'U_%4i*_med_subset.nc',start_date,runtime,i_date_s=-24, i_date_e=-14,dt_margin=3)
    vfiles = select_files(dirread_model,'V_%4i*_med_subset.nc',start_date,runtime,i_date_s=-24, i_date_e=-14,dt_margin=3)
    wfiles = select_files(dirread_model,'W_%4i*_med_subset.nc',start_date,runtime,i_date_s=-24, i_date_e=-14,dt_margin=3)   

    filenames = {'U': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': ufiles},
                    'V': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': vfiles}}
                    #'W': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': wfiles}}

    variables = {'U': 'vozocrtx',
                    'V': 'vomecrty'}
                    #'W': 'vovecrtz'}

    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}, #time_centered
                    'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}
                    #'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}

    # Setup indices for specific regions if required
    if 'ocean_region' not in model_settings.keys():
        indices = {}
    elif model_settings['ocean_region'] == 'MED': # Mediterranean
        minlat = 25
        maxlat = 50
        minlon = -5
        maxlon = 38

        test_data = xr.open_dataset(ufiles[0])
        latvals = test_data['nav_lat'].values
        lonvals = test_data['nav_lon'].values
        
        iy_min, ix_min = getclosest_ij(latvals, lonvals, minlat, minlon)
        iy_max, ix_max = getclosest_ij(latvals, lonvals, maxlat, maxlon)
        iy_min -= 1
        ix_min -= 1
        iy_max += 1
        ix_max += 1
        indices = {'lat': range(iy_min, iy_max), 'lon': range(ix_min, ix_max)}
        test_data.close()

    ## Load the fieldset
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=model_settings['allow_time_extrapolation'])

    ## Create flags for custom particle behaviour
    fieldset.add_constant('mixing_f', model_settings['mixing_f'])
    fieldset.add_constant('biofouling_f', model_settings['biofouling_f'])
    fieldset.add_constant('stokes_f', model_settings['stokes_f'])
    fieldset.add_constant('wind_f', model_settings['wind_f'])
    fieldset.add_constant('G', 9.81) # Gravitational constant [m s-1]
    
    mode_3D_f = False
    if model_settings['mode'] == '3D':
        mode_3D_f = True
    fieldset.add_constant('mode', mode_3D_f)
    
    # Add in bathymetry
    # fieldset.add_constant('z_start',0.5)        
    # variables = ('bathymetry', 'Bathymetry')
    # dimensions = {'lon': 'nav_lon', 'lat': 'nav_lat'}
    # bathymetry_mesh = os.path.join(model_settings['input_data_dir'], model_settings['bathymetry_mesh'])
    # bathymetry_field = Field.from_netcdf(bathymetry_mesh, variables, dimensions)
    # fieldset.add_field(bathymetry_field) 

  
    ## Add fields for different model settings
    if fieldset.mixing_f:
        kzfiles = select_files(dirread_model,'KZ_%4i*.nc',start_date,runtime,dt_margin=3)
        filenames_mixing = {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': kzfiles}
        variables_mixing = {'mixing_kz' : 'votkeavt'}
        dimensions_mixing = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

        mixing_fieldset = FieldSet.from_nemo(filenames_mixing,variables_mixing,dimensions_mixing)
        fieldset.add_field(mixing_fieldset.mixing_kz)    #phytoplankton primary productivity        

        


    if fieldset.biofouling_f: # or do_permanent_fouling:
        #MOi glossary: https://www.mercator-ocean.eu/wp-content/uploads/2021/11/Glossary.pdf
        # and https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-028.pdf
        
        dirread_bgc = os.path.join(model_settings['input_data_dir'], model_settings['bgc_dir'], model_settings['bgc_filename'])
        bgc_mesh = os.path.join(model_settings['input_data_dir'], model_settings['bgc_mesh'])   #mesh_mask_4th

        ppfiles = select_files(dirread_bgc,'nppv_%4i*.nc',start_date,runtime,dt_margin=8)
        phy1files = select_files(dirread_bgc,'phy_%4i*.nc',start_date,runtime,dt_margin=8)
        phy2files = select_files(dirread_bgc,'phy2_%4i*.nc',start_date,runtime,dt_margin=8)

        #TODO: cleanup these names for better readability
        filenames_bio = {'pp_phyto': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': ppfiles},
                        'bio_nanophy': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy1files},
                        'bio_diatom': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy2files}}

        variables_bio = {'pp_phyto': 'nppv', # Total Primary Production of Phyto - 'Net primary prodution of biomass expressed as arbon per unit volume in sea water' [mg m-3 day-1] or [milligrams of Carbon per cubic meter per day]
                        'bio_nanophy': 'phy', # Mole concentration of NanoPhytoplankton expressed as carbon in sea water
                        'bio_diatom': 'phy2'} # Mole concentration of Diatoms expressed as carbon in sea water

        dimensions_bio = {'pp_phyto': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_nanophy': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_diatom': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}
        
        
        fieldset.add_constant('collision_probability', model_settings['bgc_collision_probability'])
        fieldset.add_constant('K', model_settings['bgc_boltzmann_constant'])
        fieldset.add_constant('R20', model_settings['bgc_respiration_rate'])

        fieldset.add_constant('Q10', model_settings['bgc_respiration_temperature_coefficient'])
        fieldset.add_constant('Gamma', model_settings['bgc_shear_frequency'])
        fieldset.add_constant('carbon_molecular_weight', model_settings['bgc_carbon_atomic_weight'])
        fieldset.add_constant('algae_cell_volume', model_settings['bgc_algae_cell_volume'])
        fieldset.add_constant('biofilm_density', model_settings['bgc_biofilm_density'])    # density of biofilm [g m-3]
        fieldset.add_constant('algae_mortality_rate', model_settings['bgc_algae_mortality_rate'])
        fieldset.add_constant('algae_respiration_f', model_settings['bgc_algae_respiration_f'])
        bio_fieldset = FieldSet.from_nemo(filenames_bio,variables_bio,dimensions_bio)

        fieldset.add_field(bio_fieldset.pp_phyto)    #phytoplankton primary productivity 
        fieldset.add_field(bio_fieldset.bio_nanophy) #nanopyhtoplankton concentration [mmol C m-3]
        fieldset.add_field(bio_fieldset.bio_diatom)  #diatom concentration [mmol C m-3]   


    if fieldset.stokes_f:
        ## TODO: Update input_data_dir_2 once files are copied
        dirread_stokes = os.path.join(model_settings['input_data_dir_2'], model_settings['stokes_dir'], model_settings['stokes_filename'])
        wavesfiles = select_files(dirread_stokes,'%4i*.nc',start_date,runtime,dt_margin=32)

        filenames_S = {'Stokes_U': wavesfiles,
                        'Stokes_V': wavesfiles,
                        'wave_Tp': wavesfiles}

        variables_S = {'Stokes_U': 'ust',
                        'Stokes_V': 'vst',
                        'wave_Tp': 'pp1d'}

        dimensions_S = {'lat': 'latitude', # when not using the converted datasets, use 'longitude' otherwise use 'lon'
                        'lon': 'longitude',
                        'time': 'time'}
        
        fieldset_Stokes = FieldSet.from_netcdf(filenames_S, variables_S, dimensions_S, mesh='spherical')
        fieldset_Stokes.Stokes_U.units = GeographicPolar()
        fieldset_Stokes.Stokes_V.units = Geographic()
        fieldset_Stokes.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_Stokes.Stokes_U)
        fieldset.add_field(fieldset_Stokes.Stokes_V)
        fieldset.add_field(fieldset_Stokes.wave_Tp)

    if fieldset.wind_f:
        ## TODO: Update input_data_dir_2 once files are copied
        dirread_wind = os.path.join(model_settings['input_data_dir_2'], model_settings['wind_dir'], model_settings['wind_filename'])

        windfiles = select_files(dirread_wind,'%4i*.nc',start_date,runtime,dt_margin=32)

        filenames_wind = {'Wind_U': windfiles,
                        'Wind_V': windfiles}

        variables_wind = {'Wind_U': 'u10',
                        'Wind_V': 'v10'}

        dimensions_wind = {'lat': 'latitude', # when not using the converted datasets, use 'longitude' otherwise use 'lon'
                        'lon': 'longitude',
                        'time': 'time'}
        fieldset_wind = FieldSet.from_netcdf(filenames_wind, variables_wind, dimensions_wind, mesh='spherical')
        fieldset_wind.Wind_U.units = GeographicPolar()
        fieldset_wind.Wind_V.units = Geographic() 
        fieldset_wind.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_wind.Wind_U)
        fieldset.add_field(fieldset_wind.Wind_V)


    ## Apply unbeaching currents when Stokes/Wind can push particles into land cells
    if fieldset.stokes_f or fieldset.wind_f > 0:
        unbeachfiles = './data/land_current_NEMO0083_med_subset.nc'
        filenames_unbeach = {'unbeach_U': unbeachfiles, 
                            'unbeach_V': unbeachfiles}

        variables_unbeach = {'unbeach_U': 'land_current_u',
                            'unbeach_V': 'land_current_v'}   

        dimensions_unbeach = {'lat': 'lat',
                            'lon': 'lon'}
        

        fieldset_unbeach = FieldSet.from_netcdf(filenames_unbeach, variables_unbeach, dimensions_unbeach, mesh='spherical')
        fieldset_unbeach.unbeach_U.units = GeographicPolar()
        fieldset_unbeach.unbeach_V.units = Geographic()

        fieldset.add_field(fieldset_unbeach.unbeach_U)
        fieldset.add_field(fieldset_unbeach.unbeach_V)
    

    ## TODO: Write a comment for this
    fieldset.add_constant('verbose_delete',model_settings['verbose_delete'])

    return fieldset


def create_particleset(fieldset, particle_settings):
    """ Helper function to create a Parcels.ParicleSet

    Parameters
    ----------
    model_settings :
        A dictionary of model settings used to create the fieldset
    particle_settings :
        A dictionary of particle settings used to define some ....

    Returns
    -------
    fieldset
        A parcels.FieldSet object
    """    
    release_locations = particle_settings['release_locations']
    
    ## Set the longitude, latitude, depth and time of the particles
    lons = release_locations['lons']
    lats = release_locations['lats']
    depths = None
    times = None
    if 'depths' in release_locations.keys():
        depths = release_locations['depths']
    if 'times' in release_locations.keys():
        times = release_locations['times']



    ## Set particle densities 
    if type(particle_settings['particle_density']) == float:
        particle_densities = np.full(lons.shape, particle_settings['particle_density'])
    else:
        particle_densities = particle_settings['particle_density']

    ## Set particle lengths
    if type(particle_settings['particle_diameter']) == float:
        particle_diameters = np.full(lons.shape, particle_settings['particle_diameter'])
    else:
        particle_diameters = particle_settings['particle_diameter']


    ## Set wind coefficients of particles
    if 'windage_coefficient' in particle_settings.keys():
        if type(particle_settings['windage_coefficient']) == float:
            windage_coefficients = np.full(lons.shape, particle_settings['windage_coefficient'])
        else:
            windage_coefficients = particle_settings['windage_coefficient'] # Assumed to be an array of coefficients


    ## Add variables to particle based on fieldset flags
    to_write_tracer = False
    to_write_dynamic = False
    to_write_all = False
    if 'write_output_option' not in particle_settings.keys():
        pass
    elif particle_settings['write_output_option'] == 'none':
        pass
    elif particle_settings['write_output_option'] == 'tracer':
        to_write_tracer = True
    elif particle_settings['write_output_option'] == 'dynamic':
        to_write_tracer = True
        to_write_dynamic = True
    elif particle_settings['write_output_option'] == 'all':
        to_write_tracer = True
        to_write_dynamic = True
        to_write_all = True
    ## TODO: Write some catch for unimplemented options... raise error or something


    ## Add variables to particle based on fieldset flags
    variables = []
    variables.append(Variable('particle_diameter', dtype=np.float32, to_write=to_write_all))                   # Particle Diameter (assuming spherical particle) [meters] (l_pl)
    variables.append(Variable('particle_density', dtype=np.float32, to_write=to_write_all))                            # Particle Density [kg/m^3] (rho_pl)
    variables.append(Variable('settling_velocity', dtype=np.float64, initial=0., to_write=to_write_dynamic))              # Particle Sinking Velocity [m/s] (v_s)
    variables.append(Variable('seawater_density', dtype=np.float32, to_write=to_write_dynamic))
    variables.append(Variable('absolute_salinity', dtype=np.float64, to_write=to_write_dynamic))
    variables.append(Variable('windage_coefficient', dtype=np.float32, initial=0., to_write=to_write_all))
    variables.append(Variable('algae_amount', dtype=np.float64, initial=0., to_write=to_write_dynamic))

    # Create PlasticParticle class
    PlasticParticle = JITParticle
    for variable in variables:
        setattr(PlasticParticle, variable.name, variable)

    # Add kernel specific variables
    # TODO
    #if fieldset.biofouling_f:
    #if fieldset.wind_f:
        
    

    ## Create the particle set
    if not fieldset.wind_f:
        pset = ParticleSet.from_list(fieldset,
                                    PlasticParticle,
                                    lon=lons,
                                    lat=lats,
                                    time=times,
                                    depth=depths,
                                    particle_diameter=particle_diameters,
                                    particle_density=particle_densities
                                    )
    else:
        pset = ParticleSet.from_list(fieldset,
                                    PlasticParticle,
                                    lon=lons,
                                    lat=lats,
                                    time=times,
                                    depth=depths,
                                    particle_diameter=particle_diameters,
                                    particle_density=particle_densities,
                                    windage_coefficient=windage_coefficients
                                    )
    return pset

def create_particleset_from_file(fieldset, particle_settings):
    """ Helper function to create a Parcels.ParicleSet

    Parameters
    ----------
    model_settings :
        A dictionary of model settings used to create the fieldset
    particle_settings :
        A dictionary of particle settings used to define some ....

    Returns
    -------
    fieldset
        A parcels.FieldSet object
    """

    # Load release type information
    release_type = particle_settings['release_type']
    if release_type == 'coastal':
        release_file = './data/coastal_population_MPW_NEMO0083.csv'
        release_quantity_name = 'MPW_Cell'
    elif release_type == 'rivers':
        release_file = './data/river_emissions_NEMO0083.csv'
        release_quantity_name = 'Emissions'
    elif release_type == 'fisheries':
        release_file = './data/agg_data_fisheries_info.csv'
        release_quantity_name = 'fishing_hours'

    particle_locations = pd.read_csv(release_file)


    # Select specific continent/region/subregion/country/economic status if applicable:
    if 'continent' in particle_settings.keys():
        particle_locations = particle_locations[particle_locations['Continent'] == particle_settings['continent']]
    if 'region' in particle_settings.keys():
        particle_locations = particle_locations[particle_locations['Region'] == particle_settings['region']]
    if 'subregion' in particle_settings.keys():
        particle_locations = particle_locations[particle_locations['Subregion'] == particle_settings['subregion']]
    if 'country' in particle_settings.keys():
        particle_locations = particle_locations[particle_locations['Country'] == particle_settings['country']]
    if 'economicstatus' in particle_settings.keys():
        particle_locations = particle_locations[particle_locations['Economic status'] == particle_settings['economicstatus']]

    particle_locations = particle_locations.groupby(['Longitude', 'Latitude'])[release_quantity_name].agg('sum').reset_index()
    particle_locations = particle_locations[particle_locations[release_quantity_name]>0]

    release_locations = {'lons': particle_locations['Longitude'],
                        'lats':  particle_locations['Latitude']}

    particle_settings['release_locations'] = release_locations


    return create_particleset(fieldset, particle_settings)



def create_kernel(fieldset, pset):
    """_summary_

    Parameters
    ----------
    fieldset :
        A parcels.FieldSet object
    pset : _type_
        A parcels.ParticleSet object

    Returns
    -------
    kernels :
        A list of kernels used in the execution of the particle set
    """    
    kernels = []
    #kernels.append(PolyTEOS10_bsq)#pset.Kernel(PolyTEOS10_bsq)) # Set the seawater_density variable

    if fieldset.mode: # 3D mode = on
        kernels.append(AdvectionRK4_3D)#pset.Kernel(AdvectionRK4_3D))
    else:
        kernels.append(AdvectionRK4)#pset.Kernel(AdvectionRK4))
    

    if not fieldset.biofouling_f and fieldset.mode:
        kernels.append(settling_velocity)#pset.Kernel(settling_velocity))            
    elif fieldset.biofouling_f:
        kernels.append(biofouling)#pset.Kernel(biofouling))

    if fieldset.stokes_f:
        kernels.append(Stokes_drift)#pset.Kernel(Stokes_drift))
    if fieldset.wind_f:
        kernels.append(windage_drift)#pset.Kernel(windage_drift))


    if fieldset.stokes_f or fieldset.wind_f:
        kernels.append(unbeaching)


    if fieldset.mixing_f:
        kernels.append(vertical_mixing)#pset.Kernel(vertical_mixing))

    #kernels.append(checkThroughBathymetry)

    # Add statuscode kernels
    kernels.append(periodicBC)
    #if fieldset.mode:
    #    kernels.append(checkErrorThroughSurface)
    kernels.append(deleteParticle)

    return kernels

def load_default_settings():
    model_settings = {#'input_data_dir': '/storage/shared/oceanparcels/input_data/',
                    #'input_data_dir_2': '/storage/shared/oceanparcels/output_data/data_Michael/',
                    #'input_data_dir_3': '/storage/shared/oceanparcels/output_data/data_Mikael/',
                    'input_data_dir': './data/',
                    'input_data_dir_2': './data/',
                    'input_data_dir_3': './data/',
                    
                    'mode': '2D', # Options [3D, 2D]
                    
                    # Ocean model
                    'ocean_dir': '',                 # Directory of ocean model data
                    'ocean_filename': 'psy4v3r1-daily_',          # Filename style of ocean model data
                    'ocean_mesh': 'coordinates_med_subset.nc', # File location of ocean mesh
                    #'bathymetry_mesh': 'MOi/domain_ORCA0083-N006/bathymetry_ORCA12_V3.3.nc', # File location of bathymetry mesh
                    'mixing_f': False,                             # Turn on/off vertical turbulent mixing
                    
                    # Biogeochemistry model
                    'biofouling_f': False,                         # Turn on/off biofouling of particles
                    #'bgc_dir': '',                # Directory of biogeochemistry model
                    #'bgc_filename': 'biomer4v2r1-weekly_',
                    #'bgc_mesh': 'MOi/domain_ORCA025-N006/mesh_hgr_PSY4V3_deg.nc', # File location of biogeochemistry model mesh

                    #'bgc_biofilm_density': 1388.,                 # Biofilm density [g m-3]
                    #'bgc_algae_cell_volume': 2.0E-16,             # Volume of 1 algal cell [m3]
                    #'bgc_boltzmann_constant': 1.0306E-13 / (86400. ** 2.),  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)
                    #'bgc_respiration_rate': 0.1 / 86400.,         # Respiration rate, now [s-1]
                    #'bgc_respiration_temperature_coefficient': 2.13,# temperature coefficient respiration [-]
                    #'bgc_shear_frequency': 1.728E5 / 86400.,      # Shear [d-1], now [s-1]
                    #'bgc_carbon_atomic_weight': 12.,              # Atomic weight of Carbon
                    #'bgc_collision_probability': 1.,              # Collision probability [-]
                    #'bgc_algae_mortality_rate': 1.,               # TODO: Add description
                    #'bgc_algae_respiration_f': 1.,                # TODO: Add description


                    # Waves model
                    'stokes_f': True,                             # Turn on/off Stokes Drift
                    'stokes_dir': 'ERA5',                  # Directory of Stokes drift model data
                    'stokes_filename': 'ERA5_global_waves_monthly_',

                    # Wind model
                    'wind_f': True,                               # Turn on/off Windage
                    'wind_dir': 'ERA5',                     # Directory of Wind model data
                    'wind_filename': 'ERA5_global_wind_monthly_',

                    'allow_time_extrapolation': False,            # Allow extrapolation of time for fieldset
                    'verbose_delete': False,                      # Print extra information when executing delete operations
                    }

    particle_settings = {'start_date': datetime.strptime('2019-01-01-00:00:00', '%Y-%m-%d-%H:%M:%S'), # Start date of simulation
                        'runtime': timedelta(days=360),             # Runtime of simulation, use negative if releasing particles backwards in time
                        'dt_write': timedelta(days=1),             # Timestep of output
                        'dt_timestep': timedelta(minutes=20),       # Timestep of advection
                        # TODO: Could create own particle class with own sampling kernels to append after helper?
                        #'particle_class': PlasticParticle, # Particle class to use, feel free to create your own based on the base PlasticParticle class
                        }
    return model_settings, particle_settings

def load_test_settings():
    # Use the mediterranean as a test region for faster IO
    model_settings, particle_settings = load_default_settings()
    #model_settings['ocean_region'] = 'MED'
    
    return model_settings, particle_settings