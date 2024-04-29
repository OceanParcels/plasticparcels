import os
import numpy as np
import xarray as xr


import pandas as pd
from parcels import FieldSet, Field, ParticleSet, JITParticle, Variable, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.converters import Geographic, GeographicPolar
from datetime import datetime, timedelta

from PlasticParcels.kernels import PolyTEOS10_bsq, StokesDrift, WindageDrift, SettlingVelocity, Biofouling, VerticalMixing, unbeaching, periodicBC, checkErrorThroughSurface, deleteParticle, checkThroughBathymetry
from PlasticParcels.utils import select_files, getclosest_ij


# ## Is this a better way?
# class PlasticParticle(JITParticle):
#     def __init__(self):
#         self.add_variables([
#             Variable('plastic_diameter', dtype=np.float32, initial=np.nan),
#             Variable('plastic_density', dtype=np.float32, initial=np.nan),
#             Variable('windage_coefficient', dtype=np.float32, initial=0.),
#             Variable('settling_velocity', dtype=np.float64, initial=0.),
#             Variable('seawater_density', dtype=np.float32, initial=np.nan),
#             Variable('absolute_salinity', dtype=np.float64, initial=np.nan),
#             Variable('algae_amount', dtype=np.float64, initial=0.)
#             ])


## Ocean model fieldset
def create_hydrodynamic_fieldset(settings):
    """ Constructor function to create a Parcels.Fieldset for hydrodynamic data

    Parameters
    ----------
    settings :
        A dictionary of settings used to create the fieldset

    Returns
    -------
    fieldset
        A parcels.FieldSet object
    """     

    # Location of hydrodynamic data
    dirread_model = os.path.join(settings['ocean']['directory'], settings['ocean']['filename_style'])

    # Start date and runtime of the simulation
    start_date = settings['simulation']['start_date']
    runtime = int(np.ceil(settings['simulation']['runtime'].total_seconds()/86400.)) # convert to days

    # Mesh masks
    ocean_mesh = os.path.join(settings['ocean']['directory'], settings['ocean']['ocean_mesh'])      #mesh_mask

    
    # Setup input for fieldset creation
    ufiles = select_files(dirread_model,'U_%4i*.nc',start_date,runtime,dt_margin=3)
    vfiles = select_files(dirread_model,'V_%4i*.nc',start_date,runtime,dt_margin=3)
    wfiles = select_files(dirread_model,'W_%4i*.nc',start_date,runtime,dt_margin=3)   
    tfiles = select_files(dirread_model,'T_%4i*.nc',start_date,runtime,dt_margin=3)    
    sfiles = select_files(dirread_model,'S_%4i*.nc',start_date,runtime,dt_margin=3)

    filenames = {'U': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': ufiles},
                    'V': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': vfiles},
                    'W': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': wfiles},
                    'conservative_temperature': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': tfiles},
                    'absolute_salinity': {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': sfiles}}

    variables = settings['ocean']['variables']
    dimensions = settings['ocean']['dimensions']
    indices = settings['ocean']['indices']

    if settings['mode'] == '3D':
        mode_3D_f = True
    else:
        mode_3D_f = False
        indices['depth'] = range(0,2)

    
    # Load the fieldset
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions,
                                  indices=indices, allow_time_extrapolation=settings['allow_time_extrapolation'])

    ## Create flags for custom particle behaviour
    fieldset.add_constant('mixing_f', settings['mixing_f'])
    fieldset.add_constant('biofouling_f', settings['biofouling_f'])
    fieldset.add_constant('stokes_f', settings['stokes_f'])
    fieldset.add_constant('wind_f', settings['wind_f'])
    fieldset.add_constant('G', 9.81) # Gravitational constant [m s-1]
    fieldset.add_constant('mode', mode_3D_f)
    
    # Add in bathymetry
    fieldset.add_constant('z_start',0.5)        
    bathymetry_variables = ('bathymetry', 'Bathymetry')
    bathymetry_dimensions = {'lon': 'nav_lon', 'lat': 'nav_lat'}
    bathymetry_mesh = os.path.join(settings['ocean']['directory'], settings['ocean']['bathymetry_mesh'])
    bathymetry_field = Field.from_netcdf(bathymetry_mesh, bathymetry_variables, bathymetry_dimensions)
    fieldset.add_field(bathymetry_field) 

    # If vertical mixing is turned on, add in the KPP-Profile
    if fieldset.mixing_f:
        dirread_model = os.path.join(settings['ocean']['directory'], settings['ocean']['filename_style'])
        kzfiles = select_files(dirread_model,'KZ_%4i*.nc',start_date,runtime,dt_margin=3)
        filenames_mixing = {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': kzfiles}
        variables_mixing = {'mixing_kz' : 'votkeavt'}
        dimensions_mixing = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

        mixing_fieldset = FieldSet.from_nemo(filenames_mixing,variables_mixing,dimensions_mixing)
        fieldset.add_field(mixing_fieldset.mixing_kz)    #phytoplankton primary productivity        


    return fieldset


## Fieldset creation
def create_fieldset(settings):
    """ Constructor function to create a Parcels.Fieldset with all fields necessary for a PlasticParcels simulation

    Parameters
    ----------
    settings :
        A dictionary of model settings used to create the fieldset
    
    Returns
    -------
    fieldset
        A parcels.FieldSet object
    """     
    
    # First create the hydrodynamic fieldset
    fieldset = create_hydrodynamic_fieldset(settings)

    # Now add the other fields
    # Start date and runtime of the simulation
    start_date = settings['simulation']['start_date']
    runtime = int(np.ceil(settings['simulation']['runtime'].total_seconds()/86400.)) # convert to days
           


    if fieldset.biofouling_f: # or do_permanent_fouling:
        #MOi glossary: https://www.mercator-ocean.eu/wp-content/uploads/2021/11/Glossary.pdf
        # and https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-028.pdf

        # Add BGC constants to current fieldset
        for key in settings['bgc']['constants']:
            fieldset.add_constant(key, settings['bgc']['constants'][key])

        # Create a fieldset with BGC data
        dirread_bgc = os.path.join(settings['bgc']['directory'], settings['bgc']['filename_style'])
        bgc_mesh = os.path.join(settings['bgc']['directory'], settings['bgc']['bgc_mesh'])   #mesh_mask_4th

        dirread_model = os.path.join(settings['ocean']['directory'], settings['ocean']['filename_style'])
        wfiles = select_files(dirread_model,'W_%4i*.nc',start_date,runtime,dt_margin=3)         

        ppfiles = select_files(dirread_bgc,'nppv_%4i*.nc',start_date,runtime,dt_margin=8)
        phy1files = select_files(dirread_bgc,'phy_%4i*.nc',start_date,runtime,dt_margin=8)
        phy2files = select_files(dirread_bgc,'phy2_%4i*.nc',start_date,runtime,dt_margin=8)

        filenames_bio = {'pp_phyto': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': ppfiles},
                        'bio_nanophy': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy1files},
                        'bio_diatom': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy2files}}

        variables_bio = settings['bgc']['variables']
        dimensions_bio = settings['bgc']['dimensions']
        
        # Create the BGC fieldset
        bio_fieldset = FieldSet.from_nemo(filenames_bio,variables_bio,dimensions_bio)

        ## Add the fields to the main fieldset
        fieldset.add_field(bio_fieldset.pp_phyto)    #phytoplankton primary productivity 
        fieldset.add_field(bio_fieldset.bio_nanophy) #nanopyhtoplankton concentration [mmol C m-3]
        fieldset.add_field(bio_fieldset.bio_diatom)  #diatom concentration [mmol C m-3]   


    if fieldset.stokes_f:
        dirread_Stokes = os.path.join(settings['stokes']['directory'], settings['stokes']['filename_style'])
        wavesfiles = select_files(dirread_Stokes,'%4i*.nc',start_date,runtime,dt_margin=32)

        filenames_Stokes = {'Stokes_U': wavesfiles,
                        'Stokes_V': wavesfiles,
                        'wave_Tp': wavesfiles}

        variables_Stokes = settings['stokes']['variables']
        dimensions_Stokes = settings['stokes']['dimensions']
        
        fieldset_Stokes = FieldSet.from_netcdf(filenames_Stokes, variables_Stokes, dimensions_Stokes, mesh='spherical')
        fieldset_Stokes.Stokes_U.units = GeographicPolar()
        fieldset_Stokes.Stokes_V.units = Geographic()
        fieldset_Stokes.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_Stokes.Stokes_U)
        fieldset.add_field(fieldset_Stokes.Stokes_V)
        fieldset.add_field(fieldset_Stokes.wave_Tp)

    if fieldset.wind_f:
        dirread_wind = os.path.join(settings['wind']['directory'], settings['wind']['filename_style'])
        windfiles = select_files(dirread_wind,'%4i*.nc',start_date,runtime,dt_margin=32)

        filenames_wind = {'Wind_U': windfiles,
                        'Wind_V': windfiles}

        variables_wind = settings['wind']['variables']
        dimensions_wind = settings['wind']['dimensions']

        fieldset_wind = FieldSet.from_netcdf(filenames_wind, variables_wind, dimensions_wind, mesh='spherical')
        fieldset_wind.Wind_U.units = GeographicPolar()
        fieldset_wind.Wind_V.units = Geographic() 
        fieldset_wind.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_wind.Wind_U)
        fieldset.add_field(fieldset_wind.Wind_V)


    # ## Apply unbeaching currents when Stokes/Wind can push particles into land cells
    # if fieldset.stokes_f or fieldset.wind_f > 0:
    #     unbeachfiles = os.path.join(model_settings['input_data_dir_3'], '../../data/output_data/masks/land_current_NEMO0083.nc')
    #     filenames_unbeach = {'unbeach_U': unbeachfiles, 
    #                         'unbeach_V': unbeachfiles}

    #     variables_unbeach = {'unbeach_U': 'land_current_u',
    #                         'unbeach_V': 'land_current_v'}   

    #     dimensions_unbeach = {'lat': 'lat',
    #                         'lon': 'lon'}
        

    #     fieldset_unbeach = FieldSet.from_netcdf(filenames_unbeach, variables_unbeach, dimensions_unbeach, mesh='spherical')
    #     fieldset_unbeach.unbeach_U.units = GeographicPolar()
    #     fieldset_unbeach.unbeach_V.units = Geographic()

    #     fieldset.add_field(fieldset_unbeach.unbeach_U)
    #     fieldset.add_field(fieldset_unbeach.unbeach_V)
    

    fieldset.add_constant('verbose_delete',settings['verbose_delete'])

    return fieldset





def create_particleset(fieldset, particle_settings):
    """ Helper function to create a Parcels.ParticleSet

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


def create_particleset_from_map(fieldset, settings):
    """ Helper function to create a Parcels.ParicleSet from one of the available initialisation maps

    Parameters
    ----------
    fieldset :
        A Parcels.FieldSet object
    settings :
        A dictionary of model settings, simulation settings, and particle release settings
    
    Returns
    -------
    particleset
        A parcels.ParticleSet object
    """

    # Load release type information
    release_type = settings['release']['initialisation_type'] ## YOU ARE HERE, TRING TO MAKE THIS PART BETTER!

    release_quantity_names = {
        'coastal':'MPW_Cell',
        'rivers':'Emissions',
        'fisheries':'fishing_hours',
        'global_concentrations':None #Not implemented yet
    }
    release_quantity_name = release_quantity_names[release_type]

    particle_locations = pd.read_csv(settings['release_maps'][release_type])


    # Select specific continent/region/subregion/country/economic status if applicable:
    if 'continent' in settings['release'].keys():
        particle_locations = particle_locations[particle_locations['Continent'] == settings['release']['continent']]
    if 'region' in settings['release'].keys():
        particle_locations = particle_locations[particle_locations['Region'] == settings['release']['region']]
    if 'subregion' in settings['release'].keys():
        particle_locations = particle_locations[particle_locations['Subregion'] == settings['release']['subregion']]
    if 'country' in settings['release'].keys():
        particle_locations = particle_locations[particle_locations['Country'] == settings['release']['country']]
    if 'economicstatus' in settings['release'].keys():
        particle_locations = particle_locations[particle_locations['Economic status'] == settings['release']['economicstatus']]

    particle_locations = particle_locations.groupby(['Longitude', 'Latitude'])[release_quantity_name].agg('sum').reset_index()
    particle_locations = particle_locations[particle_locations[release_quantity_name]>0]

    release_locations = {'lons': particle_locations['Longitude'],
                        'lats':  particle_locations['Latitude'],
                        'plastic_amount': particle_locations[release_quantity_name]}

    ## Set the longitude, latitude, and plastic amount per particle
    lons = release_locations['lons']
    lats = release_locations['lats']
    plastic_amounts = release_locations['plastic_amount']

    # 20240429 - for simplification I have removed the option to set custom times and depths, as well as different densities, diameters, and wind_coefficients, this can be done in an alternate way
    ## TO DO:
    ## Update this to use model grid id's instead of T-points lat/lon
    ## Set particle properties
    plastic_densities = np.full(lons.shape, settings['release']['plastic_density'])
    plastic_diameters = np.full(lons.shape, settings['release']['plastic_diameter'])
    wind_coefficients = np.full(lons.shape, settings['release']['wind_coefficient'])

    # Create a plasticparticle class
    PlasticParticle = JITParticle.add_variables([
            Variable('plastic_diameter', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('plastic_density', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('wind_coefficient', dtype=np.float32, initial=0., to_write=False),
            Variable('settling_velocity', dtype=np.float64, initial=0., to_write=False),
            Variable('seawater_density', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('absolute_salinity', dtype=np.float64, initial=np.nan, to_write=False),
            Variable('algae_amount', dtype=np.float64, initial=0., to_write=False),
            Variable('plastic_amount', dtype=np.float32, initial=0., to_write=True)
            ])

    pset = ParticleSet.from_list(fieldset,
                                    PlasticParticle,
                                    lon=lons,
                                    lat=lats,
                                    plastic_diameter=plastic_diameters,
                                    plastic_density=plastic_densities,
                                    wind_coefficient=wind_coefficients,
                                    plastic_amount=plastic_amounts
                                    )


    return pset

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

    kernels.append(PolyTEOS10_bsq)#pset.Kernel(PolyTEOS10_bsq)) # Set the seawater_density variable

    if fieldset.mode: # 3D mode = on
        kernels.append(AdvectionRK4_3D)#pset.Kernel(AdvectionRK4_3D))
    else:
        kernels.append(AdvectionRK4)#pset.Kernel(AdvectionRK4))
    

    if not fieldset.biofouling_f and fieldset.mode:
        kernels.append(SettlingVelocity)#pset.Kernel(settling_velocity))            
    elif fieldset.biofouling_f and fieldset.mode: # Must be in 3D to use biofouling mode
        kernels.append(Biofouling)#pset.Kernel(biofouling))

    if fieldset.stokes_f:
        kernels.append(StokesDrift)#pset.Kernel(Stokes_drift))
    if fieldset.wind_f:
        kernels.append(WindageDrift)#pset.Kernel(windage_drift))

    if fieldset.mixing_f:
        kernels.append(VerticalMixing)#pset.Kernel(vertical_mixing))

    # ## Add the unbeaching kernel to the beginning
    # if fieldset.stokes_f or fieldset.wind_f:
    #     kernels.append(unbeaching)

    if fieldset.mode:
        kernels.append(checkThroughBathymetry)
        kernels.append(checkErrorThroughSurface)

    # Add statuscode kernels
    kernels.append(periodicBC)
    kernels.append(deleteParticle)

    return kernels

def load_default_settings():
    settings = {
                    # Model settings
                    'mode': '2D', # Options [3D, 2D]
                    'allow_time_extrapolation': False,            # Allow extrapolation of time for fieldset
                    'verbose_delete': False,                      # Print extra information when executing delete operations
                    
                    # Flags
                    'mixing_f': False,                             # Turn on/off vertical turbulent mixing
                    'biofouling_f': False,                         # Turn on/off biofouling of particles
                    'stokes_f': False,                             # Turn on/off Stokes Drift
                    'wind_f': False,                               # Turn on/off Windage

                    # Ocean model
                    'ocean' : {
                        'directory': 'data/input_data/MOi/',                 # Directory of ocean model data
                        'filename_style': 'psy4v3r1/psy4v3r1-daily_',          # Filename style of ocean model data
                        'ocean_mesh': 'domain_ORCA0083-N006/coordinates.nc', # File location of ocean mesh
                        'bathymetry_mesh': 'domain_ORCA0083-N006/bathymetry_ORCA12_V3.3.nc', # File location of bathymetry mesh
                        'variables': { #Variable names in the ocean model for velocities, temperature and salinity
                            'U': 'vozocrtx', # Variable name for the U-velocity
                            'V': 'vomecrty', # Variable name for the V-velocity
                            'W': 'vovecrtz', # Variable name for the W-velocity
                            'conservative_temperature': 'votemper', # Variable name for the temperature field
                            'absolute_salinity': 'vosaline' # Variable name for the salinity field
                            }, 
                        'dimensions':{ # The dimensions of the ocean model data, providing f-points for C-grid models
                            'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'conservative_temperature': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'},
                            'absolute_salinity': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'}
                            },
                        'indices': {}
                    },
                    
                    # Biogeochemistry model
                    'bgc' : {
                        'directory': 'data/input_data/MOi/',                # Directory of biogeochemistry model
                        'filename_style': 'biomer4v2r1/biomer4v2r1-weekly_',
                        'bgc_mesh': 'domain_ORCA025-N006/mesh_hgr_PSY4V3_deg.nc', # File location of biogeochemistry model mesh

                        'variables' : {
                             'pp_phyto': 'nppv', # Total Primary Production of Phyto - 'Net primary prodution of biomass expressed as arbon per unit volume in sea water' [mg m-3 day-1] or [milligrams of Carbon per cubic meter per day]
                            'bio_nanophy': 'phy', # Mole concentration of NanoPhytoplankton expressed as carbon in sea water
                            'bio_diatom': 'phy2' # Mole concentration of Diatoms expressed as carbon in sea water
                            },
                        'dimensions' : {
                            'pp_phyto': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_nanophy': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_diatom': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
                            },

                        'constants':{
                            'biofilm_density': 1388.,                 # Biofilm density [g m-3]
                            'algae_cell_volume': 2.0E-16,             # Volume of 1 algal cell [m3]
                            'K': 1.0306E-13 / (86400. ** 2.),  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)
                            'R20': 0.1 / 86400.,         # Respiration rate, now [s-1]
                            'Q10': 2.13,# temperature coefficient respiration [-]
                            'Gamma': 1.728E5 / 86400.,      # Shear frequency [d-1], now [s-1]
                            'carbon_molecular_weight': 12.,              # Atomic weight of Carbon
                            'collision_probability': 1.,              # Collision probability [-]
                            'algae_mortality_rate': 1.,               # TODO: Add description
                            'algae_respiration_f': 1.,                # TODO: Add description
                        }
                    },

                    # Waves model
                    'stokes' : {
                        'directory': 'data/input_data/ERA5/waves/',                  # Directory of Stokes drift model data
                        'filename_style': 'ERA5_global_waves_monthly_',
                        'variables' : {
                            'Stokes_U': 'ust',
                            'Stokes_V': 'vst',
                            'wave_Tp': 'pp1d'
                            },
                        'dimensions': {
                            'lat': 'latitude', # when not using the converted datasets, use 'longitude' otherwise use 'lon'
                            'lon': 'longitude',
                            'time': 'time'
                            }
                    },

                    # Wind model
                    'wind': {
                        'directory': 'data/input_data/ERA5/wind/',                     # Directory of Wind model data
                        'filename_style': 'ERA5_global_wind_monthly_',
                        'variables' : {
                            'Wind_U': 'u10',
                            'Wind_V': 'v10'
                            },
                        'dimensions' : {
                            'lat': 'latitude', # when not using the converted datasets, use 'longitude' otherwise use 'lon'
                            'lon': 'longitude',
                            'time': 'time'
                            }
                    },

                    'simulation': {
                        'start_date': None,
                        'runtime': None,
                        'dt_write': None,
                        'dt_timestep': None
                    },

                    'release_maps': {
                        'coastal':'data/release/generated_files/coastal_population_MPW_NEMO0083.csv',
                        'rivers':'data/release/generated_files/river_emissions_NEMO0083.csv',
                        'fisheries':'data/release/generated_files/agg_data_fisheries_info.csv',
                        'global_concentrations':None #Not implemented yet
                    }

                    }

    return settings
