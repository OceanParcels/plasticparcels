from utils import select_files
import os
#import pandas as pd # Is this required?
from parcels import FieldSet, Field, ParticleSet, JITParticle, Variable, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.converters import Geographic, GeographicPolar 
import numpy as np

from plastic_kernels import PolyTEOS10_bsq, Stokes_drift, windage_drift, settling_velocity, biofouling, permanent_biofouling

## Fieldset creation
def create_fieldset(model_settings, particle_settings):
    
    ## Directories and arrays of files
    dirread_model = os.path.join(model_settings['input_data_dir'], model_settings['ocean_dir'], model_settings['ocean_filename'])

    start_date = particle_settings['start_date']
    runtime = int(np.ceil(particle_settings['runtime'].total_seconds()/86400.)) # convert to days

    ## If running backward simulations, we need to change the start date and runtime for loading the fieldset
    if np.sign(runtime) < 0:
        start_date = start_date + particle_settings['runtime']
        runtime = -runtime
    
    
    ## Mesh masks
    ocean_mesh = os.path.join(model_settings['input_data_dir'], model_settings['ocean_mesh'])      #mesh_mask

    ## Setup input for fieldset creation
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

    variables = {'U': 'vozocrtx',
                    'V': 'vomecrty',
                    'W': 'vovecrtz',
                    'conservative_temperature': 'votemper',
                    'absolute_salinity': 'vosaline'}

    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}, #time_centered
                    'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                    'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                    'conservative_temperature': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'},
                    'absolute_salinity': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'}}

    

    ## Load the fieldset
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=model_settings['allow_time_extrapolation'])

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
    

  
    ## Add fields for different model settings
    if fieldset.mixing_f:
        kzfiles = select_files(dirread_model,'KZ_%4i*.nc',start_date,runtime,dt_margin=3)
        
        bathymetry_mesh = os.path.join(model_settings['input_data_dir'], model_settings['bathymetry_mesh'])

        filenames['mixing_kz'] = {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': kzfiles}

        variables['mixing_kz'] = 'votkeavt'

        dimensions['mixing_kz'] = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

        fieldset.add_constant('z_start',0.5)
        
        variables = ('bathymetry', 'Bathymetry')
        dimensions = {'lon': 'nav_lon', 'lat': 'nav_lat'}
        bathymetry_field = Field.from_netcdf(bathymetry_mesh, variables, dimensions)
        fieldset.add_field(bathymetry_field) 


    if fieldset.biofouling_f: # or do_permanent_fouling:
        dirread_bgc = os.path.join(model_settings['input_data_dir'], model_settings['bgc_dir'], model_settings['bgc_filename'])
        bgc_mesh = os.path.join(model_settings['input_data_dir'], model_settings['bgc_mesh'])   #mesh_mask_4th

        ppfiles = select_files(dirread_bgc,'nppv_%4i*.nc',start_date,runtime,dt_margin=8)
        phy1files = select_files(dirread_bgc,'phy_%4i*.nc',start_date,runtime,dt_margin=8)
        phy2files = select_files(dirread_bgc,'phy2_%4i*.nc',start_date,runtime,dt_margin=8)

        #TODO: cleanup these names for better readability
        filenames_bio = {'pp_phyto': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': ppfiles},
                        'bio_nanophy': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy1files},
                        'bio_diatom': {'lon': bgc_mesh, 'lat': bgc_mesh, 'depth': wfiles[0], 'data': phy2files}}

        variables_bio = {'pp_phyto': 'nppv',
                        'bio_nanophy': 'phy',
                        'bio_diatom': 'phy2'}

        dimensions_bio = {'pp_phyto': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_nanophy': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                            'bio_diatom': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}
        
        #if not fieldset.mixing_f:
        #    raise RuntimeError('Mixing must be enabled for biofouling model to work') ## We can remove this as displacement from settling velocity is now applied in the kernel
        
        fieldset.add_constant('collision_probability', model_settings['bgc_collision_probability'])
        fieldset.add_constant('K', model_settings['bgc_boltzmann_constant'])
        fieldset.add_constant('R20', model_settings['bgc_respiration_rate'])

        fieldset.add_constant('Q10', model_settings['bgc_respiration_temperature_coefficient'])
        fieldset.add_constant('Gamma', model_settings['bgc_shear_frequency'])
        fieldset.add_constant('Wt_C', model_settings['bgc_carbon_atomic_weight'])
        fieldset.add_constant('algae_cell_volume', model_settings['bgc_algae_cell_volume'])
        fieldset.add_constant('biofilm_density', model_settings['bgc_biofilm_density'])    # density of biofilm [g m-3]

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

        dimensions_S = {'lat': 'lat',
                        'lon': 'lon',
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

        filenames_wind = {'wind_U': windfiles,
                        'wind_V': windfiles}

        variables_wind = {'wind_U': 'u10',
                        'wind_V': 'v10'}

        dimensions_wind = {'lat': 'lat',
                        'lon': 'lon',
                        'time': 'time'}
        fieldset_wind = FieldSet.from_netcdf(filenames_wind, variables_wind, dimensions_wind, mesh='spherical')
        fieldset_wind.wind_U.units = GeographicPolar()
        fieldset_wind.wind_V.units = Geographic() 
        fieldset_wind.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_wind.wind_U)
        fieldset.add_field(fieldset_wind.wind_V)


    ## Apply unbeaching currents when Stokes/Wind can push particles into land cells
    if fieldset.stokes_f or fieldset.wind_f > 0:
        ## TODO: UPDATE LOCATION OF DIRECTORIES ONCE DATA HAS BEEN MOVED - &&write a readme for data locations
        unbeachfiles = os.path.join(model_settings['input_data_dir_3'], '00_data_files/land_current_NEMO0083.nc')
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
    PlasticParticle = JITParticle
    PlasticParticle.particle_diameter = Variable('particle_diameter', dtype=np.float32, to_write=to_write_all)                   # Particle Diameter (assuming spherical particle) [meters] (l_pl)
    PlasticParticle.particle_density = Variable('particle_density', dtype=np.float32, to_write=to_write_all)                            # Particle Density [kg/m^3] (rho_pl)
    PlasticParticle.settling_velocity = Variable('settling_velocity', dtype=np.float64, initial=0., to_write=to_write_dynamic)              # Particle Sinking Velocity [m/s] (v_s)
    PlasticParticle.seawater_density = Variable('seawater_density', dtype=np.float32, to_write=to_write_dynamic)
    #PlasticParticle.conservative_temperature = Variable('conservative_temperature', dtype=np.float32, to_write=to_write_tracer)
    #PlasticParticle.absolute_salinity = Variable('absolute_salinity', dtype=np.float32, to_write=to_write_tracer)

    # Add kernel specific variables
    if fieldset.biofouling_f:
        PlasticParticle.algae_amount = Variable('absolute_salinity', dtype=np.float64, to_write=to_write_dynamic)
    
    if fieldset.wind_f:
        PlasticParticle.windage_coefficient = Variable('windage_coefficient', dtype=np.float32, initial=0., to_write=to_write_all)
    

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
        
def create_kernel(fieldset, pset):
    kernels = []
    kernels.append(pset.Kernel(PolyTEOS10_bsq)) # Set the seawater_density variable

    if fieldset.mode: # 3D mode = on
        kernels.append(pset.Kernel(AdvectionRK4_3D))
    else:
        kernels.append(pset.Kernel(AdvectionRK4))
    

    if fieldset.mixing_f:
        kernels.append(pset.Kernel(vertical_mixing))
        if not fieldset.biofouling_f:
            kernels.append(pset.Kernel(settling_velocity))
    if fieldset.biofouling_f:
        kernels.append(pset.Kernel(biofouling))
    if fieldset.stokes_f:
        kernels.append(pset.Kernel(Stokes_drift))
    if fieldset.wind_f:
        kernels.append(pset.Kernel(windage_drift))

    return kernels
