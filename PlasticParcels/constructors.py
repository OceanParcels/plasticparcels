import os
import numpy as np
import xarray as xr

import pandas as pd
from parcels import FieldSet, Field, ParticleSet, JITParticle, Variable, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.converters import Geographic, GeographicPolar
from datetime import datetime, timedelta

from PlasticParcels.kernels import PolyTEOS10_bsq, StokesDrift, WindageDrift, SettlingVelocity, Biofouling, VerticalMixing, unbeaching, periodicBC, checkErrorThroughSurface, deleteParticle, checkThroughBathymetry
from PlasticParcels.utils import select_files, getclosest_ij



def create_hydrodynamic_fieldset(settings):
    """ A constructor method to create a Parcels.Fieldset from hydrodynamic model data

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
    bathymetry_variables = settings['ocean']['bathymetry_variables']
    bathymetry_dimensions = settings['ocean']['bathymetry_dimensions']
    bathymetry_mesh = os.path.join(settings['ocean']['directory'], settings['ocean']['bathymetry_mesh'])
    bathymetry_field = Field.from_netcdf(bathymetry_mesh, bathymetry_variables, bathymetry_dimensions)
    fieldset.add_field(bathymetry_field) 

    # If vertical mixing is turned on, add in the KPP-Profile
    if fieldset.mixing_f:
        dirread_model = os.path.join(settings['ocean']['directory'], settings['ocean']['filename_style'])
        kzfiles = select_files(dirread_model,'KZ_%4i*.nc',start_date,runtime,dt_margin=3)
        mixing_filenames = {'lon': ocean_mesh, 'lat': ocean_mesh, 'depth': wfiles[0], 'data': kzfiles}
        mixing_variables = settings['ocean']['vertical_mixing_variables']
        mixing_dimensions = settings['ocean']['vertical_mixing_dimensions']

        mixing_fieldset = FieldSet.from_nemo(mixing_filenames,mixing_variables,mixing_dimensions)
        fieldset.add_field(mixing_fieldset.mixing_kz)    #phytoplankton primary productivity        


    return fieldset

def create_fieldset(settings):
    """ A constructor method to create a Parcels.Fieldset with all fields necessary for a PlasticParcels simulation

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


    ## Apply unbeaching currents when Stokes/Wind can push particles into land cells
    if fieldset.stokes_f or fieldset.wind_f > 0:
        unbeachfiles = os.path.join(settings['unbeaching']['directory'], settings['unbeaching']['filename'])
        filenames_unbeach = {'unbeach_U': unbeachfiles, 
                            'unbeach_V': unbeachfiles}

        variables_unbeach = settings['unbeaching']['variables']

        dimensions_unbeach = settings['unbeaching']['dimensions']

        fieldset_unbeach = FieldSet.from_netcdf(filenames_unbeach, variables_unbeach, dimensions_unbeach, mesh='spherical')
        fieldset_unbeach.unbeach_U.units = GeographicPolar()
        fieldset_unbeach.unbeach_V.units = Geographic()

        fieldset.add_field(fieldset_unbeach.unbeach_U)
        fieldset.add_field(fieldset_unbeach.unbeach_V)
    

    fieldset.add_constant('verbose_delete',settings['verbose_delete'])

    return fieldset

def create_particleset_from_map(fieldset, settings):
    """ A constructor method to create a Parcels.ParicleSet for a PlasticParcels simulation from one of the available initialisation maps

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

     # Create a PlasticParticle class
    PlasticParticle = JITParticle
    variables = [
            Variable('plastic_diameter', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('plastic_density', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('wind_coefficient', dtype=np.float32, initial=0., to_write=False),
            Variable('settling_velocity', dtype=np.float64, initial=0., to_write=False),
            Variable('seawater_density', dtype=np.float32, initial=np.nan, to_write=False),
            Variable('absolute_salinity', dtype=np.float64, initial=np.nan, to_write=False),
            Variable('algae_amount', dtype=np.float64, initial=0., to_write=False),
            Variable('plastic_amount', dtype=np.float32, initial=0., to_write=True)
            ]
    for variable in variables:
        setattr(PlasticParticle, variable.name, variable)

    # # Create a PlasticParticle class
    # PlasticParticle = JITParticle.add_variables([
    #         Variable('plastic_diameter', dtype=np.float32, initial=np.nan, to_write=False),
    #         Variable('plastic_density', dtype=np.float32, initial=np.nan, to_write=False),
    #         Variable('wind_coefficient', dtype=np.float32, initial=0., to_write=False),
    #         Variable('settling_velocity', dtype=np.float64, initial=0., to_write=False),
    #         Variable('seawater_density', dtype=np.float32, initial=np.nan, to_write=False),
    #         Variable('absolute_salinity', dtype=np.float64, initial=np.nan, to_write=False),
    #         Variable('algae_amount', dtype=np.float64, initial=0., to_write=False),
    #         Variable('plastic_amount', dtype=np.float32, initial=0., to_write=True)
    #           ])

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

def create_kernel(fieldset):
    """ A constructor method to create a list of kernels for a PlasticParcels simulation

    Parameters
    ----------
    fieldset :
        A parcels.FieldSet object containing a range of constants to turn on/off different kernel behaviours

    Returns
    -------
    kernels :
        A list of kernels used in the execution of the particle set
    """    
    kernels = []

    kernels.append(PolyTEOS10_bsq) # To set the seawater_density variable

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

    ## Add the unbeaching kernel to the beginning
    if fieldset.stokes_f or fieldset.wind_f:
        kernels.append(unbeaching)

    if fieldset.mode:
        kernels.append(checkThroughBathymetry)
        kernels.append(checkErrorThroughSurface)

    # Add statuscode kernels
    kernels.append(periodicBC)
    kernels.append(deleteParticle)

    return kernels


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
