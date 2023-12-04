## A list of kernels used to displace particles

#TODO: Add documentation here

## @author: denes001
## @date: 2023-08-09


from parcels import Variable, JITParticle, ParcelsRandom
import math
import numpy as np

## ERA 5 https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

def Stokes_drift(particle, fieldset, time):
    """Stokes drift kernel
    
    Description
    ----------
    Using the approach in [1] assuming a Phillips wave spectrum to determine the depth dependent Stokes drift.
    Specifically, the 'Stokes drift velocity' (u_s) is computed as per Eq. (19) in [1]

    We treat the Stokes drift as a linear addition to the velocity field
        u(x,t) = u_c(x,t) + C_s * u_s(x,t)
    where u_c is the current velocity, u_s is the Stokes drift velocity, and C_s is the depth-varying decay factor.

    Parameter Requirements
    ----------
    fieldset :
        - Stokes drift field (U,V)  [m s-1]
        - Peak wave period (T_p)    [s]

    Kernel Requirements
    ----------
    Order of Operations:
        None
        
    References
    ----------
    [1] Breivik (2016) - https://doi.org/10.1016/j.ocemod.2016.01.005


    """

    # Sample the U / V components of Stokes drift
    stokes_U = fieldset.Stokes_U[time, particle.depth, particle.lat, particle.lon]
    stokes_V = fieldset.Stokes_V[time, particle.depth, particle.lat, particle.lon]
    
    # Sample the peak wave period
    T_p = fieldset.wave_Tp[time, particle.depth, particle.lat, particle.lon]

    # Only compute displacements if the peak wave period is large enough
    if T_p > 1E-14:
        # Peak wave frequency 
        omega_p = 2. * math.pi / T_p
        
        # Peak wave number
        k_p = (omega_p ** 2) / fieldset.G

        # Repeated inner term of Eq. (19) - note depth is negative in this formulation, but model depths are positive by convention
        kp_z_2 = 2. * k_p * particle.depth

        # Decay factor in Eq. (19) -- Where beta=1 for the Phillips spectrum
        decay = math.exp(-kp_z_2) - math.sqrt(math.pi * kp_z_2) * math.erfc(math.sqrt(kp_z_2))

        # Apply Eq. (19) and compute particle displacement
        particle_dlon += stokes_U * decay * particle.dt
        particle_dlat += stokes_V * decay * particle.dt        
    


### Wind related kernels ###
def windage_drift(particle, fieldset, time):
    """Leeway windage kernel
    
    Description
    ----------
    A simple windage kernel that applies a linear 'wind velocity' to the particle.
    ?????TODO: the best I can find so far is eq 1 from Impact of windage on ocean surface Lagrangian coherent structures but follow the references more
    ????? Mayber Breivik et al 2011 - Wind-induced drift of objects at sea: The leeway field method.

    We treat the windage drift as a linear addition to the velocity field
        u(x,t) = u_c(x,t) + C_w * u_w(x,t)
    where u_c is the current velocity, u_w is the wind velocity at 10m height, and C_w is the windage coefficient (usually taken to be between 0.01 and 0.05,
    depending on particle size)

    See https://doi.org/10.1007/s10652-016-9499-3 for impact of windage on LCS detection
    or https://doi.org/10.1016/j.apor.2011.01.005 for Breivik et al 2011

    Parameter Requirements
    ----------
    particle :
        - windage_coefficient - the particle windage coefficient (usually taken to be between 0.01 and 0.05, depending on particle size)
    fieldset :
        - Wind field (U,V) at 10m height above sea surface [m s-1]
    

    Kernel Requirements
    ----------
    Order of Operations:
        None          
        
    """

    # Sample the U / V components of wind
    wind_U = fieldset.wind_U[time, particle.depth, particle.lat, particle.lon]
    wind_V = fieldset.wind_V[time, particle.depth, particle.lat, particle.lon]

    # Apply windage to particles that have some exposed surface above the ocean surface
    if particle.depth < 0.5*particle.particle_diameter:
        # Compute particle displacement
        particle_dlon += particle.windage_coefficient * wind_U * particle.dt
        particle_dlat += particle.windage_coefficient * wind_V * particle.dt


## Decay is hard to justify. Use on/off kernel as above, and maybe create a windage function that is depth/diameter dependent?
# def windage_drift_decay(particle, fieldset, time):
#     """
#     Windage kernel with a sqrt decay profile
    
#     Description:
#         using
#         ?????TO DO the best I can find so far is eq 1 from Impact of windage on ocean surface Lagrangian coherent structures but follow the references more
 
#         We treat the windage drift as a linear addition to the velocity field
#             u(x,t) = u_c(x,t) + D * C_w * u_w(x,t)
#         where u_c is the current velocity, u_w is the wind velocity at 10m height, and C_w is the windage coefficient (usually taken to be between 0.01 and 0.05,
#         depending on particle size) and D is a depth-varying decay factor.

#         See https://doi.org/10.1007/s10652-016-9499-3 for impact of windage on LCS detection

#     Requirements:
#         Fieldset:
#             - Wind field (U,V) at 10m height above sea surface
#         Particle:
#             - windage_coefficient - the particle windage coefficient (usually taken to be between 0.01 and 0.05, depending on particle size)
#         Order of Operations:
#             None


#     ## TO DO: is a decay necessary? we should really just apply this at the surface
#     """

#     # U / V components of wind
#     U_wind = fieldset.wind_U[time, particle.depth, particle.lat, particle.lon]
#     V_wind = fieldset.wind_V[time, particle.depth, particle.lat, particle.lon]

#     # Where to apply decay rate
#     decay_depth = 0.2 # meters, point at which decay is 0%
    
#     # Apply to top layer only using a linear decay between 0 and 0.2m depth
#     if particle.depth > decay_depth or particle.depth < 0: # If the particle depth is above the surface, set to 0
#         pass # No change to particle position

#     else: # The particle depth is in [0, decay_depth], compute a sqrt decay
#         # Compute decay factor where decay = m * sqrt(depth) + b, passing through points (depth,decay) = (0,1) and (decay_depth,0)
#         decay = 1. - (1. / math.sqrt(decay_depth)) * math.sqrt(particle.depth)
      
#         # Compute particle displacement
#         dlon = decay * particle.windage_coefficient * U_wind * particle.dt
#         dlat = decay * particle.windage_coefficient * V_wind * particle.dt

#         # Update particle position
#         particle.lon += dlon
#         particle.lat += dlat




### Vertical diffusion related kernels ###

## TODO: Is this required? Why not just set v_s = 0 from the beginning??
def neutral_buoyancy(particle, fieldset, time):
    """Neutral buoyancy kernel

    Description
    ----------
    A kernel to treat the particle as neutrally buoyant by setting the sinking/rising velocity of a particle to 0 [m s-1]

    Parameter Requirements
    ----------
    particle :
        sinking_velocity

    Kernel Requirements
    ----------
    Order of Operations:
        Apply this kernel before any other vertical mixing kernels
    """


    # Set the particle's settling velocity to 0 [m s-1]
    particle.settling_velocity = 0.





def settling_velocity(particle, fieldset, time):
    """Settling velocity kernel
    
    Description
    ----------
    A kernel to calculate the settling velocity of a particle based on its size and density, 
    and the surrounding seawater properties (density and kinematic viscosity). Based on the
    settling velocity paper of [1] as implemented in [2].
    
    Calculation steps:
        1. Compute the seawater dynamic viscosity from Eq. (27) in [2]
        2. Compute the kinematic viscosity from Eq. (25) in [2]
        3. Compute the dimensionless particle diameter from Eq. (4) in [2], equivalent to Eq. (6) in [1]
        4. Compute the dimensionless settling velocity from Eq. (3) in [2], equivalent to Eq. (8) in [1]
        5. Compute the settling velocity of the particle from Eq. (2) in [2], equivalent to Eq. (9) in [1]


    Parameter Requirements
    ----------
    particle :
        - particle_diameter
        - particle_density
        - seawater_density Surrounding seawater density
    fieldset :
        - G Gravity constant [m s-2]
        - cons_temperature Conservative temperature field
        - abs_salinity Absolute salinity field
    
        
    Kernel Requirements:
    ----------
    Order of Operations:
        - Must run after the PolyTEOS kernel which sets the surrounding particle.seawater_density variable which this kernel relies on
        - Must run before any kernel that uses the particle.settling_velocty variable

        
    References:
    ----------
    [1] Dietrich (1982) - https://doi.org/10.1029/WR018i006p01615
    [2] Kooi et al. (2017) - https://doi.org/10.1021/acs.est.6b04702

    TODO: Add units to each variable, and the TODO's below.
    """

    # Define constants and sample fieldset variables
    g = fieldset.G  # gravitational acceleration [m s-2]
    seawater_density = particle.seawater_density  # [kg m-3]
    temperature = fieldset.conservative_temperature[time, particle.depth, particle.lat, particle.lon]
    seawater_salinity = fieldset.absolute_salinity[time, particle.depth, particle.lat, particle.lon] / 1000 # TODO: Check that we should convert here and not in the fieldset
    particle_diameter = particle.particle_diameter
    particle_density = particle.particle_density

    # Compute the kinematic viscosity of seawater 
    water_dynamic_viscosity = 4.2844E-5 + (1 / ((0.157 * (temperature + 64.993) ** 2) - 91.296)) # Eq. (26) from [2] #TODO: check this is correct, constants aren't the same... - 0.156?
    A = 1.541 + 1.998E-2 * temperature - 9.52E-5 * temperature ** 2 # Eq. (28) from [2]
    B = 7.974 - 7.561E-2 * temperature + 4.724E-4 * temperature ** 2 # Eq. (29) from [2]
    seawater_dynamic_viscosity = water_dynamic_viscosity * (1 + A * seawater_salinity + B * seawater_salinity ** 2) # Eq. (27) from [2]
    seawater_kinematic_viscosity = seawater_dynamic_viscosity / seawater_density # Eq. (25) from [2]

    # Compute the density difference of the particle
    normalised_density_difference = (particle_density - seawater_density) / seawater_density  # normalised difference in density between the plastic particle and seawater [-]
    
    # Compute the dimensionless particle diameter D_* using Eq. (4) from [2]
    dimensionless_diameter = (math.fabs(particle_density - seawater_density) * g * particle_diameter ** 3.) / (seawater_density * seawater_kinematic_viscosity ** 2.)  # [-]

    # Compute the dimensionless settling velocity w_*
    if dimensionless_diameter > 5E9: # "The boundary layer around the sphere becomes fully turbulent, causing a reduction in drag and an increase in settling velocity" - [1]
        dimensionless_velocity = 265000. # Set a maximum dimensionless settling velocity
    elif dimensionless_diameter < 0.05: # "At values of D_* less than 0.05, (9) deviates signficantly ... from Stokes' law and (8) should be used." - [1]
        dimensionless_velocity = (dimensionless_diameter ** 2.) / 5832 # Using Eq. (8) in [1]
    else:
        dimensionless_velocity = 10. ** (-3.76715 + (1.92944 * math.log10(dimensionless_diameter)) - (0.09815 * math.log10(dimensionless_diameter) ** 2.) - (
                    0.00575 * math.log10(dimensionless_diameter) ** 3.) + (0.00056 * math.log10(dimensionless_diameter) ** 4.)) # Using Eq. (9) in [1]

    # Compute the settling velocity of the particle using Eq. (5) from [1] (solving for the settling velocity)
    sign_of_density_difference = math.copysign(1., normalised_density_difference) #normalised_density_difference/math.fabs(normalised_density_difference) ## maybe use math.copysign(1,normalised_density_difference)? does copysign work in cgen?
    settling_velocity = sign_of_density_difference * (g * seawater_kinematic_viscosity * dimensionless_velocity * math.fabs(normalised_density_difference)) ** (1. / 3.)  # m s-1

    # Update the settling velocity
    particle.settling_velocity = settling_velocity   

    # Update particle depth
    particle_ddepth += particle.settling_velocity * particle.dt







def biofouling(particle, fieldset, time):
    """ Biofouling kernel

    Description
    ----------
    Kernel to compute the vertical velocity of particles due to changes in ambient algal concentrations,
    growth and death of attached algae based on Kooi et al. 2017

    model settling velocity and MEDUSA 2.0 biofilm dynamics, including modelling of the 3D mesozooplankton grazing of diatoms
    
    A kernel to calculate the settling velocity of a particle based on its size and density, and the surrounding seawater properties (density and kinematic viscosity)

    TODO: What are the assumptions here? ALSO - why would a biofouled particle use this kernel? they should use the biofouling kernel no? can remove the whole biofouling stuff from this!

    Based on the settling velocity paper of [1] Dietrich (1982) - https://doi.org/10.1029/WR018i006p01615 as implemented in [2] Kooi (2017) - https://doi.org/10.1021/acs.est.6b04702
    
    Calculation steps:
        1. Compute the seawater dynamic viscosity from Eq. (27) in [2]
        2. Compute the kinematic viscosity from Eq. (25) in [2]
        3. Compute the dimensionless particle diameter from Eq. (4) in [2], equivalent to Eq. (6) in [1]
        4. Compute the dimensionless settling velocity from Eq. (3) in [2], equivalent to Eq. (8) in [1]
        5. Compute the settling velocity of the particle from Eq. (2) in [2], equivalent to Eq. (9) in [1]




    Parameter Requirements
    ----------
    particle :
        - Particle diameter 'particle_diameter'
        - Particle density 'particle_density'
        - Surrounding seawater density 'seawater_density'
    fieldset : _type_
        - Gravity constant 'G' [m s-2]
        - Conservative temperature field 'cons_temperature'
        - Absolute salinity field 'abs_salinity'
        - Algae cell volume 'algae_cell_volume'
        - Biofilm density 'biofilm_density'


    Kernel Requirements
    ----------
        Order of Operations:

    TODO: Add units to each variable
    Hard to generalise this without something like fieldset.biofouling == True? to determine if we need to compute the biofilm component
     ----- WELL we don't since anything with biofouling would use the biofouling kernel.... 



    References:
    ----------
    [1] Dietrich (1982) - https://doi.org/10.1029/WR018i006p01615
    [2] Kooi et al. (2017) - https://doi.org/10.1021/acs.est.6b04702
    [3] Menden-Deuer, Lessard (2000) - https://doi.org/10.4319/lo.2000.45.3.0569
    """




    # Define constants and algal properties, and sample fieldset variables
    g = fieldset.G  # gravitational acceleration [m s-2]
    k = fieldset.K  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)

    algae_cell_volume = fieldset.algae_cell_volume  # Volume of 1 algal cell [m-3]    - V_a
    algae_amount = particle.algae_amount  # Number of attached algae per square meter of particle surface [number m-2]
    biofilm_density = fieldset.biofilm_density
    r20 = fieldset.R20  # respiration rate [s-1]
    q10 = fieldset.Q10  # temperature coefficient respiration [-]
    gamma = fieldset.Gamma  # shear rate [s-1]

    seawater_density = particle.seawater_density  # [kg m-3]
    temperature = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    seawater_salinity = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000 # TODO: Check that we should convert here and not in the fieldset
    particle_radius = 0.5 * particle.particle_diameter
    particle_density = particle.particle_density


    # This needs to be computed as in the settling velocity kernel? Or, at least the initial settling velocity must be set?
    initial_settling_velocity = particle.settling_velocity  # settling velocity [m s-1]

    # ------ sample fields ------
    water_dynamic_viscosity = 4.2844E-5 + (1 / ((0.157 * (temperature + 64.993) ** 2) - 91.296)) # Eq. (26) from [2] #TODO: check this is correct, constants aren't the same... - 0.156?
    A = 1.541 + 1.998E-2 * temperature - 9.52E-5 * temperature ** 2 # Eq. (28) from [2]
    B = 7.974 - 7.561E-2 * temperature + 4.724E-4 * temperature ** 2 # Eq. (29) from [2]
    seawater_dynamic_viscosity = water_dynamic_viscosity * (1 + A * seawater_salinity + B * seawater_salinity ** 2) # Eq. (27) from [2]
    seawater_kinematic_viscosity = seawater_dynamic_viscosity / seawater_density # Eq. (25) from [2]


    median_mg_carbon_per_cell = 2726e-9 # Median mg of Carbon per cell from [3]. From [2] pg7966 - "The conversion from carbon to algae cells is highly variable, ranging between 35339 to 47.76 pg per carbon cell [3]. We choose the median value, 2726 x 10^-9 mg per carbon cell."
    carbon_molecular_weight= fieldset.carbon_molecular_weight # grams C per mol of C #Wt_C
    
    
    ###---------------Growth--------------------###
    #Mole concentration of Diatoms expressed as carbon in sea water
    mol_concentration_diatoms = fieldset.bio_diatom[time, particle.depth, particle.lat, particle.lon]

    #Mole concentration of Nanophytoplankton expressed as carbon in seawater
    mol_concentration_nanophytoplankton = fieldset.bio_nanophy[time, particle.depth, particle.lat, particle.lon] #nanophytoplankton
    
    number_concentration_diatoms = mol_concentration_diatoms * (carbon_molecular_weight / median_mg_carbon_per_cell) # conversion from [mmol C m-3] to [mg C m-3] to [no. m-3]
    number_concentration_nanophytoplankton = mol_concentration_nanophytoplankton * (carbon_molecular_weight / median_mg_carbon_per_cell)
    number_concentration_total = number_concentration_diatoms + number_concentration_nanophytoplankton


    #I don't like this bit below, it does not make sense...
    if number_concentration_total < 0:
        number_concentration_total = 0.
        print('negative diat/non-diat. concentration')
    

    
    total_primary_production_of_phyto = fieldset.pp_phyto[time, particle.depth, particle.lat, particle.lon] # mg C /m3/day #pp_phyto_ 
    
    primary_production_per_cell = total_primary_production_of_phyto / number_concentration_total # primary productivity per cell, in mg C / cell / day
    
    #pp_ncell_per_cell = pp_per_cell * (1 / median_mg_carbon_per_cell) #primary productivity in terms of amount of cells per cell, in cells / cell / day
    pp_ncell_per_cell = primary_production_per_cell / median_mg_carbon_per_cell #primary productivity in terms of amount of cells per cell, in cells / cell / day
    
    if pp_ncell_per_cell < 0:
        mu_a = 0
    elif pp_ncell_per_cell > 1.85:
        mu_a = 1.85 / 86400. # maximum growth rate
    else:
        mu_a = pp_ncell_per_cell / 86400. #d-1 to s-1
    
    a_growth = mu_a * particle.a #productivity in amount of cells/m2/s
    

    # Compute the volume of the particle and potential biofilm
    particle_volume = (4. / 3.) * math.pi * particle_radius ** 3.  # volume of plastic [m3]
    particle_surface_area = 4. * math.pi * particle_radius ** 2.  # surface area of plastic particle [m2]
    algal_cell_radius = ((3. / 4.) * (algae_cell_volume / math.pi)) ** (1. / 3.)  # radius of an algal cell [m]
    biofilm_volume = algae_cell_volume * algae_amount * particle_surface_area  # volume of living biofilm [m3]
    
    total_volume = biofilm_volume + particle_volume  # volume of total (biofilm + plastic) [m3]
    total_radius = ((total_volume * (3. / (4. * math.pi))) ** (1. / 3.)) # total radius [m]    
    biofilm_thickness = total_radius - particle_radius  # biofilm thickness [m]
  
    # ------ Diffusivity -----
    total_density = (particle_volume * particle_density + biofilm_volume * biofilm_density) / total_volume  # total density [kg m-3]
    # theta_tot = 4. * math.pi * r_tot ** 2.  # surface area of total [m2]
    d_pl = k * (temperature + 273.16) / (6. * math.pi * seawater_dynamic_viscosity * total_radius)  # diffusivity of plastic particle [m2 s-1]
    d_a = k * (temperature + 273.16) / (6. * math.pi * seawater_dynamic_viscosity * algal_cell_radius)  # diffusivity of algal cells [m2 s-1]

    # ------ Encounter rates -----
    beta_abrown = 4. * math.pi * (d_pl + d_a) * (total_radius + algal_cell_radius)  # Brownian motion [m3 s-1]
    beta_ashear = 1.3 * gamma * ((total_radius + algal_cell_radius) ** 3.)  # advective shear [m3 s-1]
    beta_aset = (1. / 2.) * math.pi * total_radius ** 2. * math.fabs(initial_settling_velocity)  # differential settling [m3 s-1]
    beta_a = beta_abrown + beta_ashear + beta_aset  # collision rate [m3 s-1]

    # ------ Attached algal growth (Eq. 11 in Kooi et al. 2017) -----
    a_coll = (beta_a * number_concentration_diatoms) / particle_surface_area * fieldset.collision_probability  # [no. m-2 s-1] collisions with diatoms


    ### ------------Respiration-------------- ###
    a_resp = (q10 ** ((temperature - 20.) / 10.)) * r20 * particle.a  # [no. m-2 s-1] respiration


    a_graze = 0
    
    #--------------------Total-----------------------------
    particle.a += (a_coll + a_growth - a_resp - a_graze) * particle.dt

    dn = 2. * (r_tot)  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
    delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    dstar = (math.fabs(rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * sw_kin_visc ** 2.)  # [-]

    if dstar > 5e9:
        w_star = 265000
    elif dstar < 0.05:
        w_star = (dstar ** 2.) * 1.71E-4
    else:
        w_star = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                    0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))

    # ------ Settling velocity of particle -----
    if delta_rho > 0:  # sinks
        vs_new = (g * sw_kin_visc * w_star * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs_new = -1. * (g * sw_kin_visc * w_star * a_del_rho) ** (1. / 3.)  # m s-1
        
    
    particle.v_s = vs_new    




def PolyTEOS10_bsq(particle, fieldset, time):
    """ A seawater density kernel

    Description:
    ----------
    A kernel to calculate the seawater density surrounding a particle based on the polyTEOS10-bsq algorithm
    from Appendix A.2 of https://doi.org/10.1016/j.ocemod.2015.04.002

    Parameter Requirements
    ----------
    particle :
        seawater_density
    fieldset : _type_
        absolute_salinity
        conservative_temperature

    Kernel Requirements:
    ----------
    Order of Operations:
        This kernel must be run before any kernel that requires particle.seawater_density

    
    
    References:
     Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate
      polynomial expressions for the density and specific volume of
      seawater using the TEOS-10 standard. Ocean Modelling.
     McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003:
      Accurate and computationally efficient algorithms for potential
      temperature and density of seawater.  Journal of Atmospheric and
      Oceanic Technology, 20, 730-741.
    """
    
    Z = - particle.depth  # note: use negative depths!
    SA = fieldset.absolute_salinity[time, particle.depth, particle.lat, particle.lon]
    CT = fieldset.conservative_temperature[time, particle.depth, particle.lat, particle.lon]

    SAu = 40 * 35.16504 / 35
    CTu = 40
    Zu = 1e4
    deltaS = 32
    R000 = 8.0189615746e+02
    R100 = 8.6672408165e+02
    R200 = -1.7864682637e+03
    R300 = 2.0375295546e+03
    R400 = -1.2849161071e+03
    R500 = 4.3227585684e+02
    R600 = -6.0579916612e+01
    R010 = 2.6010145068e+01
    R110 = -6.5281885265e+01
    R210 = 8.1770425108e+01
    R310 = -5.6888046321e+01
    R410 = 1.7681814114e+01
    R510 = -1.9193502195e+00
    R020 = -3.7074170417e+01
    R120 = 6.1548258127e+01
    R220 = -6.0362551501e+01
    R320 = 2.9130021253e+01
    R420 = -5.4723692739e+00
    R030 = 2.1661789529e+01
    R130 = -3.3449108469e+01
    R230 = 1.9717078466e+01
    R330 = -3.1742946532e+00
    R040 = -8.3627885467e+00
    R140 = 1.1311538584e+01
    R240 = -5.3563304045e+00
    R050 = 5.4048723791e-01
    R150 = 4.8169980163e-01
    R060 = -1.9083568888e-01
    R001 = 1.9681925209e+01
    R101 = -4.2549998214e+01
    R201 = 5.0774768218e+01
    R301 = -3.0938076334e+01
    R401 = 6.6051753097e+00
    R011 = -1.3336301113e+01
    R111 = -4.4870114575e+00
    R211 = 5.0042598061e+00
    R311 = -6.5399043664e-01
    R021 = 6.7080479603e+00
    R121 = 3.5063081279e+00
    R221 = -1.8795372996e+00
    R031 = -2.4649669534e+00
    R131 = -5.5077101279e-01
    R041 = 5.5927935970e-01
    R002 = 2.0660924175e+00
    R102 = -4.9527603989e+00
    R202 = 2.5019633244e+00
    R012 = 2.0564311499e+00
    R112 = -2.1311365518e-01
    R022 = -1.2419983026e+00
    R003 = -2.3342758797e-02
    R103 = -1.8507636718e-02
    R013 = 3.7969820455e-01
    ss = math.sqrt((SA + deltaS) / SAu)
    tt = CT / CTu
    zz = -Z / Zu
    rz3 = R013 * tt + R103 * ss + R003
    rz2 = (R022 * tt + R112 * ss + R012) * tt + (R202 * ss + R102) * ss + R002
    rz1 = (((R041 * tt + R131 * ss + R031) * tt + (R221 * ss + R121) * ss + R021) * tt + ((R311 * ss + R211) * ss + R111) * ss + R011) * tt + (((R401 * ss + R301) * ss + R201) * ss + R101) * ss + R001
    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt + ((R330 * ss + R230) * ss + R130) * ss + R030) * tt + (((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt + ((((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt + (((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000
    particle.seawater_density = ((rz3 * zz + rz2) * zz + rz1) * zz + rz0
    #particle.sw_surface_density = rz0
    
    





















































def vertical_mixing(particle, fieldset, time):
    """
    A markov-0 kernel for vertical mixing
    
    Description:
        A simple verticle mixing kernel that uses a markov-0 process to determine the vertical displacement of a particle.
        Ross & Sharples 2004
        The deterministic component is determined using forward-difference with a given delta_z. 

    Requirements:
        Fieldset:
            ---
        Particle:
            settling_velocity - 
        Order of Operations:
            To ensure the reflecting boundary condition of the random walk component, this kernel should be performed at the very end.
            Additionally, this kernel should be performed after the rising/sinking velocity of the particle has been computed.

            

    NOTES: previously called markov_0_mixing
    """
    
    # Sample the ocean vertical eddy diffusivity field KZ
    delta_z = 0.5 # [m], used to compute the gradient of kz using a forward-difference
    kz = fieldset.mixing_kz[time, particle.depth, particle.lat, particle.lon]
    kz_delta = fieldset.mixing_kz[time, particle.depth+delta_z, particle.lat, particle.lon]
    
    # Compute the gradient of kz using a forward-difference
    dkz_dz = (kz_delta - kz) / delta_z
    
    # Compute the deterministic component of Eq. (1)
    dz_deterministic = dkz_dz * particle.dt

    # Compute the random walk component of Eq. (1)
    dz_random = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3) * math.sqrt(2 * kz)
    ## TO DO - implement the reflective boundary condition
    
    # Compute rise velocity component of Eq. (1)
    dz_wb = particle.settling_velocity * particle.dt

    # Apply Eq. (1) - computing only the difference in depth
    ddepth = dz_deterministic + dz_random  + dz_wb

    # Update particle position
    particle.depth += ddepth 

    # bathymetry_local = fieldset.bathymetry[time, fieldset.z_start, particle.lat, particle.lon]
    
    # if potential < fieldset.z_start:
    #     particle.depth = fieldset.z_start
    # elif potential > bathymetry_local:
    #     particle.depth += 0 #TO DO: keep particles 'beached'
    #     particle.hit_bottom = 1
    # elif particle.depth > 100 and potential > (bathymetry_local*0.99): # for deeper particles; since bathymetry can be quite rough (and is interpolated linearly) look at the 99% value instead
    #     particle.depth += 0
    #     particle.hit_bottom = 1
    # elif potential > 3900:
    #     particle.depth += 0
    #     particle.hit_bottom = 1
    # else:
    #     particle.depth = potential        
        


### Biofouling related kernels ###



### Sampling related kernels ###

def PolyTEOS10_bsq(particle, fieldset, time):
    '''
    # calculates density based on the polyTEOS10-bsq algorithm from Appendix A.2 of
    # https://www.sciencedirect.com/science/article/pii/S1463500315000566
    # requires fieldset.abs_salinity and fieldset.cons_temperature Fields in the fieldset
    # and a particle.density Variable in the ParticleSet
    #
    # References:
    #  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate
    #   polynomial expressions for the density and specific volume of
    #   seawater using the TEOS-10 standard. Ocean Modelling.
    #  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003:
    #   Accurate and computationally efficient algorithms for potential
    #   temperature and density of seawater.  Journal of Atmospheric and
    #   Oceanic Technology, 20, 730-741.
    '''

    Z = - particle.depth  # note: use negative depths!
    SA = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon]
    CT = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]

    SAu = 40 * 35.16504 / 35
    CTu = 40
    Zu = 1e4
    deltaS = 32
    R000 = 8.0189615746e+02
    R100 = 8.6672408165e+02
    R200 = -1.7864682637e+03
    R300 = 2.0375295546e+03
    R400 = -1.2849161071e+03
    R500 = 4.3227585684e+02
    R600 = -6.0579916612e+01
    R010 = 2.6010145068e+01
    R110 = -6.5281885265e+01
    R210 = 8.1770425108e+01
    R310 = -5.6888046321e+01
    R410 = 1.7681814114e+01
    R510 = -1.9193502195e+00
    R020 = -3.7074170417e+01
    R120 = 6.1548258127e+01
    R220 = -6.0362551501e+01
    R320 = 2.9130021253e+01
    R420 = -5.4723692739e+00
    R030 = 2.1661789529e+01
    R130 = -3.3449108469e+01
    R230 = 1.9717078466e+01
    R330 = -3.1742946532e+00
    R040 = -8.3627885467e+00
    R140 = 1.1311538584e+01
    R240 = -5.3563304045e+00
    R050 = 5.4048723791e-01
    R150 = 4.8169980163e-01
    R060 = -1.9083568888e-01
    R001 = 1.9681925209e+01
    R101 = -4.2549998214e+01
    R201 = 5.0774768218e+01
    R301 = -3.0938076334e+01
    R401 = 6.6051753097e+00
    R011 = -1.3336301113e+01
    R111 = -4.4870114575e+00
    R211 = 5.0042598061e+00
    R311 = -6.5399043664e-01
    R021 = 6.7080479603e+00
    R121 = 3.5063081279e+00
    R221 = -1.8795372996e+00
    R031 = -2.4649669534e+00
    R131 = -5.5077101279e-01
    R041 = 5.5927935970e-01
    R002 = 2.0660924175e+00
    R102 = -4.9527603989e+00
    R202 = 2.5019633244e+00
    R012 = 2.0564311499e+00
    R112 = -2.1311365518e-01
    R022 = -1.2419983026e+00
    R003 = -2.3342758797e-02
    R103 = -1.8507636718e-02
    R013 = 3.7969820455e-01
    ss = math.sqrt((SA + deltaS) / SAu)
    tt = CT / CTu
    zz = -Z / Zu
    rz3 = R013 * tt + R103 * ss + R003
    rz2 = (R022 * tt + R112 * ss + R012) * tt + (R202 * ss + R102) * ss + R002
    rz1 = (((R041 * tt + R131 * ss + R031) * tt + (R221 * ss + R121) * ss + R021) * tt + ((R311 * ss + R211) * ss + R111) * ss + R011) * tt + (((R401 * ss + R301) * ss + R201) * ss + R101) * ss + R001
    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt + ((R330 * ss + R230) * ss + R130) * ss + R030) * tt + (((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt + ((((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt + (((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000
    particle.sw_density = ((rz3 * zz + rz2) * zz + rz1) * zz + rz0
    particle.sw_surface_density = rz0
    
    









## A list of kernels used to 'recover' particles after --- status codes ---

## @author: denes001
## @date: 2023-08-09

def periodicBC(particle, fieldset, time):
    if particle.lon <= -180.:
        particle.lon += 360.
    elif particle.lon >= 180.:
        particle.lon -= 360.

def delete_particle(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    if fieldset.verbose_delete == 1:
        print('particle is deleted out of bounds at lon = ' + str(particle.lon) + ', lat =' + str(
            particle.lat) + ', depth =' + str(particle.depth))
        
def delete_particle_interp(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    if fieldset.verbose_delete == 1:
        print('particle is deleted due to an interpolation error at lon = ' + str(particle.lon) + ', lat =' + str(
            particle.lat) + ', depth =' + str(particle.depth))
    
    particle.delete()
  