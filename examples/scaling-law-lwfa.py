# This is the PICMI file for the simulation presented in Wei Lu et al
# 3D LWFA with a 200TW laser which yielded > 1GeV mono-energetic 
# electrons
#
#  This should be the only line that needs to be changed for different codes
# e.g. `from pywarpx import picmi`
#      `from fbpic import picmi`
#      `from warp import picmi`
#      'from osiris import picmi'
#
from fbpic import picmi

# Create alias fror constants
cst = picmi.constants

# Run parameters - can be in separate file
# ========================================

# Physics parameters
# ------------------

# *****************************************************************************************
# --- laser
    # Specifies a Gaussian laser distribution
    #   - name=None: Optional name of the laser
    #   - wavelength: Laser wavelength
    #   - waist: Waist of the Gaussian pulse at focus [m]
    #   - duration: Duration of the Gaussian pulse [s]
    #   - focal_position=[0,0,0]: Position of the laser focus (vector) [m]
    #   - centroid_position=[0,0,0]: Position of the laser centroid at time 0 (vector) [m]
    #   - propagation_direction=[0,0,1]: Direction of propagation (unit vector) [1]
    #   - polarization_direction=[1,0,0]: Direction of polarization (unit vector) [1]
    #   - a0: Normalized vector potential at focus
    #         Specify either a0 or E0 (E0 takes precedence).
    #   - E0: Maximum amplitude of the laser field [V/m]
    #         Specify either a0 or E0 (E0 takes precedence).
    #   - zeta: Spatial chirp at focus (in the lab frame) [m.s]
    #   - beta: Angular dispersion at focus (in the lab frame) [rad.s]
    #   - phi2: Temporal chirp at focus (in the lab frame) [s^2]
# *****************************************************************************************


import math
laser_polarization   = [0,1,0]   # Polarization angle (in rad)
propagation_direction= [1,0,0]   # direction of propagation
laser_a0             = 4.        # Normalized potential vector
laser_wavelength     = 8e-07     # Wavelength of the laser (in meters)
laser_waist          = 5e-06     # Waist of the laser (in meters)
laser_duration       = 7.89955027e-14    # Duration of the laser (in seconds)
laser_injection_loc  = 9.e-6     # Position of injection (in meters, along z)
laser_focal_distance = 100.e-6   # Focal distance from the injection (in meters)
laser_t_peak         = 30.e-15   # The time at which the laser reaches its peak
                                 # at the antenna injection location (in seconds)
centroid_position = [90.57,6.3662e-5,6.3662e-5]                           
focal_position= centroid_position# laser is focused at the start of the simulation



# 
# --- plasma
# 
plasma_density_expression = "(x < 1.0267e-4)? 0 : 1.524195e24 "
plasma_min     = [0, 0, 0]
plasma_max     = [ 2e-2,  1.27324e-4,  1.27324e-4]



# Numerics parameters
# -------------------

# --- Nb time steps
max_steps = 1000

# --- grid
nx = 4032
ny = 256
nz = 256
xmin = 0
xmax = 1.5*plasma_max[0]
ymin = plasma_min[1]
ymax = plasma_max[1]
zmin = plasma_min[2]
zmax = plasma_max[2]
moving_window_velocity = [cst.c., 0., 0]
n_macroparticle_per_cell = [2, 1, 1]

# --- geometry and solver
em_solver_method = 'Yee'  
geometry = '3D'
time_step_size = 8.409198678305278e-17 # step size in seconds


# Physics part - can be in separate file
# ======================================

# Physics components
# ------------------

# --- laser
laser = picmi.GaussianLaser(
    wavelength             = laser_wavelength,
    waist                  = laser_waist,
    duration               = laser_duration,
    focal_position         = [0., 0., laser_focal_distance + laser_injection_loc],
    centroid_position      = [0., 0., laser_injection_loc - cst.c*laser_t_peak],
    polarization_direction = [0,0, 1.],
    propagation_direction  = [0,0,1],
    a0                     = laser_a0)

# --- plasma
plasma_dist = picmi.AnalyticDistribution(
                density_expression = plasma_density_expression,
                lower_bound        = plasma_min,
                upper_bound        = plasma_max,
                fill_in            = True)
plasma = picmi.MultiSpecies(
                particle_types = ['electron'],
                names          = ['e-'],
                charge_states  = [None],
                proportions    = [1.0],
                initial_distribution=plasma_dist)
# Individual species in a `MultiSpecies` can be addressed either
# with their index (using Python indexing conventions) or with their name
# (if the user provided a name)
# Set the ionization for the species number 1 (Argon)
# and place the created electrons into the species number 2 (electron)
if picmi.codename != 'warpx':
    plasma['Argon'].activate_field_ionization(
        model           = "ADK", # Ammosov-Delone-Krainov model
        product_species = plasma['e-'])



# Numerics components
# -------------------

if geometry == '3D':
    grid = picmi.Cartesian3DGrid(
        number_of_cells           = [nx, ny, nz],
        lower_bound               = [xmin, ymin, zmin],
        upper_bound               = [xmax, ymax, zmax],
        lower_boundary_conditions = ['open','periodic', 'periodic'],
        upper_boundary_conditions = ['open','periodic', 'periodic'],
        moving_window_velocity    = moving_window_velocity)
        # Note that code-specific arguments use the code name as a prefix.
elif geometry == 'RZ':
    # In the following lists:
    # - the first element corresponds to the radial direction
    # - the second element corresponds to the longitudinal direction
    grid = picmi.CylindricalGrid(
        number_of_cells           = [nx//2, nz],
        lower_bound               = [0., zmin],
        upper_bound               = [xmax, zmax],
        lower_boundary_conditions = [ None, 'open'],
        upper_boundary_conditions = ['reflective', 'open'],
        n_azimuthal_modes         = 2,
        moving_window_zvelocity   = moving_window_velocity[-1])

smoother = picmi.BinomialSmoother( n_pass       = 4,
                                   compensation = True )
solver = picmi.ElectromagneticSolver( grid            = grid,
                                      cfl             = 1.,
                                      method          = em_solver_method,
                                      source_smoother = smoother)

# Diagnostics
# -----------
field_diag = picmi.FieldDiagnostic(name = 'diag1',
                                   grid = grid,
                                   period = 4444)
part_diag = picmi.ParticleDiagnostic(name = 'diag1',
                                     period = 4444,
                                     species = [e-])

# Simulation setup
# -----------------

# Initialize the simulation object
# Note that the time step size is obtained from the solver
sim = picmi.Simulation(solver = solver, verbose = 1)

# Inject the laser through an antenna
antenna = picmi.LaserAntenna(
                position = (0, 0, 9.e-6),
                normal_vector = (0, 0, 1.))
    
sim.add_laser(laser, injection_method = antenna)

# Add the plasma: continuously injected by the moving window
plasma_layout = picmi.GriddedLayout(
                    grid = grid,
                    n_macroparticle_per_cell = n_macroparticle_per_cell)
sim.add_species(species=plasma, layout=plasma_layout)

# Add the diagnostics
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)

# Picmi input script
# ==================

run_python_simulation = True

if run_python_simulation:
    # `sim.step` will run the code, controlling it from Python
    sim.step(max_steps)
else:
    # `write_inputs` will create an input file that can be used to run
    # with the compiled version.
    sim.set_max_step(max_steps)
    sim.write_input_file(file_name='input_script')
    # for OSIRIS, using the default input name
    # 
    sim.write_input_file(file_namee='os-stdin')



# ***********************************************************************
# ***********************************************************************
# below I am including the original osiris input deck for reference
# ***********************************************************************
# ***********************************************************************

# !----------the node configuration for this simulation----------
# node_conf 
# {
#  node_number(1:3) =  48,2,2,
#  if_periodic(1:3) = .false., .true., .true.,
# }


# !----------spatial grid----------
# grid 
# {
#   nx_p(1:3) = 4032, 256, 256,
#   coordinates = "cartesian",
# }


# !----------time step and global data dump timestep number----------
# time_step 
# {
#   dt     =   0.198,
#   ndump  =  4444, 
# }


# !----------restart information----------
# restart 
# {
#   ndump_fac = 2,  file_name  = ' ',
#   if_restart=.true.,
# }


# !----------spatial limits of the simulations----------
# !(note that this includes information about
# ! the motion of the simulation box)
# space 
# {
#   xmin(1:3) =   0.000d0 , 0.000d0 , 0.00,
#   xmax(1:3) =   806.4, 1000.0, 1000.0,
#   if_move= .true., .false., .false.,
# }


# !----------time limits ----------
# time 
# {
#   tmin = 0.0d0, tmax  = 200000.0d0 ,
# }


# !----------field solver set up----------
# el_mag_fld 
# {
# }

# !----------boundary conditions for em-fields ----------
# emf_bound 
# {
#   type(1:2,1) =  30, 30, 
#   type(1:2,2) =  1, 1,
#   type(1:2,3) =  1, 1,
# }

# !----------diagnostic for electromagnetic fields---------- 
# diag_emf 
# {
#   ndump_fac_all = 0,  
#   ndump_fac_ave = 1, 
#   n_ave(1:3)      = 2, 2,2,
#   ifdmp_efl = .true. , .true. , .true. ,
#   ifdmp_bfl = .true. , .true. , .true. ,
# }


# !----------number of particle species----------
# particles 
# {   num_species =1,   
# }


# !----------diagnostics for all particles----------
# diag_particles 
# {
#   ndump_fac = 0,
# }


# !----------information for species 1----------
# species 
# {
#   num_par_max = 2650000,
#   rqm=-1.0,
#   num_par_x(1:3) = 2, 1, 1,
#   vth(1:3) = 0.0d0 , 0.0d0 , 0.0d0 ,
#   vfl(1:3) = 0.0d0 , 0.0d0 , 0.0d0 ,
#   den_min = 0.0000001d0,
#   n_sort=20,
# }

# !----------density profile for this species----------
# !  number of points in profile along each direction
# !----------density profile for this species----------
# !  number of points in profile along each direction
# num_x 
# {  num_x = 19,   
# }

# !  actual profile
# profile 
# {
#   profile_type='channel_obj',
# }

# channel 
# {
#   c_start=806.4001,c_end=200000.0,
#   n_in=0.000875,n_out=0.0008751,n_slope=1.0,
#   center(1:2)=500.0, 500.0,
#   width=980,
# }

# !----------boundary conditions for this species----------
# spe_bound 
# {
#   type(1:2,1)=40,40,
#   type(1:2,2)=1,1,
#   type(1:2,3)=1,1,
# }

# !----------diagnostic for this species----------
# diag_species 
# {
#   ndump_fac_raw = 1,
#   ndump_fac_ene = 1,

#   ndump_fac_pha = 1,  file_name_pha = ' ',
#   ps_xmin(1:3) = 0.0, 0.0, 0.0,  ps_pmin(1:3) = -400.0, -10.0, -10.0,
#   ps_xmax(1:3) = 0.0, 0.0, 0.0, ps_pmax(1:3) = 2600.0, 10.0, 10.0,
#   ps_nx(1:3)   = 200,  128, 128,  ps_np(1:3)   =  1000,  200,  200,
  
#   if_ps_p_auto(1:3) = .true., .true., .true.,  

  
#   if_x2x1    = .true.,
#   if_x3x1    = .false.,
#   if_p1x1    = .true.,
#   if_p2x1    = .true.,
#   if_p3x1    = .true.,
#   if_x3x1    = .true.,
#   if_x3x2x1    = .true.,
#   if_p1x2    = .true.,
#   if_p2x2    = .true.,
#   if_p3x2    = .true.,
#   if_p1x3    = .false.,
#   if_p2x3    = .false.,
#   if_p3x3    = .false.,
#   if_p2p1    = .true.,
#   if_p3p1    = .true.,
#   if_p3p2    = .false.,
#   iftestpar  = .false.,
#   iforigin   = .false.,
  
#   particle_fraction = 1.0,
#   gamma_limit = 40.0d0,

# }

# !----------number of pulses----------
# pulse_sequence 
# {   num_pulses = 2,   
# }


# !----------information for pulse 1----------
# pulse 
# {
#   iflaunch = .true.,
#   wavetype=1,
#   w0(1:2)=153.0,153.0, rise=93.0,      
#   fall=93.0,      length= 0.0,
#   vosc=2.82843,     rkkp=1.0,      pol=90.0, phase=0.0,
#   start=0.50,   focus=-1.0, offset=0.0,  time=0.0,
# }
# !----------information for pulse 2----------
# pulse 
# {
#   iflaunch = .true.,
#   wavetype=1,
#   w0(1:2)=153.0,153.0, rise=93.0,      
#   fall=93.0,      length= 0.0,
#   vosc=2.82843,     rkkp=1.0,      pol=0.0, phase=90.0,
#   start=0.50,   focus=-1.0, offset=0.0,  time=0.0,
# }

# smooth 
# {
#   ifsmooth(1) = .true.,
#   smooth_level(1)=5,
#   swfj(1:3,1,1) = 1,2,1,
#   swfj(1:3,2,1) = 1,2,1,
#   swfj(1:3,3,1) = 1,2,1,
#   swfj(1:3,4,1) = 1,2,1,
#   swfj(1:3,5,1) = -5,14,-5,
#   ifsmooth(2) = .false.,
#   smooth_level(2)=2,
#   swfj(1:3,1,2) = 1,2,1,
#   swfj(1:3,2,2) = 1,2,1,
# }

# !----------diagnostic for currents----------
# diag_phy_field 
# {
#   ndump_fac_all = 0,  file_name_all = ' ',
#   ndump_fac_ave = 0, file_name_ave = ' ',
#   n_ave(1:3) = 1,1,1,
#   ifdmp_phy_field = .true. , .true. , .true. ,
# }

# smooth 
# {
#  ifsmooth(1) = .false.,
#  smooth_level(1)=2,
#  swfj(1:3,1,1) = 1,2,1,
#  swfj(1:3,2,1) = 1,2,1,
#  ifsmooth(2) = .false.,
#  smooth_level(2)=2,
#  swfj(1:3,1,2) = 1,2,1,
#  swfj(1:3,2,2) = 1,2,1,
# }

# !-----------diagnostic for charge-----------
# diag_phy_field 
# {
#   ndump_fac_all= 0,  file_name_all = ' ',
#   file_name_ave = ' ',
#   ifdmp_phy_field = .true. ,
# }

