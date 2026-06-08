module params

use consts
implicit none

!---------------------------------------------------------------------------------------------------
! internally set parameters :

real epsv                        ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
real:: cpvf = 1.                  ! switch to take into account cpv in specific hrat of air
logical:: dosubsidence = .false. ! compute subsidence from profile in lsf file
real fcorz                       ! Vertical Coriolis parameter
real coszrs                      ! Sun's zenith angle
real salt_factor                 ! correction factor for water vapor saturation over sea-water
real:: sst_mean                  ! mean sst (set ny model when dosstclimo is .true.)

real :: cfl_max = 0.70  ! maximum CFL the model is allowed to run, default for nadams = 3
                        ! and nadv_mom = 2. For nadv_mom=23 or nadv_mom3 see setparms.f90
integer:: ncycle_max = 100  ! maximum number of subcycling within dt
logical:: SLM = .false.    ! flag to run Simplified Land Model
real(8), parameter :: deg2rad = pi/180.d0 ! conversion factor from degrees to radians
real(8), parameter :: rad2deg = 180.d0/pi ! conversion factor from radians to degrees

logical:: docheck = .false. ! used for gSAM -namelists check
logical:: dosubgridcloudfraction = .false. ! cld(:,:,:) carries subgrid cloud fraction

!--------------------------------------------------------------------------------------------------
! Parameters set by PARAMETERS namelist:

logical:: LES_S = .true.  ! if true cloud threshols for stats for PBL, false - for deep clouds
real:: cwp_threshold = 0.01 ! cld wat/ice path thresh for diag. of cld fraction

logical:: dodebug = .false.   ! print some extended diagnostic on each time step (be careful with file-size!)
logical:: docheckmom = .false.  ! print info on conservation of momentum
logical:: docheckenergy = .false.  ! print info on conservation of enetgy

integer:: nadv_mom = 2  ! space-order of momentum advection (can be 2,3,and 23, which is a mix of 2 and 3)
real:: alpha_hybrid = 1. ! fraction of scheme order 2 over scheme order 3 for momentum (when nadv_mom = 23)

logical:: doimplicitdiff=.false. ! switch to implicit diffusion of momentum and scalars in z

integer:: npressure_iter = 0 !number of additional iterations of pressure solver to minimize parasite flow inside topography

real:: ug = 0.        ! Velocity of the domain's drift in x direction
real:: vg = 0.        ! Velocity of the domain's drift in y direction
real:: fcor = -999.   ! Coriolis parameter	
real:: longitude0 = 0.    ! latitude of the domain's left boundary 
real:: latitude0  = 0.    ! longitude of the domain's center 

logical:: docurrentco2 = .false. ! co2 is automatically computed depending on year
real:: nxco2 = 1.         ! factor to modify co2 concentration (e.g. nxco2 = 2 mean 2xCO2)
real:: n2ox = 0.  ! prescrined uniform n2o (g/g)
real:: ch4x = 0.  ! prescrined uniform ch4 (g/g)
real:: cfc11x = 0.  ! prescrined uniform cfc11 (g/g)
real:: cfc12x = 0.  ! prescrined uniform cfc12 (g/g)

real(8):: earth_factor = 1._8 ! radius of planet relative to real Earth i
                              ! (=radius real-Earth/radius of your planet) 
real(8):: gamma_RAVE = 1._8  ! RAVE gamma-factor (from Kuang, Blossey and Bretherton (GRL 2005))

real(8):: tabs_s =0.    ! surface temperature,K

! slab ocean model parameters. See simple_ocean.f90
real:: delta_sst = 0.        ! amplitude of sin-pattern of sst about tabs_s (ocean_type=1)
real:: depth_slab_ocean = 2. ! thickness of the slab-ocean (m)
real:: shift_sst = 0.        ! shift SST pattern in latitude, deg, positive shift to north
real:: Szero = 0.            ! mean ocean transport (W/m2) 
real:: deltaS = 0.           ! amplitude of linear variation of ocean transport (W/m2)
real:: timesimpleocean = 0.  ! time to start simple ocean
real:: sst_climo =0.         ! "climatological" SST
real:: tau_ocean =0.         ! prescribed ocean surface drag 

logical:: sfc_flx_fxd =.false. ! surface sensible flux is prescribed
logical:: sfc_tau_fxd =.false. ! surface drag is prescribed

real:: lhf_fudge = 1.    ! fudge factor for latent flux LHF over ocean. 
real:: shf_fudge = 1.    ! fudge factor for sensible flux SHF over ocean. 
real:: tau_fudge = 1.    ! fudge factor for stress TAU over ocean. 
     ! Normally these factors should be used for 
     ! tuning performance depending on the grid, based on reanalysis (e.g. ERA5)-MK, Dec2024

real::   fluxt0 =0.            ! prescribed surface sensible flux, W/m2
real::   fluxq0 =0.            ! prescribed surface latent flux, W/m2
real::   tau0   =0.            ! prescribed surface stress, m2/s2
real::   z0     = 0.035        ! default roughness length (m)

integer:: ocean_type =0        ! type of SST forcing
logical:: OCEAN =.false.       ! flag indicating that surface is water
logical:: LAND =.false.        ! flag indicating that surface is land
logical:: ISLAND=.false.       ! flag to indicate that both ocean and land are present

! specify idealized islands. Defaul is rectangular island
real:: island_x1 =  1.1        ! island boundary in units of fraction of domain in x
real:: island_x2 = -0.1        ! island boundary in units of fraction of domain in x
real:: island_y1 =  1.1        ! island boundary in units of fraction of domain in y
real:: island_y2 = -0.1        ! island boundary in units of fraction of domain in y
logical:: doroundisland  = .false. ! specify round islands at the center of the domain
real:: island_radius = 0.          ! radius (in m) of a round island

logical:: readinit = .false.       ! read 3D initial conditions from file
logical:: readlandmask = .false.   ! read landmask
logical:: readsst = .false.        ! read sst from file
logical:: readterr = .false.       ! read terrain
logical:: reado3 = .false.         ! read ozone file
logical:: readlatlon = .false.     ! read lat/lon grid from file 
logical:: readlat = .false.        ! read lat grid from file 

logical:: w3D_pressure = .true.    ! flag for vertical velocity type: 
                                   ! .true. - pressure velocity (Ps/s),
                                   ! .false. - verical velocity (m/s)

logical:: read_meters = .false.! read 3D data from file with coords in meters

! place holders for names of data files:
character(120):: initfile = ""
character(120):: landmaskfile = ""
character(120):: sstfile = ""
character(120):: terrainfile = ""
character(120):: o3file = ""
character(120):: latlonfile = ""
logical:: latlonfilebin = .false.

logical:: docap_snd_cu = .false. ! limit the magnitude of horiz. velocity in input sounding (to reduce CFL)
real:: cap_snd_cu = 0.5 ! capping Courant (CFL)
logical:: dosimplesnd = .false.! initialize the model by a single sounding with no vertical
                               ! interpolations, that is the number of the height and pressure 
                               ! levels should be the same as number of grid levels. 

real:: gmg_precision = 1.e-6  ! precision of GMG solver (the smaller the more V-cycles will be needed)

integer:: nrad = 1           ! frequency of calling the radiation routines (in time steps)
integer:: nrad_ems = 6       ! frequency of updating gaseous emiss/absorp coeffs (every nrad*nrad_ems steps) 

logical:: dodatefilename = .false.  ! insert current date to names of 2D and 3D output files

! cap winds in data (initial and nudging) 
logical:: do_cap_wind = .false. ! do it or not
real:: cap_wind = 0.            ! maximum wind speed in datasets (has to be specified if do_cap_wind = T
!------------------------------------------------------------------------------------------------
! Cyclic datasets

logical:: docyclic = .false.  ! make all time-sependent input datasets cyclic or periodic in time
integer:: cycle_period = 365    ! period of cycling dataset, days
!-------------------------------------------------------------------------------------------------
! flags and variables related to damping/nudging:

logical:: dodamping = .false.       ! do top-of-domain damping
logical:: dodamping_poles = .false. ! do additional damping around the poles
logical:: dodamping_w = .false.     ! damp vertical velocity to maintain CFL below critical
logical:: dodamping_u = .false.     ! damp horizontal velocity at domain top to maintain CFL below critical

real:: nub = 0.6           ! normalized height where top-of-domain damping begins 
real:: damping_u_cu = 0.5  ! Courant number used in damping of horizontal velocity (dodamping_pole=.true.)
real:: damping_w_cu = 0.5  ! Courant number used in damping of vertical velocity (dodamping_w=.true.)

logical:: donudging_uv = .false.    ! nudge horizontal velocities to profiles in lsf or snd files
logical:: donudging_w = .false.     ! nudge vertical velocities to nudging fields (works only with donudge3D=T)
logical:: donudging_tq = .false.    ! nudge temperature and vapor to profiles in snd file
logical:: donudging_t = .false.     ! nudge temperature to profiles in snd file
logical:: donudging_q = .false.     ! nudge vapor to profiles in snd file
logical:: donudge3D = .false.       ! do nudging to 3D fields

logical:: dospectralnudging = .false.    ! do spectral nudging
integer:: nx_spectral=-1, ny_spectral=-1 ! size of the averaging box in x and y   
integer:: nstep_spectral = 1             ! frequency of updating spectral nudging (as it is expensive)   

character(120):: nudge3D_dir=""     ! directory with data files
character(120):: nudge3D_file=""    ! name of the file-list file 
real:: nudge3D_tau = 0.             ! nudging time-scale
integer:: nudge3Dstep_start = 0     ! step when 3D nudging starts. 
integer:: nudge3Dstep_end = 9999999 ! step when 3D nudging ends.
real:: timelargescale =0.           ! time to start large-scale forcing

! nudging boundaries (between z1 and z2, where z2 > z1): 
real:: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
real:: nudging_t_z1 =-1.,  nudging_t_z2 = 1000000.
real:: nudging_q_z1 =-1.,  nudging_q_z2 = 1000000.
real:: tauls = 99999999.            ! nudging-to-large-scaler-profile time-scale
real:: tautqls = 99999999.          ! nudging-to-large-scaler-profile time-scale for scalars

logical:: dobufferzonex = .false.   ! use nudging zone at domain lateral boundaries in x
logical:: dobufferzoney = .false.   ! use nudging zone at domain lateral boundaries in y
real:: bufferzonex = 0.             ! fraction of the domain in x for nudging zone
real:: bufferzoney = 0.             ! fraction of the domain in y for nudging zone

!-------------------------------------------------------------------------------------------------------
! clouds/precipitation

logical:: docloud = .false.         ! condensation is allowed
logical:: doprecip = .false.        ! precipitation processes are allowed

!-------------------------------------------------------------------------------------------------------
! chemistry

logical:: dochem = .false.         ! chemistry occurs

!-------------------------------------------------------------------------------------------------------
! radiation:

logical:: dolongwave = .false.      ! compute longwave radiation transfer
logical:: doshortwave = .false.     ! compute solar radiation transfer
logical:: doradlat = .false.        ! compute incoming solar flux as function of latitude
logical:: doradlon = .false.        ! compute incoming solar flux as function of longitude
logical:: doradavg = .true.         ! time-avarage input fields between radiation calls
logical:: doequinox = .false.       ! force equinox conditions
logical:: doradsimple = .false.     ! use simple analytical radiation instead of radiation package
logical:: doradforcing = .false.    ! apply prescribed radiation heating rate profiles (file rad is needed)
logical:: doradhomo = .false.       ! horizontally homogenize radiation heating rates
logical:: doradhomozonal = .false.  ! zonally homogenize radiation heating rates
logical:: dosolarconstant = .false. ! specify solar constant
logical:: doperpetual = .false.     ! do perpetual solar radiation (constant and no seasonas)
logical:: doseasons = .false.       ! solar radiation depends on calendar day
! Specify solar constant and zenith angle for perpetual insolation.
! Based onn Tompkins and Graig (1998)
! Note that if doperpetual=.true. and dosolarconstant=.false.
! the insolation will be set to the daily-averaged value on day0.
real:: solar_constant = 685. ! solar constant (in W/m2)
real:: zenith_angle = 51.7   ! zenith angle (in degrees)

!-------------------------------------------------------------------------------------------------------
! surface fluxes, ocean/sst forcing

logical:: dosurface = .false.       ! surface fluxes are allowed
logical:: dosfcforcing = .false.    ! apply prescribed surface temperature and fluxes (sfc file is needed)
logical:: dosfchomo = .false.       ! horizontally homogenize surface fluxes
logical:: dossthomo = .false.       ! horizontally homogenize surface temperature
logical:: dodynamicocean = .false.  ! use actual ocean model
logical:: dooceanonly = .false.     ! run the ocean model only
logical:: doequilocean = .false.    ! compute surface temperature assuming surface flux equilibrium
logical:: doslabocean = .false.     ! use slab ocean model
logical:: dosstclimo = .false.      ! nudge SST to "climatological" value (set by sst_climo)
logical:: dossthomozonal = .false.  ! zonally homogenize surface temperature
logical:: dotc = .false.            ! modify surface transfer coefficient in oceflx.f90 for very strong winds


logical:: doupperbound = .false.    ! maintain virtical gradient from sounding at the domain top boundary
logical:: dosgs = .false.           ! use subgrid-scale parameterization
logical:: docoriolis = .false.      ! Coriolis force is allowed
logical:: docoriolisz = .false.     ! Vertical coriolis force is allowed
logical:: dometric = .true.         ! Calulate metric curviture terms (tangent) for spherical coordinates
logical:: dofplane = .true.         ! f-plane, that is Coriolis parameter is the same everywhere
logical:: dolargescale = .false.    ! apply large-scale forcing (file lsf is needed)
logical:: doensemble = .false.      ! do ensemble run by perturbing initial conditions (see nensemble and nens)
logical:: dowallx = .false.         ! put solid walls in x-domain (be careful)
logical:: dowally = .false.         ! put solid walls in y-domain (automatically used in global runs)
logical:: doregion = .false.        ! make regional run with prescribed boundary condition
logical:: docolumn = .false.        ! do single-column version of the model 
logical:: dofliplon = .false.       ! switch to degrees west (negative latitudes) when reading data
logical:: notracegases = .false.    ! if true no trace gases are used except for CO2
logical:: dotracers = .false.       ! do tracer transport
logical:: dotrsfcflux = .false.     ! do flux of tracers into the surface (deposition)
logical:: dosmoke = .false.         ! smoke cloud when docloud=.false.
logical:: doseawater = .false.      ! salty seawater or fresh (lake)
logical:: doterrain = .false.       ! do topography
logical:: dobuildings = .false.     ! compute building effects  

logical:: doseaice = .false.        ! assume that SST below 271K is a sea-ice surface
real:: seaicethickness = 1.        ! sea-ice prescrived thickness (m) when doseaice = .true.
logical:: doseaiceevol = .false.        ! parameterization of seaice thickness in on

!----------------------------------------------------------------------------------------
logical:: docup = .false.           ! call convective parameterization 
integer:: n_cup = 6                 ! frequency of calling cumulus patameterization (time steps)
!----------------------------------------------------------------------------------------------
! Run setup:

logical:: dolatlon = .false. ! do lat/lon grid
logical:: doglobal = .false. ! global model. Poles are at 90N/S. 
                             ! domain is periodic in zonal direction
                             ! dolatlon and dowally are set to .true. automatically
logical:: doglobalpresets=.false. ! set some parameters for global runs automatically
                                  ! see setparams.f90 for the presets used.
                                  ! also, it affects  settings in MICRO_SAM1MOM 
logical:: donearglobal = .false. ! global model. solid walls are determined using dlat. 
                             ! domain is periodic in zonal direction
                             ! dolatlon and dowally are set to .true. automatically
logical:: doflat = .false.   ! make cartesian coordinates even when doglobal or dolatlon=.true.
                              ! Instead the pole is surrounded by a wall at latitude figured out from the lat grid.
logical:: dohs94 = .false.   ! do Held-Suares dry test

logical:: dofixdynamics = .false.   ! 3D wind is prescribed (set in setfixdynamics.f90)
integer:: fixdynamics_type = 1 ! chose fixed wind (see setfixdynamics.f90)
logical:: donodynamics = .false.    ! wind is not evolving
logical:: donobuoyancy = .false.    ! buoyancy is not computed

!-------------------------------------------------------------------------------------------------
! Tracers transport

logical:: dopointsource = .false. ! incert sources of tracers
integer:: pointsource_start_step = 1 ! time step when sources of tracers are activated
integer:: pointsource_stop_step = 1 ! time step when sources of tracers are disactivated
integer:: tracer_source_type = 0 ! type of trace source

!----------------------------------------------------------------------------------------

integer:: perturb_type  = 0 ! type of initial noise in setperturb()

! forecasting parameters:
logical:: dofcast = .false.  ! don't change sst, LAI during run
logical:: dofcast_ice = .false.  ! don't change seaice  during run 
integer:: nensemble =0   ! the number of subensemble set of perturbations
integer:: nens = 0          ! add this number to seed in setperturb to make ensembles
!----------------------------------------------------------------------------------------

integer:: terrain_type = 0 ! type of the terrain initialization
logical:: doterrhyb = .false. ! use hybrid advection for momentum around terrain
logical:: doterrepair = .true. ! repair terrain for one point picks, holes etc

! Initial bubble parameters. Activated when perturb_type = 2
  real:: bubble_x0 = 0.
  real:: bubble_y0 = 0.
  real:: bubble_z0 = 0.
  real:: bubble_radius_hor = 0.
  real:: bubble_radius_ver = 0.
  real:: bubble_dtemp = 0.
  real:: bubble_dq = 0.

logical:: donorestart = .false. ! save no restart files

logical:: dostatcoars = .false.     ! collect coarse-grid stats (stat_coars.f90)

end module params
