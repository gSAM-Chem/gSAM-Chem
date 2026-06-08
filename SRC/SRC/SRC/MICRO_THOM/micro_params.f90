module micro_params
  !bloss: Holds microphysical parameter settings separately, so that they can be used
  !   inside both the microphysics module and the module(s) from WRF.

  logical :: doicemicro = .true. ! Turn on ice microphysical processes

  logical :: doaerosols = .false. ! Use Thompson-Eidhammer water- and ice-friendly aerosols

  real :: Nc0 = 100. ! initial/specified cloud droplet number concentration (#/cm3).

  logical:: doeiscld = .false.    ! set shallow cloud fraction for radiation based on EIS

!bloss(2025-03) start ==============
  ! option to tune snow fallspeed relationship
  !    v = a*D^b * exp(-fv*D), with units of v (m/s), D (m)
  real :: snow_fall_a_fudge = 1. ! scaling factor
  real :: snow_fall_a = 40. ! m/s
  real :: snow_fall_b = 0.55 ! exponent
  real :: snow_fall_fv = 195.
  ! Suggested options: (with v for D=10, 100, 1000um)
  ! THOM default: snow_fall_a = 40., snow_fall_b = 0.55, snow_fall_fv = 195., v = 0.07, 0.25, 0.74
  ! M2005 default: snow_fall_a = 11.72, snow_fall_b = 0.41, snow_fall_fv = 0., v = 0.1, 0.27, 0.69
  ! all data fit, Vasquez-Martin et al (2021, https://doi.org/10.5194/acp-21-7545-2021)
  !   snow_fall_a = 1.63, snow_fall_b = 0.19, snow_fall_fv = 0., v = 0.18, 0.28, 0.44
  
  ! option to tune cloud ice fallspeed relationship
  !    v = a*D^b, with units of v (m/s), D (m)
  ! NOTE v = MIN(1.2, v), so large predicted fallspeeds will not be applied.
  real :: clice_fall_a_fudge = 1. ! scaling factor
  real :: clice_fall_a = 1847.5 ! m/s
  real :: clice_fall_b = 1.0 ! exponent
  ! Suggested options: (with v for D=10, 100, 1000um)
  ! THOM default: clice_fall_a = 1847.5, clice_fall_b = 1.0, v = 0.02, 0.18, 1.85 m/s
  ! M2005 default: clice_fall_a = 700., clice_fall_b = 1.0, v = 0.01, 0.07, 0.7 m/s
  ! M2005 MK tune: clice_fall_a = 700., clice_fall_b = 0.865, v = 0.03, 0.24, 1.78 m/s
  ! spherical ice fit, Vasquez-Martin et al (2021, https://doi.org/10.5194/acp-21-7545-2021)
  !   clice_fall_a = 86617.76, clice_fall_b = 1.42, v = 0.01, 0.18, 4.76 m/s
  
  ! options to preserve cloud ice:
  logical :: do_preserve_cloud_ice = .false.
  real :: D0i_iau = 20.e-6 
  !  When true, this does four things:
  !  - disables cloud ice deposition being added to snow
  !  - cloud-ice-to-snow autoconversion only happens when Sice > 1
  !  - weakens the maximum allowed cloud ice sink in one time step to (qice/30 sec), and
  !  - zeroes out cloud-ice-to-snow autoconversion for xDi<D0i_iau
  !      and ramps up to the default value at xDi=2*D0i_iau.
  !  - cloud-ice-to-snow autoconversion threshold, units m, default = 0.1*D0s = 20.e-6 m
  !bloss(2025-03) end ==============

  ! option to allow the gamma exponent for rain, graupel and cloud ice to be specified.
  !   Note that mu=0 (the default) corresponds to an exponential distribution.
  real :: fixed_mu_r = 0.
  real :: fixed_mu_i = 0.
  real :: fixed_mu_g = 0.

  ! Fix the exponent in the gamma droplet size distribution for cloud liquid to a single value.
  logical :: dofix_mu_c = .false. 
  real :: fixed_mu_c = 10.3 ! fixed value from Geoffroy et al (2010).
  ! Full citation: Geoffroy, O., Brenguier, J.-L., and Burnet, F.:
  !   Parametric representation of the cloud droplet spectra for LES warm
  !   bulk microphysical schemes, Atmos. Chem. Phys., 10, 4835-4848,
  !   doi:10.5194/acp-10-4835-2010, 2010.

  !..Densities of rain, snow, graupel, and cloud ice. --> Needed in radiation
        REAL, PARAMETER :: rho_water = 1000.0
        REAL, PARAMETER :: rho_snow = 100.0
        REAL, PARAMETER :: rho_graupel = 500.0
        REAL, PARAMETER :: rho_cloud_ice = 890.0

        ! If using effective radii for ice in radiation (i.e., CAM radiation),
        !   this option lets you use the Thompson et al (2016) effective radii for snow and ice
        !   Otherwise, the generalized effective size (Dge) is computed and rescaled into
        !   an effective radius following eqn 3.12 in Fu (1996, JClim).
        logical :: doThompsonEtAl2016SnowEffRad = .false.

        ! Fix rain number generation from melting snow (backported from WRF V3.9)      
        logical :: BrownEtAl2017_pnr_sml_fix = .true. 

  !bloss(2018-02): Enable choice of snow moment
  !parameterizations between Field et al (2005), the default,
  !and Field et al (2007).
  logical :: doFieldEtAl2007Snow = .false.
  
  ! Field et al (2007) has two snow size distributions: tropical and mid-latitude.
  !   The size distribution is used in the computation of sedimentation.
  logical :: TropicalSnow = .true. ! if false, use mid-latitude size distribution

  logical :: do_output_process_rates = .false.
  integer, parameter ::  nproc_rates_mass_thompson_cold = 31, nproc_rates_number_thompson_cold = 17
  integer, parameter ::  nproc_rates_mass_thompson_warm = 11, nproc_rates_number_thompson_warm = 5
  integer ::  nproc_rates_mass = 1 , nproc_rates_number = 1 ! default value is one

  character*80 :: lookup_table_location = './RUNDATA/'

   logical:: dowarmcloud=.true. ! dummy parameter - don't change
   logical:: docloudfall=.false. ! dummy parameter - don't change


end module micro_params
