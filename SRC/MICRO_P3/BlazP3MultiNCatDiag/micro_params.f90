module micro_params

!bloss: move indices of different microphysical species here, so that they
!  can be accessed in both microphysics.f90 and module_mp_p3.f90.

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, incl, inr, iqr, itabs_old, iqv_old
integer, dimension(:), allocatable :: inci, iqit, iqir, iqib

! SETTINGS FOR WARM RAIN MICROPHYSICS (AND DEFAULT VALUES)

! choice of warm rain microphysics scheme
!   = 1 for Seifert-Beheng (2001)
!   = 2 for Beheng (1994)
!   = 3 for Khairoutdinov-Kogan (2000) <-- DEFAULT
integer, public :: iWarmRainScheme = 3 

! control properties of cloud drop size distribution.
logical, public :: log_predictNc = .false. ! if true, predict droplet number
real, public :: Nc0 = 200. ! specified cloud droplet number conc (CDNC, #/cm3)

logical, public :: dofix_pgam = .false. ! option to fix value of exponent in cloud water gamma distn
real, public ::    pgam_fixed = 10.3 ! Geoffroy et al (2010, doi:10.5194/acp-10-4835-2010)

! background aerosol properties for activation when log_predictNc==.true.
! This is mode 1 (nominally, accumulation mode aerosol)
real, public :: aerosol_mode1_radius = 50.e-9 ! aerosol mean size (m)
real, public :: aerosol_mode1_sigmag = 2.0 ! geometric standard deviation
real, public :: aerosol_mode1_number = 300.e6 ! aerosol number mixing ratio (kg-1)

! This is mode 2.  (Default is coarse mode with zero concentration. Could be Aitken as well.)
real, public :: aerosol_mode2_radius = 1.3e-6 ! aerosol mean size (m)
real, public :: aerosol_mode2_sigmag = 2.5 ! geometric standard deviation
real, public :: aerosol_mode2_number = 0.e6 ! aerosol number mixing ratio (kg-1)

! SETTINGS FOR ICE MICROPHYSICS
integer :: nCat =1  ! number of free ice categories
real    :: minDiam=3. !minimum ice diameter factor/fractional difference between newly nucleated and preexisting ice particles=> see the modified icecat_destination code !BG
                      !if minDiam=3., than all particles that larger than 3x preex or smaller than 1/3 of preex go into new category
logical,public  :: frzmodes = .false. !if true, each source of ice has a different ice category, => 7 currently, !BG
logical, public :: do_defeat_limits_on_ice_number_concentration = .false. !BG added
!when set to .true. this can change e.g. extent of anvil clouds, sensitivities to aerosols...

real,public :: MaxTotalIceNumber = 20.e6 ! max ice concentration in #/m3 (default = 10/cm3). !BG could go down to 5/cm3
real,public :: NumCirrusSulf = 20. !20 will effectively limit the ice nucleation for one event to 20 k /L: better limit, as no way to deplete air mass of nsulf in this code
                            !100. !100/cm3 
                            !BG max number of sulphate aerosol used for homog freezing -> taken some reasonable value for upper troposphere
                            !used this number as Liu and Penner, 2005
real,public ::  NumCirrusINP = 2.e-3 !in units per cm3! BG added some reasonable background upper tropospheric dust/het ice nuclei value for cirrus conditions
                              !BG: based on simulated upper tropospheric (200 hPa) values of dust in Tropical Western Pacific for ECHAM-HAM GCM simulations
real,public :: ramp_min = 0.1 !INP ramp parameter ramp_min*NumCirrusINP = at lower temperature level, e.g. 0.2*2.e-3


!BG done a la' Bloss
!more parameters etc to be added, if needed
!(or put here from elswhere) -> set_param.f90 for instance
   logical, public :: typeDiags_ON = .false.! logic variable, for diagnostic hydrometeor/precip rate types
   character(len=10) :: model ="WRF"    ! type of model that is used in P3_MAIN subroutine

   integer :: n_diag_2d = 1       ! the number of 2-D diagnositic variables output in P3_MAIN subroutine
   integer :: n_diag_3d = 1       ! the number of 3-D diagnositic variables output in P3_MAIN subroutine
   !lookup_file_dir = '/global/homes/g/gxlin/MMF_code/wrfv4_p3_release_092217'

   !options for coupling cloud radiative properties to information
   !  from the microphysics
   logical :: douse_reffc = .true., douse_reffi = .true.

   !bloss: Add a flag that enables the conversion of effective radius to
   !   generalized effective size (Fu, 1996) for RRTMG
   logical :: doConvertReffToDge = .true.

   !!!!!!BG do that properly with P3 lookup tables!!!
   real, parameter :: rho_cloud_ice = 500.
   real, parameter :: rho_snow      = 100.

!----------some of these currently not active!!!
!BG added for micro sensitivities:
   logical, public :: no_ice_nucleation    =.false.
   logical, public :: no_hom_ice_nucleation=.false.
   logical, public :: no_het_ice_nucleation=.false.        !no heterogeneous freezing of any kind 
   logical, public :: no_cirrus_mohler_ice_nucleation=.false. !no Mohler et al. heterogeneous freezing on "dust" at T<-35°C
   logical, public :: depo_frz_std = .true.           !use the Meyers deposition freezing parameterization (default in P3)
                                                      !no Mohler et al het freezing, no Liu-Penner Hom freezing
   logical, public :: no_lphom_ice_nucleation=.false. !no Liu-Penner Hom freezing
   logical, public :: do_meyers = .true.             !Meyers depo frz instead of Cooper 86
   logical, public :: use_preexisting_ice = .false.   !use preexisting ice in freezing by Karcher et al., 2006
   logical, public :: do_new_lp_frz =  .false.
   logical, public :: do_mesoscale_variab =  .false.  !add in-cloud mesoscale variability of Sice for new nucleation purpose
   logical, public :: dust_ramp =.false. !make a linear ramp for dust INP concentration, decreasing with temperature

   logical, public :: do_default_icnc_limit=.false.   !use P3's default ICNC upper limit of 500 IC/L
   logical, public :: dolatent=.true. !do latent heating
   !BG
   logical, public :: depofix= .false. !if true limit depo/sublim for T<220 K
   !no sublimation of ice      
   logical, public :: no_evap_i = .false.
   !no evaporation of water
   logical, public :: no_evap_w     = .false.
   !no deposition of vapor on ice
   logical, public :: no_deposition = .false.
   !no ice sedimentation
   logical, public :: do_ice_sedi   = .true.
   logical, public :: incldlath     = .false. !avg in cloud latent heating per layer, horizontally
   logical, public :: incldhomsli   = .false. !avg in cloud static liquid ice energy per layer, horizontally

  !bloss: option to use liquid saturation vapor pressure at all temperatures.   
   logical :: useLiquidESat = .false.

   !bloss: Use esat formulas for liquid and ice from Murphy & Koop (2005, QJRMS)
   !  instead of standard formulas from Flatau et al.
   logical, public :: useMurphyKoopESat = .false.
!---------------

!BG- a la Bloss in M2005 (9Jun2018): Add outputs of microphysical process rates
!   If icemicro==.true., this amount to an extra XXX outputs
!   in the statistics file.  Only 14 for warm cloud microphysics.
   logical, public :: do_output_micro_process_rates = .false.
   !currrently not in code   logical, public :: do_output_micro_process_rates_ud = .false. !separate output for updrafts and downdrafts
   logical, public :: freezing_3dout = .false. !3d output of freezing

   integer :: nmicro_proc
   integer, parameter :: nmicro_process_rates = 10 !24 ! out of 43
   !no need for that: integer, parameter :: nmicro_process_rates_warm = 14
   character(len=8), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_names = (/ &
! liquid-phase microphysical process rates:
!  (all Q process rates in kg kg-1 s-1)
!  (all N process rates in # kg-1)
 !     'qrcon    ', & ! rain condensation
 !     'qcacc    ', & ! cloud droplet accretion by rain
 !     'qcaut    ', & ! cloud droplet autoconversion to rain
  !   'ncacc    '/)!, & ! change in cloud droplet number from accretion by rain
   !  'ncautc   ', & ! change in cloud droplet number from autoconversion
   !  'ncslf    ', & ! change in cloud droplet number from self-collection
   !  'nrslf    ', & ! change in rain number from self-collection
   !  'ncnuc    ', & ! change in cloud droplet number from activation of CCN
  !    'qccon    ', & ! cloud droplet condensation
   !   'qcnuc    ', & ! activation of cloud droplets from CCN
   !   'qrevp    ', & ! rain evaporation
   !   'qcevp    ', & ! cloud droplet evaporation
   !  'nrevp    ', & ! change in rain number from evaporation
   !  'ncautr   ', & ! change in rain number from autoconversion of cloud water
! i!ce-phase microphysical process rates:
!  !(all Q process rates in kg kg-1 s-1)
!  !(all N process rates in # kg-1)
   !   'qccol     ', & ! collection of cloud water by ice
   !   'qwgrth    ', & ! wet growth rate
   !   'qidep     ', & ! vapor deposition
   !   'qrcol     ', & ! collection rain mass by ice
   !   'qinuc     ', & ! deposition/condensation freezing nuc !for T< -15 C and Sice> 5%, not really sure how
   !  'nccol     ', & ! change in cloud droplet number from collection by ice
   !  'nrcol     ', & ! change in rain number from collection by ice
   !  'ninuc     ', & ! change in ice number from deposition/cond-freezing nucleation
   !   'qisub     ', & ! sublimation of ice
   !   'qimlt     ', &!, & ! melting of ice
   !   'qcrfrz    ', &! /) ! freezing of cloud droplets (qinuc + qchetc + qcheti) and rain (qrhti+qrhetc) !BG
     !number rates
     'ninuc     ',& !deposition freezing from Meyers/ in mixed phase
     'ninuc2     ',& ! Mohler deposition freezing at cirrus level, if turned on
     'ninuc3    ',&  ! LP freezing, competition hom, het, preex at cirrus conditions
     !contact frz turned on right now
     'ncheti    ', & ! immersion freezing droplets
     'nrheti    ', & ! immersion freezing rain
     'nimul     ', & ! Hallet Mossop/ice multiplication from rime-splintering (not turned on?)
     'nimlt     ', & ! melting of ice
     'nisub     ', & ! change in ice number from sublimation
     'nislf     ', & ! change in ice number from collection within a category
     'nrchomi   '/) ! homog freezing of cloud droplets and rain, should be last!!!

   !  'ncrfrz    ', & ! number freezing of cloud droplets (qinuc + qchetc + qcheti) and rain (qrhti+qrhetc) !BG
   !  'nchetc    ', & ! contact freezing droplets
   !  'ncheti    ', & ! immersion freezing droplets
   !  'nrhetc    ', & ! contact freezing rain
   !  'nrheti    ', & ! immersion freezing rain
   !  'nrshdr    ', & ! source for rain number from collision of rain/ice above freezing and shedding
   !  'nimul     ', & ! change in Ni, ice multiplication from rime-splintering (not included in the paper)
   !  'ncshdc    ', & ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
   !  'rhorime_c ', & ! density of rime (from cloud)
   !  'rhorime_r '  /)! density of rime (from rain)
 !only if you have multiple ice categories (i.e. nCat>1)
 !nicol ! change of N due to ice-ice collision between categories
 !qicol ! change of q due to ice-ice collision between categories

   character(len=80), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_longnames = (/ &
     !number rates in units output per timestep per m-3
     'ninuc depo frz mixed phase',& !deposition freezing from Meyers/ in mixed phase
     'ninuc2 Mohler depo frz at cirrus     ',& ! Mohler deposition freezing at cirrus level, if turned on
     'ninuc3 LP cirrus freezing    ',&  ! LP freezing, competition hom, het, preex at cirrus conditions !contact frz turned on right now
     'ncheti immersion frz droplets    ', & ! immersion freezing droplets
     'nrheti immersion frz rain   ', & ! immersion freezing rain
     'nimul rime-splintering     ', & ! Hallet Mossop/ice multiplication from rime-splintering (not turned on?)
     'nimlt melting    ', & ! melting of ice
     'nisub sublimation    ', & ! change in ice number from sublimation
     'nislf self-collection    ', & ! change in ice number from collection within a category
     'nrchomi hom frz of droplets and rain  '/) ! homog freezing of cloud droplets and rain, should be last!!!
 end module micro_params
