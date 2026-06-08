module micro_params

!bloss: move indices of different microphysical species here, so that they
!  can be accessed in both microphysics.f90 and module_mp_p3.f90.

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, incl, inr, iqr
integer, dimension(:), allocatable :: inci, iqit, iqir, iqib

! these non-advected fields hold the value of qv and tabs from the last time step.
integer :: iqv_old, itabs_old 

! SETTINGS FOR WARM RAIN MICROPHYSICS (AND DEFAULT VALUES)

! choice of warm rain microphysics scheme
!   = 1 for Seifert-Beheng (2001)
!   = 2 for Beheng (1994)
!   = 3 for Khairoutdinov-Kogan (2000) <-- DEFAULT
integer, public :: iWarmRainScheme = 3 

! control properties of cloud drop size distribution.
logical, public :: log_predictNc = .false. ! if true, predict droplet number
real, public :: Nc0 = 100. ! specified cloud droplet number conc (CDNC, #/cm3)

logical, public :: dofix_pgam = .false. ! option to fix value of exponent in cloud water gamma distn
real, public ::    pgam_fixed = 10.3 ! Geoffroy et al (2010, doi:10.5194/acp-10-4835-2010)

logical, public :: doSAMsvp = .true. ! Use svp calculation from SRC/sat.f90

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

! choice of ice nucleation scheme
!   = 1 for P3 Default (Cooper/Myers for T < -15C)
!   = 2 for Cooper/Myers for -15C < T < -37C, Shi et al (2015) below -37C
integer, public :: iIceNucleationScheme = 1

! general options for ice formation
logical, public :: do_meyers = .true.             !Meyers depo frz instead of Cooper 86
real,public :: MaxTotalIceNumber = 20.e6 ! max ice concentration in #/m3 (default = 10/cm3). !BG could go down to 5/cm3
logical, public :: do_defeat_limits_on_ice_number_concentration = .false. !BG added
!when set to .true. this can change e.g. extent of anvil clouds, sensitivities to aerosols...

!  If using iIceNucleationScheme==2, these further options are available:
!    Note that do_new_lp_frz==T and use_preexisting_ice==T by default.
logical, public :: do_new_lp_frz = .true. ! Freezing at cold temperatures follows Liu and Penner as implemented in Shi et al
logical, public :: use_preexisting_ice = .true.   !use preexisting ice in freezing by Karcher et al., 2006
logical, public :: do_mesoscale_variab =  .false.  !add in-cloud mesoscale variability of Sice for new nucleation purpose

logical, public :: enable_cirrus_mohler_ice_nucleation=.false. !Sensitivity study to use Mohler et al (2006) freezing.
logical, public :: enable_lphom_ice_nucleation=.false. !Sensitivity study to use only homogeneous freezing from Liu & Penner.

real,public :: NumCirrusSulf = 20. ! Max number of sulphate aerosol used for homog freezing [#/cm3]
real,public :: NumCirrusINP = 2.e-3 ! Max Number of dust/heterogeneous INP at cirrus temperatures [#/cm3]
                                    !BG: based on upper tropospheric dust in ECHAM-HAM over West Pacific warm pool
logical, public :: dust_ramp =.false. ! option for dust concentration to ramp down for T < -40 C
real,public :: ramp_min = 0.1  ! when dust_ramp==T, dust concentration below -70C is ramp_min*NumCirrusINP

real    :: minDiam=3. !minimum ice diameter factor/fractional difference between newly nucleated and preexisting ice particles=> see the modified icecat_destination code !BG
                      !if minDiam=3., than all particles that larger than 3x preex or smaller than 1/3 of preex go into new category
logical,public  :: frzmodes = .false. !if true, each source of ice has a different ice category, => 7 currently, !BG

!BG done a la' Bloss
!more parameters etc to be added, if needed
!(or put here from elswhere) -> set_param.f90 for instance
   logical, public :: typeDiags_ON = .true.! logic variable, for diagnostic hydrometeor/precip rate types
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

   !bloss: Option to split ice optics into cloud ice and snow
   !  and sum their IWC and Projected area across ice categories if nCat>1.
   logical :: doSplitIceOptics = .false. !.true.
   integer :: iSplitIceOption = 3 ! Different definitions of snow.
   
   !!!!!!BG do that properly with P3 lookup tables!!!
   real, parameter :: rho_cloud_ice = 500.
   real, parameter :: rho_snow      = 100.

   logical, public :: depo_frz_std = .true.           !use the Meyers deposition freezing parameterization (default in P3)

   logical, public :: do_default_icnc_limit=.false.   !use P3's default ICNC upper limit of 500 IC/L

   !----------some of these currently not active!!!
!BG added for micro sensitivities:
   logical, public :: no_ice_nucleation    =.false.
   logical, public :: no_hom_ice_nucleation=.false.
   logical, public :: no_het_ice_nucleation=.false.        !no heterogeneous freezing of any kind 

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
   logical, public :: do_output_3d_micro_process_rates = .false.
   !currrently not in code   logical, public :: do_output_micro_process_rates_ud = .false. !separate output for updrafts and downdrafts
   logical, public :: freezing_3dout = .false. !3d output of freezing

   integer :: nmicro_proc
   integer, parameter :: nmicro_process_rates = 17 !24 ! out of 43
   !no need for that: integer, parameter :: nmicro_process_rates_warm = 14
   character(len=8), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_names = (/ &
! liquid-phase microphysical process rates:
!  (all Q process rates in kg kg-1 s-1)
!  (all N process rates in # kg-1)
     'qidep     ', & ! vapor deposition
     'qisub     ', & ! sublimation of ice
     'qimlt     ', & ! melting of ice
     'qcrfrz    ', & ! freezing of cloud droplets (qinuc + qchetc + qcheti) and rain (qrhti+qrhetc) !BG
     'qinuc     ', & ! ice mass generated by mixed-phase/Cooper/Myers ice nucleation
     'qinuc2     ', & ! ice mass generated by Mohler (het) freezing of cirrus
     'qinuc3     ', & ! ice mass generated by heterogeneous ice nucleation 
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
 !only if you have multiple ice categories (i.e. nCat>1)
 !nicol ! change of N due to ice-ice collision between categories
 !qicol ! change of q due to ice-ice collision between categories

   character(len=80), dimension(nmicro_process_rates), parameter, public :: &
        micro_process_rate_longnames = (/ &
! liquid-phase microphysical process rates:
!  (all Q process rates in kg kg-1 s-1)
!  (all N process rates in # kg-1)
     'qidep   , vapor deposition                                                      ', &
     'qisub   , sublimation of ice                                                    ', &
     'qimlt   , melting of ice                                                        ', &
     'qcrfrz  , freezing of cloud droplets and rain drops                             ', &
     'qinuc   , ice mass generated by mixed-phase/Cooper/Myers ice nucleation         ', & ! 
     'qinuc2  , ice mass generated by Mohler (het) freezing of cirrus                 ', & ! 
     'qinuc3  , ice mass generated by hetero/homogeneous ice nucleation (Liu & Penner)', & ! 
     'ninuc   , depo frz mixed phase',& !deposition freezing from Meyers/ in mixed phase
     'ninuc2  , Mohler depo frz at cirrus     ',& ! Mohler deposition freezing at cirrus level, if turned on
     'ninuc3  , LP cirrus freezing    ',&  ! LP freezing, competition hom, het, preex at cirrus conditions !contact frz turned on right now
     'ncheti  , immersion frz droplets    ', & ! immersion freezing droplets
     'nrheti  , immersion frz rain   ', & ! immersion freezing rain
     'nimul   , rime-splintering     ', & ! Hallet Mossop/ice multiplication from rime-splintering (not turned on?)
     'nimlt   , melting    ', & ! melting of ice
     'nisub   , sublimation    ', & ! change in ice number from sublimation
     'nislf   , self-collection    ', & ! change in ice number from collection within a category
     'nrchomi , hom frz of droplets and rain  '/) ! homog freezing of cloud droplets and rain, should be last!!!
 !only if you have multiple ice categories (i.e. nCat>1)
 !nicol ! change of N due to ice-ice collision between categories
 !qicol ! change of q due to ice-ice collision between categories

   logical:: dowarmcloud=.true. ! dummy parameter - don't change
   logical:: docloudfall=.false. ! dummy parameter - don't change
   logical:: doeiscld = .false. ! set shallow cloud fraction for radiation based on EIS
   logical:: do_scale_dependence_of_autoconv = .false. ! autoconversion in KK2000 depends on grid spacing
   real deltax ! dx*mu - actual grid spacing in x.
   real:: icefall_fudge=1.0 ! ice fall terminal velocity's fudge factor

   ! Ice size dependent fallspeed adjustment.  Piecewise linear based on
   !   qice = Nice*(4/3)*pi*rho_ice*Rice^3
   ! where qice and Nice are mass/number mixing ratios and Rice is ice mass radius.
   !   Rice = ( (qice/Nice) * ( 3 / (4*pi*rho_ice) ) ) ^ (1/3)
   logical :: doSizeDependentIceFallFudge = .false.
   real:: icefall_fudge_small=1. ! terminal velocity fudge factor for Rice<=25 um
   real:: icefall_fudge_medium=1. ! terminal velocity fudge factor for Rice = 75 um 
   real:: icefall_fudge_large=1. ! terminal velocity fudge factor for Rice >=150 um

   real :: IceNucleiRadius = 1.e-6 ! Radius of Newly Nucleated Ice Crystals, m

 end module micro_params
