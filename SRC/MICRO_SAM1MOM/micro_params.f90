module micro_params

use grid, only: nx,ny,nzm

implicit none

!  Microphysics stuff:

! Densities of hydrometeors

real, parameter :: rhor = 1000. ! Density of water, kg/m3
real, parameter :: rhos = 100.  ! Density of snow, kg/m3
real, parameter :: rhog = 400.  ! Density of graupel, kg/m3
!real, parameter :: rhog = 917.  ! hail - Lin 1983    

! Default temperatures limits for various hydrometeors (changed by namelist MICRO_SAM1MOM)

real :: tbgmin = 253.16    ! Minimum temperature for cloud water., K
real :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
real :: tprmin = 268.16    ! Minimum temperature for rain, K
real :: tprmax = 283.16    ! Maximum temperature for snow+graupel, K
real :: tgrmin = 223.16    ! Minimum temperature for graupel, K
real :: tgrmax = 283.16    ! Maximum temperature for graupel, K

! Terminal velocity coefficients

real, parameter :: a_rain = 842. ! Coeff.for rain term vel 
real, parameter :: b_rain = 0.8  ! Fall speed exponent for rain
real, parameter :: a_snow = 4.84 ! Coeff.for snow term vel
real, parameter :: b_snow = 0.25 ! Fall speed exponent for snow
!real, parameter :: a_grau = 40.7! Krueger (1994) ! Coef. for graupel term vel
real, parameter :: a_grau = 94.5 ! Lin (1983) (rhog=400)
!real, parameter :: a_grau = 127.94! Lin (1983) (rhog=917)
real, parameter :: b_grau = 0.5  ! Fall speed exponent for graupel

! Autoconversion

logical:: doKKauto = .false. ! switch to use KK (2000) autoconversion formula 
logical :: do_scale_dependence_of_autoconv = .false. ! make resolution-dependent autoconversion rate
logical:: donomicro = .false. ! supress cloud microphysics 
logical:: dowarmcloud = .false. ! if true, no ice microphysics 
logical:: donograupel = .false. ! if true, no graupel in microphysics 
logical:: doeiscld = .false.    ! set shallow cloud fraction for radiation based on EIS
logical:: docloudfall = .false. ! compute sedimentation of cloud wat
logical:: dovticeifs = .false. ! use constant ice terminal velocity similar to IFS
logical:: do_dependence_N0r_on_rainrate = .false.  ! change intersept parameter N0 for rain
                                                   ! depending on local rain - MK 2025
logical:: doprecipdiss = .false.    ! add dissipative heating due to precipitation - MK 2005

real:: auto_fudge=1.0 !  fudge factor to change autoconversion in KK &
                      !  (need do_scale_dependence_of_autoconv = .false.)
real:: icefall_fudge=1.0 ! ice fall terminal velocity's fudge factor (used in ice_fall.f90)
real:: Nc_land = 300.    ! drop concentration over land (used when doKKauto = .true.)
real:: Nc_ocn = 50.      ! drop concentration over ocean (used when doKKauto = .true.)
real:: sigmag = 1.5      ! relative cloud drop dispersion (for sedimentation)
real :: qcw0 = 1.e-3      ! Threshold for water autoconversion, g/g  
real :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
real :: alphaelq = 1.e-3  ! autoconversion of cloud water rate coef
real :: betaelq = 1.e-3   ! autoconversion of cloud ice rate coef

! Accretion

logical:: doKKaccr = .false. ! switch to use KK (2000) accretion formula for rain
real, parameter :: erccoef = 1.0   ! Rain/Cloud water collection efficiency
real, parameter :: esccoef = 1.0   ! Snow/Cloud water collection efficiency
real, parameter :: esicoef = 0.1   ! Snow/cloud ice collection efficiency
real, parameter :: egccoef = 1.0   ! Graupel/Cloud water collection efficiency
real, parameter :: egicoef = 0.1   ! Graupel/Cloud ice collection efficiency

! Interseption parameters for exponential size spectra

real, parameter :: nzeror = 8.e6   ! Intercept coeff. for rain  
real, parameter :: nzeros = 3.e6   ! Intersept coeff. for snow
real, parameter :: nzerog = 4.e6   ! Intersept coeff. for graupel
!real, parameter :: nzerog = 4.e4   ! hail - Lin 1993 

real, parameter :: qp_threshold = 1.e-12 ! minimal rain/snow water content


! Misc. microphysics variables

real*4 gam3       ! Gamma function of 3
real*4 gams1      ! Gamma function of (3 + b_snow)
real*4 gams2      ! Gamma function of (5 + b_snow)/2
real*4 gams3      ! Gamma function of (4 + b_snow)
real*4 gamg1      ! Gamma function of (3 + b_grau)
real*4 gamg2      ! Gamma function of (5 + b_grau)/2
real*4 gamg3      ! Gamma function of (4 + b_grau)
real*4 gamr1      ! Gamma function of (3 + b_rain)
real*4 gamr2      ! Gamma function of (5 + b_rain)/2
real*4 gamr3      ! Gamma function of (4 + b_rain)
      
real accrsc(nzm),accrsi(nzm),accrrc(nzm)
real accrgc(nzm),accrgi(nzm)
real evaps1(nx,ny,nzm),evaps2(nx,ny,nzm),evapr1(nx,ny,nzm),evapr2(nx,ny,nzm)
real evapg1(nx,ny,nzm),evapg2(nx,ny,nzm)
real:: nzeror_factor(nx,ny,nzm) = 1. ! intersept parameter adjustment factor - MK 2025
            
real a_bg, a_pr, a_gr 


!bloss: These parameters are used in MICRO_M2005 cloud optics routines.
real, parameter :: rho_snow = rhos, rho_water = rhor, rho_cloud_ice = 917.

end module micro_params
