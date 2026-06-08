module micro_params

use grid, only: nzm

implicit none

!  Microphysics stuff:

! Densities of hydrometeors

real, parameter :: rhor = 1000. ! Density of water, kg/m3
real, parameter :: rhog = 400.  ! Density of graupel, kg/m3
!real, parameter :: rhog = 917.  ! hail - Lin 1983    

! Temperatures limits for various hydrometeors

real, parameter :: thfreez = 233.16   ! homogeneous freezing temperature
real, parameter :: tprmin = 268.16    ! Minimum temperature for rain, K
real, parameter :: tprmax = 273.16    ! Maximum temperature for snow+graupel, K

! Terminal velocity coefficients

real, parameter :: a_rain = 842. ! Coeff.for rain term vel 
real, parameter :: b_rain = 0.8  ! Fall speed exponent for rain
!real, parameter :: a_grau = 40.7! Krueger (1994) ! Coef. for graupel term vel
real, parameter :: a_grau = 94.5 ! Lin (1983) (rhog=400)
!real, parameter :: a_grau = 127.94! Lin (1983) (rhog=917)
real, parameter :: b_grau = 0.5  ! Fall speed exponent for graupel

! Autoconversion

real, parameter :: Nc0 = 200.  ! cm-3, droplet concentration for autoconversion
real :: qcw0 = 1.e-3      ! Threshold for water autoconversion, g/g
real :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
real, parameter :: alphaelq = 1.e-3  ! autoconversion of cloud water rate coef
real, parameter :: betaelq = 1.e-3   ! autoconversion of cloud ice rate coef

real:: icefall_fudge=0.3 ! ice fall terminal velocity's fudge factor (used in ice_fall.f90)

! Accretion

real, parameter :: erccoef = 1.0   ! Rain/Cloud water collection efficiency
real, parameter :: eiccoef = 1.0   ! Snow/Cloud water collection efficiency
real, parameter :: egccoef = 1.0   ! Graupel/Cloud water collection efficiency
real, parameter :: egicoef = 0.1   ! Graupel/Cloud ice collection efficiency

! Interseption parameters for exponential size spectra

real, parameter :: nzeror = 8.e6   ! Intercept coeff. for rain  
real, parameter :: nzerog = 4.e6   ! Intersept coeff. for graupel
!real, parameter :: nzerog = 4.e4   ! hail - Lin 1993 

! coefficient relating ice/snow diameter to mass
real, parameter :: a_ice = 0.095
real, parameter :: b_ice = 2.1
real, parameter :: d_ice_min = 5. ! minimum ice/snow diameter (micron)

real accric(nzm),accrrc(nzm),coefice(nzm),accrgc(nzm)
real evapi1(nzm),evapi2(nzm),evapr1(nzm),evapr2(nzm),evapg1(nzm),evapg2(nzm)

real niz(nzm) ! ice crystal concentration as a function of temperature
real av_ice(nzm), bv_ice(nzm) ! ice terminal velocity coeffs
real gami(nzm) ! gamma parameter in gamma-distribution for ice/snow
real qi_min(nzm) ! minimum allowed ice/snow mixing ratio (based on d_ice_min)
            

real, parameter :: qp_threshold = 1.e-12 ! minimal rain/snow water content
real a_pr, a_gr 

logical:: docloudfall=.false. ! dummy parameter - don't change
logical:: dowarmcloud=.true. ! dummy parameter - don't change

end module micro_params
