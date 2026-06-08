module consts

implicit none

!   Constants:

real(8), parameter :: cp = 1004.64d0             ! Specific heat of air, J/kg/K
real(8), parameter :: cpv = 1870.0             ! Specific heat of water vapor, J/kg/K
real(8), parameter :: ggr = 9.79764d0          ! Gravity acceleration, m/s2
real(8), parameter :: lcond = 2.501d+06     ! Latent heat of condensation, J/kg
real(8), parameter :: lfus = 0.337d+06      ! Latent heat of fusion, J/kg
real(8), parameter :: lsub = 2.834d+06      ! Latent heat of sublimation, J/kg
real(8), parameter :: rv = 461.5d0              ! Gas constant for water vapor, J/kg/K
real(8), parameter :: rgas = 287.04d0            ! Gas constant for dry air, J/kg/K
real(8), parameter :: diffelq = 2.21d-05     ! Diffusivity of water vapor, m2/s
real(8), parameter :: therco = 2.40d-02      ! Thermal conductivity of air, J/m/s/K
real(8), parameter :: muelq = 1.717d-05      ! Dynamic viscosity of air
!real(8), parameter :: cpw= 4187.d0            ! specific heat of fresh water
real(8), parameter :: cpw = 3991.86795711963d0 ! specific hat of seawater

real(8), parameter :: fac_cond = lcond/cp
real(8), parameter :: fac_fus = lfus/cp
real(8), parameter :: fac_sub = lsub/cp

real(8), parameter :: rad_earth = 6371229.d0   ! radius of the Earth (to make 40000 km circumference)
real(8), parameter :: pi = acos(-1._8) !3.141592653589793d0
real(8), parameter :: sigmaSB = 5.670373d-8 ! Stefan-Boltzmann constant, Wm-2K-4
real(8), parameter :: emis_water = 0.98



end module consts
