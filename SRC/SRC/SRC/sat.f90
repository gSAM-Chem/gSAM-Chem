! Saturation vapor pressure and mixing ratio. 
! Based on ECMWF IFS Documentation - Cy47r3
!   Part IV: Physical Processes, section 12.7

! ===========================
!  Saturation vapor pressures

real function esatw(t) ! returns hPa
  implicit none
  real t	! temperature (K)
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 17.502 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =  32.19 ! T1 in K
  ! equation 12.9 in IFS documentation Cy47r3
  esatw = e0 * exp( a1 * (T-T0) / (T-T1) )
end function esatw
        
real function esati(t) ! returns hPa
  implicit none
  real t	! temperature (K)
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 22.587 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =   -0.7 ! T1 in K
  ! equation 12.10 in IFS documentation Cy47r3
  esati = e0 * exp( a1 * (T-T0) / (T-T1) )
end function esati
        
! ===============================
!  Saturation vapor mixing ratios

real function qsatw(t,p) ! returns kg H2O / kg dry air
  implicit none
  real t	! temperature (K)
  real p	! pressure    (mb)
  real esat
  real, external :: esatw
  esat = esatw(t)
  qsatw = 0.622 * esat/max(esat,p-esat)
end function qsatw
        
real function qsati(t,p) ! returns kg H2O / kg dry air
  implicit none
  real t	! temperature (K)
  real p	! pressure    (mb)
  real esat
  real, external :: esati
  esat=esati(t)
  qsati=0.622 * esat/max(esat,p-esat)
end function qsati
        
! ===============================
!  d(esat)/dT

! This formulation is inconsistent with definition of esatw - MK
!real function dtesatw(t)
!  use consts, only: lcond, rv
!  implicit none
!  real t	! temperature (K)
!  real esat	! saturation vapor pressure (hPa)
!  real, external :: esatw
!  esat = esatw(t)
!  ! equation 12.8 in IFS documentation Cy47r3
!  dtesatw = lcond*esat / (rv*t*t)
!end function dtesatw
        
! this formulation is found by taking a true derivative of esatw - MK
real function dtesatw(t)
  implicit none
  real t        ! temperature (K)
  real esat     ! saturation vapor pressure (hPa)
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 17.502 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =  32.19 ! T1 in K
  real, external :: esatw
  esat = esatw(t)
  dtesatw = esat * a1 * (T0-T1)/((T-T1)*(T-T1))
end function dtesatw

        
!real function dtesati(t)
!  use consts, only: lsub, rv
!  implicit none
!  real t	! temperature (K)
!  real esat	! saturation vapor pressure (hPa)
!  real, external :: esati
!  
!  esat = esati(t)
!
!  ! equation 12.8 in IFS documentation Cy47r3
!  dtesati = lsub*esat / (rv*t*t)
!end function dtesati
        
! this formulation is found by taking a true derivative of esati - MK
real function dtesati(t)
  implicit none
  real t        ! temperature (K)
  real esat     ! saturation vapor pressure (hPa)
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 22.587 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =   -0.7 ! T1 in K
  real, external :: esati
  esat = esati(t)
  dtesati = esat * a1 * (T0-T1)/((T-T1)*(T-T1))
end function dtesati

! ===============================
!  d(qsat)/dT

! this implementation produces large error for warm T, warmer than 310K,
! although it may bot affect the convergence of condensation compytation much.  -MK
!real function dtqsatw(t,p)
!  implicit none
!  real t	! temperature (K)
!  real p	! pressure    (mb)
!  real, external :: dtesatw
!  dtqsatw=0.622*dtesatw(t)/p
!end function dtqsatw

! this formulation is found by taking a true derivative of dqsatw - MK
real function dtqsatw(t,p)
  implicit none
  real t        ! temperature (K)
  real p        ! pressure    (mb)
  real, external :: esatw
  real esat, dtesatw
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 17.502 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =  32.19 ! T1 in K
  esat = esatw(t)
  dtesatw = esat * a1 * (T0-T1)/((T-T1)*(T-T1))
  dtqsatw=0.622*dtesatw/(p-esat)*(1.+esat/(p-esat))
end function dtqsatw

        
!real function dtqsati(t,p)
!  implicit none
!  real t	! temperature (K)
!  real p	! pressure    (mb)
!  real, external :: dtesati
!  dtqsati=0.622*dtesati(t)/p
!end function dtqsati
      
! this formulation is found by taking a true derivative of dqsati - MK
real function dtqsati(t,p)
  implicit none
  real t        ! temperature (K)
  real p        ! pressure    (mb)
  real, external :: esati
  real esat,dtesati
  real, parameter :: e0 = 6.1121 ! esat(T0) in hPa
  real, parameter :: a1 = 22.587 ! coeficient, dimensionless
  real, parameter :: T0 = 273.16 ! T0 in K
  real, parameter :: T1 =   -0.7 ! T1 in K
  esat = esati(t)
  dtesati = esat * a1 * (T0-T1)/((T-T1)*(T-T1))
  dtqsati=0.622*dtesati/(p-esat)*(1.+esat/(p-esat))
end function dtqsati

! ===============================
!  dew point

real function tdew(t,q,p)
  ! compute the dew-point temperqature given temperature (K), 
  ! specific humidity (kg/kg), and pressure (mb)
  ! MK 2021
  implicit none
  real, intent(in) :: t
  real, intent(in) :: q
  real, intent(in) :: p
  real qsatw, dtqsatw
  integer niter, niter_max
  real t1
  t1 = t ! initial guess
  niter_max = 3
  tdew =1000.
  niter=0
  do while(abs(tdew-t1).gt.0.1.and.niter.lt.niter_max)
    niter = niter+1
    tdew = t1 - (qsatw(t1,p)-q)/dtqsatw(t1,p)
    t1 = tdew
  end do
end function tdew

