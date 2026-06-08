subroutine albedo_ocn(lchnk,ncol, landmask, coszrs, asdir, aldir, asdif, aldif )
  !-----------------------------------------------------------------------
  ! Computes surface albedos over ocean for Slab Ocean Model (SOM)

  ! and the surface (added by Marat Khairoutdinov)

  ! Two spectral surface albedos for direct (dir) and diffuse (dif)
  ! incident radiation are calculated. The spectral intervals are:
  !   s (shortwave)  = 0.2-0.7 micro-meters
  !   l (longwave)   = 0.7-5.0 micro-meters
  !
  ! Uses knowledge of surface type to specify albedo, as follows:
  !
  ! Ocean           Uses solar zenith angle to compute albedo for direct
  !                 radiation; diffuse radiation values constant; albedo
  !                 independent of spectral interval and other physical
  !                 factors such as ocean surface wind speed.
  !
  ! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
  ! Approximation for Solar Radiation in the NCAR Community Climate Model,
  ! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
  !
  use shr_kind_mod, only: r4 => shr_kind_r4
  use ppgrid
  implicit none
  !------------------------------Arguments--------------------------------
  !
  ! Input arguments

  !
  integer lchnk,ncol

  real(r4) coszrs    ! Cosine solar zenith angle
  real snowm ! snow mass (mm)
  real landmask	! landmask: 0 - ocean, 1 - land
  real icemask	! landmask: 1 - ice, 0 - no ice
  real albedovis_v, albedovis_s  ! visible albedo of vegetation and soil
  real albedonir_v, albedonir_s  ! near IR albedo of land
  real soilw  ! soil wetness
  real phi_1, phi_2, LAI
  !
  ! Output arguments
  !
  real(r4) asdir     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r4) aldir     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r4) asdif     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r4) aldif     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
  !
  !---------------------------Local variables-----------------------------
  !
  real(r4), parameter :: adif = 0.06
  real ka, explai, wetfactor
  !
  !
  ! Initialize all ocean surface albedos to zero
  !
  asdir = 0.
  aldir = 0.
  asdif = 0.
  aldif = 0.

  if (coszrs.le.0.) return

  if (landmask.eq.0.) then
     ! Ice-free ocean albedos function of solar zenith angle only, and
     ! independent of spectral interval (Briegleb et al 1986):
     !   aldir = .026/(coszrs**1.7 + .065) + .15*(coszrs - 0.10)*(coszrs - 0.50)*(coszrs - 1.00)
     ! From IFS Doc: Taylor et al. (1996)
        aldir = .037/(1.1*coszrs**1.4 + .15)
        asdir = aldir
        aldif = adif
        asdif = adif
  else ! seaice
        aldir = 0.33
        asdir = 0.78
        aldif = 0.33
        asdif = 0.78
  endif
      
  return
end subroutine albedo_ocn

