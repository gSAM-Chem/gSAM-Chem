module cam_rad_parameterizations
  use parkind, only: kind_rb !bloss(120709): use RRTM real kind here.
  use shr_orb_mod
  !
  ! CAM 3.0 parameterizations related to radiation
  !
  implicit none
  private
  
  real(kind_rb), parameter :: tmelt = 273.16  ! specify melting temperature
  !       Tabulated values of re(T) in the temperature interval
  !       180 K -- 274 K; hexagonal columns assumed
  !
  real(kind_rb), parameter :: iceSizeTableMinTemp = 180. 
  real(kind_rb), dimension(95), parameter :: &
    retab = (/ 5.92779, 6.26422, 6.61973, 6.99539, 7.39234,	&
        7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,	&
        10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
        15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,	&
        20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,	&
        27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
        31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
        34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,	&
        38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,	&
        42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,	&
        50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,	&
        65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,	&
        93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
        124.954, 130.630, 136.457, 142.446, 148.608, 154.956,	&
        161.503, 168.262, 175.248, 182.473, 189.952, 197.699,	&
        205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)


  public :: computeRe_Liquid, computeRe_Ice, albedo, albedo_slm

contains
  !-----------------------------------------------------------------------
  elemental real(kind_rb) function computeRe_Liquid(temperature, landfrac, icefrac, snowh) &
     result(rel)
    real(kind_rb),           intent(in) :: temperature, landfrac
    real(kind_rb), optional, intent(in) :: icefrac, snowh  ! Snow depth over land, water equivalent (m)

    real(kind_rb), parameter ::  rliqland  =  8.0, & ! liquid drop size if over land
                        rliqocean = 14.0, & ! liquid drop size if over ocean
                        rliqice   = 14.0    ! liquid drop size if over sea ice
 
    ! jrm Reworked effective radius algorithm
    ! Start with temperature-dependent value appropriate for continental air
    rel = rliqland + (rliqocean - rliqland) * min(1.0_kind_rb, max(0.0_kind_rb, (tmelt - temperature) * 0.05))
    
    if(present(snowh)) & ! Modify for snow depth over land
      rel = rel + (rliqocean - rel) * min(1.0_kind_rb, max(0.0_kind_rb, snowh*10.))
      
    ! Ramp between polluted value over land to clean value over ocean.
    rel = rel + (rliqocean-rel) * min(1.0_kind_rb, max(0.0_kind_rb, 1.0_kind_rb - landfrac))
    
    if(present(icefrac)) & ! Ramp between the resultant value and a sea ice value in the presence of ice.
      rel = rel + (rliqice-rel) * min(1.0_kind_rb, max(0.0_kind_rb, icefrac))

  end function computeRe_Liquid
  !-----------------------------------------------------------------------
  elemental real(kind_rb) function computeRe_Ice(temperature) result(rei) 
    real(kind_rb),           intent(in) :: temperature
  
    real(kind_rb)    :: fraction
    integer :: index
    !
    !
    index = int(temperature - (iceSizeTableMinTemp - 1.))
    index = min(max(index, 1), 94)
    fraction = temperature - int(temperature)
    rei = retab(index) * (1. - fraction) + retab(index+1) * fraction
  end function computeRe_Ice
  !-----------------------------------------------------------------------
  subroutine albedo_slm(nx, landmask, icemask, coszrs, snowm, albedovis_v, albedonir_v,  &
                    albedovis_s, albedonir_s, soilw, phi_1, phi_2, LAI, IMPERV, &
                    asdir, aldir, asdif, aldif)
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean, ice and land 
  
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    use slm_params, only : albvis_snow, albnir_snow, albvis_ice, albnir_ice
    integer, intent(in) :: nx ! number of grid points
    real(kind_rb), dimension(:), intent( in) :: coszrs   ! Cosine of solar zenith angle
    real(kind_rb), dimension(:), intent(out) :: asdir, & ! Srf alb for direct rad   0.2-0.7 micro-ms
                                       aldir, & ! Srf alb for direct rad   0.7-5.0 micro-ms
                                       asdif, & ! Srf alb for diffuse rad  0.2-0.7 micro-ms
                                       aldif    ! Srf alb for diffuse rad  0.7-5.0 micro-ms
    real, dimension(:), intent( in) :: snowm ! snow mass (mm)
    real, dimension(:), intent( in) :: landmask ! landmask: 0 - ocean, 1 - land
    real, dimension(:), intent( in) :: icemask  ! landmask: 1 - ice, 0 - no ice
    real, dimension(:), intent( in) :: albedovis_v, albedovis_s  ! visible albedo of vegetation and soil
    real, dimension(:), intent( in) :: albedonir_v, albedonir_s  ! near IR albedo of land
    real, dimension(:), intent( in) :: soilw  ! soil wetness
    real, dimension(:), intent( in) :: phi_1, phi_2, LAI, IMPERV
    real(kind_rb), parameter :: adif = 0.06
    !real(kind_rb), parameter :: adif = 0.07 !AAW 12/5/17 RCEMIP albedo value
    integer i
    real explai, explai1, wetfactor

  !
  ! Initialize all ocean surface albedos to zero
  !
  asdir(:) = 0.
  aldir(:) = 0.
  asdif(:) = 0.
  aldif(:) = 0.

  do i = 1,nx

  if (coszrs(i).le.0.) then
    aldir(i) = 0.
    asdir(i) = 0.
    aldif(i) = 0.
    asdif(i) = 0.
    cycle
  end if

  if (landmask(i).eq.0.) then
    if(icemask(i).eq.0) then
     ! Ice-free ocean albedos function of solar zenith angle only, and
     ! independent of spectral interval (Briegleb et al 1986):
     !   aldir(i) = .026/(coszrs(i)**1.7 + .065) +  &
     !               .15*(coszrs(i) - 0.10)*(coszrs(i) - 0.50)*(coszrs(i) - 1.00)
     ! From IFS Doc: Taylor et al. (1996)
        aldir(i) = .037/(1.1*coszrs(i)**1.4 + .15)
        asdir(i) = aldir(i)
        aldif(i) = adif
        asdif(i) = adif
    else
     ! albedo of sea-ice
     ! based on values close to Ebert and Curry (JGR 1993)
        aldir(i) = 0.33
        asdir(i) = 0.78
        aldif(i) = 0.33
        asdif(i) = 0.78
    endif

  else ! land

    if(icemask(i).eq.0) then
     if(snowm(i).eq.0.) then
        wetfactor = 1.-0.5*soilw(i)
        explai = exp(-(phi_1(i)/max(0.01,coszrs(i))+phi_2(i))*LAI(i))
        asdir(i) = albedovis_v(i)*(1.-explai)+albedovis_s(i)*wetfactor*explai
        aldir(i) = albedonir_v(i)*(1.-explai)+albedonir_s(i)*wetfactor*explai
        asdir(i) = (1.-IMPERV(i))*asdir(i) +IMPERV(i)*0.2*explai
        aldir(i) = (1.-IMPERV(i))*aldir(i) +IMPERV(i)*0.4*explai
        explai1 =explai
        ! assume diffuse albedo is just albedo for coszrs=1 (Marat K):
        explai = exp(-(phi_1(i)+phi_2(i))*LAI(i))
        asdif(i) = albedovis_v(i)*(1.-explai)+albedovis_s(i)*wetfactor*explai
        aldif(i) = albedonir_v(i)*(1.-explai)+albedonir_s(i)*wetfactor*explai
        asdif(i) = (1.-IMPERV(i))*asdif(i) +IMPERV(i)*0.2*explai
        aldif(i) = (1.-IMPERV(i))*aldif(i) +IMPERV(i)*0.4*explai
     else ! albedo of vegetated land covered with snow
        explai = exp(-(phi_1(i)/max(0.01,coszrs(i))+phi_2(i))*LAI(i))
! albedos of snow are for dry snow in Ebert abd Curry (1993)
        asdir(i) = albedovis_v(i)*(1.-explai)+albvis_snow*explai
        aldir(i) = albedonir_v(i)*(1.-explai)+albnir_snow*explai
        ! assume diffuse albedo is just albedo for coszrs=1 (Marat K):
        explai = exp(-(phi_1(i)+phi_2(i))*LAI(i))
        asdif(i) = albedovis_v(i)*(1.-explai)+albvis_snow*explai
        aldif(i) = albedonir_v(i)*(1.-explai)+albnir_snow*explai
     end if
!     print*,LAI(i),coszrs(i),explai1, explai,asdir(i),aldir(i),asdif(i),aldif(i)
!     stop
    else
     ! albedo of glacier
        aldir(i) = albnir_ice
        asdir(i) = albvis_ice
        aldif(i) = albnir_ice
        asdif(i) = albvis_ice
    endif

  endif

  end do ! i

  end subroutine albedo_slm

  !-----------------------------------------------------------------------
  subroutine albedo(ocean, coszrs, ts, asdir, aldir, asdif, aldif)
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean
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

    logical,            intent( in) :: ocean
    real(kind_rb), dimension(:), intent( in) :: coszrs   ! Cosine of solar zenith angle
    real(kind_rb), dimension(:), intent( in) :: ts   ! surface temperature
    real(kind_rb), dimension(:), intent(out) :: asdir, & ! Srf alb for direct rad   0.2-0.7 micro-ms
                                       aldir, & ! Srf alb for direct rad   0.7-5.0 micro-ms
                                       asdif, & ! Srf alb for diffuse rad  0.2-0.7 micro-ms
                                       aldif    ! Srf alb for diffuse rad  0.7-5.0 micro-ms
    !real(kind_rb), parameter :: adif = 0.06
    real(kind_rb), parameter :: adif = 0.07 !AAW 12/5/17 RCEMIP albedo value
    !-----------------------------------------------------------------------
    if (ocean) then
      !
      !
      ! Ice-free ocean albedos function of solar zenith angle only, and
      ! independent of spectral interval:
      !
      where(coszrs <= 0)
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        where(ts > 271.)
          aldir(:) = ( .026 / (coszrs(:)**1.7 + .065)) + &
             (.15*(coszrs(:) - 0.10) * (coszrs(:) - 0.50) * (coszrs(:) - 1.00) )
          asdir(:) = aldir(:)
          aldif(:) = adif
          asdif(:) = adif
        elsewhere
         ! albedo of sea-ice/snow surface
           aldir(:) = 0.45
           asdir(:) = 0.75
           aldif(:) = 0.45
           asdif(:) = 0.75
        end where
      end where
    else ! land
      where(coszrs <= 0)
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        ! Albedos for land type I (Briegleb)
        asdir(:) = 1.4 * 0.06 / ( 1. + 0.8 * coszrs(:))
        asdif(:) = 1.2 * 0.06
        aldir(:) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(:))
        aldif(:) = 1.2 * 0.24
      end where
    endif
  end subroutine albedo


end module cam_rad_parameterizations
