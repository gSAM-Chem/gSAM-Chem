!---------------------------------------------------------------------------------------
! vars.f90
! Coded by Marat Khairoutdinov
! Contains prognostic and diagnostic arrays nessesary to run the dynamical core.
! Also defines all the 2D-output arrays, statistical profiles, nudging arrays, and forcing arrays.
!---------------------------------------------------------------------------------------

module vars

use grid

implicit none
!---------------------------------------------------------------------------------------
! prognostic variables:

real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! liquid/ice water static energy 

!--------------------------------------------------------------------
! diagnostic variables:

real p      (0:nx, (1-YES3D):ny, nzm, 3)    ! perturbation pi-pressure (pressure/rho)
real pp     (0:nx, (1-YES3D):ny, nzm)       ! pressure (hPa) ! actual 3D pressure
real tabs   (nx, ny, nzm)                ! temperature
real qv     (nx, ny, nzm)                ! water vapor
real qcl    (nx, ny, nzm)                ! liquid water  (condensate)
real qpl    (nx, ny, nzm)                ! liquid water  (precipitation)
real qci    (nx, ny, nzm)                ! ice water  (condensate)
real qpi    (nx, ny, nzm)                ! ice water  (precipitation)
real qpg    (nx, ny, nzm)                ! graupel water  (part of qpi)
real cld    (nx, ny, nzm)                ! cloud fraction (from microphysics)
        
!--------------------------------------------------------------------
! time-tendencies for prognostic variables

real dudt   (nx, ny, nzm, 3)
real dvdt   (nx, ny, nzm, 3)
real dwdt   (nx, ny, nz,  3)

! for forward-in-time diffusion:
real dudtd  (nx, ny, nzm)
real dvdtd  (nx, ny, nzm)
real dwdtd  (nx, ny, nz)

!----------------------------------------------------------------
! Temporary storage array:

real misc(nx, ny, nz), misc1(nx, ny, nz)
! storage of velocity on previous time step to avoid extra boundary exchenge:
real u1   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v1   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w1   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real(8), allocatable ::  pfy(:,:,:)  ! the previous-step solution of GMG pressure solver
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny), raf (nx,ny), rafq(nx,ny)
real eis(nx,ny) ! estimated inversion strength
real precsfc(nx,ny) ! surface precip. rate averaged over statistics sampling period
real precinst(nx,ny) ! surface precip. rate (instantaneous)
real precinstsoil(nx,ny) ! surface precip. rate on soil (instantaneous)
real precsfcsnow(nx,ny) ! surface snow rate averaged over statistics sampling period

real shf_ocean(nx,ny), lhf_ocean(nx,ny), shf_all(nx,ny)  ! diagnostics 
real shf_land(nx,ny), lhf_land(nx,ny), lhf_all(nx,ny), evp_all(nx,ny) ! diagnostics 
real shf_top(nx,ny), lhf_top(nx,ny)
                
integer k200,k500,k850,k700,k950  ! vertical indexes of arrays that correspond to given pressure.
                             ! For example, k200 is an index closest to 200 mb pressure level.
!-----------------------------------------------------------------
! Some usel diagnostics (for stat file):

real:: t2m(nx,ny)=0., q2m(nx,ny)=0., u10m(nx,ny)=0., v10m(nx,ny)=0., soilflux(nx,ny)=0.
!-----------------------------------------------------------------
! horizotally averaged profiles to compute local perturbations:

real   t0(nzm), q0(nzm), qv0(nzm), tabs0(nzm), tl0(nzm), &
       tv0(nzm), u0(nzm), v0(nzm), &
       tg0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), &
       t01(nzm), q01(nzm), qp0(nzm), qn0(nzm)
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

real  sst1obs, lhobs, shobs, tauobs, t2mobs, q2mobs, wind10mobs
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:
 
real   prespot(nzm)  ! (1000./pres)**R/cp
real   prespoti(nz)  ! (1000./presi)**R/cp
real   rho(nzm)	  ! air density at pressure levels,kg/m3 
real   rhow(nz)   ! air density at vertical velocity levels,kg/m3
real   bet(nzm)	  ! = ggr/tv0
real   gamaz(nzm) ! ggr/cp*z
real   wsub(nz)   ! Large-scale subsidence velocity,m/s
real   qtend(nzm) ! Large-scale tendency for total water
real   ttend(nzm) ! Large-scale tendency for temp.
real:: taudamp(nzm)=0. ! spange-layer damping coefficient (see damping.f90)

!----------------------------------------------------------------
! ocean model interface arrays:

real(8) shf_ocn(nx,ny)   ! sensible heat flux (W/m2)
real(8) lhf_ocn(nx,ny)   ! latent heat flux (W/m2)
real(8) taux_ocn(nx,ny)    ! stress in x (N/m2)
real(8) tauy_ocn(nx,ny)    ! stress in y (N/m2)
real(8) prec_ocn(nx,ny)    ! precipitation (kg/s/m2)
real(8) roff_ocn(nx,ny)   ! river runoff (kg/s/m2)
real(8) lw_ocn(nx,ny)      ! net LW flux (W/m2)
real(8) sw_ocn(nx,ny)      ! net SW flux (W/m2)
real(8) u_ocn(nx,ny)      ! ocean current in x (m/s)
real(8) v_ocn(nx,ny)      ! ocean current in y (m/s)
!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing units
integer nsnd	! number of observed soundings
integer nzlsf	! number of large-scale forcing profiles
integer nzrfc	! number of radiative forcing profiles
integer nzsnd	! number of observed soundings

real, allocatable :: dqls(:,:) ! Large-scale tendency for total water
real, allocatable :: dtls(:,:) ! Large-scale tendency for temp.
real, allocatable :: ugls(:,:) ! Large-scale wind in X-direction
real, allocatable :: vgls(:,:) ! Large-scale wind in Y-direction
real, allocatable :: wgls(:,:) ! Large-scale subsidence velocity,m/s
real, allocatable :: pres0ls(:)! Surface pressure, mb
real, allocatable :: zls(:,:)  ! Height
real, allocatable :: pls(:,:)  ! Pressure
real, allocatable :: dayls(:)  ! Large-scale forcing arrays time (days) 
real, allocatable :: dtrfc(:,:)! Radiative tendency for pot. temp.
real, allocatable :: dayrfc(:) ! Radiative forcing arrays time (days) 
real, allocatable :: prfc(:,:) ! Pressure/Height
real, allocatable :: sstsfc(:) ! SSTs
real, allocatable :: shsfc(:)   ! Sensible heat flux,W/m2
real, allocatable :: lhsfc(:)  ! Latent heat flux,W/m2
real, allocatable :: tausfc(:) ! Surface drag,m2/s2
real, allocatable :: t2msfc(:) ! Temperature at 2m, K
real, allocatable :: td2msfc(:) ! Dewpoint temperature at 2m, K
real, allocatable :: wind10msfc(:) ! Wind speed at 10m, m/s
real, allocatable :: daysfc(:) ! Surface forcing arrays time (days) 
real, allocatable :: usnd(:,:) ! Observed zonal wind
real, allocatable :: vsnd(:,:) ! Observed meriod wind
real, allocatable :: tsnd(:,:) ! Observed Abs. temperature
real, allocatable :: qsnd(:,:) ! Observed Moisture
real, allocatable :: zsnd(:,:) ! Height
real, allocatable :: psnd(:,:) ! Pressure
real, allocatable :: daysnd(:) ! number of sounding samples

real sstobs(nx,ny,2) ! prescribed SSTs
real:: daysstobs(2) = 0.  ! time arrays (days)

! 3D fields for nudging:
real, allocatable :: uobs(:,:,:,:) 
real, allocatable :: vobs(:,:,:,:)
real, allocatable :: wobs(:,:,:,:)
real, allocatable :: tobs(:,:,:,:)
real, allocatable :: qobs(:,:,:,:)
real:: dayfld3Dobs(2) = 0.  ! time arrays (days)
real, allocatable :: tau_nudge3D(:,:) ! inverse nudging time for scalars and w
real, allocatable :: tau_nudge3Du(:,:) ! inverse nudging time for u
real, allocatable :: tau_nudge3Dv(:,:) ! inverse nudging time for v
real, allocatable :: tend_in_u(:,:,:)  ! spectral nudging tendency for u
real, allocatable :: tend_in_v(:,:,:)  ! spectral nudging tendency for v

 
real:: u_w(dimy1_u:dimy2_u,nzm) = 0. ! west x-boundanry velocities for region
real:: u_e(dimy1_u:dimy2_u,nzm) = 0. ! east x-boundanry velocities for region
real:: v_s(dimx1_v:dimx2_v,nzm) = 0. ! south y-boundanry velocities for region
real:: v_n(dimx1_v:dimx2_v,nzm) = 0. ! north y-boundanry velocities for region

!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real(8)  fcory(0:ny)          !  Coriolis parameter xy-distribution
real(8)  fcorzy(1:ny)         !  z-Coriolis parameter xy-distribution
real(8)  latitude(nx,0:ny+1)  ! latitude (degrees)
real(8)  longitude(nx,ny)     ! longitude(degrees)
real(8)  sstxy(nx,ny)	  !  skin surface temperature (over land) and SST  (perturbation from t00)
real sstxy0(nx,ny)	  !  prescribed sst
real:: elevation(nx,ny) = 0.     ! actual (not grid) surface elevation, m
real:: elevationg(nx,ny) = 0.    ! surface elevation truncated to grid levelas, m
integer landmask(nx,ny)  ! land mask, 1 for land, 0 for ocean
real seaice_h(nx,ny)     ! sea-ice thickness. Activated for skin-ocean model and doseaiceevol=T) 

real, parameter :: t00 = 273.15   ! reference temperature for sstxy()
!---------------------------------------------------------------------
! 2D output arrays:
real prec_xy(nx,ny)      ! precip. rate 
real preca_xy(nx,ny)     ! mean precip. rate 
real precac_xy(nx,ny)    ! accumulated precip. rate 
real precs_xy(nx,ny)     ! frozen precip. rate 
real precsa_xy(nx,ny)    ! frozen mean precip. rate 
real precsac_xy(nx,ny)   ! frozen accumulated precip. rate 
real taux_xy(nx,ny)      ! stress in x
real tauy_xy(nx,ny)      ! stress in y
real taux_top_xy(nx,ny)  ! stress in x at top
real tauy_top_xy(nx,ny)  ! stress in y at top
real shf_xy(nx,ny)       ! sensible heat flux 
real lhf_xy(nx,ny)       ! latent heat flux
real shf_top_xy(nx,ny)   ! sensible heat flux at top
real lhf_top_xy(nx,ny)   ! latent heat flux at top
real lwnt_xy(nx,ny)      ! net lw at TOA
real swnt_xy(nx,ny)      ! net sw at TOA
real lwntc_xy(nx,ny)     ! clear-sky mean net lw at TOA
real swntc_xy(nx,ny)     ! clear-sky mean net sw at TOA
real lwns_xy(nx,ny)      ! net lw at SFC
real swns_xy(nx,ny)      ! net sw at SFC
real lwnsc_xy(nx,ny)     ! clear-sky mean net lw at SFC
real swnsc_xy(nx,ny)     ! clear-sky mean net sw at SFC
real lwds_xy(nx,ny)      ! downward lw at SFC
real swds_xy(nx,ny)      ! downward sw at SFC
real lwdsc_xy(nx,ny)     ! clear sky downward lw at SFC
real swdsc_xy(nx,ny)     ! clear sky downward sw at SFC
real lwnta_xy(nx,ny)     ! mean net lw at TOA
real lwntac_xy(nx,ny)    ! accumulated net lw at TOA 
real swnta_xy(nx,ny)     ! mean net sw at TOA
real swntac_xy(nx,ny)    ! accumulated net sw at TOA
real lwntca_xy(nx,ny)    ! clear-sky mean net lw at TOA
real swntca_xy(nx,ny)    ! clear-sky mean net sw at TOA
real lwnsa_xy(nx,ny)     ! mean net lw at SFC
real lwnsac_xy(nx,ny)    ! accumulated net lw at SFC
real swnsa_xy(nx,ny)     ! mean net sw at SFC
real swnsac_xy(nx,ny)    ! accumulated net sw at SFC
real lwnsca_xy(nx,ny)    ! clear-sky mean net lw at SFC
real swnsca_xy(nx,ny)    ! clear-sky mean net sw at SFC
real lwdsa_xy(nx,ny)     ! mean downward lw at SFC
real lwdsac_xy(nx,ny)    ! accumulated downward lw at SFC
real swdsa_xy(nx,ny)     ! mean downward sw at SFC
real swdsac_xy(nx,ny)    ! accumulated downward sw at SFC
real lwdsca_xy(nx,ny)    ! clear sky downward lw at SFC
real swdsca_xy(nx,ny)    ! clear sky downward sw at SFC
real solin_xy(nx,ny)     ! solar TOA insolation
real solina_xy(nx,ny)    ! solar TOA insolation
real pw_xy(nx,ny)        ! precipitable water
real pwobs_xy(nx,ny)     ! observed precipitable water (when nudged to obs)
real pws_xy(nx,ny)       ! saturation precipitable water
real fse_xy(nx,ny)       ! frozen moist static energy
real ta_xy(nx,ny)        ! mass-weighted column temperature
real cw_xy(nx,ny)        ! cloud water path
real iw_xy(nx,ny)        ! cloud ice path
real rw_xy(nx,ny)        ! rain path
real sw_xy(nx,ny)        ! snow path
real gw_xy(nx,ny)        ! grouple path
real cld_xy(nx,ny)       ! cloud frequency
real eis_xy(nx,ny)       ! Estimated Inversion Strength (EIS)
real u200_xy(nx,ny)      ! u-wind at 200 mb
real uobs200_xy(nx,ny)   ! observed u-wind at 200 mb
real usfc_xy(nx,ny)      ! u-wind at at first level above surface
real v200_xy(nx,ny)      ! v-wind at 200 mb
real vobs200_xy(nx,ny)   ! observed v-wind at 200 mb
real vsfc_xy(nx,ny)      ! v-wind at first level above surface
real gustsfc_xy(nx,ny)   ! maximum wind at at first level above surface
real u10m_xy(nx,ny)      ! 10-m u-wind 
real v10m_xy(nx,ny)      ! 10-m u-wind 
real u10ma_xy(nx,ny)     ! average 10-m u-wind 
real v10ma_xy(nx,ny)     ! average 10-m v-wind 
real gust10m_xy(nx,ny)   ! maximum 10-m wind 
real t2m_xy(nx,ny)       ! 2-m temperature 
real q2m_xy(nx,ny)       ! 2-m humidity 
real w500_xy(nx,ny)      ! w at 500 mb
real phi500_xy(nx,ny)    ! Perturbation of 500 mb geopotential, m
real ZdBZ_xy(nx,ny)      ! Composite radar reflectivity, dBZ
real omega200_xy(nx,ny)  ! omega at 200 mb
real omega500_xy(nx,ny)  ! omega at 500 mb
real omega700_xy(nx,ny)  ! omega at 700 mb
real omega850_xy(nx,ny)  ! omega at 850 mb
real rh200_xy(nx,ny)     ! relative humidity at 200 mb
real rh500_xy(nx,ny)     ! relative humidity at 500 mb
real rh700_xy(nx,ny)     ! relative humidity at 700 mb
real rh850_xy(nx,ny)     ! relative humidity at 850 mb
real qocean_xy(nx,ny)    ! ocean cooling in W/m2
real sst_xy(nx,ny)       ! surface temperature
real tsfc_xy(nx,ny)      ! temperature  near the surface
real qsfc_xy(nx,ny)      ! vapor  near the surface
real sst_min_xy(nx,ny)   ! minimum surface temperature
real sst_max_xy(nx,ny)   ! maximum surface temperature
real cflz_max_xy(nx,ny)  ! maximum cfl in vertical
real cflh_max_xy(nx,ny)  ! maximum cfl in horizonal
real cfl_max_xy(nx,ny)   ! maximum cfl in horizonal
real zcflz_max_xy(nx,ny) ! height of maximum cfl in vertical
real zcflh_max_xy(nx,ny) ! height of maximum cfl in horizonal
real zcfl_max_xy(nx,ny)  ! height of maximum cfl in horizonal
real soil_wet_xy(nx,ny)  ! soil wetness
real snow_depth_xy(nx,ny)! snow depth
real snow_melt_xy(nx,ny) ! snow melt rate
real grnd_xy(nx,ny)      ! in-ground flux (W/m2)
real tsoil_xy(nx,ny)     ! soil temperature at 30-cm depth
real ra_xy(nx,ny)        ! aerodynamic resistance
real rc_xy(nx,ny)        ! undercanopy aerodynamic resistance
real albvis_xy(nx,ny)    ! diffuse visible albedo
real albnir_xy(nx,ny)    ! diffuse near-infrared albedo
real alb_xy(nx,ny)       ! diagnostic field of SWNT with cos(zenith angle)=1 everywhere
real albc_xy(nx,ny)      ! diagnostic field of SWNT with cos(zenith angle)=1 everywhere (clear-sky)
real vor200_xy(nx,ny)    ! vorticity 200 mb
real vor850_xy(nx,ny)    ! vorticity 850 mb

!----------------------------------------------------------------------
!	quantities sampled for statitistics purposesi (*.stat file):

real &
    twle(nz), twsb(nz), precflux(nz), prectot, wvp, cwp, iwp, rwp, gwp, swp, wvpobs, &
    uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
    radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
    radqrlw(nz), radqrsw(nz), radqrclw(nz), radqrcsw(nz), &
    w_max, u_max, s_acld, s_acldcold, s_ar, s_arthr, s_sst, s_pw, s_pwobs, &
    s_acldl, s_acldm, s_acldh,  ncmn, nrmn, z_inv, z_cb, z_ct, z_cbmn, z_ctmn, &
    z2_inv, z2_cb, z2_ct, cwpmean, cwp2, precmean, snowmean, prec2, precmax, nrainy, ncloudy, &
    s_acldisccp, s_acldlisccp, s_acldmisccp, s_acldhisccp, s_ptopisccp, &
    s_acldmodis, s_acldlmodis, s_acldmmodis, s_acldhmodis, s_ptopmodis, &
    s_acldmisr, s_ztopmisr, s_relmodis, s_reimodis, s_lwpmodis, s_iwpmodis, &
    s_tbisccp, s_tbclrisccp, s_acldliqmodis, s_acldicemodis, &
    s_cldtauisccp,s_cldtaumodis,s_cldtaulmodis,s_cldtauimodis,s_cldalbisccp, &
    s_flns,s_flnt,s_flntoa,s_flnsc,s_flntoac,s_flds,s_fsns, s_fldsc, s_fsdsc, &
    s_fsnt,s_fsntoa,s_fsnsc,s_fsntoac,s_fsds,s_solin, & 
    s_lhf, s_shf, s_lhfc, s_shfc, s_lhfs, s_shfs, s_soilflux, s_tsoil, s_tcanopy, &
    s_t2m, s_q2m, s_u10m, s_v10m, s_snowd, s_snowt, &
    tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz),tkelebuoy(nz), &
    t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
    q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
    twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
    qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
    momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
    momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tdiss(nz), tbuoy(nz), tpdiss(nz), tlatqi(nz),qifall(nz),qpfall(nz), &
    tlatqc(nz),qcfall(nz), fluxwallt(nz,4), fluxwallq(nz,4)


! energy conservation diagnostics:
 
  real(8) total_water_before, total_water_after
  real(8) total_water_evap, total_water_prec, total_water_ls, total_water_adv

!===========================================================================
! UW ADDITIONS (By Peter Blossey)

! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
     icondavg_satdn, icondavg_satup, icondavg_env
real, allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
real, allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
character(LEN=8), allocatable :: condavgname(:) ! array of short names
character(LEN=25), allocatable :: condavglongname(:) ! array of long names

real   qlsvadv(nzm) ! Large-scale vertical advection tendency for total water
real   tlsvadv(nzm) ! Large-scale vertical advection tendency for temperature
real   ulsvadv(nzm) ! Large-scale vertical advection tendency for zonal velocity
real   vlsvadv(nzm) ! Large-scale vertical advection tendency for meridional velocity

real   qnudge(nzm) ! Nudging of horiz.-averaged total water profile
real   tnudge(nzm) ! Nudging of horiz.-averaged temperature profile
real   unudge(nzm) ! Nudging of horiz.-averaged zonal velocity
real   vnudge(nzm) ! Nudging of horiz.-averaged meridional velocity

real   qstor(nzm) ! Storage of horiz.-averaged total water profile
real   tstor(nzm) ! Storage of horiz.-averaged temperature profile
real   ustor(nzm) ! Storage of horiz.-averaged zonal velocity
real   vstor(nzm) ! Storage of horiz.-averaged meridional velocity

real   utendcor(nzm) ! coriolis acceleration of zonal velocity
real   vtendcor(nzm) ! coriolis acceleration of meridional velocity

real   u850_xy(nx,ny) ! zonal velocity at 850 mb
real   v850_xy(nx,ny) ! meridional velocity at 850 mb
real   uobs850_xy(nx,ny) ! observed zonal velocity at 850 mb (when nudged to obs)
real   vobs850_xy(nx,ny) ! observed meridional velocity at 850 mb

! Surface pressure
real psfc_xy(nx,ny) ! pressure (in millibar) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
real swvp_xy(nx,ny)  ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
real cloudtopheight(nx,ny), echotopheight(nx,ny), cloudtoptemp(nx,ny), &
     cloudcover(nx,ny)

!  Use these variables to make sure the run does not exceed nelapse_minutes
!  The time is checked during stepout each time a full set of restart files
!  are written
integer :: nelapse_minutes = -99999, oldstep
real(8) :: time_init_seconds, time_old_seconds


! END UW ADDITIONS
!===========================================================================


!===========================================================================
! UW ADDITIONS (by James Yoon)

real interactive_soil_wetness(nx, ny) ! soil wetness as an interactive array

! END JY ADDITIONS
!===========================================================================



end module vars
