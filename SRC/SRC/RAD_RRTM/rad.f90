module rad
    ! -------------------------------------------------------------------------- 
    !
    ! Interface to RRTM longwave and shortwave radiation code.
    !   Robert Pincus, November 2007
    !  
    ! Modified by Peter Blossey, July 2009.
    !   - interfaced to RRTMG LW v4.8 and SW v3.8.
    !   - reversed indices in tracesini to match bottom-up RRTM indexing.
    !   - fixed issue w/doperpetual caused by eccf=0.
    !   - added extra layer in calls to RRTM to improve heating rates at
    !        model top.  Required changes to inatm in rrtmg_lw_rad.f90
    !        and inatm_sw in rrtmg_sw_rad.f90.
    !   - fixed bug that would zero out cloudFrac if LWP>0 and IWP==0.
    !   - changed definition of interfaceT(:,1) to SST.
    !   - added o3, etc. profiles to restart file.  Only call tracesini if nrestart==0.
    !
    ! Modified by Peter Blossey, July 2009.
    !   - update to RRTMG LW v4.84
    !   - add extra layer to model top within rad_driver, rather than inatm_*w.
    !   - use instantaneous fields for radiation computation, rather than time-averaged.
    !
    ! Modified by Peter Blossey, May 2014.
    !   - update to RRTMG SW v3.9 which include fix in lookup tables for cloud
    !      liquid optical properties.  See AER's description of the update here:
    !        http://rtweb.aer.com/rrtmg_sw_whats_new.html
    !
    ! Modified by Peter Blossey (with help from Robert Pincus), Sept 2015.
    !   - Moved all "use" statements to top of module.  This seems to be
    !      better fortran90 coding practice.
    !   - Coupled radiation more tightly with instrument simulators in 
    !      SRC/SIMULATORS/ by providing tau_067 and emis_105 with values from
    !      closest RRTMG bands.
    !   - Added tighter microphysics-coupling for M2005 and Thompson microphysics.
    !      MICRO_M2005 uses CAM5 cloud optics lookup tables for cloud liquid, cloud ice
    !      and snow.  MICRO_WRF (Thompson microphysics) uses the RRTMG lookup tables
    !      for cloud liquid, cloud ice and snow.  Both routines use mass and other cloud
    !      properties (mostly effective radius or generalized effective diameter for ice)
    !      to determine the optical properties.  See the new modules for computing
    !      cloud properties in m2005_cloud_optics.f90 and thompson_cloud_optics.f90.
    !      These new cloud optics routines are enabled by default.  They can be disabled
    !      in the microphysics namelists with 
    !          dorrtm_cloud_optics_from_effrad_LegacyOption = .true.
    !
    ! Modified by Peter Blossey (with help from Robert Pincus), Feb 2016.
    !   - providing separate optical depths for liquid, cloud ice and
    !      snow to better couple with the MODIS simulator in SRC/SIMULATORS/.
    !
    ! Modified by Marat Khairoutdinov, May 2023
    ! - added call to albedo_slm to compute albedoes based on the input from the 
    !   Simplified Land Model (SLM) as well as included module slm_vars 
    ! - implemented the radiation computation for the case of terrain by shifting the 
    !   radiation input profiles to start above the local topography height and adding ghost levels
    !   above the domain top to preseve the total number of levels for radiation calls
    ! -------------------------------------------------------------------------- 

    ! -------------------------------------------------------------------------- 
  use parkind, only : kind_rb, kind_im 
  use params, only : docheck, ggr, rgas, coszrs, cp, ocean,  &
       doshortwave, dolongwave, doradhomo, doradhomozonal,     &
       doseasons, doperpetual, dosolarconstant, doequinox, &
       solar_constant, zenith_angle, nxco2, notracegases, nrad, SLM, &
       n2ox, ch4x, cfc11x, cfc12x, reado3, o3file, dosubgridcloudfraction, docurrentco2
  use grid, only : nx, ny, nz, nzm, compute_reffc, compute_reffi, &
       icycle, dtn, dt, nstop, nstep, nstat, nrestart, day, day0, year, &
       pres, presi, z, dz, adz, dompi, masterproc, nsubdomains,          &
       dostatis, dostatisrad, nelapse, nrestart_skip, case, rundatadir, &
       doisccp, domodis, domisr, nrestart_steps,&
       restart_sep, caseid, case_restart, caseid_restart, rank, wgtxys, wgtxyt, wgt, wgtw
  use vars, only : landmask, t, tabs, qv, qcl, qci, cld, sstxy, rho, t00, tabs0, & !bloss: cld(:,:,:) for partial cloudiness
       latitude, longitude,                         &
                                ! Domain-average diagnostic fields
       radlwup, radlwdn, radswup, radswdn, radqrlw, radqrsw, radqrclw, radqrcsw, &
                                ! 2D diagnostics
       lwns_xy, lwnt_xy, lwds_xy, swns_xy, swnt_xy, swds_xy,  &
       lwdsc_xy, swdsc_xy, &
       lwnsa_xy, lwnta_xy, lwdsa_xy, swnsa_xy, swnta_xy, swdsa_xy, solin_xy, &      
       lwnsc_xy, lwntc_xy, swnsc_xy, swntc_xy, &
                                ! 1D diagnostics
       s_flns, s_fsns, s_flnt, s_flntoa, s_fsnt, s_fsntoa, &
       s_flnsc, s_fsnsc, s_flntoac, s_fsntoac, s_solin, &
       s_fsds, s_flds, s_fsdsc, s_fldsc
   use terrain, only: terra, terraw, k_terra
   use slm_vars, only: albedovis_v, albedonir_v, albedovis_s, albedonir_s, &
                      phi_1, phi_2, snow_mass, LAI, IMPERV, icemask, soilw, landtype

   use consts, only: emis_water

  !
  ! Radiation solvers
  !
  use rrlw_ncpar, only: cpdair, maxAbsorberNameLength, &
       status, getAbsorberIndex
  use rrtmg_sw_init, only: rrtmg_sw_ini
  use rrtmg_lw_init, only: rrtmg_lw_ini
  use rrtmg_sw_rad, only : rrtmg_sw
  use rrtmg_lw_rad, only : rrtmg_lw
  use rrtmg_sw_rad_nomcica, only : rrtmg_sw_nomcica => rrtmg_sw
  use rrtmg_lw_rad_nomcica, only : rrtmg_lw_nomcica => rrtmg_lw
  use mcica_subcol_gen_lw, only: mcica_subcol_lw
  use mcica_subcol_gen_sw, only: mcica_subcol_sw
  use rrtmg_lw_cldprop, &
                     only: cldprop
  use rrtmg_sw_cldprop, & 
                    only : cldprop_sw
  use parrrtm,      only : nbndlw, ngptlw ! Number of LW bands and g-points
  use parrrsw,      only : nbndsw, ngptsw, naerec, jpband ! Number of SW bands and g-points
  use cam_rad_parameterizations, only : &
       computeRe_Liquid, computeRe_Ice, albedo, albedo_slm
  use microphysics, only : micro_scheme_name, reffc, reffi, &
                           CloudIceMassMixingRatio, &
                           dorrtm_cloud_optics_from_effrad_LegacyOption
  use m2005_cloud_optics, &
                    only : m2005_cloud_optics_init, compute_m2005_cloud_optics
  use thompson_cloud_optics, &
                    only : thompson_cloud_optics_init, compute_thompson_cloud_optics

  use shr_orb_mod, only: shr_orb_params, shr_orb_decl, shr_orb_cosz

  implicit none

  real, dimension(nx, ny):: lwnsxy, swnsxy, & ! Long- and short-wave radiative heating (W/m2)  
      swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, swdsxy, coszrsxy, &
      swusvisdxy, swusnirdxy
  character(5), parameter :: RAD_NAME="RRTMG_MCICA"

  !private

  ! Public procedures
  public :: rad_driver, write_rad
  
  ! Public data
  public :: qrad, lwnsxy, swnsxy, &
       do_output_clearsky_heating_profiles, &
       tau_067, emis_105, rad_reffc, rad_reffi, &  ! For instrument simulators
       tau_067_cldliq, tau_067_cldice, tau_067_snow ! for MODIS simulator: wants individual phases


  real, dimension(nx, ny, nzm) :: qrad = 0. ! Radiative heating rate (K/s) 
  real, dimension(nx, ny, nzm) :: tau_067, emis_105, &  ! Optical thickness at 0.67 microns, emissivity at 10.5 microns, for instrument simulators
                                  tau_067_cldliq, tau_067_cldice, tau_067_snow, & ! for MODIS simulator: wants individual phases
                                  rad_reffc, rad_reffi  ! Particle sizes assumed for radiation calculation when microphysics doesn't provide them

  logical, parameter :: do_output_clearsky_heating_profiles = .true.

  logical :: do_partial_cloudiness_in_radiation = .false.

  real ozone  (nx, ny, nzm+1) ! ozone 3D field (g/g) if read from 3D initfile

  !
  ! Constants
  !
  real, parameter :: Pi = 3.14159265358979312
  real(kind = kind_rb), parameter :: scon = 1367. ! solar constant 
  integer, parameter :: iyear = 1999
  !
  ! Molecular weights (taken from CAM shrc_const_mod.F90 and physconst.F90)
  !
  real, parameter :: mwdry =  28.966, & ! molecular weight dry air
                     mwco2 =  44.,    & ! molecular weight co2
                     mwh2o =  18.016, & ! molecular weight h2o
                     mwn2o =  44.,    & ! molecular weight n2o
                     mwch4 =  16.,    & ! molecular weight ch4
                     mwf11 = 136.,    & ! molecular weight cfc11
                     mwf12 = 120.,    & ! molecular weight cfc12
                     mwo3  =  48.       ! ozone, strangely missing
  ! mixingRatioMass = mol_weight/mol_weight_air * mixingRatioVolume
  
  !
  ! Global storage
  !
  logical :: initialized = .false., use_m2005_cloud_optics = .false., &
       use_thompson_cloud_optics = .false., have_cloud_optics = .false.
  real(KIND=kind_rb) :: land_frac = 1.
  
  real, dimension(nx, ny) :: &
     lwDownSurface, lwDownSurfaceClearSky, lwUpSurface, lwUpSurfaceClearSky, & 
                                           lwUpToa,     lwUpToaClearSky,     &
     lwDownTom,                            lwUpTom,  &
     swDownSurface, swDownSurfaceClearSky, swUpSurface, swUpSurfaceClearSky, & 
     swDownToa,                            swUpToa,     swUpToaClearSky,     &
     swDownTom,                            swUpTom,  &
     insolation_TOA, &
     visDownSurface, visDownSurfaceDiffuse, nirDownSurface, nirDownSurfaceDiffuse, &
     CosZenithAngle

  real(kind = kind_rb), dimension(nx, nz+1)  :: swUp, swDown

  !bloss(072009): changed from mass mixing ratios to volume mixing ratios
  !                 because we're now using rrtmg_lw.nc sounding for trace gases.
  ! Profiles of gas volume mixing ratios 
  !bloss(120409): add level to account for trace gases above model top.
  real(kind_rb), dimension(nzm+1) :: o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4
  
  integer :: nradsteps ! current number of steps done before
				       !   calling radiation routines

  real, dimension(nx, ny) ::  p_factor ! perpetual-sun factor
   
  !
  ! Earth's orbital characteristics
  !   Calculated in shr_orb_mod which is called by rad_driver
  !
  real(kind_rb), save ::  eccf,  & ! eccentricity factor (1./earth-sun dist^2)
                 eccen, & ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
                 obliq, & ! Earth's obliquity angle (deg) (-90 to +90) (typically 22-26)
                 mvelp, & ! Earth's moving vernal equinox at perhelion (deg)(0 to 360.0)
                 !
                 ! Orbital information after processed by orbit_params
                 !
                 obliqr, &  ! Earth's obliquity in radians
                 lambm0, &  ! Mean longitude of perihelion at the
                            ! vernal equinox (radians)
                 mvelpp     ! Earth's moving vernal equinox longitude
                            ! of perihelion plus pi (radians)
!
   logical:: initrad = .true. ! dummy variable for compatibility with BUILDINGS (later set to .false.)

   ! seed for random number routine used in MCICA subcolumn cloud generator
   integer, save :: permuteseed = 199
   integer :: irng = 1 ! choice of Mersenne Twister for random number generator 

contains 
  ! ----------------------------------------------------------------------------
  subroutine rad_driver 

    implicit none
    
    ! -------------------------------------------------------------------------- 

    ! Local variables 
    !
    ! Input and output variables for RRTMG SW and LW
    !   RRTM specifies the kind of real variables in 
    !   Only one column dimension is allowed parkind
    !   RRTM is indexed from bottom to top (like the CRM) 
    !
    !bloss: add layer to top to improve top-of-model heating rates.
    real(kind = kind_rb), dimension(nx, nzm+1) ::     &
        layerP,     layerT, layerMass,         & ! layer mass is for convenience
        h2ovmr,   o3vmr,    co2vmr,   ch4vmr, n2ovmr,  & ! Volume mixing ratios for H2O, O3, CH4, N20, CFCs
        o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
        swHeatingRate, swHeatingRateClearSky,  &
        lwHeatingRate, lwHeatingRateClearSky, &
        duflx_dt, duflxc_dt, &
        layerM   ! introduced instead of rho*dz*adz for topography case implementation - MK
        
    ! Arrays for cloud optical properties or the physical properties needed by the RRTM internal parmaertizations
    real(kind = kind_rb), dimension(nx, nzm+1) ::     &
         LWP, IWP, liqRe, iceRe, cloudFrac             ! liquid/ice water path (g/m2) and size (microns)
    real(kind = kind_rb), dimension(nbndlw, nx, nz+1) :: cloudTauLW = 0. 
    real(kind = kind_rb), dimension(nbndsw, nx, nz+1) :: cloudTauSW = 0., cloudSsaSW = 0., &
                                                         cloudAsmSW = 0., cloudForSW = 0., &
                                                         cloudTauSW_cldliq = 0., &
                                                         cloudTauSW_cldice = 0., &
                                                         cloudTauSW_snow = 0.

    ! arrays for simulations with partial cloudiness that use MCICA subcolumn generators
    real(kind = kind_rb), dimension(:,:,:), allocatable :: LWPmcl, IWPmcl, cloudFracmcl ! spread across g-points
    real(kind = kind_rb), dimension(:,:), allocatable :: liqRemcl, iceRemcl ! Eff Radius not spread across g-points
    real(kind = kind_rb), dimension(:,:,:), allocatable :: cloudTauLWmcl 
    real(kind = kind_rb), dimension(:,:,:), allocatable :: cloudTauSWmcl, &
         cloudSsaSWmcl, cloudAsmSWmcl, cloudForSWmcl
!!$    real(kind = kind_rb), dimension(ngptlw, nx, nzm+1) :: LWPmcl, IWPmcl, cloudFracmcl ! spread across g-points
!!$    real(kind = kind_rb), dimension(nx, nzm+1) :: liqRemcl, iceRemcl ! Eff Radius not spread across g-points
!!$    real(kind = kind_rb), dimension(ngptlw, nx, nz+1) :: cloudTauLWmcl 
!!$    real(kind = kind_rb), dimension(ngptsw, nx, nz+1) :: cloudTauSWmcl, cloudSsaSWmcl, &
!!$                                                         cloudAsmSWmcl, cloudForSWmcl

    ! arrays for aerosol radiative properties
    real(kind = kind_rb), dimension(nx, nzm+1, nbndlw) :: dummyTauAerosolLW = 0. 
    real(kind = kind_rb), dimension(nx, nzm+1, nbndsw) :: dummyAerosolProps = 0. 
    real(kind = kind_rb), dimension(nx, nzm+1, naerec) :: dummyAerosolProps2 = 0. 

    ! Arguments to RRTMG cloud optical depth routines
    real(kind = kind_rb), dimension(nbndlw, nzm+1) :: prpLWIn
    real(kind = kind_rb), dimension(nzm+1, nbndlw) :: tauLWOut
    real(kind = kind_rb), dimension(nbndsw, nzm+1) :: prpSWIn
    real(kind = kind_rb), dimension(nzm+1, jpband) :: tauSWOut, scaled1, scaled2, scaled3 
    integer                                        :: ncbands  
                                
    !bloss: add layer to top to improve top-of-model heating rates.
    real(kind = kind_rb), dimension(nx, nz+1)  :: interfaceP, interfaceT,        &
                                swUpClearSky, swDownClearSky, &
                                lwUp, lwDown, lwUpClearSky, lwDownClearSky, &
                                visDown, visdDown, nirDown, nirdDown
                                
    real(kind = kind_rb), dimension(nx) :: surfaceT, solarZenithAngleCos
                                      ! Surface direct/diffuse albedos for 0.2-0.7 (s) 
                                      !   and 0.7-5.0 (l) microns
    real(kind = kind_rb), dimension(nx) ::  asdir, asdif, aldir, aldif 

                                
    real(kind = kind_rb), dimension(nx, nbndlw) :: surfaceEmissivity = emis_water ! redefined below
    integer :: overlap_partial_cloudiness = 2 ! maximum/random overlap 
    integer :: overlap_no_partial_cloudiness = 1 ! maximum/random overlap 
    integer :: idrv ! option to have rrtm compute d(OLR)/dTABS for both full and clear sky
    
    integer :: lat, lon, i, j, k, m, ierr, tmp_overlap
    integer n_ozone_call
    real(kind = kind_rb) :: dayForSW, delta
    !
    ! 8 byte reals, I guess, used by MPI
    !
    !bloss: add layer to top to improve top-of-model heating rates.
    real(kind = 8),    dimension(nzm+1) :: radHeatingProfile, tempProfile
    
    real, external :: qsatw, qsati

    !bloss: extra arrays for handling liquid-only and ice-only cloud optical depth
    !   computation for MODIS simulator
    real(kind = kind_rb), dimension(nx, nzm+1) ::     &
        dummyWP, dummyRe, cloudFrac_liq, cloudFrac_ice ! liquid/ice water path (g/m2) and size (microns)
    real coef, tmpz(ny,nzm)
    ! ----------------------------------------------------------------------------

    if(docheck.or.icycle == 1) then   ! Skip subcycles (i.e. when icycle /= 1) 
    
      !------------------------------------------------------
      ! Initialize if necessary 
      !
      if(.not. initialized)  then
          call initialize_radiation
          do k=1,nzm+1
            ozone(:,:,k) = o3(k)
          end do
      end if

      nradsteps = nradsteps+1
    
      ! MK: ovewrite the trace gas data and ozone data is for the case if 3D ozone data are read from a file

      n_ozone_call = nint(86400./dt)
      if(reado3.and.(nstep.eq.1.or.mod(nstep,n_ozone_call).eq.0)) then

        call readlevels(trim(o3file),ozone)
        call task_barrier()
        ozone(:,:,nzm+1) = ozone(:,:,nzm)
        ozone(:,:,:) = 0.6034 * ozone(:,:,:) ! convert concentration from mixing ratio to vmr
        if(n2ox.ne.0.) n2o(:) = n2ox
        if(ch4x.ne.0.) ch4(:) = ch4x
        if(cfc11x.ne.0.) cfc11(:) = cfc11x
        if(cfc12x.ne.0.) cfc12(:) = cfc12x
        if(masterproc) then
          print*,'  z    o3   n2o   ch4   cfc11   cfc12  cfc22'
          i=nx/2
          j=ny/2
          do k=1,nzm
               write(6,'(6g12.4)') z(k),ozone(i,j,k),n2o(k),ch4(k),cfc11(k),cfc12(k)
          end do
        end if
        call fminmax_print('o3:',ozone,1,nx,1,ny,nzm)
        if(minval(ozone).le.0.) then
          print*,'rank=',rank,'ozone concentration should be positive:',minval(ozone),minloc(ozone)
          call task_abort()
        end if

        if(docheck) return
      end if
      !----------------------------------------------------
      ! Update radiation variables if the time is due
      !
    
      !kzm Oct.14, 03 changed == nrad to >= nrad to handle the
      ! case when a smaller nrad is used in restart  
      if(nstep == 1 .or. nradsteps >= nrad) then 
        !
        ! Initialize 1D diagnostics
        !
        radlwup(:) = 0.; radlwdn(:) = 0.
        radswup(:) = 0.; radswdn(:) = 0.
        radqrlw(:) = 0.; radqrsw(:) = 0.
        radqrclw(:) = 0.; radqrcsw(:) = 0.
        qrad(:, :, :) = 0. 

        ! Compute assumed cloud particle sizes if microphysics doesn't provide them 
        if(compute_reffc) then
          rad_reffc(1:nx,1:ny,1:nzm) = reffc(1:nx,1:ny,1:nzm)
        else
          rad_reffc(1:nx,1:ny,1:nzm) = &
                        MERGE(computeRe_Liquid(REAL(tabs(1:nx,1:ny,1:nzm),KIND=kind_rb), land_frac), &
                              0._kind_rb,                                         &
                              qcl(1:nx,1:ny,1:nzm) > 0._kind_rb )
        end if 
        if(compute_reffi) then
          rad_reffi(1:nx,1:ny,1:nzm) = reffi(1:nx,1:ny,1:nzm)
        else
          rad_reffi(1:nx,1:ny,1:nzm) = &
                        MERGE(computeRe_Ice(REAL(tabs(1:nx,1:ny,1:nzm),KIND=kind_rb)), &
                              0._kind_rb,                       &
                              qci(1:nx,1:ny,1:nzm) > 0._kind_rb )
        end if 

        if(do_partial_cloudiness_in_radiation) then
           ALLOCATE(LWPmcl(ngptlw, nx, nzm+1), IWPmcl(ngptlw, nx, nzm+1), &
                CloudFracmcl(ngptlw, nx, nzm+1), &
                liqRemcl(nx, nzm+1), iceRemcl(nx, nzm+1), & ! Note: Eff Radius not spread across g-points
                cloudTauLWmcl(ngptlw, nx, nz+1), cloudTauSWmcl(ngptsw, nx, nz+1), &
                cloudSsaSWmcl(ngptsw, nx, nz+1), cloudAsmSWmcl(ngptsw, nx, nz+1), &
                cloudForSWmcl(ngptsw, nx, nz+1), STAT=ierr)
           if(ierr.ne.0) then
              write(*,*) 'Error: Cannot allocate McICA arrays in RAD_RRTM/rad.f90'
              call task_abort()
           end if
        end if

        ! The radiation code takes a 1D vector of columns, so we loop 
        !   over the y direction
        !

        do lat = 1, ny 

          do lon = 1, nx
            m = k_terra(lon,lat)-1
	  ! set trace gas concentrations.  Assumed to be uniform in the horizontal.
            o3vmr   (lon, 1:nzm+1-m) =  ozone(lon, lat, m+1:nzm+1)
            co2vmr  (lon, 1:nzm+1-m) = co2  (m+1:nzm+1) 
            ch4vmr  (lon, 1:nzm+1-m) = ch4  (m+1:nzm+1) 
            n2ovmr  (lon, 1:nzm+1-m) = n2o  (m+1:nzm+1) 
            o2vmr   (lon, 1:nzm+1-m) = o2   (m+1:nzm+1) 
            cfc11vmr(lon, 1:nzm+1-m) = cfc11(m+1:nzm+1) 
            cfc12vmr(lon, 1:nzm+1-m) = cfc12(m+1:nzm+1) 
            cfc22vmr(lon, 1:nzm+1-m) = cfc22(m+1:nzm+1) 
            ccl4vmr (lon, 1:nzm+1-m) = ccl4 (m+1:nzm+1) 
            o3vmr   (lon, nzm+2-m:nzm+1) = ozone(lon, lat, nzm+1)                               
            co2vmr  (lon, nzm+2-m:nzm+1) = co2  (nzm+1)                               
            ch4vmr  (lon, nzm+2-m:nzm+1) = ch4  (nzm+1)                               
            n2ovmr  (lon, nzm+2-m:nzm+1) = n2o  (nzm+1)                               
            o2vmr   (lon, nzm+2-m:nzm+1) = o2   (nzm+1)                               
            cfc11vmr(lon, nzm+2-m:nzm+1) = cfc11(nzm+1)                               
            cfc12vmr(lon, nzm+2-m:nzm+1) = cfc12(nzm+1)                               
            cfc22vmr(lon, nzm+2-m:nzm+1) = cfc22(nzm+1)                               
            ccl4vmr (lon, nzm+2-m:nzm+1) = ccl4 (nzm+1)                               

            !
            ! Fill out 2D arrays needed by RRTMG 
            !
            layerT(lon, 1:nzm-m) = tabs(lon, lat, m+1:nzm) 
            layerT(lon, nzm-m+1:nzm+1) = tabs(lon, lat, nzm)  

            interfaceT(lon, 2:nzm+1) = (layerT(lon, 1:nzm) + layerT(lon, 2:nzm+1)) / 2. 
            !
            ! Extrapolate temperature at top and bottom interfaces
            !   from lapse rate within the layer
            !
            interfaceT(lon, 1)  = sstxy(lon, lat) + t00 
            !bloss(120709): second option for interfaceT(lon,1):
            !bloss  interfaceT(lon, 1)  = layerT(lon, 1)   + (layerT(lon, 1)   - interfaceT(lon, 2))   
            interfaceT(lon, nzm+2) = 2.*layerT(lon, nzm+1) - interfaceT(lon, nzm+1)

            layerP(lon, 1:nzm-m) = pres (m+1:nzm) 
            layerM(lon, 1:nzm-m) = rho (m+1:nzm)*adz(m+1:nzm)*dz
            if(m.gt.0) then
              coef = exp(-ggr*100./(rgas*tabs(lon, lat, nzm)))
              do k=nzm-m+1,nzm
               layerP(lon, k) = layerP(lon, k-1)*coef 
               layerM(lon, k) = layerM(lon, k-1)*coef 
              end do
            end if
            interfaceP(lon, 1:nz-m) = presi(m+1:nz) 
            if(m.gt.0) then
              do k=nz-m+1,nz
                interfaceP(lon, k) = interfaceP(lon, k-1)*coef
              end do
            end if
        
            ! add layer to top, top pressure <= 0.01 Pa.
            layerP(lon, nzm+1) = 0.5*interfaceP(lon, nz) 
            interfaceP(lon, nz+1) = MIN( 1.e-4_kind_rb, 0.25*layerP(lon,nzm+1)) 
            ! correct so that the interface pressure was in between layer pressures       
            if(m.gt.0) then
              do k=nz-m,nz
                interfaceP(lon, k) = (layerP(lon, k) + layerP(lon, k-1)) / 2.
              end do
              layerP(lon, nzm+1) = 0.5*interfaceP(lon, nz) 
              interfaceP(lon, nz+1) = MIN( 1.e-4_kind_rb, 0.25*layerP(lon,nzm+1)) 
            end if

            ! Convert hPa to Pa in layer mass calculation (kg/m2) 
            layerMass(lon, 1:nzm+1) = &
               100.*(interfaceP(lon,1:nz) - interfaceP(lon,2:nz+1))/ ggr
              
            ! avoid re-sorting rad_reff* based on terrain
            ! fill liqRe and iceRe directly
            liqRe(lon,1:nzm-m) = rad_reffc(lon,lat,m+1:nzm)
            iceRe(lon,1:nzm-m) = rad_reffi(lon,lat,m+1:nzm)

            ! fill layers above model grid with top value
            liqRe(lon,nzm-m+1:nzm) = rad_reffc(lon,lat,nzm)
            iceRe(lon,nzm-m+1:nzm) = rad_reffi(lon,lat,nzm)

            lwHeatingRate(:, :) = 0.; swHeatingRate(:, :) = 0. 

            ! ---------------------------------------------------
            !
            ! Compute cloud IWP/LWP and particle sizes - convert from kg to g
            !
            LWP(lon, 1:nzm-m) = qcl(lon, lat, m+1:nzm) * 1.e3 * layerMass(lon, 1:nzm-m)
            LWP(lon, nzm-m+1:nzm+1) = 0. ! zero out extra layer

            if(ALLOCATED(CloudIceMassMixingRatio)) then
               IWP(lon, 1:nzm-m) = CloudIceMassMixingRatio(lon, lat, m+1:nzm) &
                                    * 1.e3 * layerMass(lon, 1:nzm-m) ! For P3 Ice w/doeiscld==F
            else
               IWP(lon, 1:nzm-m) = qci(lon, lat, m+1:nzm) * 1.e3 * layerMass(lon, 1:nzm-m) 
            end if
            IWP(lon, nzm-m+1:nzm+1) = 0. ! zero out extra layer

            if(do_partial_cloudiness_in_radiation) then
               cloudFrac(lon, 1:nzm-m) = cld(lon, lat, m+1:nzm)
               cloudFrac(lon, nzm-m+1:nzm+1) = 0. ! zero out extra layer
            end if
          end do ! lon
          !
          ! No concentrations for these gases
          ! 
          cfc22vmr(:, :) = 0. 
          ccl4vmr(:, :)  = 0.

          if(.not.do_partial_cloudiness_in_radiation) then
             ! define cloud fraction as zero or one based on presence of liquid and/or ice
             cloudFrac(:,:) = MERGE(1., 0., LWP(:,:)>0. .or. IWP(:,:)>0.)
             tmp_overlap = overlap_no_partial_cloudiness
          else
             tmp_overlap = overlap_partial_cloudiness
          end if

          if(have_cloud_optics) then
             if(use_m2005_cloud_optics) then
              call compute_m2005_cloud_optics(nx, nzm, lat, layerMass, cloudFrac, &
                  cloudTauLW, cloudTauSW, cloudSsaSW, cloudAsmSW, cloudForSW, &
                  cloudTauSW_cldliq, cloudTauSW_cldice, cloudTauSW_snow )
            elseif(use_thompson_cloud_optics) then
              call compute_thompson_cloud_optics(nx, nzm, lat, layerMass, cloudFrac, &
                  cloudTauLW, cloudTauSW, cloudSsaSW, cloudAsmSW, cloudForSW, &
                  cloudTauSW_cldliq, cloudTauSW_cldice, cloudTauSW_snow )
            end if
            !
            ! Normally simulators are run only when the sun is up,
            !    but in case someone decides to use nighttime values...
            !
            ! TODO(2024-03): For now, scale the optical depth headed for the simulators by
            !    the cloud fraction.  In the future, we could have the simulators
            !    work on the MCICA realizations of the possible cloud fraction profiles.
            !
            if(doisccp .or. domodis .or. domisr) then
              ! band 9 is 625 - 778 nm, needed is 670 nm
              tau_067 (1:nx,lat,1:nzm) = cloudTauSW(9,1:nx,1:nzm)*CloudFrac(1:nx,1:nzm)
              tau_067_cldliq (1:nx,lat,1:nzm) = cloudTauSW_cldliq(9,1:nx,1:nzm)*CloudFrac(1:nx,1:nzm)
              tau_067_cldice (1:nx,lat,1:nzm) = cloudTauSW_cldice(9,1:nx,1:nzm)*CloudFrac(1:nx,1:nzm)
              tau_067_snow (1:nx,lat,1:nzm) = cloudTauSW_snow(9,1:nx,1:nzm)*CloudFrac(1:nx,1:nzm)
              ! band 6 is 820 - 980 cm-1, we need 10.5 micron
              emis_105(1:nx,lat,1:nzm) = 1. - exp(-cloudTauLW(6,1:nx,1:nzm)*CloudFrac(1:nx,1:nzm))
            end if

         else ! if(have_cloud_optics) 

!DEBUG ONLY            if(dothirtypercentcloud) cloudFrac(:,:) = 0.3*cloudFrac(:,:)
          
            if(do_partial_cloudiness_in_radiation) then
               ! convert grid-mean LWP and IWP to in-cloud values
               LWP(:,:) = LWP(:,:) / (1.e-6 + CloudFrac(:,:))
               IWP(:,:) = IWP(:,:) / (1.e-6 + CloudFrac(:,:))
            end if
            !
            ! Limit particle sizes to range allowed by RRTMG parameterizations, add top layer
            !
            where(LWP(:,1:nzm) > 0.) & 
              liqRe(:,1:nzm) = max(2.5_kind_rb, min( 60._kind_rb,liqRe(:,1:nzm)))
            where(IWP(:,1:nzm) > 0.) & 
                 iceRe(:,1:nzm) = max(5.0_kind_rb, min(140._kind_rb,iceRe(:,1:nzm)))
                 
            liqRe(:,nzm+1) = 0._kind_rb
            iceRe(:,nzm+1) = 0._kind_rb
            
            if(doisccp .or. domodis .or. domisr) then 
              !
              ! Compute cloud optical depths directly so we can provide to instrument simulators
              !   Ice particle size should be "generalized effective size" from Fu et al. 1998
              !   doi:10.1175/1520-0442(1998)011<2223:AAPOTI>2.0.CO;2
              !   This would normally require some conversion, I guess
              !
              prpLWIn = 0.; prpSWIn = 0. 
              do i = 1, nx
                call cldprop   (nzm+1, 2, 3, 1, cloudFrac(i,:), prpLWIn, &
                                IWP(i,:), LWP(i,:), iceRe(i,:), liqRe(i,:), ncbands, tauLWOut)
                ! Last three output arguments from cldprop_sw are *delta-scaled* optical properties - 
                !   RRTM needs unscaled variables, so we need to provide physical quantities to RRTMG,
                ! which will call cldprop_sw again
                call cldprop_sw(nzm+1, 2, 3, 1, cloudFrac(i,:), &
                                prpSWIn, prpSWIn, prpSWIn, prpSWIn, IWP(i,:), LWP(i,:), iceRe(i,:), liqRe(i,:), &
                                tauSWOut, scaled1, scaled2, scaled3)
                tau_067 (i,lat,1:nzm) =           tauSWOut(1:nzm,24)*CloudFrac(i,1:nzm) ! RRTMG SW bands run from 16 to 29 (parrrsw.f90); we want 9th of these
                                                                     ! band 9 is 625 - 778 nm, needed is 670 nm
                emis_105(i,lat,1:nzm) = 1. - exp(-tauLWOut(1:nzm,6)*CloudFrac(i,1:nzm)) ! band 6 is 820 - 980 cm-1, we need 10.5 micron 
              end do
           end if

            if(domodis) then 
              !
              ! Compute separate cloud optical depths for liquid and ice clouds for input to 
              !   MODIS simulator, which wants these things separately.
              !
              cloudFrac_liq(:,:) = MERGE(CloudFrac(:,:), 0._kind_rb, LWP(:,:)>0.)
              cloudFrac_ice(:,:) = MERGE(CloudFrac(:,:), 0._kind_rb, IWP(:,:)>0.)
              prpLWIn = 0.; prpSWIn = 0.; dummyRe = 0.; dummyWP = 0.
              do i = 1, nx
                ! See above comment.  We want unscaled optical depth from cloud liquid at 670nm
                call cldprop_sw(nzm+1, 2, 3, 1, cloudFrac_liq(i,:), &
                                prpSWIn, prpSWIn, prpSWIn, prpSWIn, dummyWP(i,:), LWP(i,:), dummyRe(i,:), liqRe(i,:), &
                                tauSWOut, scaled1, scaled2, scaled3)
                tau_067_cldliq (i,lat,1:nzm) = tauSWOut(1:nzm,24)*CloudFrac(i,1:nzm) ! RRTMG SW band number 9 (625 - 778 nm), needed is 670 nm

                ! Same for cloud ice
                call cldprop_sw(nzm+1, 2, 3, 1, cloudFrac_ice(i,:), &
                                prpSWIn, prpSWIn, prpSWIn, prpSWIn, IWP(i,:), dummyWP(i,:), iceRe(i,:), dummyRe(i,:), &
                                tauSWOut, scaled1, scaled2, scaled3)
                tau_067_cldice (i,lat,1:nzm) = tauSWOut(1:nzm,24)*CloudFrac(i,1:nzm) ! RRTMG SW band number 9 (625 - 778 nm), needed is 670 nm

                tau_067_snow (i,lat,1:nzm) = 0.! snow is not radiatively active here.
              end do
            end if

         end if ! if(have_cloud_optics)
          ! ---------------------------------------------------
          
          !
          ! Volume mixing fractions for gases.
          !bloss(072009): Note that o3, etc. are now in ppmv and don't need conversions.
          !
          do lon = 1, nx
            m = k_terra(lon,lat)-1
            h2ovmr(lon, 1:nzm-m)   = mwdry/mwh2o * qv(lon, lat, 1+m:nzm) 
            h2ovmr(lon, nzm-m+1:nzm+1)   = h2ovmr(lon, nzm-m) ! extrapolate above model top
          end do

          ! ---------------------------------------------------------------------------------
          if (dolongwave) then

            surfaceT(:) = sstxy(1:nx, lat) + t00

            do lon = 1, nx
             if(landtype(lon,lat).eq.0) then
               surfaceEmissivity(lon,:) = emis_water  ! over water
             else
               surfaceEmissivity(lon,:) = 1. ! emissivity over land is already baked in in sstxy - MK
             end if
            end do
            
            idrv = 0
            duflx_dt(:,:) = 0.
            duflxc_dt(:,:) = 0.

            if(do_partial_cloudiness_in_radiation) then
               call mcica_subcol_lw(1, nx, nzm+1, overlap_partial_cloudiness, permuteseed, irng, layerP, &
                    cloudFrac, IWP, LWP, iceRe, liqRe, CloudTauLW, &
                    CloudFracmcl, IWPmcl, LWPmcl, iceRemcl, liqRemcl, cloudTauLWmcl)

!debug              if(lat.eq.1) then
!debug                 do i = 1,nx
!debug                    write(*,*) 'i = ', i
!debug                    do k = 1,nzm
!debug                       if(cloudFrac(i,k).gt.0.) then
!debug                          write(*,992) k, cloudFrac(i,k), (cloudFracmcl(m,i,k),m=1,10)
!debug                          write(*,993) k, lwp(i,k), (lwpmcl(m,i,k),m=1,10)
!debug                          write(*,994) k, iwp(i,k), (iwpmcl(m,i,k),m=1,10)
!debug  992                     format('k = ',i3,'cld = ', f10.4, ', *mcl = ', 10f10.4)
!debug  993                     format('k = ',i3,'lwp = ', e10.2, ', *mcl = ', 10e10.2)
!debug  994                     format('k = ',i3,'iwp = ', e10.2, ', *mcl = ', 10e10.2)
!debug                       end if
!debug                    end do
!debug                    write(*,*)
!debug                 end do
!debug              end if

               if(have_cloud_optics) then
                  call rrtmg_lw (nx, nzm+1, overlap_partial_cloudiness, idrv,            &
                       layerP, interfaceP, layerT, interfaceT, surfaceT, &
                       h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                       cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, surfaceEmissivity,  &
                       0, 0, 0, cloudFracmcl, &
                       CloudTauLWmcl, IWPmcl, LWPmcl, iceRemcl, liqRemcl, &
                       dummyTauAerosolLW, &
                       lwUp,lwDown, lwHeatingRate, lwUpClearSky, lwDownClearSky, lwHeatingRateClearSky, &
                       duflx_dt, duflxc_dt)
               else
                  call rrtmg_lw (nx, nzm+1, overlap_partial_cloudiness, idrv,            & 
                       layerP, interfaceP, layerT, interfaceT, surfaceT, &
                       h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                       cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, surfaceEmissivity,  &
                       2, 3, 1, cloudFracmcl, &
                       CloudTauLWmcl, IWPmcl, LWPmcl, iceRemcl, liqRemcl, &
                       dummyTauAerosolLW, &
                       lwUp,lwDown, lwHeatingRate, lwUpClearSky, lwDownClearSky, lwHeatingRateClearSky, &
                       duflx_dt, duflxc_dt)
               end if

            else ! do_partial_cloudiness_in_radiation == .false.

               if(have_cloud_optics) then
                  call rrtmg_lw_nomcica (nx, nzm+1, overlap_no_partial_cloudiness, idrv,            &
                       layerP, interfaceP, layerT, interfaceT, surfaceT, &
                       h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                       cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, surfaceEmissivity,  &
                       0, 0, 0, cloudFrac, &
                       CloudTauLW, IWP, LWP, iceRe, liqRe, &
                       dummyTauAerosolLW, &
                       lwUp,lwDown, lwHeatingRate, lwUpClearSky, lwDownClearSky, lwHeatingRateClearSky, &
                       duflx_dt, duflxc_dt)
               else
                  call rrtmg_lw_nomcica (nx, nzm+1, overlap_no_partial_cloudiness, idrv,            & 
                       layerP, interfaceP, layerT, interfaceT, surfaceT, &
                       h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                       cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, surfaceEmissivity,  &
                       2, 3, 1, cloudFrac, &
                       CloudTauLW, IWP, LWP, iceRe, liqRe, &
                       dummyTauAerosolLW, &
                       lwUp,lwDown, lwHeatingRate, lwUpClearSky, lwDownClearSky, lwHeatingRateClearSky, &
                       duflx_dt, duflxc_dt)
               end if

            end if ! if(do_partial_cloudiness_in_radiation)
            
            !bloss: Recompute heating rate using layer density.
            !       This will provide better energy conservation since other
            !       flux difference terms in energy budget are computed this way.
            lwHeatingRate(:,1:nzm) = &
                 (lwUp(:,1:nzm) - lwUp(:,2:nz) + lwDown(:,2:nz) - lwDown(:,1:nzm)) &
             ! incorrect for the case of topography - MK
             !    /spread(cp*rho(1:nzm)*dz*adz(1:nzm), dim=1, ncopies=nx )
             ! Correct way - MK
                /(cp*layerM(:,1:nzm))
            lwHeatingRateClearSky(:,1:nzm) = &
                 (lwUpClearSky(:,1:nzm) - lwUpClearSky(:,2:nz) &
                  + lwDownClearSky(:,2:nz) - lwDownClearSky(:,1:nzm)) &
            !     /spread(cp*rho(1:nzm)*dz*adz(1:nzm), dim=1, ncopies=nx )
                /(cp*layerM(:,1:nzm))
            !
            ! Add fluxes to average-average diagnostics
            !
            do lon = 1,nx
              m = k_terra(lon,lat)-1
              do k=1+m,nz
               radlwup(k) = radlwup(k) + lwUp(lon, k-m)*wgtw(lat,k)*terraw(lon,lat,k)
               radlwdn(k) = radlwdn(k) + lwDown(lon, k-m)*wgtw(lat,k)*terraw(lon,lat,k)
              end do
              do k=1+m,nzm
               radqrlw(k) = radqrlw(k) + lwHeatingRate(lon, k-m)*wgt(lat,k)*terra(lon,lat,k)
               radqrclw(k) = radqrclw(k) + lwHeatingRateClearSky(lon, k-m)*wgt(lat,k)*terra(lon,lat,k)
               qrad(lon, lat, k) = lwHeatingRate(lon, k-m)
              end do
            end do ! lon
            !
            ! 2D diagnostic fields
            !
            lwDownSurface        (:, lat) = lwDown(:, 1)
            lwDownSurfaceClearSky(:, lat) = lwDownClearSky(:, 1)
            lwUpSurface          (:, lat) = lwUp(:, 1)
            lwUpSurfaceClearSky  (:, lat) = lwUpClearSky(:, 1)
            lwUpToa              (:, lat) = lwUp(:, nz+1) !bloss: nz+1 --> TOA
            lwUpTom              (:, lat) = lwUp(:, nz) 
            lwDownTom            (:, lat) = lwDown(:, nz) 
            lwUpToaClearSky      (:, lat) = lwUpClearSky(:, nz+1)

            !bloss: Increment random seed
            permuteseed = 199 + MOD(permuteseed+ngptlw,199)

         end if ! dolongwave
         ! ---------------------------------------------------------------------------------
         
         ! ---------------------------------------------------------------------------------
          if(doshortwave) then

            ! Solar insolation depends on several choices
            !
            !---------------
            if(doseasons) then 
              ! The diurnal cycle of insolation will vary
              ! according to time of year of the current day.
              dayForSW = day 
            else
              ! The diurnal cycle of insolation from the calendar
              ! day on which the simulation starts (day0) will be
              ! repeated throughout the simulation.
              ! add option to run perpetual equinox - MK
              if(doequinox) then
                dayForSW = 80. + day - int(day) 
              else
                dayForSW = int(day0) + day - int(day) 
              end if
            end if 
            !---------------
            if(doperpetual) then
               if (dosolarconstant) then
                  ! fix solar constant and zenith angle as specified
                  ! in prm file.
                  solarZenithAngleCos(:) = cos(zenith_angle * pi/180.)
                  eccf = solar_constant/(1367.)
               else
                  ! perpetual sun (no diurnal cycle) - Modeled after Tompkins
                  solarZenithAngleCos(:) = cos(zenith_angle * pi/180.)
                   ! Adjst solar constant by mean value 
                  eccf = sum(p_factor(:, lat)/solarZenithAngleCos(:)) / real(nx) 
               end if
            else
               call shr_orb_decl (dayForSW, eccen, mvelpp, lambm0, obliqr, delta, eccf)
               solarZenithAngleCos(:) =  &
                 zenith(dayForSW, real(pi * latitude(:, lat)/180., kind_rb), &
                                  real(pi * longitude(:, lat)/180., kind_rb) )
            end if
            !---------------
            ! coszrs is found in params.f90 and used in the isccp simulator
            coszrs = max(0._kind_rb, solarZenithAngleCos(1))
            
            !
            ! We only call the shortwave if the sun is above the horizon. 
            !   We assume that the domain is small enough that the entire 
            !   thing is either lit or shaded
            !
          ! MK  if(all(solarZenithAngleCos(:) >= tiny(solarZenithAngleCos))) then 

              if(lat.eq.1.AND.masterproc) print *, "Let's do some shortwave" 
              if(SLM) then
                 call albedo_slm(nx, real(landmask(1:nx,lat)), real(icemask(1:nx,lat)), &
                          solarZenithAngleCos(:), real(snow_mass(1:nx,lat)), &
                          real(albedovis_v(1:nx,lat)), real(albedonir_v(1:nx,lat)),  &
                          real(albedovis_s(1:nx,lat)), real(albedonir_s(1:nx,lat)), &
                          real(soilw(1,1:nx,lat)),  real(phi_1(1:nx,lat)), real(phi_2(1:nx,lat)), &
                          real(LAI(1:nx,lat)), real(IMPERV(1:nx,lat)), & 
                          asdir(:), aldir(:), asdif(:), aldif(:))
              else
                 ! original albedo of standard SAM without SLM
                 call albedo(ocean, solarZenithAngleCos(:), surfaceT, &
                          asdir(:), aldir(:), asdif(:), aldif(:)           )
              end if
              if(lat.eq.1.AND.masterproc) then
                print *, "Range of zenith angles", minval(solarZenithAngleCos), maxval(solarZenithAngleCos)
                print *, "Range of surface albedo (asdir)", minval(asdir), maxval(asdir)
                print *, "Range of surface albedo (aldir)", minval(aldir), maxval(aldir)
                print *, "Range of surface albedo (asdif)", minval(asdif), maxval(asdif)
                print *, "Range of surface albedo (aldif)", minval(aldif), maxval(aldif)
                print *, "Range of cloud fraction", minval(cloudFrac), maxval(cloudFrac)

             end if
             
              if(do_partial_cloudiness_in_radiation) then

                 call mcica_subcol_sw(1, nx, nzm+1, overlap_partial_cloudiness, permuteseed, irng, layerP, &
                      cloudFrac, IWP, LWP, iceRe, liqRe, &
                      cloudTauSW, cloudSsaSW, cloudAsmSW, cloudForSW, &
                      CloudFracmcl, IWPmcl, LWPmcl, iceRemcl, liqRemcl, &
                      cloudTauSWmcl, cloudSsaSWmcl, cloudAsmSWmcl, cloudForSWmcl)

                 if(have_cloud_optics) then
                    call rrtmg_sw(nx, nzm+1, overlap_partial_cloudiness,                     &
                         layerP, interfaceP, layerT, interfaceT, surfaceT, &
                         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
                         asdir, asdif, aldir, aldif, &
                         solarZenithAngleCos, eccf, 0, scon,   &
                         0, 0, 0, cloudFracmcl, &
                         cloudTauSWmcl, cloudSsaSWmcl, cloudAsmSWmcl, cloudForSWmcl, &
                         IWPmcl, LWPmcl, iceRemcl, liqRemcl,  &
                         dummyAerosolProps, dummyAerosolProps, dummyAerosolProps, dummyAerosolProps2, &
                         swUp, swDown, swHeatingRate, swUpClearSky, swDownClearSky, swHeatingRateClearSky, &
                         visDown, visdDown, nirDown, nirdDown)
                 else 
                    call rrtmg_sw(nx, nzm+1, overlap_partial_cloudiness,                     & 
                         layerP, interfaceP, layerT, interfaceT, surfaceT, &
                         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
                         asdir, asdif, aldir, aldif, &
                         solarZenithAngleCos, eccf, 0, scon,   &
                         2, 3, 1, cloudFracmcl, &
                         cloudTauSWmcl, cloudSsaSWmcl, cloudAsmSWmcl, cloudForSWmcl, &
                         IWPmcl, LWPmcl, iceRemcl, liqRemcl,  &
                         dummyAerosolProps, dummyAerosolProps, dummyAerosolProps, dummyAerosolProps2, &
                         swUp, swDown, swHeatingRate, swUpClearSky, swDownClearSky, swHeatingRateClearSky, &
                         visDown, visdDown, nirDown, nirdDown)
                 end if

              else ! do_partial_cloudiness_in_radiation == .false.

                 if(have_cloud_optics) then
                    call rrtmg_sw_nomcica(nx, nzm+1, overlap_no_partial_cloudiness,                     &
                         layerP, interfaceP, layerT, interfaceT, surfaceT, &
                         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
                         asdir, asdif, aldir, aldif, &
                         solarZenithAngleCos, eccf, 0, scon,   &
                         0, 0, 0, cloudFrac, &
                         cloudTauSW, cloudSsaSW, cloudAsmSW, cloudForSW, &
                         IWP, LWP, iceRe, liqRe,  &
                         dummyAerosolProps, dummyAerosolProps, dummyAerosolProps, dummyAerosolProps2, &
                         swUp, swDown, swHeatingRate, swUpClearSky, swDownClearSky, swHeatingRateClearSky, &
                         visDown, visdDown, nirDown, nirdDown)
                 else 
                    call rrtmg_sw_nomcica(nx, nzm+1, overlap_no_partial_cloudiness,                     & 
                         layerP, interfaceP, layerT, interfaceT, surfaceT, &
                         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
                         asdir, asdif, aldir, aldif, &
                         solarZenithAngleCos, eccf, 0, scon,   &
                         2, 3, 1, cloudFrac, &
                         cloudTauSW, cloudSsaSW, cloudAsmSW, cloudForSW, &
                         IWP, LWP, iceRe, liqRe,  &
                         dummyAerosolProps, dummyAerosolProps, dummyAerosolProps, dummyAerosolProps2, &
                         swUp, swDown, swHeatingRate, swUpClearSky, swDownClearSky, swHeatingRateClearSky, &
                         visDown, visdDown, nirDown, nirdDown)
                 end if

              end if ! if(do_partial_cloudiness_in_radiation)
            
              !bloss: Recompute heating rate using layer density.
              !       This will provide better energy conservation since other
              !       flux difference terms in energy budget are computed this way.
              swHeatingRate(:,1:nzm) = &
                   (swUp(:,1:nzm) - swUp(:,2:nz) + swDown(:,2:nz) - swDown(:,1:nzm)) &
                 ! incorrect in the case of topography -MK
                 !  /spread(cp*rho(1:nzm)*dz*adz(1:nzm), dim=1, ncopies=nx )
                 ! correct way:  -MK
                   /(cp*layerM(:,1:nzm))
              swHeatingRateClearSky(:,1:nzm) = &
                   (swUpClearSky(:,1:nzm) - swUpClearSky(:,2:nz) &
                   + swDownClearSky(:,2:nz) - swDownClearSky(:,1:nzm)) &
                 !  /spread(cp*rho(1:nzm)*dz*adz(1:nzm), dim=1, ncopies=nx )
                   /(cp*layerM(:,1:nzm))
  
              !
              ! Add fluxes to average-average diagnostics
              !
              do lon=1,nx
               m = k_terra(lon,lat)-1
               do k=1+m,nz
                radswup(k) = radswup(k) + swUp(lon, k-m)*wgtw(lat,k)*terraw(lon,lat,k)
                radswdn(k) = radswdn(k) + swDown(lon, k-m)*wgtw(lat,k)*terraw(lon,lat,k)
               end do
               do k=1+m,nzm
                radqrsw(k) = radqrsw(k) + swHeatingRate(lon, k-m)*wgt(lat,k)*terra(lon,lat,k)
                radqrcsw(k) = radqrcsw(k) + swHeatingRateClearSky(lon, k-m)*wgt(lat,k)*terra(lon,lat,k)
                qrad(lon, lat, k) = qrad(lon, lat, k) + swHeatingRate(lon, k-m)
               end do
              end do ! lon
              !
              ! 2D diagnostic fields
              !
              swDownSurface        (:, lat) = swDown(:, 1)
              swDownSurfaceClearSky(:, lat) = swDownClearSky(:, 1)
              swUpSurface          (:, lat) = swUp(:, 1)
              swUpSurfaceClearSky  (:, lat) = swUpClearSky(:, 1)
              swDownToa            (:, lat) = swDown(:, nz+1) !bloss: nz+1 --> TOA
              swDownTom            (:, lat) = swDown(:, nz) 
              swUpToa              (:, lat) = swUp(:, nz+1)
              swUpTom              (:, lat) = swUp(:, nz)
              swUpToaClearSky      (:, lat) = swUpClearSky(:, nz+1) 
              insolation_TOA       (:, lat) = swDown(:, nz+1)
              visDownSurface       (:, lat) = visDown(:, 1)
              visDownSurfaceDiffuse(:, lat) = visdDown(:, 1)
              nirDownSurface       (:, lat) = nirDown(:, 1)
              nirDownSurfaceDiffuse(:, lat) = nirdDown(:, 1)
              CosZenithAngle(:,lat) = solarZenithAngleCos(:)
          !MK  else
          !    !
          !    ! 2D diagnostic fields - nighttime values
          !    !
          !    swDownSurface        (:, lat) = 0.0
          !    swDownSurfaceClearSky(:, lat) = 0.0
          !    swUpSurface          (:, lat) = 0.0
          !    swUpSurfaceClearSky  (:, lat) = 0.0
          !    swDownToa            (:, lat) = 0.0
          !    swDownTom            (:, lat) = 0.0
          !    swUpToa              (:, lat) = 0.0
          !    swUpTom              (:, lat) = 0.0
          !    swUpToaClearSky      (:, lat) = 0.0
          !    insolation_TOA       (:, lat) = 0.0
          !    visDownSurface       (:, lat) = 0.0
          !    visDownSurfaceDiffuse(:, lat) = 0.0
          !    nirDownSurface       (:, lat) = 0.0
          !    nirDownSurfaceDiffuse(:, lat) = 0.0
          !    CosZenithAngle(:,lat) = solarZenithAngleCos(:)
          !  end if 

          !bloss: Increment random seed
          permuteseed = 199 + MOD(permuteseed+ngptsw,199)

         end if ! do shortwave
       ! ---------------------------------------------------------------------------------
        end do ! Loop over y dimension
! ---------------------------------------------------------------------------------

        if(do_partial_cloudiness_in_radiation) then
           DEALLOCATE(LWPmcl, IWPmcl, cloudFracmcl, liqRemcl, iceRemcl, &
                cloudTauLWmcl, cloudTauSWmcl, cloudSsaSWmcl, cloudAsmSWmcl, &
                cloudForSWmcl, STAT=ierr)
           if(ierr.ne.0) then
              write(*,*) 'Error: Cannot deallocate McICA arrays in RAD_RRTM/rad.f90'
              call task_abort()
           end if
        end if

        !
        ! 2D diagnostics
        !
        nradsteps = 0 ! re-initialize nradsteps
        
        if(masterproc) then 
          if(doshortwave) then 
            if(doperpetual) then
              print *,'radiation: perpetual sun, solin=', sum(swDownToa(:, :)) / float(nx*ny)
            else
              print *,'radiation: coszrs=', coszrs,&
                      ' solin=', sum(swDownToa(:, :)) / float(nx*ny)
            end if
          end if
          if(dolongwave) print *,'longwave radiation is called'
        end if
        
        if(doradhomozonal) then
          call mean_x_3D(qrad,tmpz)
            do k=1,nzm
             do j=1,ny
               qrad(:,j,k) = tmpz(j,k)
             end do
            end do
         end if

        
        if(doradhomo) then    
          !
          ! Homogenize radiation if desired
          !
      !    fortran's sum() function should be used with caution. It can be quite inaccurate. -MK
      !    Besides, the averaging should be waighted by latitude.
      !    radHeatingProfile(1:nzm) = sum(sum(qrad(:, :, :), dim = 1), dim = 1) / (nx * ny) 

           do k=1,nzm
            radHeatingProfile(k) = 0.
            do j=1,ny
             do i=1,nx
              radHeatingProfile(k) = radHeatingProfile(k) + qrad(i,j,k)*wgt(j,k)
             end do
            end do
            radHeatingProfile(k) = radHeatingProfile(k) /real(nx*ny)
           end do

          
          !
          ! Homogenize across the entire domain
          !
          if(dompi) then
            tempProfile(1:nzm) = radHeatingProfile(1:nzm)
            call task_sum_real8(radHeatingProfile, tempProfile, nzm)
            radHeatingProfile(1:nzm) = tempProfile(1:nzm) / real(nsubdomains)
          end if 
          
          qrad(:, 1, :) = spread(radHeatingProfile(1:nzm), dim = 1, ncopies = nx)
          if(ny > 1) &
            qrad(:, :, :) = spread(qrad(:, 1, :), dim = 2, ncopies = ny)
        end if
      end if ! nradsteps >= nrad
      
      !------------------------------------------------------------------------
      !
      ! Update 2d diagnostic fields 
      !
      ! Net surface and toa fluxes
      !
      ! First two for ocean evolution
      lwnsxy(:, :) = lwUpSurface(:, :) - lwDownSurface(:, :)  ! Net LW 
      swnsxy(:, :) = swDownSurface(:, :) - swUpSurface(:, :)  ! Net SW 
      lwdsxy(:, :) = lwDownSurface(:, :)  ! LW downwards
      swdsxy(:, :) = swDownSurface(:, :)   ! SW downwards

      swdsvisxy(:,:) = visDownSurface(:,:)
      swdsvisdxy(:,:) = visDownSurfaceDiffuse(:,:)
      swdsnirxy(:,:) = nirDownSurface(:,:)
      swdsnirdxy(:,:) = nirDownSurfaceDiffuse(:,:)
      coszrsxy(:,:) = CosZenithAngle(:,:)
     
      initrad = .false. ! needed for some other modules that use surface radiation fluxes
 
      lwns_xy(:, :) = lwns_xy(:, :) + &
                          lwUpSurface(:, :) - lwDownSurface(:, :)  ! Net LW upwards
      swns_xy(:, :) = swns_xy(:, :) + &
                          swDownSurface(:, :) - swUpSurface(:, :)  ! New SW downwards
      lwnsa_xy(:, :) = lwnsa_xy(:, :) + &
                          lwUpSurface(:, :) - lwDownSurface(:, :)  ! Net LW upwards
      swnsa_xy(:, :) = swnsa_xy(:, :) + &
                          swDownSurface(:, :) - swUpSurface(:, :)  ! New SW downwards
      lwds_xy(:, :) = lwds_xy(:, :) - lwDownSurface(:, :)  ! Net LW upwards
      lwdsc_xy(:, :) = lwdsc_xy(:, :) - lwDownSurfaceClearSky(:, :)  ! Clear sky LW upwards
      swds_xy(:, :) = swds_xy(:, :) + swDownSurface(:, :)  ! New SW downwards
      swdsc_xy(:, :) = swdsc_xy(:, :) + swDownSurfaceClearSky(:, :)  ! Clear sky SW downwards
      lwdsa_xy(:, :) = lwdsa_xy(:, :) - lwDownSurface(:, :)  ! Net LW upwards
      swdsa_xy(:, :) = swdsa_xy(:, :) + swDownSurface(:, :)  ! New SW downwards
      ! Net LW at Toa is upwards 
      lwnt_xy(:, :) = lwnt_xy(:, :) + lwUpToa(:, :) 
      lwnta_xy(:, :) = lwnta_xy(:, :) + lwUpToa(:, :) 
      ! Net SW at Toa  
      swnt_xy(:, :) = swnt_xy(:, :) + swDownToa(:, :) - swUpToa(:, :) 
      swnta_xy(:, :) = swnta_xy(:, :) + swDownToa(:, :) - swUpToa(:, :) 
     
      !
      ! Net surface and toa clear sky fluxes
      ! 
      lwnsc_xy(:, :) = lwnsc_xy(:, :) + &
                           lwUpSurfaceClearSky(:, :) - lwDownSurfaceClearSky(:, :) 
      swnsc_xy(:, :) = swnsc_xy(:, :) + &
                           swDownSurfaceClearSky(:, :) - swUpSurfaceClearSky(:, :)
      lwntc_xy(:, :) = lwntc_xy(:, :) +                           lwUpToaClearSky(:, :) 
      swntc_xy(:, :) = swntc_xy(:, :) + swDownToa(:, :) - swUpToaClearSky(:, :)
      
      ! TOA Insolation
      solin_xy(:, :) = solin_xy(:, :) + swDownToa(:, :) 
      
      
      !------------------------------------------------------------------------
      !
      ! Update 1D diagnostics
      ! 
    
      if(dostatisrad) then
        s_flns = s_flns + sum((lwUpSurface(:, :) - lwDownSurface(:, :))*wgtxys(:,:)) ! lwnsxy
        s_fsns = s_fsns + sum((swDownSurface(:, :) - swUpSurface(:, :))*wgtxys(:,:)) ! swnsxy
        s_flntoa = s_flntoa + sum(lwUpToa(:, :)*wgtxyt(:,:))                       ! lwntxy
        s_flnt = s_flnt + sum((lwUpTom(:, :) - lwDownTom(:, :))*wgtxyt(:,:))         ! lwntmxy
        s_fsntoa = s_fsntoa + sum((swDownToa(:, :) - swUpToa(:, :))*wgtxyt(:,:))         ! swntxy
        s_fsnt = s_fsnt + sum((swDownTom(:, :) - swUpTom(:, :))*wgtxyt(:,:))         ! swntxy
        s_flnsc = s_flnsc + &
          sum((lwUpSurfaceClearSky(:, :) - lwDownSurfaceClearSky(:, :))*wgtxys(:,:)) ! lwnscxy
        s_fsnsc = s_fsnsc + &
          sum((swDownSurfaceClearSky(:, :) - swUpSurfaceClearSky(:, :))*wgtxys(:,:)) ! swnscxy 
        s_flntoac = s_flntoac + &
          sum(lwUpToaClearSky(:, :)*wgtxyt(:,:))                                   ! lwntcxy 
        s_fsntoac = s_fsntoac + &
          sum((swDownToa(:, :) - swUpToaClearSky(:, :))*wgtxyt(:,:))                 ! swntcxy
        s_solin = s_solin + sum(swDownToa(:, :)*wgtxyt(:,:))                       ! solinxy 
        ! 
        ! I think the next two are supposed to be downwelling fluxes at the surface
        !
        s_fsds = s_fsds + sum(swDownSurface(:, :)*wgtxys(:,:)) 
        s_flds = s_flds + sum(lwDownSurface(:, :)*wgtxys(:,:)) 
        s_fsdsc = s_fsdsc + sum(swDownSurfaceClearSky(:, :)*wgtxys(:,:)) 
        s_fldsc = s_fldsc + sum(lwDownSurfaceClearSky(:, :)*wgtxys(:,:)) 
      end if ! if(dostatis)
       
      if(mod(nstep,nrestart_steps*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) &
                 call write_rad() ! write radiation restart file
  
    end if ! if icycle == 1

    !
    ! Add radiative heating to liquid ice static energy variable
    !
    t(1:nx, 1:ny, 1:nzm) = t(1:nx, 1:ny, 1:nzm) + qrad(:, :, :) * dtn
  end subroutine rad_driver
  ! ----------------------------------------------------------------------------
  subroutine initialize_radiation

    implicit none

    real(KIND=kind_rb) :: cpdair

    !bloss  subroutine shr_orb_params
    !bloss  inputs:  iyear, log_print
    !bloss  ouptuts: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
    
    call shr_orb_params(iyear, eccen, obliq, mvelp, obliqr, lambm0, mvelpp, .false.)
 
    ! sets up initial mixing ratios of trace gases.
    if(nrestart.eq.0) call tracesini()
 
    if(nrestart == 0) then
      qrad    (:, :, :) = 0.
      radlwup(:) = 0.
      radlwdn(:) = 0.
      radswup(:) = 0.
      radswdn(:) = 0.
      radqrlw(:) = 0.
      radqrsw(:) = 0.
      nradsteps = 0
    else
       call read_rad()
    endif
 
    if(doperpetual) then
      ! perpetual sun (no diurnal cycle)
      p_factor(:, :) = perpetual_factor(real(day0, kind_rb), &
                                        real(latitude(1:nx, 1:ny), kind_rb), &
                                        real(longitude(1:nx, 1:ny), kind_rb) )
    end if
 
    cpdair = cp
    call rrtmg_sw_ini(cpdair)
    call rrtmg_lw_ini(cpdair)
    
    if(trim(micro_scheme_name()) == 'm2005' .and. & 
       (compute_reffc .or. compute_reffi) .and. &
       (.NOT.dorrtm_cloud_optics_from_effrad_LegacyOption)) then
       call m2005_cloud_optics_init
       use_m2005_cloud_optics = .true.
       have_cloud_optics = .true.
    end if

    if( (trim(micro_scheme_name()) == 'thompson') .and. &
       (compute_reffc .or. compute_reffi) .and. &
       (.NOT.dorrtm_cloud_optics_from_effrad_LegacyOption)) then
       call thompson_cloud_optics_init
       use_thompson_cloud_optics = .true.
       have_cloud_optics = .true.
    end if

    if(dosubgridcloudfraction) then
       do_partial_cloudiness_in_radiation = .true.
    end if

    land_frac = MERGE(0., 1., ocean)
    initialized = .true.
  end subroutine initialize_radiation
  ! ----------------------------------------------------------------------------
  !
  ! Trace gas profiles
  !
  ! ----------------------------------------------------------------------------
  subroutine tracesini()
    use netcdf

    implicit none
    !
    ! Initialize trace gaz vertical profiles
    !   The files read from the top down 
    !
    !bloss(072009): Get trace gas profiles from rrtmg_lw.nc, the data
    !                 file provided with RRTMG.  These are indexed from
    !                 bottom to top and are in ppmv, so that no conversion
    !                 is needed for use with RRTMG.
    !
    integer k, m, ierr
    real :: godp ! gravity over delta pressure
    real :: plow, pupp
    real(kind=kind_rb), dimension(nzm+1) :: tmp_pres ! level pressure with extra level at TOA
    real(kind=kind_rb), dimension(nz+1) :: tmp_presi ! interface pressure with extra level at TOA
    integer(kind=kind_im) :: ncid, varID, dimIDab, dimIDp

    integer(kind=kind_im) :: Nab, nPress, ab
    real(kind=kind_rb) :: wgtlow, wgtupp, pmid
    real(kind=kind_rb), allocatable, dimension(:) :: pMLS
    real(kind=kind_rb), allocatable, dimension(:,:) :: trace, trace_in
    character(LEN=nf90_max_name) :: tmpName

    integer, parameter :: nTraceGases = 9
    real(kind=kind_rb), dimension(nzm+1) :: tmpTrace
    real(kind=kind_rb), dimension(nz+1,nTraceGases) :: trpath
    character(len = maxAbsorberNameLength), dimension(nTraceGases), parameter :: &
         TraceGasNameOrder = (/        &
     				'O3   ',  &
     				'CO2  ',  &
     				'CH4  ',  &
     				'N2O  ',  & 
     				'O2   ',  &
     				'CFC11',  &
     				'CFC12',  &
     				'CFC22',  &
     				'CCL4 '  /)
	
    real factor
    ! ---------------------------------

    !bloss: RRTMG radiation orders levels from surface to top of model.
    !       This routine was originally written for CCM/CAM radiation
    !       which orders levels from top down.  As a result, we need to
    !       reverse the ordering here to make things work for RRTMG.

!DD add code to read from master only and bcast trace profiles
  if(masterproc) then

    ! Read profiles from rrtmg data file.
    status(:)   = nf90_NoErr
    status(1)   = nf90_open(trim(rundatadir)//'/rrtmg_lw.nc',nf90_nowrite,ncid)
	
    status(2)   = nf90_inq_dimid(ncid,"Pressure",dimIDp)
    status(3)   = nf90_inquire_dimension(ncid, dimIDp, tmpName, nPress)

    status(4)   = nf90_inq_dimid(ncid,"Absorber",dimIDab)
    status(5)   = nf90_inquire_dimension(ncid, dimIDab, tmpName, Nab)

    allocate(pMLS(nPress), trace(nTraceGases,nPress), trace_in(Nab,nPress), STAT=ierr)
    pMLS = 0.
    trace = 0.
    if(ierr.ne.0) then
      write(*,*) 'ERROR: could not declare arrays in tracesini'
      call task_abort()
    end if

    status(6)   = nf90_inq_varid(ncid,"Pressure",varID)
    status(7)   = nf90_get_var(ncid, varID, pMLS)

    status(8)   = nf90_inq_varid(ncid,"AbsorberAmountMLS",varID)
    status(9) = nf90_get_var(ncid, varID, trace_in)

    do m = 1,nTraceGases
      call getAbsorberIndex(TRIM(tracegasNameOrder(m)),ab)
      trace(m,1:nPress) = trace_in(ab,1:nPress)
      where (trace(m,:)>2.)
        trace(m,:) = 0.
      end where
    end do

    if(MAXVAL(ABS(status(1:8+nTraceGases))).ne.0) then
      write(*,*) 'Error in reading trace gas sounding from'//trim(rundatadir)//'/rrtmg_lw.nc'
      call task_abort()
    end if

    !bloss(120409): copy level, interface pressures into local variable.
    tmp_pres(1:nzm) = pres(1:nzm)
    tmp_presi(1:nz) = presi(1:nz)

    ! Add level at top of atmosphere (top interface <= 0.01 Pa)
    tmp_pres(nz) = 0.5*presi(nz)
    tmp_presi(nz+1) = MIN(1.e-4_kind_rb, 0.25*tmp_pres(nz))

    !bloss: modify routine to compute trace gas paths from surface to
    ! top of atmosphere.  Then, interpolate these paths onto the
    !  interface pressure levels of the model grid, with an extra level
    !  between the model top and the top of atmosphere.  Differencing these
    !   paths and dividing by dp/g will give the mean mass concentration
    !    in that level.
    !
    !  This procedure has the advantage that the total trace gas path
    !   will be invariant to changes in the vertical grid.

    ! trace gas paths at surface are zero.
    trpath(1,:) = 0.

    do k = 2,nz+1
      ! start with trace path at interface below.
      trpath(k,:) = trpath(k-1,:)

      ! if pressure greater than sounding, assume concentration at bottom.
      if (tmp_presi(k-1).gt.pMLS(1)) then
        trpath(k,:) = trpath(k,:) &
             + (tmp_presi(k-1) - MAX(tmp_presi(k),pMLS(1)))/ggr & ! dp/g
             *trace(:,1)                                 ! *tr
      end if

      do m = 2,nPress
        ! limit pMLS(m:m-1) so that they are within the model level
        !  tmp_presi(k-1:k).
        plow = MIN(tmp_presi(k-1),MAX(tmp_presi(k),pMLS(m-1)))
        pupp = MIN(tmp_presi(k-1),MAX(tmp_presi(k),pMLS(m)))

        if(plow.gt.pupp) then
          pmid = 0.5*(plow+pupp)

          wgtlow = (pmid-pMLS(m))/(pMLS(m-1)-pMLS(m))
          wgtupp = (pMLS(m-1)-pmid)/(pMLS(m-1)-pMLS(m))
!!$          write(*,*) pMLS(m-1),pmid,pMLS(m),wgtlow,wgtupp
          ! include this level of the sounding in the trace gas path
          trpath(k,:) = trpath(k,:) &
               + (plow - pupp)/ggr*(wgtlow*trace(:,m-1)+wgtupp*trace(:,m)) ! dp/g*tr
        end if
      end do

      ! if pressure is off top of trace gas sounding, assume
      !  concentration at top
      if (tmp_presi(k).lt.pMLS(nPress)) then
        trpath(k,:) = trpath(k,:) &
             + (MIN(tmp_presi(k-1),pMLS(nPress)) - tmp_presi(k))/ggr & ! dp/g
             *trace(:,nPress)                               ! *tr
      end if

    end do

    if(notracegases) then
     factor=0.01  ! don't make compte zeros as it may blow up the code
    else
     factor=1.
    end if
      

    do m = 1,nTraceGases
      do k = 1,nzm+1
        godp = ggr/(tmp_presi(k) - tmp_presi(k+1))
        tmpTrace(k) = (trpath(k+1,m) - trpath(k,m))*godp
      end do
      if(TRIM(TraceGasNameOrder(m))=='O3') then
        o3(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='CO2') then
        if(docurrentco2.and.year.gt.2059) then
          co2(1:nzm+1) = 1.e-6*co2_ppm(year)
        else
          co2(1:nzm+1) = tmpTrace(1:nzm+1)*nxco2
        end if
      elseif(TRIM(TraceGasNameOrder(m))=='CH4') then
        ch4(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='N2O') then
        n2o(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='O2') then
        o2(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='CFC11') then
        cfc11(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='CFC12') then
        cfc12(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='CFC22') then
        cfc22(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      elseif(TRIM(TraceGasNameOrder(m))=='CCL4') then
        ccl4(1:nzm+1) = tmpTrace(1:nzm+1)*factor
      end if
    
    end do

    print*,'RRTMG rrtmg_lw.nc trace gas profile: number of levels=',nPress
    print*,'gas traces vertical profiles (ppmv):'
    print*,' p (hPa) ', ('       ',TraceGasNameOrder(m),m=1,nTraceGases)
    do k=1,nzm+1
      write(*,999) tmp_pres(k),o3(k),co2(k),ch4(k),n2o(k),o2(k), &
           cfc11(k),cfc12(k), cfc22(k),ccl4(k)
999   format(f8.2,12e12.4)
    end do
    print*,'done...'

    deallocate(pMLS, trace, STAT=ierr)
  endif
    
  if(dompi) then
      call task_bcast_real8(0,o3,nzm+1)
      call task_bcast_real8(0,co2,nzm+1)
      call task_bcast_real8(0,ch4,nzm+1)
      call task_bcast_real8(0,n2o,nzm+1)
      call task_bcast_real8(0,o2,nzm+1)
      call task_bcast_real8(0,cfc11,nzm+1)
      call task_bcast_real8(0,cfc12,nzm+1)
      call task_bcast_real8(0,cfc22,nzm+1)
      call task_bcast_real8(0,ccl4,nzm+1)
  end if

  end subroutine tracesini
  ! ----------------------------------------------------------------------------
  !
  ! Astronomy-related procedures 
  ! 
  ! ----------------------------------------------------------------------------
  elemental real(kind_rb) function zenith(calday, clat, clon)
     implicit none
     real(kind_rb), intent(in ) :: calday, & ! Calendar day, including fraction
                                   clat,   & ! Current centered latitude (radians)
                                   clon      ! Centered longitude (radians)

     real(kind_rb)     :: delta, & ! Solar declination angle in radians
                          eccf
     integer  :: i     ! Position loop index

     call shr_orb_decl (calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
     !
     ! Compute local cosine solar zenith angle
     !
     zenith = shr_orb_cosz(calday, clat, clon, delta)
  end function zenith
  ! ----------------------------------------------------------------------------
  elemental real(kind_rb) function perpetual_factor(day, lat, lon)
    implicit none
    real(kind_rb), intent(in) :: day, lat, lon ! Day (without fraction); centered lat/lon (degrees) 
    real(kind_rb)     :: delta, & ! Solar declination angle in radians
                         eccf
    
    !  estimate the factor to multiply the solar constant
    !  so that the sun hanging perpetually right above
    !  the head (zenith angle=0) would produce the same
    !  total input the TOA as the sun subgect to diurnal cycle.
    !  coded by Marat Khairoutdinov, 2004
    
    ! Local:
    real(kind_rb) :: tmp
    real(kind_rb) :: dttime, dayy 
    real(kind_rb) :: coszrs, ttime
    real(kind_rb) :: clat, clon
    
    if(doequinox) then
      dayy = 80.
    else
      dayy = day
    end if
    ttime = dayy
    dttime = dt*float(nrad)/86400.
    tmp = 0.
    
    clat = pi * lat/180.
    clon = pi * lon/180.
    do while (ttime.lt.dayy+1.)
      call shr_orb_decl (ttime, eccen, mvelpp, lambm0, obliqr, delta, eccf)
       coszrs = zenith(ttime, clat, clon)
       tmp = tmp + min(dttime, dayy+1. - ttime)*max(0._kind_rb, eccf * coszrs)
       ttime = ttime+dttime
    end do
    
    perpetual_factor = tmp
    
  end function perpetual_factor
  ! ----------------------------------------------------------------------------
  !
  ! Writing and reading binary restart files
  !
  ! ----------------------------------------------------------------------------
  subroutine write_rad()
    integer :: irank, ii

    !bloss: added a bunch of statistics-related stuff to the restart file
    !         to nicely handle the rare case when nrad exceeds nstat and 
    !         the model restarts with mod(nstep,nrad)~=0.  This would cause
    !         many of the radiation statistics to be zero before the next
    !         multiple of nrad.

    if(masterproc) print*,'Writting radiation restart file...'

    if(restart_sep) then
      open(56, file = trim(constructRestartFileName(case, caseId, rank)), &
           status='unknown',form='unformatted')
      write(56) nsubdomains
	  write(56) nradsteps, qrad, ozone, radlwup, radlwdn, radswup, radswdn, &
        radqrlw, radqrsw, radqrclw, radqrcsw, &
        lwDownSurface, lwDownSurfaceClearSky, lwUpSurface, lwUpSurfaceClearSky, & 
        lwUpToa,     lwUpToaClearSky,     &
        swDownSurface, swDownSurfaceClearSky, swUpSurface, swUpSurfaceClearSky, & 
        swDownToa,                            swUpToa,     swUpToaClearSky,     &
        insolation_TOA, &
          o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
        visDownSurface, visDownSurfaceDiffuse, nirDownSurface, nirDownSurfaceDiffuse, &
        CosZenithAngle 
      close(56)
    else
      do irank = 0, nsubdomains-1
        call task_barrier()
        if(irank == rank) then
          open(56, file = trim(constructRestartFileName(case, caseId, nSubdomains)), &
               status='unknown',form='unformatted')
          if(masterproc) then
            write(56) nsubdomains
          else
            read (56)
            do ii=0,irank-1 ! skip records
              read(56)
            end do
          end if
          write(56) nradsteps, qrad , ozone, radlwup, radlwdn, radswup, radswdn, &
               radqrlw, radqrsw, radqrclw, radqrcsw, &
               lwDownSurface, lwDownSurfaceClearSky, lwUpSurface, lwUpSurfaceClearSky, & 
               lwUpToa,     lwUpToaClearSky,     &
               swDownSurface, swDownSurfaceClearSky, swUpSurface, swUpSurfaceClearSky, & 
               swDownToa,                            swUpToa,     swUpToaClearSky,     &
               insolation_TOA, &
          o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
          visDownSurface, visDownSurfaceDiffuse, nirDownSurface, nirDownSurfaceDiffuse, &
          CosZenithAngle
          close(56)
        end if
      end do

    end if ! restart_sep

	if(masterproc) print *,'Saved radiation restart file. nstep=',nstep
    call task_barrier()
  end subroutine write_rad
  ! ----------------------------------------------------------------------------
  subroutine read_rad()
    integer ::  irank, ii

    if(masterproc) print*,'Reading radiation restart file...'

    if(restart_sep) then
    
      if(nrestart.ne.2) then
        open(56, file = trim(constructRestartFileName(case, caseid, rank)), &
             status='unknown',form='unformatted')
      else
        open(56, file = trim(constructRestartFileName(case_restart, caseid_restart, rank)), &
             status='unknown',form='unformatted')
      end if
      read (56)
      read(56) nradsteps, qrad, ozone, radlwup, radlwdn, radswup, radswdn, &
        radqrlw, radqrsw, radqrclw, radqrcsw, &
        lwDownSurface, lwDownSurfaceClearSky, lwUpSurface, lwUpSurfaceClearSky, & 
        lwUpToa,     lwUpToaClearSky,     &
        swDownSurface, swDownSurfaceClearSky, swUpSurface, swUpSurfaceClearSky, & 
        swDownToa,                            swUpToa,     swUpToaClearSky,     &
        insolation_TOA, &
          o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
        visDownSurface, visDownSurfaceDiffuse, nirDownSurface, nirDownSurfaceDiffuse, &
        CosZenithAngle
      close(56)
      
    else
    
      do irank=0,nsubdomains-1
        call task_barrier()
        if(irank == rank) then
          if(nrestart.ne.2) then
            open(56, file = trim(constructRestartFileName(case, caseId, nSubdomains)), &
                 status='unknown',form='unformatted')
          else
            open(56, file = trim(constructRestartFileName(case, caseId_restart, nSubdomains)), &
                 status='unknown',form='unformatted')
          end if
          read (56)
          do ii=0,irank-1 ! skip records
             read(56)
          end do
          read(56) nradsteps, qrad, ozone, radlwup, radlwdn, radswup, radswdn, &
               radqrlw, radqrsw, radqrclw, radqrcsw, &
               lwDownSurface, lwDownSurfaceClearSky, lwUpSurface, lwUpSurfaceClearSky, & 
               lwUpToa,     lwUpToaClearSky,     &
               swDownSurface, swDownSurfaceClearSky, swUpSurface, swUpSurfaceClearSky, & 
               swDownToa,                            swUpToa,     swUpToaClearSky,     &
               insolation_TOA, &
          o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4, &
          visDownSurface, visDownSurfaceDiffuse, nirDownSurface, nirDownSurfaceDiffuse, &
          CosZenithAngle
          close(56)
        end if
      end do
      
    end if ! restart_sep
    
    if(rank == nsubdomains-1) then
         print *,'Case:',caseid
         print *,'Restart radiation at step:',nstep
         print *,'Time:',nstep*dt
    endif
    
    call task_barrier()
  end subroutine read_rad      
  
  ! ----------------------------------------------------------------------------
  function constructRestartFileName(case, caseid, index) result(name) 
    character(len = *), intent(in) :: case, caseid
    integer,            intent(in) :: index
    character(len=256) :: name
    
    character(len=4) :: indexChar

    integer, external :: lenstr

    write(indexChar,'(i4)') index

    name = './RESTART/' // trim(case) //'_'// trim(caseid) //'_'// &
              indexChar(5-lenstr(indexChar):4) //'_restart_rad.bin'
!bloss              trim(indexChar) //'_restart_rad.bin'

  end function constructRestartFileName
  ! ----------------------------------------------------------------------------

!==============================================================
! Formila fit to co2 observed record since 1959.
! Based on period from 1959 to 2025.
! MK 2025

real function co2_ppm(year)
implicit none
integer, intent(in):: year
real tmp
tmp = year
co2_ppm = 1.496911e-4 * tmp**3 - 8.766766e-1 * tmp**2 + 1.712398e3 * tmp - 1.115248e6
end function co2_ppm


 
end module rad
