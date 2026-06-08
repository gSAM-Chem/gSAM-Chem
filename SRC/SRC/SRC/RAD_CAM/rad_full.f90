subroutine rad_full()

  ! Interface to the longwave and shortwave radiation code from the
  ! NCAR Community Atmosphere Model (CAM3.0).
  !
  ! Originally written as interface to CCM3 radiation code by Marat
  !     Khairoutdinov
  ! Adapted to CAM3.0 radiation code by Peter Blossey, August 2004.
  !
  use rad
  use ppgrid
  use vars
  use params
  use shr_orb_mod, only: shr_orb_params
  use radae,        only: radaeini, initialize_radbuffer
  use pkg_cldoptics, only: cldefr, cldems
  use aer_optics, only: aer_optics_initialize
  use microphysics, only: reffc, reffi, &
         SnowMassMixingRatio, reffs, dosnow_radiatively_active, doallice_radiatively_active
  use terrain, only: k_terra, terra, terraw
  use slm_vars, only: albedovis_v,albedonir_v,albedovis_s,albedonir_s, &
                      phi_1,phi_2,snow_mass,LAI,IMPERV,icemask,soilw
  use consts, only: emis_water
  use buildings, only: doshadows, shadow_mask
  use stat_coars

  implicit none

  ! Local space:

  real(r4) pmid(pcols,pver)	! Level pressure (Pa)
  real(r4) pint(pcols,pverp)	! Model interface pressure (Pa)
  real(r4) massl(pcols,pver)	! Level mass (g/m2)
  real(r4) pmidrd(pcols,pver)	! Level pressure (dynes/cm2)
  real(r4) pintrd(pcols,pverp)	! Model interface pressure (dynes/cm2)
  real(r4) pmln(pcols,pver)	! Natural Log of pmid
  real(r4) piln(pcols,pverp)	! Natural Log of pint
  real(r4) tlayer(pcols,pver)	! Temperature
  real(r4) qlayer(pcols,pver)	! Specific humidity
  real(r4) cldr(pcols,pverp)	! Fractional cloud cover
  real(r4) cliqwp(pcols,pver)	! Cloud liquid water path
  real(r4) cicewp(pcols,pver)	! Cloud ice water path

  real(r4) fice(pcols,pver)	! Fractional ice content within cloud
  real(r4) rel(pcols,pver)	! Liquid effective drop radius (micron)
  real(r4) rei(pcols,pver)	! Ice effective drop size
  real(r4) o3vmr(pcols,pver)	! Ozone volume mixing ratio
  real(r4) o3mmr(pcols,pver)	! Ozone mass mixing ratio

  !bloss(2016-02-09): add variables to handle radiatively active snow
  real(r4) SnowWaterPath(pcols,pver)	! Snow water path
  real(r4) re_snow(pcols,pver)	! Snow effective particle size
  real(r4) CloudTauLW(pcols,pver)	! cloud optical depth in longwave

  real(r4) qlayer1(pcols,pver)	! Specific humidity
  real(r4) cld1(pcols,pverp)	! Fractional cloud cover
  real(r4) cliqwp1(pcols,pver)	! Cloud liquid water path
  real(r4) cicewp1(pcols,pver)	! Cloud ice water path
  real(r4) fice1(pcols,pver)	! Fractional ice content within cloud
  real(r4) rel1(pcols,pver)	! Liquid effective drop radius (micron)
  real(r4) rei1(pcols,pver)	! Ice effective drop size
  real(r4) SnowWaterPath1(pcols,pver)	! Snow water path
  real(r4) re_snow1(pcols,pver)	! Snow effective particle size

  integer lchnk              ! chunk identifier
  integer ncol               ! number of atmospheric columns
  integer nmxrgn(pcols)      ! Number of maximally overlapped regions

  real(r4) emis(pcols,pver)     ! cloud emissivity (fraction)
  real(r4) landfrac(pcols)      ! Land fraction (seems deprecated)
  real(r4) icefrac(pcols)       ! Ice fraction
  real(r4) psurface(pcols)      ! Surface pressure
  real(r4) player(pcols,pver)   ! Midpoint pressures
  real(r4) landm(pcols)         ! Land fraction
  real(r4) snowm(pcols)         ! snow depth, water equivalent (meters)

  real(r4) pmxrgn(pcols,pverp)  ! Maximum values of pmid for each

  real(r4) qrl(pcols,pver)	! Longwave heating rate (K/s)
  real(r4) qrs(pcols,pver)	! Shortwave heating rate (K/s)

  real(r4) fnl(pcols,pverp)	! Net Longwave Flux at interfaces
  real(r4) fns(pcols,pverp)	! Net Shortwave Flux at interfaces
  real(r4) fcnl(pcols,pverp)	! Net Clearsky Longwave Flux at interfaces
  real(r4) fcns(pcols,pverp)	! Net Clearsky Shortwave Flux at interfaces
  real(r4) flu(pcols,pverp)	! Longwave upward flux
  real(r4) fld(pcols,pverp)	! Longwave downward flux
  real(r4) fsu(pcols,pverp)	! Shortwave upward flux
  real(r4) fsd(pcols,pverp)	! Shortwave downward flux

  !	aerosols:

  real(r4) rh(pcols,pver)		! relative humidity for aerorsol 
  real(r4) aer_mass(pcols,pver,naer_all)  ! aerosol mass mixing ratio
  integer, parameter :: nspint = 19 ! # spctrl intrvls in solar spectrum
  integer, parameter :: naer_groups= 7 ! # aerosol grp for opt diagnostcs
  real(r4) aertau(nspint,naer_groups) ! Aerosol column optical depth
  real(r4) aerssa(nspint,naer_groups) ! Aero col-avg single scattering albedo
  real(r4) aerasm(nspint,naer_groups) ! Aerosol col-avg asymmetry parameter
  real(r4) aerfwd(nspint,naer_groups) ! Aerosol col-avg forward scattering

  !       Diagnostics:

  ! Longwave radiation
  real(r4) flns(pcols)          ! Surface cooling flux
  real(r4) flnt(pcols)          ! Net outgoing flux
  real(r4) flnsc(pcols)         ! Clear sky surface cooing
  real(r4) flntc(pcols)         ! Net clear sky outgoing flux
  real(r4) flwds(pcols)         ! Down longwave flux at surface
  real(r4) flwdsc(pcols)        ! Clear-sky down longwave flux at surface

  !bloss: New in CAM3.0.
  real(r4) flut(pcols)          ! Upward flux at top of model
  real(r4) flutc(pcols)         ! Upward clear-sky flux at top of model

  ! Shortwave radiation
  real(r4) solin(pcols)        ! Incident solar flux
  real(r4) fsns(pcols)         ! Surface absorbed solar flux
  real(r4) fsnt(pcols)         ! Flux Shortwave Downwelling Top-of-Model
  real(r4) fsntoa(pcols)      ! Total column absorbed solar flux
  real(r4) fsds(pcols)         ! Flux Shortwave Downwelling Surface
  real(r4) fsdsc(pcols)        ! Clearsky Flux Shortwave Downwelling Surface

  real(r4) fsnsc(pcols)        ! Clear sky surface absorbed solar flux
  real(r4) fsntc(pcols)        ! Clear sky total column absorbed solar flx
  real(r4) fsntoac(pcols)      ! Clear sky total column absorbed solar flx
  real(r4) sols(pcols)         ! Direct solar rad incident on surface (< 0.7)
  real(r4) soll(pcols)         ! Direct solar rad incident on surface (>= 0.7)
  real(r4) solsd(pcols)        ! Diffuse solar rad incident on surface (< 0.7)
  real(r4) solld(pcols)        ! Diffuse solar rad incident on surface (>= 0.7)
  real(r4) solsud(pcols)       ! Up Diffuse solar rad incident on surface (< 0.7)
  real(r4) sollud(pcols)       ! Up Diffuse solar rad incident on surface (>= 0.7)
  real(r4) fsnirtoa(pcols)     ! Near-IR flux absorbed at toa
  real(r4) fsnrtoac(pcols)     ! Clear sky near-IR flux absorbed at toa
  real(r4) fsnrtoaq(pcols)     ! Near-IR flux absorbed at toa >= 0.7 microns

  real(r4) frc_day(pcols)      ! = 1 for daylight, =0 for night columns
  real(r4) coszrs_in(pcols)    ! cosine of solar zenith angle
  real(r4) coszrs_in1(pcols)    ! cosine of solar zenith angle

  real(r4) asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r4) aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r4) asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r4) aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms

  real(r4) lwupsfc(pcols)   ! Longwave up flux in CGS units

  real(r4) qtot
  real(r4) dayy
  integer i,j,k,m,tmp_count,nrad_call,n_ozone_call
  real(r4) coef,factor,tmp(1),www1
  real(8) qradz(nzm),buffer(nzm)
  real perpetual_factor
  real(r4) clat(pcols),clon(pcols)
  real(r4) pii
  real(r4) tmp_ggr, tmp_cp, tmp_eps, tmp_ste, tmp_pst, eccf1
  integer nx1, ny1, nz1, shad
  logical dosolardia,solardia ! do diagnostic of visible TOA everywhere for 2D output
  logical fopened
  real tmpz(ny,nzm)
  real co2_ppm
!-----------------------------------------------------------------------------------

  if(.not.docheck.and.icycle.ne.1) goto 999  ! ugly way to handle the subcycles. add rad heating.

  nrad_call = nrad_ems*nrad
  n_ozone_call = nint(86400./dt)
  pii = atan2(0.,-1.)

  ncol = 1 ! compute one column of radiation at a time.

  if(collect_coars) call coars_fld(t(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm), &
                                         fld_tend(:,:,:,2))
  dosolardia = nstep.eq.1.or.mod(nstep,nsaveM).eq.0.and.nstep.ge.nsaveMstart &
                                   .and.nstep.le.nsaveMend
  !-------------------------------------------------------
  ! Initialize some stuff
  !


  if(initrad) then

     if(masterproc) print*,'Initializing radiation (RAD_CAM) ...'

     ! Marat: remove all the ndiv's for this version of Global SAM 
     begchunk = 1
     endchunk = nx*ny

     !bloss  subroutine initialize_radbuffer
     !bloss  inputs:  none
     !bloss  ouptuts: none (allocates and initializes abs/ems arrays)
     call initialize_radbuffer()

     !bloss  subroutine shr_orb_params
     !bloss  inputs:  iyear, log_print
     !bloss  ouptuts: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
     call shr_orb_params(iyear    , eccen  , obliq , mvelp     ,     &
           &               obliqr   , lambm0 , mvelpp, .false.)

     !bloss  subroutine radaeini
     !bloss  inputs:  pstdx (=1013250 dynes/cm2), mwdry (mwair) and mwco2.
     !bloss  ouptuts: none (sets up lookup tables for abs/ems computat.)
     call radaeini( 1.013250e6_r4, mwdry, mwco2 )

     !bloss  subroutine aer_optics_initialize
     !bloss  inputs:  none
     !bloss  ouptuts: none (sets up lookup tables for aerosol properties)
     call aer_optics_initialize()

     if(.not.allocated(tabs_rad)) then
       allocate(tabs_rad(nx,ny,nzm))
       allocate(qv_rad(nx,ny,nzm))
       allocate(qc_rad(nx,ny,nzm))
       allocate(qi_rad(nx,ny,nzm))
       allocate(qs_rad(nx,ny,nzm))
       allocate(cld_rad(nx,ny,nzm))
       allocate(rel_rad(nx,ny,nzm))
       allocate(rei_rad(nx,ny,nzm))
       allocate(res_rad(nx,ny,nzm))
     end if

     if(nrestart.eq.0) then

        do k=1,nzm
           do j=1,ny
              do i=1,nx
	         tabs_rad(i,j,k)=0.
	         qv_rad(i,j,k)=0.
	         qc_rad(i,j,k)=0.
	         qi_rad(i,j,k)=0.
	         cld_rad(i,j,k)=0.
	         rel_rad(i,j,k)=25.
	         rei_rad(i,j,k)=25.
	         qrad(i,j,k)=0.
              end do
           end do
        end do
        nradsteps=0	  
        do k=1,nz
           radlwup(k) = 0.
           radlwdn(k) = 0.
           radswup(k) = 0.
           radswdn(k) = 0.
           radqrlw(k) = 0.
           radqrsw(k) = 0.
           radqrclw(k) = 0.
           radqrcsw(k) = 0.
        end do

        if(compute_reffc) rel_rad(:,:,:) = 0.
        if(compute_reffi) rei_rad(:,:,:) = 0.

        if(dosnow_radiatively_active) then
          qs_rad(:,:,:) = 0.
          res_rad(:,:,:) = 0.
        end if

     else

        call read_rad()

     endif

     if(doperpetual) then
           ! perpetual sun (no diurnal cycle)
           do j=1,ny
              do i=1,nx
                 p_factor(i,j) = perpetual_factor(day0, latitude(i,j))
              end do
           end do
     end if

  endif

  if(initrad.or.mod(nstep,n_ozone_call).eq.0) then

   ! initialize or update  mixing ratios of trace gases.

     if(reado3) then
        n2o(:) = 0.
        ch4(:) = 0.
        cfc11(:) = 0.
        cfc12(:) = 0.
        call readlevels(trim(o3file),ozone)
        call task_barrier()
     else
        call tracesini()
     end if
     if(n2ox.ne.0.) n2o(:) = n2ox
     if(ch4x.ne.0.) ch4(:) = ch4x
     if(cfc11x.ne.0.) cfc11(:) = cfc11x
     if(cfc12x.ne.0.) cfc12(:) = cfc12x
     if(masterproc) then
       print*,'  z    o3   n2o   ch4   cfc11   cfc12'
       i=nx/(1+COL1)
       j=ny/(1+YES3D)
       do k=1,nzm
             write(6,'(6g12.4)') z(k),ozone(i,j,k),n2o(nz-k),ch4(nz-k),cfc11(nz-k),cfc12(nz-k)
       end do
      end if
      call fminmax_print('o3:',ozone,1,nx,1,ny,nzm)
      if(minval(ozone).le.0.) then
        print*,'rank=',rank,'ozone concentration should be positive:',minval(ozone),minloc(ozone) 
        call task_abort()
      end if
   end if

  !bloss  subroutine radini
  !bloss  inputs:  ggr, cp, epislo (=0.622), stebol (=5.67e-8), pstd
  !bloss  outputs: none (although it initializes constants, computes
  !                       ozone path lengths).
  tmp_ggr = ggr
  tmp_cp  = cp
  tmp_eps = mwh2o/mwdry
  tmp_ste = 5.67e-8_r4
  tmp_pst = 1.013250e6_r4
  call radini(tmp_ggr, tmp_cp, tmp_eps, tmp_ste, tmp_pst)
  
  !bloss  initialize co2 mass mixing ratio
  if(docurrentco2.and.year.gt.1959) then
   co2vmr =  1.e-6*co2_ppm(year)
  else
   co2vmr = 3.670e-4_r4 * nxco2
  end if
  co2mmr = co2vmr*rmwco2 ! rmwco2: ratio of mw of co2 to that of dry air
  if ((nstep.eq.1).and.(icycle.eq.1)) then
     if (masterproc) write(*,*) 'CO2 MMR = ', co2mmr
     if (masterproc) write(*,*) 'CO2 VMR = ', co2vmr
  end if

  if(docheck) return ! return if it is a simple check run (-nemalists option)

  ! Initialize aerosol mass mixing ratio to zero.
  ! TODO: come up with scheme to input aerosol concentrations 
  ! similar to the current scheme for trace gases.
  aer_mass = 0.

  !------------------------------------------------------
  !  Accumulate thermodynamical fields over nrad steps 
  !

  if(doradavg) then

   do k=1,nzm
      do j=1,ny
         do i=1,nx
            tabs_rad(i,j,k)=tabs_rad(i,j,k)+tabs(i,j,k)
            qv_rad(i,j,k)=qv_rad(i,j,k)+qv(i,j,k)
            qc_rad(i,j,k)=qc_rad(i,j,k)+qcl(i,j,k)
         end do
      end do
   end do
 
   if(doallice_radiatively_active) then
     do k=1,nzm
       do j=1,ny
         do i=1,nx
           qi_rad(i,j,k)=qi_rad(i,j,k)+qci(i,j,k)+qpi(i,j,k)
         end do
       end do
     end do
   else
     ! accumulate cloud liquid and cloud ice mass mixing ratios
     do k=1,nzm
       do j=1,ny
         do i=1,nx
           qi_rad(i,j,k)=qi_rad(i,j,k)+qci(i,j,k)
         end do
       end do
     end do
   end if

   if(dosnow_radiatively_active) then
     do k=1,nzm
       do j=1,ny
         do i=1,nx
           qs_rad(i,j,k)=qs_rad(i,j,k)+SnowMassMixingRatio(i,j,k)
           if(cld(i,j,k).gt.0..or.SnowMassMixingRatio(i,j,k).gt.0.) cld_rad(i,j,k) = cld_rad(i,j,k)+1.
         end do
       end do
     end do
   else
     do k=1,nzm
       do j=1,ny
         do i=1,nx
           cld_rad(i,j,k) = cld_rad(i,j,k)+cld(i,j,k)
         end do
       end do
     end do
   end if
   ! Accumulate effective radius by weighting it with mass
   if(compute_reffc) then
     rel_rad(1:nx,1:ny,1:nzm) = rel_rad(1:nx,1:ny,1:nzm) + reffc(1:nx,1:ny,1:nzm) * qcl(1:nx,1:ny,1:nzm) 
   end if
   if(compute_reffi) then
     rei_rad(1:nx,1:ny,1:nzm) = rei_rad(1:nx,1:ny,1:nzm) + reffi(1:nx,1:ny,1:nzm) * qci(1:nx,1:ny,1:nzm) 
   end if
   if(dosnow_radiatively_active) then
     res_rad(1:nx,1:ny,1:nzm) = res_rad(1:nx,1:ny,1:nzm) + reffs(1:nx,1:ny,1:nzm) * SnowMassMixingRatio(1:nx,1:ny,1:nzm) 
   end if
  
  end if ! doradavg

  nradsteps=nradsteps+1
  !----------------------------------------------------
  ! Update radiation variables if the time is due
  !

  !kzm Oct.14, 03 changed .eq.nrad to .ge.nrad to handle the
  ! case when a smaller nrad is used in restart  

  if(initrad.or.nradsteps.ge.nrad) then 
     ! Compute radiation fields for averaged thermodynamic fields

    !-----------------------------------------------------------
    ! Check if it is time to compute gas absortion coefficients for
    ! longwave radiation. 

     if(initrad.or.mod(nstep,nrad_call).eq.0) then
        doabsems = .true. 
        initrad = .false.
        if(masterproc) print*,'gaseous abs/emis coeffs recomputed. nrad_call = ',nrad_call
     else
        doabsems = .false. 
     end if

     do k=1,nz
        radlwup(k) = 0.
        radlwdn(k) = 0.
        radswup(k) = 0.
        radswdn(k) = 0.
        radqrlw(k) = 0.
        radqrsw(k) = 0.
        radqrclw(k) = 0.
        radqrcsw(k) = 0.
     end do

     coef=1./float(nradsteps)

     if(doradavg) then
 
       if(compute_reffc) then
         do k=1,nzm
            do j=1,ny
               do i=1,nx
                   rel_rad(i,j,k) = max(2.5,min(60.,rel_rad(i,j,k)/(1.e-8+qc_rad(i,j,k))))
               end do
            end do
         end do
       end if
       if(compute_reffi) then
         do k=1,nzm
            do j=1,ny
               do i=1,nx
                   rei_rad(i,j,k) = max(5.,min(250.,rei_rad(i,j,k)/(1.e-8+qi_rad(i,j,k))))
               end do
            end do
         end do
       end if
       do k=1,nzm
          do j=1,ny
             do i=1,nx
               tabs_rad(i,j,k)=tabs_rad(i,j,k)*coef
               qv_rad(i,j,k)=qv_rad(i,j,k)*coef
               qc_rad(i,j,k)=qc_rad(i,j,k)*coef
               qi_rad(i,j,k)=qi_rad(i,j,k)*coef
               cld_rad(i,j,k)=cld_rad(i,j,k)*coef
             end do
          end do
       end do

       if(dosnow_radiatively_active) then
         !bloss(2016-02-09): Compute average snow effective radius and mass mixing ratio.
         !  Note that the average effective radius is mass-weighted.
         do k=1,nzm
           do j=1,ny
             do i=1,nx
               res_rad(i,j,k) = max(5.,min(250.,res_rad(i,j,k)/(1.e-8+qs_rad(i,j,k))))
             end do
           end do
         end do
         qs_rad(:,:,:)=qs_rad(:,:,:)*coef
       end if

     else

! use instantaneous fields (except water vapor) for radiation calculation:

       do k=1,nzm
          do j=1,ny
            do i=1,nx
              tabs_rad(i,j,k)=tabs(i,j,k)
              qv_rad(i,j,k)=qv(i,j,k)
              qc_rad(i,j,k)=qcl(i,j,k)
           end do
          end do
       end do

       if(doallice_radiatively_active) then
         do k=1,nzm
           do j=1,ny
             do i=1,nx
               qi_rad(i,j,k)=qci(i,j,k)+qpi(i,j,k)
             end do
           end do
         end do
       else
         do k=1,nzm
           do j=1,ny
             do i=1,nx
               qi_rad(i,j,k)=qci(i,j,k)
             end do
           end do
         end do
       end if
       if(dosnow_radiatively_active) then
         do k=1,nzm
           do j=1,ny
             do i=1,nx
               qs_rad(i,j,k)=SnowMassMixingRatio(i,j,k)
               if(cld(i,j,k).gt.0.or.SnowMassMixingRatio(i,j,k).gt.0.) cld_rad(i,j,k) = 1.
             end do
           end do
         end do
       else
         do k=1,nzm
           do j=1,ny
             do i=1,nx
               cld_rad(i,j,k) = cld(i,j,k)
             end do
           end do
         end do
       end if
       ! Accumulate effective radius by weighting it with mass
       if(compute_reffc) then
         rel_rad(1:nx,1:ny,1:nzm) = max(2.5,min(60.,reffc(1:nx,1:ny,1:nzm)))
       end if
       if(compute_reffi) then
         rei_rad(1:nx,1:ny,1:nzm) = max(5.,min(250.,reffi(1:nx,1:ny,1:nzm)))
       end if
       if(dosnow_radiatively_active) then
         res_rad(1:nx,1:ny,1:nzm) = max(5.,min(250.,reffs(1:nx,1:ny,1:nzm)))
       end if

     end if  ! doradavg

!     if(dodebug) then
!        call fminmax_print('tabs_rad:',tabs_rad(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!        call fminmax_print('qv_rad:',qv_rad(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!        call fminmax_print('qi_rad:',qi_rad(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!        call fminmax_print('qc_rad:',qc_rad(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!        call fminmax_print('cld_rad:',cld_rad(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!        call fminmax_print('ozone:',ozone(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!     end if

     do j=1,ny
! parallelized i loop so that it could work for 2D case -Marat
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,m,lchnk,pint,o3vmr,piln,tlayer,qlayer,qtot,pmid,pmln,cldr,cliqwp,fice,massl,rh,qrl,qrs, &
!$OMP              eccf1,flu,fld,fsu,fsd,o3mmr,pmidrd,pintrd,lwupsfc,emis,nmxrgn,pmxrgn, &
!$OMP              rel,rei,re_snow,SnowWaterPath,cloudTauLW,flns,flnt,flnsc,flntc,flwds,flwdsc, &
!$OMP              flut,flutc,fnl,fcnl,dayy,coszrs_in,clat,clon,coszrs,asdir,aldir,asdif,aldif, &
!$OMP              cicewp,frc_day,fsnt,fsntc,fsntoa,solin,player,psurface,landm,icefrac,snowm,landfrac, &
!$OMP              fsntoac,fsnirtoa,fsnrtoac,fsnrtoaq,fsns,fsnsc,fsdsc,fsds,sols,soll, &
!$OMP              solsd,solld,aertau,aerssa,aerasm,aerfwd,fns,fcns)
        do i=1,nx
           lchnk = i+nx*(j-1)
           qrl = 0.
           qrs = 0.
           k=nz
           pint(1,nz)= presi(k_terra(i,j))*100.
           piln(1,nz) = log(pint(1,nz))
           do m=k_terra(i,j),nzm
              k=k-1
              o3vmr(:,k)=0.6034*ozone(i,j,m)
              pint(1,k)= presi(m+1)*100.
              piln(1,k) = log(pint(1,k))
              tlayer(1,k)=tabs_rad(i,j,m)
              qtot = (qc_rad(i,j,m)+qi_rad(i,j,m))
              qlayer(1,k)=max(1.e-7_r4,qv_rad(i,j,m))
              if(qtot.gt.0.) then
                 cldr(1,k) = min(0.99_r4,cld_rad(i,j,m))
                 cliqwp(1,k) = qtot/cldr(1,k)  ! make it "in-cloud"
                 fice(1,k) = qi_rad(i,j,m)/qtot
              else
                 if(dosnow_radiatively_active) then
                   cldr(1,k) = min(0.99_r4,cld_rad(i,j,m))
                 else
                   cldr(1,k) = 0.
                 end if
                 cliqwp(1,k) = 0.
                 fice(1,k) = 0.
              endif
              if(dosolardia) then
                qtot = (qcl(i,j,m)+qci(i,j,m))
                qlayer1(1,k)=max(1.e-7_r4,qv(i,j,m))
                if(qtot.gt.0.) then
                   cliqwp1(1,k) = qtot
                   fice1(1,k) = qci(i,j,m)/qtot
                   cld1(1,k) = min(0.99_r4,cld(i,j,m))
                else
                   if(dosnow_radiatively_active) then
                     if(SnowMassMixingRatio(i,j,m).gt.0.) then
                       cld1(1,k) = min(0.99_r4,cld(i,j,m))
                     end if
                   else
                     cld1(1,k) = 0.
                   end if
                   cliqwp1(1,k) = 0.
                   fice1(1,k) = 0.
                endif
              else
                 cld1(1,k) = 0.
                 cliqwp1(1,k) = 0.
                 fice1(1,k) = 0.
              end if
           end do
           m = k-1
           ! reconstruct profiles above domain top assuming isothermal atmospehere
           ! and space steps of 100 m.
           do k=m,1,-1
            pint(1,k) = pint(1,k+1)*exp(-ggr*100./(rgas*tabs_rad(i,j,nzm)))
            piln(1,k) = log(pint(1,k))
            o3vmr(1,k)=0.6034*ozone(i,j,nzm)
            tlayer(1,k)=tabs_rad(i,j,nzm)
            qlayer(1,k)=max(1.e-7_r4,qv_rad(i,j,nzm))
            cldr(1,k) = 0.
            cliqwp(1,k) = 0.
            fice(1,k) = 0.
            qlayer1(1,k)=max(1.e-7_r4,qv(i,j,nzm))
            cld1(1,k) = 0.
            cliqwp1(1,k) = 0.
            fice1(1,k) = 0.
            qrl(:,k)=0.
            qrs(:,k)=0.
           end do
           do k=1,nzm
             pmid(1,k) = 0.5*(pint(1,k)+pint(1,k+1))
             pmln(1,k) = log(pmid(1,k))
             massl(:,k)=1000.*(pint(:,k+1)-pint(:,k))/ggr
             cliqwp(1,k) = cliqwp(1,k)*massl(1,k)
             cliqwp1(1,k) = cliqwp1(1,k)*massl(1,k)
             rh(1,k)=0.
           end do

        ! Default reff computation:
           !bloss  subroutine cldefr
           !bloss  inputs:  lchnk, ncol, landfrac, icefrac, pres0, 
           !                pmid, landm, icrfrac, snowh
           !bloss  outputs: rel, rei (liq/ice effective radii)
           player = pmid
           psurface = 100.*pres0
           if(landmask(i,j).eq.0) then
            icefrac = 0.
            snowm = 0.
            landm = 0.
            landfrac = 0.
           else
              icefrac = 0.
              snowm = snow_mass(i,j)
              landfrac = 1.
              landm = 1.
           end if
           call cldefr(lchnk,ncol,landfrac,tlayer,rel,rei, & ! CAM3 interface
              psurface,player,landm,icefrac, snowm)
            if(compute_reffc) then
              k=nz
              do m=k_terra(i,j),nzm
                k=k-1
                rel(1,k) = rel_rad(i,j,m)
              end do
              rel(1,1:k-1) = rel(1,k)
            else
              k=nz
              do m=k_terra(i,j),nzm
                k=k-1
                rel_rad(i,j,m) = rel(1,k)
              end do
            end if
            if(compute_reffi) then
              k=nz
              do m=k_terra(i,j),nzm
                k=k-1
                rei(1,k) = rei_rad(i,j,m)
              end do
              rei(1,1:k-1) = rei(1,k)
            else
              k=nz
              do m=k_terra(i,j),nzm
                k=k-1
                rei_rad(i,j,m) = rei(1,k)
              end do
            end if
            if(dosolardia) then
              if(compute_reffc) then
                k=nz
                do m=k_terra(i,j),nzm
                  k=k-1
                  rel1(1,k) = max(2.5,min(60.,reffc(i,j,m)))
                end do
                rel1(1,1:k-1) = rel1(1,k)
              else
                rel1(1,:) = rel(1,:)
              end if
              if(compute_reffi) then
                k=nz
                do m=k_terra(i,j),nzm
                  k=k-1
                  rei1(1,k) = max(5.,min(250.,reffi(i,j,m)))
                end do
                rei1(1,1:k-1) = rei1(1,k)
              else
                rei1(1,:) = rei(1,:)
              end if
            end if

           !bloss  subroutine cldems
           !bloss  inputs:  lchnk, ncol, cliqwp, rei, fice
           !bloss  outputs: emis (cloud emissivity)
           call cldems(lchnk,ncol,cliqwp,fice,rei,emis)

           if(dosnow_radiatively_active) then
             !bloss(2016-02-08): Make snow radiatively active in longwave

             ! get the combined longwave optical depth for cloud liquid and cloud ice
             cloudTauLW(:,:) = - log( 1. - emis(:,:) )

             ! get effective radius and layer water path for snow.
             k=nz
             do m=k_terra(i,j),nzm
                k=k-1
                re_snow(1,k) = res_rad(i,j,m)
                SnowWaterPath(1,k) = qs_rad(i,j,m)*massl(1,k)
             end do
             re_snow(1,1:k-1) = re_snow(1,k)
             SnowWaterPath(1,1:k-1) = 0.
             if(dosolardia) then
               k=nz
               do m=k_terra(i,j),nzm
                 k=k-1
                 re_snow1(1,k) = max(5.,min(250.,reffs(i,j,m)))
                 SnowWaterPath1(1,k) = SnowMassMixingRatio(i,j,m)*massl(1,k)
               end do
               re_snow1(1,1:k-1) = re_snow1(1,k)
               SnowWaterPath1(1,1:k-1) = 0.
             end if

             ! taken from cldems
             !   ice absorption coefficient is kabsi = 0.005 + 1. / rei
             !   LW optical depth is 1.66*kabsi*iwp
             cloudTauLW(:,:) = cloudTauLW(:,:) &
                  + 1.66 * ( 0.005 + 1. / re_snow(:,:) ) * SnowWaterPath(:,:) 

             ! re-compute emissivity
             emis(:,:) = 1. - exp( - cloudTauLW(:,:) )

             ! Note that cloud fraction already accounts for presence of snow
           else
             re_snow(:,:) = 0.
             SnowWaterPath(:,:) = 0.
             cloudTauLW(:,:) = 0.
             re_snow1(:,:) = 0.
             SnowWaterPath1(:,:) = 0.
           end if

           cicewp = fice*cliqwp
           cliqwp = cliqwp - cicewp
           solardia = dosolardia
           if(solardia) then
            cicewp1 = fice1*cliqwp1
            cliqwp1 = cliqwp1 - cicewp1
           end if


           !bloss  subroutine radinp
           !bloss  inputs:  lchnk, ncol, pmid, pint, o3vmr
           !bloss  outputs: pmidrd, pintrd, eccf, o3mmr
           call radinp(lchnk,ncol,pmid,pint,o3vmr,pmidrd, &
                pintrd,eccf,o3mmr)
           eccf1 = eccf ! needed for opemMP

           if(landmask(i,j).eq.0) then
             lwupsfc(1) = emis_water*stebol*(sstxy(i,j)+t00)**4 ! CGS units
           else
             lwupsfc(1) = stebol*(sstxy(i,j)+t00)**4 ! CGS units
           end if

           if(dolongwave) then

              ! bloss: Set number of maximally overlapped regions 
              !        and maximum pressure so that only a single 
              !        region will be computed.
              nmxrgn = 1
              pmxrgn = 1.2e6

              !bloss  subroutine radclwmx
              !bloss  inputs:  lchnk, ncol, nmxrgn, pmxrgn, lwupsfc,
              !                tlayer, qlayer, o3vmr, pmid, pint, pmln, piln,
              !                n2o, ch4, cfc11, cfc12, cldr, emis, aer_mass
              !bloss  outputs: qrl,   ! Longwave heating rate
              !                flns,  ! Surface cooling flux
              !                flnt,  ! Net outgoing flux
              !                flut,  ! Upward flux at top of model
              !                flnsc, ! Clear sky surface cooing
              !                flntc, ! Net clear sky outgoing flux
              !                flutc, ! Upward clear-sky flux at top of model
              !                flwds, ! Down longwave flux at surface
              !                flwdsc,! Clear-sky down longwave flux at surface
              !                fcnl,  ! clear sky net flux at interfaces
              !                fnl    ! net flux at interfaces 
              call radclwmx(lchnk   ,ncol    ,                   &
                   lwupsfc ,tlayer  ,qlayer  ,o3vmr   , &
                   pmidrd  ,pintrd  ,pmln    ,piln    ,          &
                   n2o     ,ch4     ,cfc11   ,cfc12   , &
                   cldr     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                   flns    ,flnt    ,flnsc   ,flntc   ,flwds   , flwdsc, &
                   flut    ,flutc   , &
                   aer_mass,fnl     ,fcnl    ,flu    ,fld)

           !   if(rank.eq.135.and.i.eq.16.and.j.eq.7) then
           !     print*,'lwupsfc=',lwupsfc,'tlayer=',tlayer,'qlayer=',qlayer, &
           !     'o3vmr=',o3vmr,'pmidrd=',pmidrd,'pintrd=',pintrd,'pmln=',pmln,'piln=',piln, &
           !     'cldr=',cldr,'emis=',emis,'pmxrgn=',pmxrgn,'nmxrgn=',nmxrgn,'qrl=',qrl
           !   end if
              ! convert radiative heating from units of J/kg/s to K/s
              qrl = qrl/cp
              !
              ! change toa/surface fluxes from cgs to mks units
              !
              flnt     = 1.e-3*flnt
              flntc    = 1.e-3*flntc
              flns     = 1.e-3*flns
              flnsc    = 1.e-3*flnsc
              flwds    = 1.e-3*flwds
              flwdsc    = 1.e-3*flwdsc
              flut     = 1.e-3*flut
              flutc    = 1.e-3*flutc
           endif


              if (doseasons) then
                 ! The diurnal cycle of insolation will vary
                 ! according to time of year of the current day.
                 dayy = day + 0.5*nrad*dt/86400.
              else
                 ! The diurnal cycle of insolation from the calendar
                 ! day on which the simulation starts (day0) will be
                 ! repeated throughout the simulation.
                 if(doequinox) then
                   dayy = 80. + day - int(day) + 0.5*nrad*dt/86400.
                 else
                   dayy = int(day0) + day - int(day) + 0.5*nrad*dt/86400.
                 end if
              end if
              if(doperpetual) then
                 if (dosolarconstant) then
                    ! fix solar constant and zenith angle as specified
                    ! in prm file.
                    coszrs_in(1) = cos(zenith_angle*pii/180.)
                    eccf1 = solar_constant/(1367.)
                 else
                    ! perpetual sun (no diurnal cycle) - Modeled after Tompkins
                    coszrs_in(1) = cos(zenith_angle*pie/180.)
                    eccf1 = p_factor(i,j)/coszrs_in(1) ! Adjst solar constant
                 end if
              else
                 !bloss  subroutine zenith
                 !bloss  inputs:  dayy, latitude, longitude, ncol
                 !bloss  outputs: coszrs  ! Cosine solar zenith angle
                 clat(1) = pie*latitude(i,j)/180.
                 clon(1) = pie*longitude(i,j)/180.
                 coszrs_in(1) = coszrs
                 call zenith(dayy,clat,clon,coszrs_in,ncol)
              end if

         333  continue
              if(doshortwave.or.solardia) then

              if(solardia) then
                 coszrs_in1(1) = 1. ! for ALBZ statistics
                 shad = 0
              else
                 coszrs_in1(1) = coszrs_in(1) 
                 if(doshadows) then
                   if(i+j.eq.2) call buildings_shadows()
                   shad = shadow_mask(i,j)
                 else 
                   shad = 0.
                 end if
              end if

	      coszrs = coszrs_in1(1) ! needed for the isccp simulator


              !bloss  subroutine albedo
              !bloss  inputs: OCEAN (land/ocean flag), coszrs_in
              !bloss  outputs: 
              !     asdir  ! Srf alb for direct rad   0.2-0.7 micro-ms
              !     aldir  ! Srf alb for direct rad   0.7-5.0 micro-ms
              !     asdif  ! Srf alb for diffuse rad  0.2-0.7 micro-ms
              !     aldif  ! Srf alb for diffuse rad  0.7-5.0 micro-ms
              if(SLM) then
               call albedo(1,1,real(landmask(i,j)),real(icemask(i,j)),coszrs_in1(1),real(snow_mass(i,j)), &
                 real(albedovis_v(i,j)),real(albedonir_v(i,j)), &
                 real(albedovis_s(i,j)),real(albedonir_s(i,j)), real(soilw(1,i,j)), &
                 real(phi_1(i,j)),real(phi_2(i,j)),real(LAI(i,j)),real(IMPERV(i,j)), shad, &
                 asdir(1),aldir(1),asdif(1),aldif(1))
              else
               call albedo_ocn(1,1,real(landmask(i,j)),coszrs_in1(1),asdir(1),aldir(1),asdif(1),aldif(1))
              end if

              ! bloss: Set number of maximally overlapped regions 
              !        and maximum pressure so that only a single 
              !        region will be computed.
              nmxrgn = 1
              pmxrgn = 1.2e6

              !bloss: set up day fraction.
              frc_day(1) = 0.
              if (coszrs_in1(1).gt.0.) frc_day(1) = 1.

              !bloss  subroutine radcswmx
              !bloss  inputs:  
              !     lchnk             ! chunk identifier
              !     ncol              ! number of atmospheric columns
              !     pmid     ! Level pressure
              !     pint     ! Interface pressure
              !     qlayer   ! Specific humidity (h2o mass mix ratio)
              !     o3mmr    ! Ozone mass mixing ratio
              !     aer_mass   ! Aerosol mass mixing ratio
              !     rh       ! Relative humidity (fraction)
              !     cld      ! Fractional cloud cover
              !     cicewp   ! in-cloud cloud ice water path
              !     cliqwp   ! in-cloud cloud liquid water path
              !     csnowp   ! in-cloud snow water path -- bloss(2016-02-09)
              !     rel      ! Liquid effective drop size (microns)
              !     rei      ! Ice effective drop size (microns)
              !     res      ! snow effective particle size (microns) -- bloss (2016-02-09)
              !     eccf     ! Eccentricity factor (1./earth-sun dist^2)
              !     coszrs_in! Cosine solar zenith angle
              !     asdir    ! 0.2-0.7 micro-meter srfc alb: direct rad
              !     aldir    ! 0.7-5.0 micro-meter srfc alb: direct rad
              !     asdif    ! 0.2-0.7 micro-meter srfc alb: diffuse rad
              !     aldif    ! 0.7-5.0 micro-meter srfc alb: diffuse rad
              !     scon     ! solar constant
              !bloss  in/outputs: 
              !     pmxrgn   ! Maximum values of pressure for each
              !              !    maximally overlapped region. 
              !     nmxrgn   ! Number of maximally overlapped regions
              !bloss  outputs: 
              !     solin     ! Incident solar flux
              !     qrs       ! Solar heating rate
              !     fsns      ! Surface absorbed solar flux
              !     fsnt      ! Total column absorbed solar flux
              !     fsntoa    ! Net solar flux at TOA
              !     fsds      ! Flux shortwave downwelling surface
              !     fsnsc     ! Clear sky surface absorbed solar flux
              !     fsdsc     ! Clear sky surface downwelling solar flux
              !     fsntc     ! Clear sky total column absorbed solar flx
              !     fsntoac   ! Clear sky net solar flx at TOA
              !     sols      ! Direct solar rad on surface (< 0.7)
              !     soll      ! Direct solar rad on surface (>= 0.7)
              !     solsd     ! Diffuse solar rad on surface (< 0.7)
              !     solld     ! Diffuse solar rad on surface (>= 0.7)
              !     fsnirtoa  ! Near-IR flux absorbed at toa
              !     fsnrtoac  ! Clear sky near-IR flux absorbed at toa
              !     fsnrtoaq  ! Net near-IR flux at toa >= 0.7 microns
              !     frc_day   ! = 1 for daylight, =0 for night columns
              !     aertau    ! Aerosol column optical depth
              !     aerssa    ! Aerosol column avg. single scattering albedo
              !     aerasm    ! Aerosol column averaged asymmetry parameter
              !     aerfwd    ! Aerosol column averaged forward scattering
              !     fns       ! net flux at interfaces
              !     fcns      ! net clear-sky flux at interfaces
              !     fsu       ! upward shortwave flux at interfaces
              !     fsd       ! downward shortwave flux at interfaces
              if(solardia) then
                call radcswmx(i,j,rank,lchnk   ,ncol    ,                   &
                   pintrd  ,pmid    ,qlayer1  ,rh      ,o3mmr   , &
                   aer_mass  ,cld1     ,cicewp1  ,cliqwp1  ,SnowWaterPath1  ,rel1     , &
                   rei1     ,re_snow1     ,eccf1    ,coszrs_in1,scon    ,solin   , &
                   asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                   pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                   fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                   fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                   solsd   ,solld   ,frc_day ,                   &
                   aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                   fcns    ,fsu     ,fsd     ,solsud   ,sollud)
              else
                call radcswmx(i,j,rank,lchnk   ,ncol    ,                   &
                   pintrd  ,pmid    ,qlayer  ,rh      ,o3mmr   , &
                   aer_mass  ,cldr     ,cicewp  ,cliqwp  ,SnowWaterPath  ,rel     , &
                   rei     ,re_snow     ,eccf1    ,coszrs_in1,scon    ,solin   , &
                   asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                   pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                   fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                   fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                   solsd   ,solld   ,frc_day ,                   &
                   aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                   fcns    ,fsu     ,fsd     ,solsud   ,sollud)
                k=nz
                do m=k_terra(i,j),nzm
                  k=k-1
                  rel_rad(i,j,m) = rel(1,k)
                  rei_rad(i,j,m) = rei(1,k)
                end do

              end if
              ! convert radiative heating from units of J/kg/s to K/s
              qrs = qrs/cp
              !
              ! change toa/surface fluxes from cgs to mks units
              !
              fsnt     = 1.e-3*fsnt
              fsntc    = 1.e-3*fsntc
              fsntoa   = 1.e-3*fsntoa
              fsntoac  = 1.e-3*fsntoac
              fsnirtoa = 1.e-3*fsnirtoa
              fsnrtoac = 1.e-3*fsnrtoac
              fsnrtoaq = 1.e-3*fsnrtoaq
              fsns     = 1.e-3*fsns
              fsnsc    = 1.e-3*fsnsc
              fsds     = 1.e-3*fsds
              fsdsc    = 1.e-3*fsdsc
              solin    = 1.e-3*solin
           endif
           if(solardia) then
             alb_xy(i,j) = 1.-fsntoa(1)/solin(1)
             albc_xy(i,j) = 1.-fsntoac(1)/solin(1)
             solardia = .false.
          !   if(alb_xy(i,j).lt.0.) then
          !    print*,'>>>snowpath:',SnowWaterPath1  ,'rel1:',rel1     , &
          !         'rei1:',rei1     ,'res:',re_snow1,'fsntoa:',fsntoa(1),'solin:',solin(1),'cld1:',cld1, &
          !         'cicewp1:',cicewp1  ,'cliqwp1:',cliqwp1   
          !   end if
             goto 333
           end if
           
           !
           ! Satellite simulator diagnostics using time-averaged values
           !
           if(doisccp .or. domodis .or. domisr) then 
             tau_067_cldliq (i,j,nzm:1:-1) = compute_tau_l(cliqwp(1,1:nzm), rel(1,1:nzm))
             tau_067_cldice (i,j,nzm:1:-1) = compute_tau_i(cicewp(1,1:nzm), rei(1,1:nzm)) 
             if(dosnow_radiatively_active) then
               tau_067_snow (i,j,nzm:1:-1) = compute_tau_i(SnowWaterPath(1,1:nzm), re_snow(1,1:nzm)) 
             else
               tau_067_snow (i,j,:) = 0. 
             end if
             tau_067 (i,j,:) = tau_067_cldliq (i,j,:) + tau_067_cldice (i,j,:) &
                  + tau_067_snow (i,j,:)
             emis_105(i,j,1:nzm) = emis(1,nzm:1:-1)
           end if 
           
           qrad(i,j,:)=0. 
           k=nz
           do m=k_terra(i,j),nzm
              k=k-1
              www1 = wgtw(j,m)*terraw(i,j,m)
              qrad(i,j,m)=qrl(1,k)+qrs(1,k)
              radlwup(m)=radlwup(m)+flu(1,k)*1.e-3*www1
              radlwdn(m)=radlwdn(m)+fld(1,k)*1.e-3*www1
              radqrlw(m)=radqrlw(m)+qrl(1,k)*www1
              radswup(m)=radswup(m)+fsu(1,k)*1.e-3*www1
              radswdn(m)=radswdn(m)+fsd(1,k)*1.e-3*www1
              radqrsw(m)=radqrsw(m)+qrs(1,k)*www1
              !bloss: clearsky heating rates
              radqrclw(m)=radqrclw(m)+(fcnl(1,k+1)-fcnl(1,k))/massl(1,k)/cp*www1
              radqrcsw(m)=radqrcsw(m)-(fcns(1,k+1)-fcns(1,k))/massl(1,k)/cp*www1
           enddo

           lwnsxy(i,j) = flns(1)
           lwntxy(i,j) = flut(1)
           lwntmxy(i,j) = flnt(1)
           lwnscxy(i,j) = flnsc(1)
           lwntcxy(i,j) = flntc(1)
           lwdsxy(i,j) = flwds(1)
           lwdscxy(i,j) = flwdsc(1)

           swnsxy(i,j) = fsns(1)
           swntxy(i,j) = fsntoa(1)
           swntmxy(i,j) = fsnt(1)
           swnscxy(i,j) = fsnsc(1)
           swntcxy(i,j) = fsntoac(1)
           swdsxy(i,j) = fsds(1)
           swdscxy(i,j) = fsdsc(1)
           swdsvisxy(i,j) = sols(1)
           swdsnirxy(i,j) = soll(1)
           swdsvisdxy(i,j) = solsd(1)
           swdsnirdxy(i,j) = solld(1)
           swusvisdxy(i,j) = solsud(1)
           swusnirdxy(i,j) = sollud(1)
           solinxy(i,j) = solin(1)
           coszrsxy(i,j) = coszrs_in1(1)
           albvisxy(i,j) = asdif(1)
           albnirxy(i,j) = aldif(1)

        end do
!$OMP END PARALLEL DO 
     end do

!    compute solar and ir fluxes on buildiings walls and on surfaces around:
     if(dobuildings) call buildings_rad_walls()

     ! MODIS simulator diagnostics
     if(domodis) then 
       rad_reffc(1:nx,1:ny,1:nzm) = rel_rad(1:nx,1:ny,1:nzm)
       rad_reffi(1:nx,1:ny,1:nzm) = rei_rad(1:nx,1:ny,1:nzm) 
     end if

     tabs_rad(:,:,:)=0.
     qv_rad(:,:,:)=0.
     qc_rad(:,:,:)=0.
     qi_rad(:,:,:)=0.
     cld_rad(:,:,:)=0.
     if(compute_reffc) then
       rel_rad(:,:,:) = 0. 
     end if
     if(compute_reffi) then
       rei_rad(:,:,:) = 0. 
     end if
     if(nstep.ne.1) nradsteps=0
     
     if(dosnow_radiatively_active) then
       !bloss(2016-02-12): Zero out accumulated snow variables.
       qs_rad(:,:,:) = 0.
       res_rad(:,:,:) = 0.
     end if
     
     if(masterproc.and.doshortwave.and..not.doperpetual) &
          print*,'radiation: coszrs=',coszrs_in1(1),' solin=',solinxy(nx,ny)
     if(masterproc.and.doshortwave.and.doperpetual) &
          print*,'radiation: perpetual sun, solin=',solinxy(1,1)
     if(masterproc.and..not.doshortwave) &
          print*,'longwave radiation is called'

   ! Homogenize radiation:

     if(doradhomozonal) then
      call mean_x_3D(qrad,tmpz)
      do k=1,nzm
         do j=1,ny
           qrad(:,j,k) = tmpz(j,k) 
         end do
      end do
     end if

     if(doradhomo) then    

        factor = 1./dble(nx*ny)
        do k=1,nzm
          qradz(k) = 0.
           do j=1,ny
             do i=1,nx
              qradz(k) = qradz(k) + qrad(i,j,k)*wgt(j,k)
             end do
           end do
           qradz(k) = qradz(k) * factor
           buffer(k) = qradz(k)
        end do

        factor = 1./float(nsubdomains)
        if(dompi) call task_sum_real8(qradz,buffer,nzm)

        do k=1,nzm
           qradz(k)=buffer(k)*factor
           do j=1,ny
             do i=1,nx
               qrad(i,j,k) = qradz(k) 
             end do
           end do
        end do

     end if


  endif ! (nradsteps.eq.nrad) 

!--------------------------------------------------------
! Prepare statistics:
  
  do j=1,ny
     do i=1,nx
        ! Net surface and toa fluxes
        lwnt_xy(i,j) = lwnt_xy(i,j) + lwntxy(i,j) 
        swnt_xy(i,j) = swnt_xy(i,j) + swntxy(i,j)
        lwntc_xy(i,j) = lwntc_xy(i,j) + lwntcxy(i,j)
        swntc_xy(i,j) = swntc_xy(i,j) + swntcxy(i,j)
        lwns_xy(i,j) = lwns_xy(i,j) + lwnsxy(i,j) 
        swns_xy(i,j) = swns_xy(i,j) + swnsxy(i,j)
        lwnsc_xy(i,j) = lwnsc_xy(i,j) + lwnscxy(i,j)
        swnsc_xy(i,j) = swnsc_xy(i,j) + swnscxy(i,j)
        lwds_xy(i,j) = lwds_xy(i,j) + lwdsxy(i,j) 
        swds_xy(i,j) = swds_xy(i,j) + swdsxy(i,j)
        lwdsc_xy(i,j) = lwdsc_xy(i,j) + lwdscxy(i,j) 
        swdsc_xy(i,j) = swdsc_xy(i,j) + swdscxy(i,j)
        lwnta_xy(i,j) = lwnta_xy(i,j) + lwntxy(i,j) 
        swnta_xy(i,j) = swnta_xy(i,j) + swntxy(i,j)
        lwntca_xy(i,j) = lwntca_xy(i,j) + lwntcxy(i,j) 
        swntca_xy(i,j) = swntca_xy(i,j) + swntcxy(i,j)
        lwnsa_xy(i,j) = lwnsa_xy(i,j) + lwnsxy(i,j) 
        swnsa_xy(i,j) = swnsa_xy(i,j) + swnsxy(i,j)
        lwnsca_xy(i,j) = lwnsca_xy(i,j) + lwnscxy(i,j) 
        swnsca_xy(i,j) = swnsca_xy(i,j) + swnscxy(i,j)
        lwdsa_xy(i,j) = lwdsa_xy(i,j) + lwdsxy(i,j) 
        swdsa_xy(i,j) = swdsa_xy(i,j) + swdsxy(i,j)
        lwdsca_xy(i,j) = lwdsca_xy(i,j) + lwdscxy(i,j) 
        swdsca_xy(i,j) = swdsca_xy(i,j) + swdscxy(i,j)
        ! TOA Insolation
        solin_xy(i,j) = solin_xy(i,j) + solinxy(i,j)
        solina_xy(i,j) = solina_xy(i,j) + solinxy(i,j)
        albvis_xy(i,j) = albvis_xy(i,j) + albvisxy(i,j)
        albnir_xy(i,j) = albnir_xy(i,j) + albnirxy(i,j)
     end do
  end do

!----------------------------------------------------------------
  if(dostatisrad) then

     do j=1,ny
        do i=1,nx
           www1=wgty(j)
           s_flns = s_flns + lwnsxy(i,j)*www1 
           s_flnsc = s_flnsc + lwnscxy(i,j)*www1 
           s_flnt = s_flnt + lwntmxy(i,j)*www1 
           s_flntoa = s_flntoa + lwntxy(i,j)*www1 
           s_flntoac = s_flntoac + lwntcxy(i,j)*www1 
           s_fsnt = s_fsnt + swntmxy(i,j)*www1 
           s_fsntoa = s_fsntoa + swntxy(i,j)*www1 
           s_fsntoac = s_fsntoac + swntcxy(i,j)*www1 
           s_fsns = s_fsns + swnsxy(i,j)*www1 
           s_fsnsc = s_fsnsc + swnscxy(i,j)*www1 
           s_fsds = s_fsds + swdsxy(i,j)*www1 
           s_fsdsc = s_fsdsc + swdscxy(i,j)*www1 
           s_flds = s_flds + lwdsxy(i,j)*www1 
           s_fldsc = s_fldsc + lwdscxy(i,j)*www1 
           s_solin = s_solin + solinxy(i,j)*www1 
        end do
     end do
  end if
!----------------------------------------------------------------
!  Write the radiation-restart file:
  
  
  if(mod(nstep,nrestart_steps*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) then
  
     if(.not.donorestart) call write_rad() ! write radiation restart file
  
  endif

!-------------------------------------------------------
! Update the temperature field:	


999 continue

  if(save2Dradac) then
    do j=1,ny
     do i=1,nx
        lwntac_xy(i,j) = lwntac_xy(i,j) + lwntxy(i,j) * dtn
        swntac_xy(i,j) = swntac_xy(i,j) + swntxy(i,j) * dtn
        lwnsac_xy(i,j) = lwnsac_xy(i,j) + lwnsxy(i,j) * dtn
        swnsac_xy(i,j) = swnsac_xy(i,j) + swnsxy(i,j) * dtn
        lwdsac_xy(i,j) = lwdsac_xy(i,j) + lwdsxy(i,j) * dtn
        swdsac_xy(i,j) = swdsac_xy(i,j) + swdsxy(i,j) * dtn
     end do
    end do
  end if

  do k=1,nzm
     do j=1,ny
        do i=1,nx
	   t(i,j,k)=t(i,j,k)+qrad(i,j,k)*dtn
        end do
     end do
  end do

if(collect_coars) call coars_tend(t(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),2)

contains
  ! -------------------------------------------------------------------------------
  elemental function compute_tau_l(lwp, re_l)
    real(r4), intent(in) :: lwp, re_l
    real(r4)             :: compute_tau_l
    !
    ! Diagnose optical thickness (nominally at 0.67 microns) from water clouds
    !
    ! This version comes from radcswmx.f90 in the CAM3 radiation package, which 
    !   expects L/IWP in g/m2 and particle sizes in microns
  
    real(r4), parameter :: abarl = 2.817e-02, bbarl = 1.305, & 
                       abari = 3.448e-03, bbari = 2.431

    compute_tau_l = 0.
    if(re_l > 0.) compute_tau_l= (abarl + bbarl/re_l) * lwp  
  
  end function compute_tau_l

  ! -------------------------------------------------------------------------------
  elemental function compute_tau_i(iwp, re_i)
    real(r4), intent(in) :: iwp, re_i
    real(r4)             :: compute_tau_i
    !
    ! Diagnose optical thickness (nominally at 0.67 microns) from ice clouds
    !
    ! This version comes from radcswmx.f90 in the CAM3 radiation package, which 
    !   expects L/IWP in g/m2 and particle sizes in microns
  
    real(r4), parameter :: abarl = 2.817e-02, bbarl = 1.305, & 
                       abari = 3.448e-03, bbari = 2.431

    compute_tau_i = 0.
    if(re_i > 0.) compute_tau_i = (abari + bbari/re_i) * iwp 
  
  end function compute_tau_i
  ! -------------------------------------------------------------------------------
end subroutine rad_full


real function perpetual_factor(day, lat)
use shr_kind_mod, only: r4 => shr_kind_r4
use grid, ONLY: dt
use ppgrid, only: pie, eccf
use params, only: doequinox, nrad
implicit none

!  estimate the factor to multiply the solar constant
!  so that the sun would produce the same
!  total input the TOA as the sun subgect to diurnal cycle.
!  depends only on latitude as averaging is done over one day
!
! Input arguments
!
real(8), intent(in) ::  day             ! Calendar day, without fraction
real(8), intent(in) ::  lat                ! Current centered latitude (degrees)
real(8) lon          

! Local:
real(r4) :: tmp
real(r4) :: ttime
real(r4) :: coszrs(1) 
real(r4) :: clat(1), clon(1)
real :: dttime,iday 

if(doequinox) then
 iday = 80.
else
 iday = int(day)
end if
ttime = iday
dttime = dt*float(nrad)/86400.
tmp = 0.
pie = acos(-1.)
lon = 0.

do while (ttime.lt.iday+1.)
  clat = pie*lat/180.
  clon = pie*lon/180.
  call zenith(ttime, clat, clon, coszrs, 1)
  tmp = tmp + min(dttime,day+1-ttime)*max(0._r4,eccf*coszrs(1))
  ttime = ttime+dttime
end do

perpetual_factor = tmp

end function perpetual_factor


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
end

