	
subroutine setparm()
	
!       initialize parameters:

use vars
!use micro_params
use params
use microphysics, only: micro_setparm
use chemistry, only: chem_setparm
use sgs, only: sgs_setparm
use movies, only : irecc
use instrument_diagnostics, only: zero_instr_diag
use terrain, only: alpha_min
use cup, only: cup_init
implicit none
	
integer icondavg, ierr 
real tmp

!--------------------------------------

NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, & 
 dolongwave, doshortwave, dosgs, dz, doconstdz, &
 docoriolis, docoriolisz, dosurface, dolargescale, doradforcing, &
 fluxt0,fluxq0,tau0,tabs_s,z0,nelapse, dt, dx, dy,  &
 fcor, ug, vg, nstop, caseid, case_restart,caseid_restart, &
 nstat, nstatfrq, nprint, nrestart, doradsimple, &
 nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
 donudging_uv, donudging_tq, donudging_t, donudging_q, tauls,tautqls,&
 nudging_uv_z1, nudging_uv_z2, nudging_t_z1, nudging_t_z2, &
 nudging_q_z1, nudging_q_z2, dofplane, timelargescale, longitude0, latitude0, day0, nrad, &
 OCEAN,LAND,ISLAND,SFC_FLX_FXD,SFC_TAU_FXD, doensemble, nensemble, dowallx, dowally, &
 nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, docolumn, save2Dbin, save2Davg, save3Dbin, &
 save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
 doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
 domodis, domisr, dodynamicocean, ocean_type, delta_sst, &
 depth_slab_ocean, Szero, deltaS, timesimpleocean, &
 dosolarconstant, solar_constant, zenith_angle, rundatadir, dotracers, perturb_type, &
 doSAMconditionals, dosatupdnconditionals, doscamiopdata, iopfile, dozero_out_day0, &
 nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
 nmovie, nmoviestart, nmovieend, nrestart_skip, bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
 bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, dossthomo, &
 rad3Dout, nxco2, dosimfilesout, notracegases, &
 doradlat, doradlon, ncycle_max, doseawater, doterrain, dobufferzonex, &
 terrain_type, dobufferzoney, island_x1, island_x2, island_y1, island_y2, &
 doroundisland, island_radius, bufferzonex, bufferzoney, &
 sst_climo, tau_ocean, dosstclimo, dolatlon, dodamping_poles, dodebug, &
 doglobal, readinit, landmaskfile, sstfile, initfile, terrainfile, &
 readsst, readterr, readlandmask, earth_factor, dodamping_w, &
 reado3, o3file, n2ox, ch4x, cfc11x, cfc12x, &
 nsave2DZ, nsave2DZstart, nsave2DZend, save2DZsep, dogzip2DZ, dogzip2DL, &
 nsave2DL, nsave2DLstart, nsave2DLend, save2DLbin, save2DLavg, save2DLsep, dogzip2DL, &
 dopointsource, pointsource_start_step, year0, nens, date0, dofcast, save2Drada, &
 readlatlon, latlonfile, doslabocean, ncallocean, damping_u_cu, damping_w_cu, dohs94, &
 nub, readlat, gmg_precision, doimplicitdiff, dodamping_u, LES_S, &
 dooceanonly, nstatcoars, nstatcoarsstart, nstatcoarsend, savecoarsbin, &
 donudge3D,nudge3D_dir,nudge3D_file, nudge3D_tau, dofcast_ice, nfiles3D, nfilescoars, &
 docap_snd_cu, cap_snd_cu, tracer_source_type,  dofixdynamics,fixdynamics_type, &
 donodynamics, doflat, doradavg, dosaferestart, gamma_RAVE, & 
 nadv_mom, alpha_hybrid, doequilocean, dodatefilename, donobuoyancy, dostatcoars, &
 donorestart, doequinox, pointsource_stop_step, latlonfilebin, &
 nsaveM, nsaveMstart, nsaveMend, saveMbin, saveMsep, saveMavg, dogzipM, save2Dradac, &
 doterrhyb, doterrepair, dometric, npressure_iter, doregion, donudging_w, &
 docheckmom, docheckenergy, dotrsfcflux, dlon, dlat, dofliplon, &
 dobuildings, dossthomozonal, doseaice, seaicethickness, doradhomozonal, &
 nrad_ems, doseaiceevol, docup, n_cup, dotc, save3Dnetcdf, save2Dnetcdf, save2DZnetcdf, &
 saveMnetcdf, save2DLnetcdf, donearglobal, w3D_pressure, docyclic, cycle_period, &
 read_meters, nrestart_steps, doglobalpresets, dosimplesnd, cwp_threshold, &
 shift_sst, lhf_fudge, shf_fudge, nelapse_minutes, tau_fudge, snow2Dout, &
 nudge3Dstep_start, nudge3Dstep_end, dospectralnudging, nx_spectral, ny_spectral, &
 nstep_spectral, do_cap_wind, cap_wind, dochem

!----------------------------------
!  Read namelist variables from the standard input:
!------------

open(55,file='./CASES/'//trim(case)//'/prm', status='old',form='formatted') 
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        rewind(55)
        read (55,PARAMETERS)
        call task_abort()
end if
close(55)

! write namelist values out to file for documentation
if(masterproc) then
      open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted')
      write (55,nml=PARAMETERS)
      write(55,*) 
      close(55)
end if

!--------------------------------------
! First check if the domains sizes are legit, that is divisible by 2,3, and 5 only:

  if(mod(nx_gl,2).ne.0.or.RUN3D.and.mod(ny_gl,2).ne.0) then
   if(masterproc) print*,'nx_gl and ny_gl (if 3D) should be even numbers. Exit...'
   call task_abort()
  end if

  if(.not.Productof235(nx_gl).or.RUN3D.and.(.not.Productof235(ny_gl))) then
   if(masterproc) print*,'nx_gl and ny_gl should be divisible only by 2,3,5. Exit...'
   call task_abort()
  end if


!------------------------------------
!  Set parameters 

if(nrestart.eq.2.and.trim(caseid).eq.trim(caseid_restart)) then
  if(masterproc) print*,'branch name',trim(caseid_restart),' is the same as main runs name', trim(caseid) 
  stop
end if

if(masterproc) then
  print*,'nrestart=',nrestart
  print*,'Case: ',trim(case)
  print*,'Caseid:',trim(caseid)
  if(nrestart.eq.2) print*,'branch from Case: ',trim(case_restart)
  if(nrestart.eq.2) print*,'branch from Caseid: ',trim(caseid_restart)
  print*
end if

if(mod(nstop,nstat).ne.0.or.nelapse.ne.999999999.and.mod(nelapse,nstat).ne.0) then
    if(masterproc) print*, 'nstop or nelapse should be divisible by nstat. exit...'
    call task_abort()
end if
if(mod(nstop,nsave2D).ne.0.or.nelapse.ne.999999999.and.mod(nelapse,nsave2D).ne.0) then
    if(masterproc) print*, 'nstop or nelapse should be divisible by nsave2D. exit...'
    call task_abort()
end if
if(mod(nstop,nsave2DL).ne.0.or.nelapse.ne.999999999.and.mod(nelapse,nsave2DL).ne.0) then
    if(masterproc) print*, 'nstop or nelapse should be divisible by nsave2DL. exit...'
    call task_abort()
end if
if(mod(nstop,nsaveM).ne.0.or.nelapse.ne.999999999.and.mod(nelapse,nsaveM).ne.0) then
    if(masterproc) print*, 'nstop or nelapse should be divisible by nsaveM. exit...'
    call task_abort()
end if
! set the frequency of restart
if(nrestart_steps.eq.1) nrestart_steps = max(nstat,nsave2D,nsaveM,nsave2DL)
if(masterproc) print*,'nrestart_steps=',nrestart_steps
if(mod(nrestart_steps,nstat).ne.0.or. &
   mod(nrestart_steps,nsave2D).ne.0.or. &
   mod(nrestart_steps,nsaveM).ne.0.or. &
   mod(nrestart_steps,nsave2DL).ne.0) then
    if(masterproc) print*, 'inconsistent choice for nstat,nsave2D,nsaveM,nsave2DL'
    if(masterproc) print*, 'and nrestart_steps.'
    if(masterproc) print*,'nstat=',nstat
    if(masterproc) print*,'nsave2D=',nsave2D
    if(masterproc) print*,'nsaveM=',nsaveM
    if(masterproc) print*,'nsave2DL=',nsave2DL
    if(masterproc) print*,'nrestart_steps=',nrestart_steps
    if(masterproc) print*, 'Consult the User Guide, Exit...'
    call task_abort()
end if

if(save3Dnetcdf) nfiles3D = 1

if(day0.ne.0..and.date0.ne.0) then
 if(masterproc) print*, 'day0 and date0 cannot be both set. exitting...'
 call task_abort()
end if

if(date0.ne.0) then
   if(doequinox) then
    day0 = 80.   ! perpetual equinox run
    date0 = 0321000000
   end if
   dodatefilename = .true.
   call date_int2char(date0,datechar)
   if(masterproc) print*,'date0=',date0
   if(masterproc) print*,'datechar=',datechar
   call dayofyear_from_date(datechar,day0)
   if(masterproc) print*,'day0=',day0
   year0 = date0/10000000000_8
   if(masterproc) print*,'year0=',year0
else
   if(doequinox) then
    day0 = 80.5   ! perpetual equinox run
   end if
!   call date_from_dayofyear(day0,year0,datechar)
!   if(masterproc) print*,'datechar=',datechar
!   call dayofyear_from_date(datechar,day0)
   if(masterproc) print*,'day0=',day0
end if

! Allow only special cases for separate output:

if(save2Davg)  save2Drada = .false.

if(docolumn.and.nx_gl.ne.1) then
 if(masterproc) print*,'docolumn = T, but nx_gl.ne.1. Exit...'
 call task_abort()
end if

if(sizeof(gmg_precision).gt.4) gmg_precision = gmg_precision * 1.e-3

if(donearglobal) dolatlon = .true.

if(masterproc) print*,'read external data in meter: read_meters=',read_meters

if(doglobal) then
   if(.not.RUN3D) then
      if(masterproc)print*,'Error: when doglobal, domain should be 3D.  exit...'
      call task_abort()
   end if
   if(nrestart.ne.2.and.(dx.ne.0..or.dy.ne.0.)) then
      if(masterproc)print*,'Error: when doglobal, dx or dy should not be set. exit...'
      call task_abort()
   end if
   dolatlon = .true.
   dowallx = .false.
   dowally = .true.
   dodamping = .true.
   dodamping_poles = .true.
   latitude0 = 0.
   longitude0 = 0.
   dlon = 360._8/nx_gl
   dx = dlon*deg2rad*rad_earth/earth_factor
   if(.not.(readlatlon.or.readlat)) then
     dlat = 180._8/ny_gl
     dy = dlat*deg2rad*rad_earth/earth_factor
   end if
else if(dolatlon) then
   if(readlat.and.dlat.ne.0.) then
    if(masterproc) print*,'if readlat=T, dolatlon=T, dlat cannot be set by namelist! STOP.'
    call task_abort()
   end if
   if(readlatlon.and.(dlat.ne.0.or.dlon.ne.0.)) then
    if(masterproc) print*,'if readlatlon=T, dolatlon=T, dlat or dlon  cannot be set by namelist! STOP.'
    call task_abort()
   end if
   dowally = .true.
   if(dlat.eq.0..and..not.(readlat.or.readlatlon)) then
    if(masterproc) print*,'when dolatlon=.true. dlat needs to be set by namelist or', &
           'latitudes need to be read from file (readlat=T or readlatlon=T) STOP'
    call task_abort()
   end if
   if(.not.donearglobal.and.dlon.eq.0..and..not.readlatlon) then
    if(masterproc) print*,'when dolatlon=.true. dlon needs to be set by namelist or', &
           'longitudes need to be read from file (readlatlon=T) STOP'
    call task_abort()
   end if
  ! it needs not be equal to grid spacing in y at the equator. - MK
   if(.not.(readlatlon.or.readlat)) then
     if(donearglobal) then
        if(dlat*ny_gl.gt.180.) then
          if(masterproc) print*,'dlat is too large for near-global run: dlat*ny_gl=',dlat*ny_gl
          call task_abort()
        end if 
        dlon = 360./nx_gl
     end if
     dy = dlat*deg2rad*rad_earth/earth_factor
    ! not that dx is always the grid spacings (in m) 
    ! at the equator when dolatlon = .true.. - MK
     if(masterproc) print*,'dlat =',dlat
     if(masterproc) print*,'dy =',dy
     if(.not.readlatlon) then
      dx = dlon*deg2rad*rad_earth/earth_factor
      if(masterproc) print*,'dlon =',dlon
      if(masterproc) print*,'dx =',dx
     end if
   end if
else
   if(dx.eq.0..or.dy.eq.0.) then
    if(masterproc) print*,'dx and dy are not set. STOP'
    call task_abort()
   end if
end if

if(doglobalpresets.and..not.doglobal) then
    if(masterproc) print*,'doglobalpresents can be set only for global runs'
    call task_abort()
end if
if(doglobalpresets) then
   doimplicitdiff = .true.
   dodamping_w = .true.
   dodamping_u = .true.
   damping_u_cu = 0.25
   damping_w_cu = 0.3
   doseawater = .true.
end if

if(masterproc) then
 print*
 if(dolatlon) then
  print*,'Latitude-Longitude grid'
 else
  print*,'Cartesian grid'
 end if
 print*
end if

if(readlatlon.and.readlat) then
  if(masterproc)print*,'Error: both readlatlon and readlat are true. exit...'
    call task_abort()
end if

if(readterr) terrain_type = 0
if(dobuildings.and..not.doterrain) then
  if(masterproc) print*,'buildings=T  requires doterrain=T'
  call task_abort()
end if

if(doregion) then
 if(.not.donudge3D) then
   if(masterproc) print*,'doregion cannot be true without 3D nudging (donudge3D=.true.) Exit.'
   call task_abort()
 end if
 dowallx = .true.
 dowally = .true.
 dobufferzonex = .true.
 dobufferzoney = .true.
 if(bufferzonex.eq.0.) bufferzonex = 0.1
 if(bufferzoney.eq.0.) bufferzoney = 0.1
 if(dospectralnudging.and.(nx_spectral.le.0.or.ny_spectral.le.0)) then
   if(masterproc) print*,' spectral nudging requires setting nx_spectral and ny_spectral. exit.'
   call task_abort()
 end if
else
 if(dospectralnudging) then
   if(masterproc) print*,' spectral nudging requires doregion=.true. exit'
   call task_abort()
 end if
end if

if(do_cap_wind.and.cap_wind.eq.0.) then
   if(masterproc) print*,'cap_wind should be set when do_cap_wind=.true. exit'
   call task_abort()
end if

if(.not.doterrain) npressure_iter = 0

if(RUN2D) dy=dx

if(RUN2D.and.YES3D.eq.1) then
  if(masterproc)print*,'Error: 2D run and YES3D is set to 1. Exitting...'
  call task_abort()
endif

if(RUN3D.and.YES3D.eq.0) then
  if(masterproc)print*,'Error: 3D run and YES3D is set to 0. Exitting...'
  call task_abort()
endif

if(docoriolis.and..not.dofplane.and..not.dowally) then
   if(masterproc) print*,'incompatible parameters: dofplane=',dofplane,' and dowally=',dowally
   call task_abort()
end if
if((doradlat.or.doradlon).and..not.dowally) then
   if(masterproc) print*,'incompatible parameters: doradlat, doradlone=',doradlat,doradlon,' and dowally=',dowally
   call task_abort()
end if

dowally = dowally.and.RUN3D ! dowally is set only for 3D runs

if(doequilocean.and.doslabocean) then
 if(masterproc) print*, 'doequilocean and doslabocean cannot be .true. at the same time'
 call task_abort()
end if

if(ny.eq.1) dy=dx

if(doterrhyb.and.nadv_mom.ne.23) then
  if(masterproc) print*,'doterrhyb=T requires that adv_mom was set to 23. Exit...'
  call task_abort()
end if

if(nadv_mom.eq.3) then
   nadams = 2
   cfl_max =  0.55
   alpha_hybrid = 0.
else if(nadv_mom.eq.23) then
   nadams = 3
   if(doterrhyb) then
      tmp = min(alpha_hybrid,alpha_min)
   else
      tmp = alpha_hybrid
   end if
   cfl_max = 0.393433 + 0.235064*tmp + 0.0919493*tmp**2
end if

if(masterproc) then
   print*,'nadams=',nadams
   print*,'nadv_mom=',nadv_mom
   if(nadv_mom.eq.23) print*,'alpha_hybrid = ',alpha_hybrid
   print*,'cfl_max = ',cfl_max
end if


notopened2D = .true.
notopened2DL = .true.
notopened2DZ = .true.
notopened3D = .true.
notopenedmom = .true.

call zero_instr_diag() ! initialize instruments output 
call sgs_setparm() ! read in SGS options from prm file.
call micro_setparm() ! read in microphysical options from prm file.
call chem_setparm()  ! read in chemistry options from prm

if(dosmoke) then
  epsv = 0.
  cpvf = 0.
  docloud = .false.
  doprecip = .false.
  dolongwave = .false.
  doshortwave = .false.
else    
  epsv=0.61
  cpvf = 1.
endif   

if(navgmom_x.lt.0.or.navgmom_y.lt.0) then  
   nstatmom        = 1
   nstatmomstart    = 999999999
   nstatmomend      = 999999999
end if
if(ncoars_x.lt.0.or.ncoars_y.lt.0) then
   nstatcoars        = 1
   nstatcoarsstart    = 999999999
   nstatcoarsend      = 999999999
end if


if(doseawater) then
  salt_factor = 0.981
else
  salt_factor = 1.
end if

if(dohs94) perturb_type = 10

SLM = dosurface.and.(LAND.or.ISLAND).and..not.(SFC_FLX_FXD.or.dosfcforcing)

if(tabs_s.ne.0.) sstxy(:,:) = tabs_s - t00

dtfactor = 1. ! just initialize it

if(docup) call cup_init() ! initialize convective parameterization arrays

if(tautqls.eq.99999999.) tautqls = tauls
          
        !===============================================================
        ! UW ADDITION

!bloss: set up conditional averages
ncondavg = 1 ! always output CLD conditional average
if(doSAMconditionals) ncondavg = ncondavg + 2
if(dosatupdnconditionals) ncondavg = ncondavg + 3
if(allocated(condavg_factor)) then ! avoid double allocation when nrestart=2
  DEALLOCATE(condavg_factor,condavg_mask,condavgname,condavglongname)
end if
ALLOCATE(condavg_factor(nzm,ncondavg), & ! replaces old cloud_factor, core_factor
     condavg_mask(nx,ny,nzm,ncondavg), & ! nx x ny x nzm indicator arrays
     condavgname(ncondavg), & ! short names (e.g. CLD, COR, SATUP)
     condavglongname(ncondavg), & ! long names (e.g. cloud, core, saturated updraft)
     STAT=ierr)
if(ierr.ne.0) then
     write(*,*) '**************************************************************************'
     write(*,*) 'ERROR: Could not allocate arrays for conditional statistics in setparm.f90'
     call task_abort()
end if
        
! indicators that can be used to tell whether a particular average
!   is present.  If >0, these give the index into the condavg arrays
!   where this particular conditional average appears.
icondavg_cld = -1
icondavg_cor = -1
icondavg_cordn = -1
icondavg_satup = -1
icondavg_satdn= -1
icondavg_env = -1

icondavg = 0
icondavg = icondavg + 1
condavgname(icondavg) = 'CLD'
condavglongname(icondavg) = 'cloud'
icondavg_cld = icondavg

if(doSAMconditionals) then
   icondavg = icondavg + 1
   condavgname(icondavg) = 'COR'
   condavglongname(icondavg) = 'core'
   icondavg_cor = icondavg

   icondavg = icondavg + 1
   condavgname(icondavg) = 'CDN'
   condavglongname(icondavg) = 'downdraft core'
   icondavg_cordn = icondavg
end if
           
if(dosatupdnconditionals) then
   icondavg = icondavg + 1
   condavgname(icondavg) = 'SUP'
   condavglongname(icondavg) = 'saturated updrafts'
   icondavg_satup = icondavg

   icondavg = icondavg + 1
   condavgname(icondavg) = 'SDN'
   condavglongname(icondavg) = 'saturated downdrafts'
   icondavg_satdn = icondavg

   icondavg = icondavg + 1
   condavgname(icondavg) = 'ENV'
   condavglongname(icondavg) = 'unsaturated environment'
   icondavg_env = icondavg
end if
           
! END UW ADDITIONS
!===============================================================

irecc = 1


contains

FUNCTION ProductOf235(n) RESULT(res)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    LOGICAL :: res
    INTEGER :: temp
    if(n == 1) then
     res = .true.  ! Default result is false
     return
    else
     res = .false.  ! Default result is false
    end if
    temp = n
    DO WHILE (MOD(temp, 2) == 0)
        temp = temp / 2
    END DO
    DO WHILE (MOD(temp, 3) == 0)
        temp = temp / 3
    END DO
    DO WHILE (MOD(temp, 5) == 0)
        temp = temp / 5
    END DO

    IF (temp == 1) THEN
        res = .TRUE.  ! Number is a product of only 2, 3, and 5
    END IF
END FUNCTION ProductOf235


end
