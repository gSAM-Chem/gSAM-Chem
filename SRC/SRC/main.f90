!---------------------------------------------------------------------------------------
! main.f90
! Coded by Marat Khairoutdinov
! Main program that initializes everything and contains a main time loop calling all the
! components of the model.
!---------------------------------------------------------------------------------------

program sam

use vars
use hbuffer
use microphysics
use sgs
use chemistry, only: chem_init, chem_proc
use tracers
use movies, only: init_movies
use cup, only: cup_tend
use terrain

implicit none

integer nn, nstatsteps ! misc local variables

!-------------------------------------------------------------------
! Initilize MPI. determine the rank of the current task and of the neighbour's ranks

call task_init() 

!------------------------------------------------------------------
! These calls set up the walltime computation to help simulations stop
!   before their queue time runs out.  Setting nelapse_minutes to the
!   specified walltime in the job submission file should make this work.
!   Job scripts that pay attention to the ReadyToRestart file can
!   sometimes allow jobs to be automatically re-submitted if nstop has
!   not been reached.
   call get_system_time_real8(time_init_seconds)
   time_old_seconds = time_init_seconds
   if(masterproc) call write_ready_to_restart_file('FALSE')
   oldstep = 0
   if(masterproc) write(*,*) 'Initial System Time = ', time_init_seconds
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()
!-------------------------------------------------------------------
! Initialize timing statistics structures. 

call t_initializef ()
call t_startf ('initialize')

!------------------------------------------------------------------

call init()     ! initialize some simple stuff
call setparm()  ! read namelist PARAMETERS and set many parameters and constants

if(dooceanonly) goto 111 ! ugly goto switch to run ocean-only model

!------------------------------------------------------------------
! Initialize or restart from the previous run:

if(nrestart.eq.0) then  ! new run
   day=day0 
   dtn = dt
   call setgrid() ! initialize vertical grid structure and terrain arrays
   call setsurface() ! set surface mask, sst, terrain, SLM, etc
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then  ! restart from previous run
   call read_all()
   call setgrid() 
   call setsurface()
   call diagnose()
   call buildings_init()
   call sgs_init()
   call micro_init()  ! initialize microphysics
   if (dochem) call chem_init()   ! initialize chemistry

elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setparm() ! overwrite parameters
   call setgrid() 
   call setsurface()
   call diagnose()
   call buildings_init()
   call sgs_init()
   call micro_init()  !initialize microphysics
   if (dochem) call chem_init()   ! initialize chemistry

else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif
call stat_2Dinit()
call stat_Minit()
call init_movies()
call tracers_init() ! initialize tracers
call setforcing()
if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()

!------------------------------------------------------------------
total_water_before = total_water()
total_water_after = total_water()

nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
if(nrestart.eq.0.and.doterrain) call terrain_fill(1)
if(dodynamicocean) call ocn_initialize()

call stepout(-1)

dts = 0.
dtst = 0.

call t_stopf ('initialize')
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------
do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  call kurant()

  time = time + dtn
  if(nstep.eq.0.or.time.ge.(nstep+1)*dt) then
    nstep = nstep + 1
    nstep_run = nstep_run + 1
    icycle = 1
    nelapse = nelapse - 1
    dtfactor = 1.
  else
    icycle = icycle + 1 
    dtfactor = 0.
    ncycle = icycle
  end if
  day = day0 + nstep*dt/86400.
  dtp = dtp + dtn
  dtpa = dtpa + dtn
  dts = dts + dtn
  dtstat = dtstat + dtn
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  total_water_before = total_water()
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.
  total_water_adv = 0.

     dt3(na) = dtn

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
        dtst = dtst + dtn
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()

!---------------------------------------------
!  	initialize stuff: 

     call zero()

!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!-----------------------------------------------------------------------------------
!       cumulus and SGS cloudiness parameterizations:

     if(docup) call cup_tend()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!-----------------------------------------------------------
!       Buoyancy term:
     
     call buoyancy()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) call radiation()     

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!-----------------------------------------------------------
!       Coriolis force:
     
     call coriolis()
 
!----------------------------------------------------------
!    Building physics
 
     if(dobuildings) call buildings_wall_proc()
!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)
!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

!---------------------------------------------------------
!    update provisional velocities.
 
     call adamsA()

!----------------------------------------------------------
!    Damping of momentum and implicit vertical diffusion:

     if(dosgs.and.doimplicitdiff) then

       call diffuse_damping_mom_z()

     else

       call damping()

     end if

!---------------------------------------------------------
!    update provisional velocities by imcoplete pressure gradient
 
     call adamsB()

!---------------------------------------------------------
!    Poisson equation for pressure
!    compute non-divergent velocities on n+1 level

     call pressure()

!----------------------------------------------------------
!     Update wind and scalar boundaries for scalar advection:

     call boundaries(1)
     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

     if(.not.docup.and.docloud.or.dosmoke) then

       call micro_proc()

     end if
!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!    Chemical processes (minus emissions)

      if(dochem) call chem_proc() 

!-----------------------------------------------------------
!    fill-in scalars inside terrain:

      if(doterrain) call terrain_fill(1)

!-----------------------------------------------------------
!   call ocean model

    if(dodynamicocean) call dyn_ocean() 

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

!----------------------------------------------------------

! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

  if(icycle.gt.1) cycle

  total_water_after = total_water()
          
!----------------------------------------------------------
!  collect statistics, write restart files, print running data, etc.

   call stepout(nstatsteps)
  
!----------------------------------------------------------
!  print timing statistics table:
   if(mod(nstep,1000).eq.0) call t_prf(rank)
!----------------------------------------------------------
end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------
! outside the main time loop
!----------------------------------------------------------
!
! run ocean model only

111 if(dooceanonly) then

  call setgrid()
  
  call ocn_initialize()

  call ocn_stepout(1)

  call t_stopf ('initialize')

  do while(nstep.lt.nstop.and.nelapse.gt.0)

    nstep = nstep + 1
    time = time + dt
    day = day0 + nstep*dt/86400.
    nelapse = nelapse - 1

    call ocn(real(dt*ncallocean,4))

  end do

end if
!----------------------------------------------------------
!  End the run.
!----------------------------------------------------------
call t_prf(rank)

if(masterproc) print*,"Run has ended normally!"
call task_barrier() ! make sure every MPI task reached this point before ending.
call task_stop()

end program sam
