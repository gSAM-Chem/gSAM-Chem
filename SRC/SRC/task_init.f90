subroutine task_init
		
!   Check things, initialize multitasking:

use grid
use microphysics, only: micro_scheme_name
use rad, only: RAD_NAME
use advection, only: ADV_SCHEME
use slm_vars, only: slm_setparm
use sgs, only: sgs_init
use params, only: SLM, docheck,dolongwave,doshortwave
use tracers, only: tracers_init
implicit none

integer itasks,ntasks
integer narg
character(10) :: arg
integer, external :: lenstr


! print simple information on the executable 
! when it is run from command line with -info flag:

call get_command_argument(1,arg)
if(LEN_TRIM(arg).ne.0) then
 if(TRIM(arg).ne."-info".and.TRIM(arg).ne."-namelists") then
   print*,'uknown argument ',TRIM(arg)
   stop
 end if
 if(TRIM(arg).eq."-info") then
   print*,'gSAM version: ',trim(version)
   print*,'microphysics: ',micro_scheme_name()
   print*,'radiation: ',RAD_NAME
   print*,'scalar advection: ',ADV_SCHEME
   print*,'DOMAIN SIZE:'
   print*,'nx_gl =',nx_gl
   print*,'ny_gl =',ny_gl
   print*,'nzm =',nzm
   print*,'nsubdomains_x=',nsubdomains_x
   print*,'nsubdomains_y=',nsubdomains_y
   print*,'number of MPI tasks:',nsubdomains_x*nsubdomains_y
   stop
 end if
 if(TRIM(arg).eq."-namelists") then
   open(8,file='./CaseName',status='old',form='formatted')
   read(8,'(a)') case
   masterproc = .true.
   docheck = .true.
   call t_initializef ()
   call setparm()
   day = day0
   if(SLM) call slm_setparm()
   call setgrid()
   call setsurface() ! set surface mask, sst, terrain, SLM, etc
   call setdata() ! initialize all variables
!   call sgs_init()
   call tracers_init() ! initialize tracers
   call setforcing()
   call printout()
   if(dolongwave.or.doshortwave) call radiation()
   call nudging()
   print*
   print*,'All checks passed. All systems go!'
   print*,'Dont forget to delete the old *.stat file if it is a new run.'
   stop
 end if
end if


if(YES3D .ne. 1 .and. YES3D .ne. 0) then
  print*,'YES3D is not 1 or 0. STOP'
  stop
endif

if(YES3D .eq. 1 .and. ny_gl .lt. 4) then
  print*,'ny_gl is too small for a 3D case.STOP'
  stop
endif

if(YES3D .eq. 0 .and. ny_gl .ne. 1) then
  print*,'ny_gl should be 1 for a 2D case. STOP'
  stop
endif

if(nsubdomains.eq.1) then

  rank =0
  ntasks = 1
  dompi = .false.

else

  dompi = .true.

  call task_start(rank, ntasks)

!  call systemf('hostname')  ! not MT-safe for multithread and parallel program

  call task_barrier()

  call task_ranks()
        
  call task_barrier()

end if ! nsubdomains.eq.1

do itasks=0,nsubdomains-1
   if(itasks.eq.rank) then
    open(8,file='./CaseName',status='old',form='formatted')
    read(8,'(a)') case
    close (8)
   endif
end do

masterproc = rank.eq.0

if(masterproc) print *,'number of MPI tasks:',ntasks
	

end
