#include "fppmacros"

subroutine readsurface(filename,field,dayobs)

! read surface data 

use grid
use params, only: docyclic, cycle_period
implicit none
character*(*), intent(in) ::  filename
real, intent(out) ::  field(nx,ny,2)   ! two time saples if available
real, intent(out) :: dayobs(2) ! days corresponding to field samples
integer i,nn,it,jt,nobs,nx1,ny1
real(4), allocatable ::  fld(:,:)
real(4), allocatable :: days(:)
real coef, day1
if(masterproc) print*,filename
call task_rank_to_index(rank,it,jt)
if(masterproc) print*,'reading data from ',trim(filename)
open(11,file=trim(filename),status='old',form='unformatted',  &
        BUFFEREDYES ACTION='READ')
read(11) nobs
if(masterproc) print*,'nobs=',nobs
read(11) nx1
if(masterproc) print*,'nx (obs)=',nx1
read(11) ny1
if(masterproc) print*,'ny (obs)=',ny1
if(nx1.ne.nx_gl.or.ny1.ne.ny_gl) then
  if(masterproc) print*,'dimensions of data in file '//trim(filename)//' is' // &
              'different from model domain size: nx_gl=',nx_gl, &
              'ny_gl=',ny_gl,'  Quitting...'
  call task_abort()
end if
allocate (days(1:nobs))
read(11) days(1:nobs)
if(masterproc) print*,'days=',days
if(docyclic.and.(days(nobs).ge.cycle_period.or.days(1).lt.0.)) then
  if(masterproc) print*,'error reading file ',filename, &
         'docyclic=',docyclic,'time in should not exceed cycle_period=',cycle_period
  call task_abort
end if
if(days(1).gt.365.) then
    days(1:nobs) = days(1:nobs) - (year0-1)*365
end if
if(docyclic) then
  day1=day-int(day)+mod(int(day),cycle_period)
else
  day1 = day 
end if
if(masterproc) print*,'day1=',day1
nn=1
if(.not.docyclic.and.nobs.gt.1.and.(day1.gt.days(nobs).or.day1.lt.days(1))) then
 if(masterproc) print*,'error reading file ',filename, 'nobs=',nobs,' day1=',day1, &
               ' days(1)=',days(1),' days(nobs)=',days(nobs),'docyclic=',docyclic
 call task_abort
end if
allocate(fld(1:nx_gl,1:ny_gl))
if(nobs.gt.1.and.docyclic.and.(day1.le.days(1).or.day1.gt.days(nobs))) then
  read(11) fld(1:nx_gl,1:ny_gl)
  field(1:nx,1:ny,2) = fld(1+it:nx+it,1+jt:ny+jt)
  do i=1,nobs-2
    read(11)
  end do
  read(11) fld(1:nx_gl,1:ny_gl)
  field(1:nx,1:ny,1) = fld(1+it:nx+it,1+jt:ny+jt)
  if(day1.le.days(1)) then
    dayobs(2) = days(1)+day-day1
    dayobs(1) = days(nobs)+day-day1-cycle_period
  else
    dayobs(2) = days(1)+day-day1+cycle_period
    dayobs(1) = days(nobs)+day-day1
  end if
else
  do i=1,nobs-1
   if(day1.gt.days(i)) then
     nn=i
   endif
  end do
  do i=1,nn-1
     read(11)
  end do
  field = 0.
  read(11) fld(1:nx_gl,1:ny_gl)
  field(1:nx,1:ny,1) = fld(1+it:nx+it,1+jt:ny+jt)
  dayobs(1) = days(nn)+day-day1
  dayobs(2) = 0.
  if(nobs.gt.1) then
   read(11) fld(1:nx_gl,1:ny_gl)
   field(1:nx,1:ny,2) = fld(1+it:nx+it,1+jt:ny+jt)
   dayobs(2) = days(nn+1)+day-day1
  end if
end if
close (11)
if(masterproc) print*,'min/max:',minval(fld),maxval(fld)
deallocate(fld)
if(nobs.gt.1.and.masterproc) print*,'dayobs(1)=',dayobs(1),'dayobs(2)=',dayobs(2)

deallocate(days)

call task_barrier()

end

