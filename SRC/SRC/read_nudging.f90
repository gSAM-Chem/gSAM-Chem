#include "fppmacros"

subroutine read_nudging(filedir,listfilename,u,v,w,t,q,dayobs)

! read 3D data from nudging files

!Changes:
! BUG: 12/14/25 - reading w changed from zi levels to z levels - MK
! BUG: 12/14/25 - fixed interpolation weights for w (adz in wrong order)  - MK

use grid
use vars, only: rhow
use consts, only: ggr
use params, only: doregion, w3D_pressure, read_meters
use terrain, only: terrau, terrav, terraw
use params, only: docyclic, cycle_period
implicit none
character*(*), intent(in) ::  listfilename
character*(*), intent(in) ::  filedir
real, intent(out) ::  u(nx+1,ny,nzm,2)   
real, intent(out) ::  v(nx,ny+1,nzm,2)   
real, intent(out) ::  w(nx,ny,nz,2)   
real, intent(out) ::  t(nx,ny,nzm,2)   
real, intent(out) ::  q(nx,ny,nzm,2)   
real, intent(out) :: dayobs(2) ! days corresponding to field samples
character(120) filename
integer i,j,k,nn,nobs,nx1,ny1,nz1
real(4), allocatable :: days(:)
real coef, day1
real(8), allocatable :: xx(:), yy(:), xx_u(:), yy_v(:)
allocate(xx(nx_gl),yy(ny_gl),xx_u(nx_gl+1),yy_v(ny_gl+1))
if(read_meters) then
  xx(:) = x_gl(:)
  yy(:) = y_gl(:)
  xx_u(:) = xu_gl(:)
  yy_v(:) = yv_gl(:)
else
  xx(:) = lon_gl(:)
  yy(:) = lat_gl(:)
  xx_u(:) = lonu_gl(:)
  yy_v(:) = latv_gl(:)
end if
if(masterproc) print*,'reading data from ',trim(listfilename)
if(masterproc) print*,'vertical velocity is assumed pressure velocity? w3D_pressure=',w3D_pressure
open(11,file=trim(listfilename),status='old',form='formatted')
read(11,*) nobs
if(masterproc) print*,'nobs=',nobs
allocate (days(1:nobs))
do i=1,nobs
  read(11,*) days(i)
  if(days(1).gt.365.) then
    days(i) = days(i) - (year0-1)*365
  end if
end do
if(docyclic.and.(days(nobs).gt.cycle_period.or.days(1).lt.0.)) then
  print*,'error reading file ',listfilename,'docyclic=',docyclic,'time in should not exceed 365 days'
  call task_abort
end if
if(nobs.gt.1.and.masterproc) print*,'range of days',minval(days),maxval(days), &
                             ' interval: ',days(2)-days(1)
if(docyclic) then
  day1=day-int(day)+mod(int(day),cycle_period)
else
  day1 = day
end if
if(masterproc) print*,'day1=',day1
nn=1
if(.not.docyclic.and.nobs.gt.1.and.(day1.gt.days(nobs).or.day1.lt.days(1))) then
 if(masterproc) print*,'error reading file ',listfilename, 'nobs=',nobs,' day1=',day1, &
               ' days(1)=',days(1),' days(nobs)=',days(nobs),'docyclic=',docyclic
 call task_abort
end if
if(nobs.gt.1.and.docyclic.and.(day1.le.days(1).or.day1.gt.days(nobs))) then
  filename = ""
  read(11,'(a)') filename
  open(12,file=trim(filedir)//"/"//trim(filename),form="unformatted", &
          BUFFEREDYES ACTION='READ')
  if(masterproc) print*,"opened "//trim(filedir)//"/"//trim(filename)
  read(12) nx1,ny1,nz1
  call read_field3D (12,u(1:nx+1,1:ny,1:nzm,2),nx1,ny1,nz1,xx_u,yy,z,nx+1,ny,nzm)
  call read_field3D (12,v(1:nx,1:ny+1,1:nzm,2),nx1,ny1,nz1,xx,yy_v,z,nx,ny+1,nzm)
  call read_field3D (12,w(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call read_field3D (12,t(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call read_field3D (12,q(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  close(12)
  u(1:nx+1,1:ny,1:nzm,2) = u(1:nx+1,1:ny,1:nzm,2) * terrau(1:nx+1,1:ny,1:nzm)
  v(1:nx,1:ny+1,1:nzm,2) = v(1:nx,1:ny+1,1:nzm,2) * terrav(1:nx,1:ny+1,1:nzm)
  if(w3D_pressure) then
   do k=nzm,2,-1
    w(1:nx,1:ny,k,2) = -(adz(k-1)*w(1:nx,1:ny,k,2)+adz(k)*w(1:nx,1:ny,k-1,2)) &
              /(adz(k)+adz(k-1))/(rhow(k)*ggr)! convert from Ps/s to m/s
   end do
  end if
  w(1:nx,1:ny,1,2) = 0.
  w(1:nx,1:ny,nz,2) = 0.
  w(1:nx,1:ny,1:nzm,2) = w(1:nx,1:ny,1:nzm,2) * terraw(1:nx,1:ny,1:nzm)
  if(doregion) then
      if(masterproc) print*,'before don-divergence correction (1):'
      call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
      call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
      call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
      call nondivergent(u(1:nx+1,1:ny,1:nzm,2),v(1:nx,1:ny+1,1:nzm,2),w(1:nx,1:ny,1:nz,2))
      call check_div(u,v,w)
  end if
  call fminmax_print('t:',t(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
  call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,2),1,nx+1,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,2),1,nx,1,ny+1,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
  call fminmax_print('q:',q(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
  do i=1,nobs-2
    read(11,*)
  end do
  filename = ""
  read(11,'(a)') filename
  open(12,file=trim(filedir)//"/"//trim(filename),form="unformatted", &
          BUFFEREDYES ACTION='READ')
  if(masterproc) print*,"opened "//trim(filedir)//"/"//trim(filename)
  read(12) nx1,ny1,nz1
  call read_field3D (12,u(1:nx+1,1:ny,1:nzm,1),nx1,ny1,nz1,xx_u,yy,z,nx+1,ny,nzm)
  call read_field3D (12,v(1:nx,1:ny+1,1:nzm,1),nx1,ny1,nz1,xx,yy_v,z,nx,ny+1,nzm)
  call read_field3D (12,w(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call read_field3D (12,t(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call read_field3D (12,q(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  close(12)
  u(1:nx+1,1:ny,1:nzm,1) = u(1:nx+1,1:ny,1:nzm,1) * terrau(1:nx+1,1:ny,1:nzm)
  v(1:nx,1:ny+1,1:nzm,1) = v(1:nx,1:ny+1,1:nzm,1) * terrav(1:nx,1:ny+1,1:nzm)
  if(w3D_pressure) then
   do k=nzm,2,-1
    w(1:nx,1:ny,k,1) = -(adz(k-1)*w(1:nx,1:ny,k,1)+adz(k)*w(1:nx,1:ny,k-1,1)) &
              /(adz(k)+adz(k-1))/(rhow(k)*ggr)! convert from Ps/s to m/s
   end do
  end if
  w(1:nx,1:ny,1,1) = 0.
  w(1:nx,1:ny,nz,1) = 0.
  w(1:nx,1:ny,1:nzm,1) = w(1:nx,1:ny,1:nzm,1) * terraw(1:nx,1:ny,1:nzm)
  if(doregion) then
      if(masterproc) print*,'before don-divergence correction (2):'
      call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
      call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
      call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
      call nondivergent(u(1:nx+1,1:ny,1:nzm,1),v(1:nx,1:ny+1,1:nzm,1),w(1:nx,1:ny,1:nz,1))
      call check_div(u,v,w)
  end if
  call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
  call fminmax_print('t:',t(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
  call fminmax_print('q:',q(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
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
     read(11,*)
  end do
  if(dayobs(1).eq.0..and.dayobs(2).eq.0) then
   filename = ""
   read(11,'(a)') filename
   if(masterproc) print*,'nn=',nn,trim(filename)
   open(12,file=trim(filedir)//"/"//trim(filename),form="unformatted", &
           BUFFEREDYES ACTION='READ')
   if(masterproc) print*,"opened "//trim(filedir)//"/"//trim(filename)
   read(12) nx1,ny1,nz1
   call read_field3D (12,u(1:nx+1,1:ny,1:nzm,1),nx1,ny1,nz1,xx_u,yy,z,nx+1,ny,nzm)
   call read_field3D (12,v(1:nx,1:ny+1,1:nzm,1),nx1,ny1,nz1,xx,yy_v,z,nx,ny+1,nzm)
   call read_field3D (12,w(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   call read_field3D (12,t(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   call read_field3D (12,q(1:nx,1:ny,1:nzm,1),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   close(12)
   u(1:nx+1,1:ny,1:nzm,1) = u(1:nx+1,1:ny,1:nzm,1) * terrau(1:nx+1,1:ny,1:nzm)
   v(1:nx,1:ny+1,1:nzm,1) = v(1:nx,1:ny+1,1:nzm,1) * terrav(1:nx,1:ny+1,1:nzm)
   if(w3D_pressure) then
    do k=nzm,2,-1
     w(1:nx,1:ny,k,1) = -(adz(k-1)*w(1:nx,1:ny,k,1)+adz(k)*w(1:nx,1:ny,k-1,1)) &
               /(adz(k)+adz(k-1))/(rhow(k)*ggr)! convert from Ps/s to m/s
    end do
   end if
   w(1:nx,1:ny,1,1) = 0.
   w(1:nx,1:ny,nz,1) = 0.
   w(1:nx,1:ny,1:nzm,1) = w(1:nx,1:ny,1:nzm,1) * terraw(1:nx,1:ny,1:nzm)
   if(doregion) then
      if(masterproc) print*,'before don-divergence correction (3):'
      call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
      call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
      call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
      call nondivergent(u(1:nx+1,1:ny,1:nzm,1),v(1:nx,1:ny+1,1:nzm,1),w(1:nx,1:ny,1:nz,1))
      call check_div(u,v,w)
   end if
  else
   read(11,*) 
   if(nobs.gt.1) then
     u(1:nx+1,1:ny,1:nzm,1) = u(1:nx+1,1:ny,1:nzm,2)
     v(1:nx,1:ny+1,1:nzm,1) = v(1:nx,1:ny+1,1:nzm,2)
     w(1:nx,1:ny,1:nzm,1) = w(1:nx,1:ny,1:nzm,2)
     t(1:nx,1:ny,1:nzm,1) = t(1:nx,1:ny,1:nzm,2)
     q(1:nx,1:ny,1:nzm,1) = q(1:nx,1:ny,1:nzm,2)
   end if
  end if
  call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
  call fminmax_print('t:',t(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
  call fminmax_print('q:',q(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
  dayobs(1) = days(nn)
  dayobs(2) = 0.
  if(nobs.gt.1) then
   filename = ""
   read(11,'(a)') filename
   open(12,file=trim(filedir)//"/"//trim(filename),form="unformatted", &
           BUFFEREDYES ACTION='READ')
   if(masterproc) print*,"opened "//trim(filedir)//"/"//trim(filename)
   read(12) nx1,ny1,nz1
   call read_field3D (12,u(1:nx+1,1:ny,1:nzm,2),nx1,ny1,nz1,xx_u,yy,z,nx+1,ny,nzm)
   call read_field3D (12,v(1:nx,1:ny+1,1:nzm,2),nx1,ny1,nz1,xx,yy_v,z,nx,ny+1,nzm)
   call read_field3D (12,w(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   call read_field3D (12,t(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   call read_field3D (12,q(1:nx,1:ny,1:nzm,2),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
   close(12)
   u(1:nx+1,1:ny,1:nzm,2) = u(1:nx+1,1:ny,1:nzm,2) * terrau(1:nx+1,1:ny,1:nzm)
   v(1:nx,1:ny+1,1:nzm,2) = v(1:nx,1:ny+1,1:nzm,2) * terrav(1:nx,1:ny+1,1:nzm)
   if(w3D_pressure) then
    do k=nzm,2,-1
      w(1:nx,1:ny,k,2) = -(adz(k-1)*w(1:nx,1:ny,k,2)+adz(k)*w(1:nx,1:ny,k-1,2)) &
                /(adz(k)+adz(k-1))/(rhow(k)*ggr)! convert from Ps/s to m/s
    end do
   end if
   w(1:nx,1:ny,1,2) = 0.
   w(1:nx,1:ny,nz,2) = 0.
   w(1:nx,1:ny,1:nzm,2) = w(1:nx,1:ny,1:nzm,2) * terraw(1:nx,1:ny,1:nzm)
   if(doregion) then
      if(masterproc) print*,'before don-divergence correction (4):'
      call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,1),1,nx+1,1,ny,nzm)
      call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,1),1,nx,1,ny+1,nzm)
      call fminmax_print('w:',w(1:nx,1:ny,1:nzm,1),1,nx,1,ny,nzm)
      call nondivergent(u(1:nx+1,1:ny,1:nzm,2),v(1:nx,1:ny+1,1:nzm,2),w(1:nx,1:ny,1:nz,2))
      call check_div(u,v,w)
   end if
   call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm,2),1,nx+1,1,ny,nzm)
   call fminmax_print('v:',v(1:nx,1:ny+1,1:nzm,2),1,nx,1,ny+1,nzm)
   call fminmax_print('w:',w(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
   call fminmax_print('t:',t(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
   call fminmax_print('q:',q(1:nx,1:ny,1:nzm,2),1,nx,1,ny,nzm)
   dayobs(2) = days(nn+1)
  end if
end if
close (11)
deallocate(days)
deallocate(xx, yy, xx_u, yy_v)
end




