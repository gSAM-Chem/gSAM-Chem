! shallow-water/flood over topography model -MK 02/2025

module runoff

implicit none

real, allocatable ::  h(:,:)   ! water depth, m
real, allocatable ::  h0(:,:)  ! terrain heigh, m

real, allocatable ::  terraul(:,:),terraur(:,:) ! flags for flow over terrain
real, allocatable ::  terravl(:,:),terravr(:,:)

real, allocatable ::  fhx(:,:), fhy(:,:) ! fluxes

CONTAINS

subroutine move_water()

use grid
use params, only: dowallx, dowally, doregion
use mpi_stuff, only: dompi, masterproc
use consts, only: ggr
use vars, only: elevation,landmask
use slm_vars, only: mws
real, parameter :: eps = 1.e-9
real, parameter :: hmin = 0.05
real d, u, v, rdx, rdy
integer i,j

!------------------------------------------------------------
!  Initialize arrays if needed
if(.not.allocated(h0)) then

 allocate (h(-1:nx+2,-1:ny+2))
 allocate (h0(-1:nx+2,-1:ny+2))
 allocate (fhx(nx+1,ny))
 allocate (fhy(nx,ny+1))
 allocate (terraul(0:nx+1,0:ny+1),terraur(0:nx+1,0:ny+1))
 allocate (terravl(0:nx+1,0:ny+1),terravr(0:nx+1,0:ny+1))

 h0(1:nx,1:ny) = elevation(1:nx,1:ny)
 if(dompi) then
   call task_exchange(h0,-1,nx+2,-1,ny+2,1,2,2,2,2)
 else
   call bound_exchange(h0,-1,nx+2,-1,ny+2,1,2,2,2,2)
 end if
 if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
   if(doregion) then
    do i=-1,0
     h0(i,:) = h0(1,:)
    end do
   end if
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   if(doregion) then
    do i=nx+1,nx+2
     h0(i,:) = h0(nx,:)
    end do
   end if
  end if
 end if
 if(RUN3D.and.dowally) then
  if(rank.lt.nsubdomains_x) then
   if(doregion) then
    do j=-1,0
     h0(:,j) = h0(:,1)
    end do
   end if
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
   if(doregion) then
    do j=ny+1,ny+2
     h0(:,j) = h0(:,ny)
    end do
   end if
  end if
 end if

end if

! End initializtion block
!------------------------------------------------------------

h(1:nx,1:ny) = mws(1:nx,1:ny)*0.001  ! convert to m from mm

if(dompi) then
  call task_exchange(h,-1,nx+2,-1,ny+2,1,2,2,2,2)
else
  call bound_exchange(h,-1,nx+2,-1,ny+2,1,2,2,2,2)
end if
if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
   if(doregion) then
    do i=-1,0
     h(i,:) = h(1,:)
    end do
   end if
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   if(doregion) then
    do i=nx+1,nx+2
     h(i,:) = h(nx,:)
    end do
   end if
  end if
end if
if(RUN3D.and.dowally) then
  if(rank.lt.nsubdomains_x) then
   if(doregion) then
    do j=-1,0
     h(:,j) = h(:,1)
    end do
   end if
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
   if(doregion) then
    do j=ny+1,ny+2
     h(:,j) = h(:,ny)
    end do
   end if
  end if
end if


do j=0,ny+1
 do i=0,nx+1
   terraur(i,j) = merge(1.,0.,h0(i,j)+h(i,j).gt.h0(i+1,j))
   terraul(i,j) = merge(1.,0.,h0(i,j)+h(i,j).gt.h0(i-1,j))
   terravr(i,j) = merge(1.,0.,h0(i,j)+h(i,j).gt.h0(i,j+1))
   terravl(i,j) = merge(1.,0.,h0(i,j)+h(i,j).gt.h0(i,j-1))
 end do
end do

rdx = dtn/dx
do j=1,ny
   do i=1,nx+1
     if(terraul(i,j).gt.0..or.terraur(i-1,j).gt.0.) then
       d = h(i-1,j) + h0(i-1,j) - h(i,j) - h0(i,j)
       u = sign(sqrt(real(ggr) * abs(d)),d)
       u = max(-0.25,min(0.25,u*rdx))
       fhx(i,j) = min(0.,u)*h(i,j)*terraul(i,j)+max(0.,u)*h(i-1,j)*terraur(i-1,j)
     else
       fhx(i,j) = 0.
     end if
   end do
end do
rdy = dtn/dy
do j=1,ny+1
   do i=1,nx
     if(terravl(i,j).gt.0..or.terravr(i,j-1).gt.0.) then
       d = h(i,j-1) + h0(i,j-1) - h(i,j) - h0(i,j)
       v = sign(sqrt(real(ggr) * abs(d)),d)
       v = max(-0.25,min(0.25,v*rdy))
       fhy(i,j) = min(0.,v)*h(i,j)*terravl(i,j)+max(0.,v)*h(i,j-1)*terravr(i,j-1)
     else
       fhy(i,j) = 0.
     end if
  end do
end do
do j=1,ny
   do i=1,nx
     mws(i,j) = mws(i,j) - 1000.*imu(j)*((fhx(i+1,j)-fhx(i,j))+(fhy(i,j+1)*muv(j+1)-fhy(i,j)*muv(j)))
   end do
end do
!if(any(fhy(nx/2,1:ny).gt.0.)) then
!write(*,'(8g12.5)') h0(nx/2,1:ny)
!print*
!write(*,'(8g12.5)') h(nx/2,1:ny)
!print*
!write(*,'(8g12.5)') fhy(nx/2,1:ny)
!print*
!write(*,'(8g12.5)') mws(nx/2,1:ny)
!!stop
!end if
!if(any(fhx(:,ny/2).gt.0.)) then
!write(*,'(8g12.5)') h0(1:nx,ny/2)
!print*
!write(*,'(8g12.5)') h(1:nx,ny/2)
!print*
!write(*,'(8g12.5)') fhx(1:nx,ny/2)
!print*
!write(*,'(8g12.5)') mws(1:nx,ny/2)
!!stop
!end if



end subroutine move_water
end module runoff
