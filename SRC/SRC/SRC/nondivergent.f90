! make input 3D velocity field non-divergent (globally, not locally!)
! MK Dec 2022

subroutine nondivergent(u, v, w)

use grid
use terrain
use vars, only : p, rho, rhow
use params, only : doregion, dolatlon
implicit none

! velcity components. Pay attention to dimensions...

real, intent(inout) :: u(nx+1, ny, nzm)
real, intent(inout) :: v(nx, ny+YES3D, nzm)
real, intent(inout) :: w(nx, ny, nz)

! Local:

real(8) buf1(2), buf2(2), veladj, dxy, coef
real(8) rdx,rdy,rdz,rup,rdn
integer i,j,k,ib,jb,kb,ic,jc,kc,it,jt


if(.not.doregion) then
  if(masterproc) print*,'nondivergent() can be called only when doregion = .true.! Exitting...'
  call task_abort
end if

call t_startf ('nondivergent')

! adjust lateral boundary velocities so there is no net mass flux in/out domain

call task_rank_to_index(rank,it,jt)

buf1(:) = 0.
dxy = dx/dy
if(it.eq.0) then
  do k=1,nzm
    do j=1,ny
     buf1(1) = buf1(1)+u(1,j,k)*rho(k)*adz(k)*ady(j)
     buf1(2) = buf1(2)+rho(k)*adz(k)*ady(j)*terrau(1,j,k)
    end do
   end do
end if
if(it+nx.eq.nx_gl) then
   do k=1,nzm
    do j=1,ny
     buf1(1) = buf1(1)-u(nx+1,j,k)*rho(k)*adz(k)*ady(j)
     buf1(2) = buf1(2)+rho(k)*adz(k)*ady(j)*terrau(nx+1,j,k)
    end do
   end do
end if
if(jt.eq.0) then
   do k=1,nzm
    do i=1,nx
     buf1(1) = buf1(1)+v(i,1,k)*muv(1)*rho(k)*adz(k)*dxy
     buf1(2) = buf1(2)+muv(1)*rho(k)*adz(k)*dxy*terrav(i,1,k)
    end do
   end do
end if
 if(jt+ny.eq.ny_gl) then
   do k=1,nzm
    do i=1,nx
     buf1(1) = buf1(1)-v(i,ny+YES3D,k)*muv(ny+YES3D)*rho(k)*adz(k)*dxy
     buf1(2) = buf1(2)+muv(ny+YES3D)*rho(k)*adz(k)*dxy*terrav(i,ny+YES3D,k)
    end do
   end do
end if
call task_barrier()
if(dompi) then
  call task_sum_real8(buf1,buf2,2)
  veladj = buf2(1)/buf2(2)
  call task_barrier()
else
  veladj = buf1(1)/buf1(2)
end if
if(masterproc) print*,'veladj=',veladj
!if(it.eq.0) u(1,1:ny,1:nzm) = u(1,1:ny,1:nzm)-veladj*terrau(1,1:ny,1:nzm)
!if(it+nx.eq.nx_gl) u(nx+1,1:ny,1:nzm) = u(nx+1,1:ny,1:nzm)+veladj*terrau(nx+1,1:ny,1:nzm)
!if(jt.eq.0) v(1:nx,1,1:nzm) = v(1:nx,1,1:nzm)-veladj*terrav(1:nx,1,1:nzm)
!if(jt+ny.eq.ny_gl) v(1:nx,ny+YES3D,1:nzm) = v(1:nx,ny+YES3D,1:nzm)+veladj*terrav(1:nx,ny+YES3D,1:nzm)


do i=1,nx+1
 coef = (2._8*(i+it)-real(nx_gl,8)-2._8)/real(nx_gl,8)
 u(i,1:ny,1:nzm) = u(i,1:ny,1:nzm)+veladj*coef*terrau(i,1:ny,1:nzm)
end do
do j=1,ny+YES3D
 coef = (2._8*(j+jt)-real(ny_gl,8)-2._8)/real(ny_gl,8)
 v(1:nx,j,1:nzm) = v(1:nx,j,1:nzm)+veladj*coef*terrav(1:nx,j,1:nzm)
end do

call t_stopf('nondivergent')

end subroutine nondivergent

