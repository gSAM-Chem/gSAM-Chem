! compute minimum and maximum divergence of the wind and print the result

subroutine check_div(u,v,w)

use grid
use vars, only: rho, rhow
implicit none

real, intent(inout) :: u(nx+1, ny, nzm)
real, intent(inout) :: v(nx, ny+YES3D, nzm)
real, intent(inout) :: w(nx, ny, nz)

real div(nx,ny,nzm), divmax(1), divmin(1), divmax1(1), divmin1(1)
real rdx, rdy, rdz, coef
integer i,j,k,jc,ic

 if(ny.ne.1) then
  do k=1,nzm
   coef = rho(k)*adz(k)*dz
   rdz = 1./coef
   do j=1,ny
     jc = j+1*YES3D
     rdx = imu(j)/dx
     rdy = imu(j)/(dy*ady(j))
     do i=1,nx
      ic = i+1
      div(i,j,k) = (u(ic,j,k)-u(i,j,k))*rdx + (muv(jc)*v(i,jc,k)-muv(j)*v(i,j,k))*rdy + &
                   (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     end do
   end do
  end do
 else
   j = 1
   rdx = 1./dx
   do k=1,nzm
     coef = rho(k)*adz(k)*dz
     rdz = 1./coef
      do i=1,nx-1
      ic = i+1
       div(i,j,k) = (u(ic,j,k)-u(i,j,k))*rdx +(w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
      end do
   end do
 endif
 divmax(1) = maxval(div)
 divmin(1) = minval(div)

  if(dompi) then
     call task_max_real(divmax(1),divmax1(1),1)
     call task_min_real(divmin(1),divmin1(1),1)
     divmax(1) = divmax1(1)
     divmin(1) = divmin1(1)
  end if

if(masterproc) write(6,*) 'div:',divmax(1),divmin(1)

end subroutine check_div

