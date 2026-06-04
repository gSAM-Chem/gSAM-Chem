
subroutine zero
	
use vars
use microphysics, only: mkwsb
use sgs, only: sgswsb
use tracers, only: trwsb
use stat_coars, only: coars_init

implicit none
	
integer i,j,k,it,jt
real(8) coef

call t_startf ('zero')
	
dudt(:,:,:,na) = 0.
dvdt(:,:,:,na) = 0.
dwdt(:,:,:,na) = 0.
misc(:,:,:) = 0.

call coars_init()

! compute modified advection velocities:

do k=1,nzm
  do j=dimy1_u,dimy2_u
   coef = rho(k)*adz(k)*ady(j)/dx
   do i=dimx1_u,dimx2_u
     u1(i,j,k) = coef*u(i,j,k)
   end do
  end do
end do
do k=1,nzm
  do j=dimy1_v,dimy2_v
   coef = rho(k)*muv(j)*adz(k)/dy
   do i=dimx1_v,dimx2_v
     v1(i,j,k) = coef*v(i,j,k)
   end do
  end do
end do
do k=1,nz
  do j=dimy1_w,dimy2_w
   coef = rhow(k)*ady(j)*mu(j)/dz
   do i=dimx1_w,dimx2_w
     w1(i,j,k) = coef*w(i,j,k)
   end do
  end do
end do

twsb(:) = 0.
uwsb(:) = 0.
vwsb(:) = 0.
sgswsb(:,:) = 0.
mkwsb(:,:) = 0.
trwsb(:,:) = 0.

call t_stopf ('zero')

end
