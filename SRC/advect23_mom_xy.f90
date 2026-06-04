subroutine advect23_mom_xy		
	
!  momentum tendency due to hybrid 2nd-order centered and 3rd-upstream-biased horizontal advection

use vars
use params, only: gamma_RAVE, dowally, alpha_hybrid
use terrain, only: alphah

implicit none
	
real fu(0:nx,1-YES3D:ny) 
real fv(0:nx,1-YES3D:ny)
real fw(0:nx,1-YES3D:ny)
real d12, d4, uuu, wg
integer i, j, k

d12 = 1. / 12.
d4 = 1. / 4.

if(RUN3D) then

!--------------------------------
! Advection of U in x:

do k = 1,nzm	
 do j = 1, ny	
  do i = 0, nx 
   uuu = u1(i+1,j,k)+u1(i,j,k) 
   wg = alphah(i,j,k)
   if(uuu.ge.0) then
      fu(i,j) = 2.*u(i+1,j,k)+5.*u(i,j,k)-u(i-1,j,k)
   else
      fu(i,j) = 2.*u(i,j,k)+5.*u(i+1,j,k)-u(i+2,j,k) 
   end if
   fu(i,j) = uuu*(d12*(1.-wg)*fu(i,j)+d4*wg*(u(i,j,k)+u(i+1,j,k)))
  end do 
  do i = 1, nx	  
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - igu(j,k)*(fu(i,j)-fu(i-1,j))
  end do 
 end do 
end do

!--------------------------------
! Advection of V in x:

do k = 1,nzm
 do j = 1, ny
  do i = 0, nx
   uuu = u1(i+1,j,k)+u1(i+1,j-1,k)
   wg = 0.25*(alphah(i,j,k)+alphah(i,j-1,k)+alphah(i+1,j,k)+alphah(i+1,j-1,k))
   if(uuu.ge.0.) then
      fv(i,j) = 2.*v(i+1,j,k)+5.*v(i,j,k)-v(i-1,j,k)
   else
      fv(i,j) = 2.*v(i,j,k)+5.*v(i+1,j,k)-v(i+2,j,k)
   end if
   fv(i,j) = uuu*(d12*(1.-wg)*fv(i,j)+d4*wg*(v(i,j,k)+v(i+1,j,k)))
  end do
  do i = 1, nx
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - igv(j,k)*(fv(i,j)-fv(i-1,j))
  end do
 end do
end do

!--------------------------------
! Advection of W in x:

do k = 2,nzm
 do j = 1, ny
  do i = 0, nx
   uuu = u1(i+1,j,k)+u1(i+1,j,k-1)
   wg = 0.25*(alphah(i,j,k)+alphah(i,j,k-1)+alphah(i+1,j,k)+alphah(i+1,j,k-1))
   if(uuu.ge.0) then
      fw(i,j) = 2.*w(i+1,j,k)+5.*w(i,j,k)-w(i-1,j,k)
   else
      fw(i,j) = 2.*w(i,j,k)+5.*w(i+1,j,k)-w(i+2,j,k)
   end if
   fw(i,j) = uuu*(d12*(1.-wg)*fw(i,j)+d4*wg*(w(i,j,k)+w(i+1,j,k)))
  end do
  do i = 1, nx
    dwdt(i,j,k,na) = dwdt(i,j,k,na) - igw(j,k)*(fw(i,j)-fw(i-1,j))*gamma_RAVE**2
  end do
 end do
end do

!--------------------------------
! Advection of U in y:

do k=1,nzm
 do j = 0, ny 
  do i = 1, nx
   uuu = v1(i,j+1,k)+v1(i-1,j+1,k) 
   wg = 0.25*(alphah(i,j,k)+alphah(i,j+1,k)+alphah(i-1,j,k)+alphah(i-1,j+1,k)) 
   if(uuu.ge.0.) then
      fu(i,j) = 2.*u(i,j+1,k)+5.*u(i,j,k)-u(i,j-1,k)
   else
      fu(i,j) = 2.*u(i,j,k)+5.*u(i,j+1,k)-u(i,j+2,k)
   end if
   fu(i,j) = uuu*(d12*(1.-wg)*fu(i,j)+d4*wg*(u(i,j,k)+u(i,j+1,k)))
  end do
 end do 
 do j = 1,ny	
  do i = 1, nx
   dudt(i,j,k,na) = dudt(i,j,k,na) - igu(j,k)*(fu(i,j) - fu(i,j-1))
  end do
 end do 
end do ! k

!--------------------------------
! Advection of V in y:

do k=1,nzm
 do j = 0, ny
  do i = 1, nx
   uuu = v1(i,j+1,k)+v1(i,j,k)
   wg = alphah(i,j,k)
   if(uuu.ge.0.) then
      fv(i,j) = 2.*v(i,j+1,k)+5.*v(i,j,k)-v(i,j-1,k)
   else
      fv(i,j) = 2.*v(i,j,k)+5.*v(i,j+1,k)-v(i,j+2,k)
   end if
   fv(i,j) = uuu*(d12*(1.-wg)*fv(i,j)+d4*wg*(v(i,j,k)+v(i,j+1,k)))
  end do
 end do

 do j = 1,ny
  do i = 1, nx
   dvdt(i,j,k,na) = dvdt(i,j,k,na) - igv(j,k)*(fv(i,j) - fv(i,j-1))
  end do
 end do
end do

!--------------------------------
! Advection of W in y:

do k=2,nzm
 do j = 0, ny
  do i = 1, nx
   uuu = v1(i,j+1,k)+v1(i,j+1,k-1)
   wg = 0.25*(alphah(i,j,k)+alphah(i,j,k-1)+alphah(i,j+1,k)+alphah(i,j+1,k-1))
   if(uuu.ge.0.) then
      fw(i,j) = 2.*w(i,j+1,k)+5.*w(i,j,k)-w(i,j-1,k)
   else
      fw(i,j) = 2.*w(i,j,k)+5.*w(i,j+1,k)-w(i,j+2,k)
   end if
   fw(i,j) = uuu*(d12*(1.-wg)*fw(i,j)+d4*wg*(w(i,j,k)+w(i,j+1,k)))
  end do
 end do
 do j = 1,ny
  do i = 1, nx
   dwdt(i,j,k,na)= dwdt(i,j,k,na)- igw(j,k)*(fw(i,j)-fw(i,j-1))*gamma_RAVE**2
  end do
 end do
end do ! k

!=======================================

! 2D model:

else

j=1

do k = 1,nzm
 do i = 0, nx
   uuu = u1(i+1,j,k)+u1(i,j,k)
   wg = alphah(i,j,k)
   if(uuu.ge.0.) then
      fu(i,j) = 2.*u(i+1,j,k)+5.*u(i,j,k)-u(i-1,j,k)
   else
      fu(i,j) = 2.*u(i,j,k)+5.*u(i+1,j,k)-u(i+2,j,k)
   end if
   fu(i,j) = uuu*(d12*(1.-wg)*fu(i,j)+d4*wg*(u(i,j,k)+u(i+1,j,k)))
 end do
 do i = 1, nx
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - igu(j,k)*(fu(i,j)-fu(i-1,j))
 end do

 do i = 0, nx
   uuu = u1(i+1,j,k)+u1(i,j,k)
   wg = alphah(i,j,k)
   if(uuu.ge.0.) then
      fv(i,j) = 2.*v(i+1,j,k)+5.*v(i,j,k)-v(i-1,j,k)
   else
      fv(i,j) = 2.*v(i,j,k)+5.*v(i+1,j,k)-v(i+2,j,k)
   end if
   fv(i,j) = uuu*(d12*(1.-wg)*fv(i,j)+d4*wg*(v(i,j,k)+v(i+1,j,k)))
 end do
 do i = 1, nx
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - igv(j,k)*(fv(i,j)-fv(i-1,j))
 end do
end do

do k = 2,nzm
 do i = 0, nx
   uuu = u1(i+1,j,k)+u1(i+1,j,k-1)
   wg = 0.25*(alphah(i,j,k)+alphah(i,j,k-1)+alphah(i,j+1,k)+alphah(i,j+1,k-1))
   if(uuu.ge.0.) then
      fw(i,j) = 2.*w(i+1,j,k)+5.*w(i,j,k)-w(i-1,j,k)
   else
      fw(i,j) = 2.*w(i,j,k)+5.*w(i+1,j,k)-w(i+2,j,k)
   end if
   fw(i,j) = uuu*(d12*(1.-wg)*fw(i,j)+d4*wg*(w(i,j,k)+w(i+1,j,k)))
 end do
 do i = 1, nx
    dwdt(i,j,k,na) = dwdt(i,j,k,na) - igw(j,k)*(fw(i,j)-fw(i-1,j))
 end do
end do


endif

end subroutine advect23_mom_xy

