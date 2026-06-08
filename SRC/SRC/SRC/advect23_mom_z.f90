
subroutine advect23_mom_z

!  momentum tendency due to hybrid 2nd central and  3rd-biased vertical advection

use vars
use terrain
use params, only: gamma_RAVE, alpha_hybrid

implicit none


real fuz(nx,ny,nz),fvz(nx,ny,nz),fwz(nx,ny,nzm)
integer i, ib, j, jb, k, kc, kb, ka, kd
real d12, d4, www, wg

d12=1./12.
d4=1./4.

do j=1,ny
 do i=1,nx
  fuz(i,j,1) = 0.
  fvz(i,j,1) = 0.
  fuz(i,j,nz) = 0.
  fvz(i,j,nz) = 0.
  fwz(i,j,1) = 0.
  fwz(i,j,nzm) = 0.
 end do
end do

uwle(1) = 0.
vwle(1) = 0.	 

if(RUN3D) then

do k=2,nzm
 kc = min(nzm,k+1)
 kb = k-1
 ka = max(1,k-2)
 uwle(k) = 0.
 vwle(k) = 0.
 do j=1,ny
  jb = j-1
  do i=1,nx
   ib = i-1
   www = w1(i,j,k)+w1(i-1,j,k) 
   wg = 0.25*(alphah(i,j,k)+alphah(ib,j,k)+alphah(i,j,kc)+alphah(ib,j,kc))
   if(www.ge.0.) then
     fuz(i,j,k) = 2.*u(i,j,k)+5.*u(i,j,kb)-u(i,j,ka)
   else
     fuz(i,j,k) = 2.*u(i,j,kb)+5.*u(i,j,k)-u(i,j,kc)
   end if
   fuz(i,j,k) = www*(d12*(1.-wg)*fuz(i,j,k)+d4*wg*(u(i,j,k)+u(i,j,kb)))   
   www = w1(i,j,k)+w1(i,jb,k) 
   wg = 0.25*(alphah(i,j,k)+alphah(i,jb,k)+alphah(i,j,kc)+alphah(i,jb,kc))
   if(www.ge.0.) then
     fvz(i,j,k) = 2.*v(i,j,k)+5.*v(i,j,kb)-v(i,j,ka)
   else
     fvz(i,j,k) = 2.*v(i,j,kb)+5.*v(i,j,k)-v(i,j,kc)  
   end if
   fvz(i,j,k) = www*(d12*(1.-wg)*fvz(i,j,k)+d4*wg*(v(i,j,k)+v(i,j,kb)))
   uwle(k) = uwle(k)+fuz(i,j,k)*wgtu(j,k)*terrau(i,j,k)
   vwle(k) = vwle(k)+fvz(i,j,k)*wgtv(j,k)*terrav(i,j,k)
  end do
 end do
end do

else

j=1

do k=2,nzm
 kc = min(nzm,k+1)
 kb = k-1
 ka = max(1,k-2)
 uwle(k) = 0.
 vwle(k) = 0.
  do i=1,nx
   www = w1(i,j,k)+w1(i-1,j,k) 
   wg = 0.25*(alphah(i,j,k)+alphah(ib,j,k)+alphah(i,j,kc)+alphah(ib,j,kc))
   if(www.ge.0.) then
     fuz(i,j,k) = 2.*u(i,j,k)+5.*u(i,j,kb)-u(i,j,ka)
     fvz(i,j,k) = 2.*v(i,j,k)+5.*v(i,j,kb)-v(i,j,ka)
   else
     fuz(i,j,k) = 2.*u(i,j,kb)+5.*u(i,j,k)-u(i,j,kc)
     fvz(i,j,k) = 2.*v(i,j,kb)+5.*v(i,j,k)-v(i,j,kc)      
   end if
   fuz(i,j,k) = www*(d12*(1.-wg)*fuz(i,j,k)+d4*wg*(u(i,j,k)+u(i,j,kb)))
   fvz(i,j,k) = www*(d12*(1.-wg)*fvz(i,j,k)+d4*wg*(v(i,j,k)+v(i,j,kb)))
   uwle(k) = uwle(k)+fuz(i,j,k)
   vwle(k) = vwle(k)+fvz(i,j,k)
  end do
end do


endif
	
do k=1,nzm
 kc = k+1
 do j=1,ny
  do i=1,nx
   dudt(i,j,k,na)=dudt(i,j,k,na)-igu(j,k)*(fuz(i,j,kc)-fuz(i,j,k))
   dvdt(i,j,k,na)=dvdt(i,j,k,na)-igv(j,k)*(fvz(i,j,kc)-fvz(i,j,k))
  end do
 end do
end do

do k=1,nzm
 kb = max(1,k-1)
 kc = k+1
 kd = min(nz,k+2)
 do j=1,ny
  do i=1,nx
   www=w1(i,j,kc)+w1(i,j,k)
   wg = alphah(i,j,k)
   if(www.ge.0.) then
     fwz(i,j,k)=2.*w(i,j,kc)+5.*w(i,j,k)-w(i,j,kb)
   else
     fwz(i,j,k)=2.*w(i,j,k)+5.*w(i,j,kc)-w(i,j,kd)
   end if
   fwz(i,j,k) = www*(d12*(1.-wg)*fwz(i,j,k)+d4*wg*(w(i,j,kc)+w(i,j,k)))
  end do
 end do
end do

do k=2,nzm
 kb=k-1
 do j=1,ny
  do i=1,nx
   dwdt(i,j,k,na)=dwdt(i,j,k,na)-igw(j,k)*(fwz(i,j,k)-fwz(i,j,kb))*gamma_RAVE**2
  end do
 end do
end do ! k
dwdt(:,:,1,na) = 0.
dwdt(:,:,nz,na) = 0.

end subroutine advect23_mom_z

