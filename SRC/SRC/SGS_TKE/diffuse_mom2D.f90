
subroutine diffuse_mom2D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk
use params, only: docolumn, dowallx, doimplicitdiff
use terrain
implicit none

real rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251

integer i,j,k,ic,ib,kc,kcu
real tkx, tkz, rhoi, iadzw, iadz, tk0
real fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

rdx2=1./dx/dx
rdx25=0.25*rdx2

j=1

tk0 = 1./64./dt ! filter 2*delta noise

if(.not.docolumn) then

do k=1,nzm

 kc=k+1
 kcu=min(kc,nzm)
 rdx21=rdx2 
 rdx251=rdx25
 
   do i=0,nx
    ic=i+1
    tkx=max(rdx21*tk(i,j,k),tk0)*terra(i,j,k)
    fu(i,j,k)=-tkx*(u(ic,j,k)-u(i,j,k))*min(terrau(ic,j,k),terrau(i,j,k))
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))*min(terrav(ic,j,k),terrav(i,j,k))
    tkx=(rdx251*(tk(i,j,k)*terra(i,j,k)+tk(ic,j,k)*terra(ic,j,k) &
               +tk(i,j,kcu)*terra(i,j,kcu)+tk(ic,j,kcu)*terra(ic,j,kcu))+tk0) 	&
            *4./(terra(i,j,k)+terra(ic,j,k)+terra(i,j,kcu)+terra(ic,j,kcu)+1.e-10)
    fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*min(terraw(ic,j,kc),terraw(i,j,kc)))
   end do 
   do i=1,nx
    ib=i-1
    dudtd(i,j,k) = -(fu(i,j,k)-fu(ib,j,k))
    dvdtd(i,j,k) = -(fv(i,j,k)-fv(ib,j,k))
    dwdtd(i,j,kc) = -(fw(i,j,k)-fw(ib,j,k))
   end do  

end do 

end if ! .not.docolumn

if(doimplicitdiff) return

!-------------------------
rdz=1./dz

do k=1,nzm-1
 kc=k+1
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2=rdz*rdz
 rdz25=0.25*rdz2
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk(i,j,k)*terra(i,j,k)
    fw(i,j,kc)=-tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz*min(terraw(i,j,kc),terraw(i,j,k))
    tkz=rdz25*(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
             *4./(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))*iadzw)*rhow(kc) 	
    fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))*iadzw*rhow(kc)
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
end do

do i=1,nx
 k=k_terra(i,j)
 tkz=rdz2*tk(i,j,nzm)
 fw(i,j,k) = 0.
 fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)
 fu(i,j,k)=fluxbu(i,j) * rdz * rhow(k)
 fv(i,j,k)=fluxbv(i,j) * rdz * rhow(k)
 fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
 fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
 uwsb(k) = uwsb(k) + fu(i,j,k)
 vwsb(k) = vwsb(k) + fv(i,j,k)
end do
	 

do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do i=1,nx
    dudtd(i,j,k) = dudtd(i,j,k) - (fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdtd(i,j,k) = dvdtd(i,j,k) - (fv(i,j,kc)-fv(i,j,k))*rhoi
  end do
end do ! k

do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do i=1,nx	 
    dwdtd(i,j,k) = dwdtd(i,j,k) - (fw(i,j,k+1)-fw(i,j,k))*rhoi
  end do
end do ! k


end subroutine diffuse_mom2D


