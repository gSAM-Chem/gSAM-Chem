
subroutine diffuse_mom2D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_z
use params, only: docolumn, dowallx
use terrain
implicit none

real rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
real dxz,dzx

integer i,j,k,ic,ib,kc,kcu
real tkx, tkz, rhoi, iadzw, iadz
real fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

rdx2=1./dx/dx
rdx25=0.25*rdx2

dxz=dx/dz

j=1

if(.not.docolumn) then

if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
         v(0,j,k) = v(1,j,k)
         w(0,j,k) = w(1,j,k)
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
         v(nx+1,j,k) = v(nx,j,k)
         w(nx+1,j,k) = w(nx,j,k)
    end do
  end if

end if



do k=1,nzm

 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 rdx21=rdx2 * grdf_x(k)
 rdx251=rdx25 * grdf_x(k)
 
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk(i,j,k)*terra(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*min(terrau(ic,j,k),terrau(i,j,k))
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))*min(terrav(ic,j,k),terrav(i,j,k))
    tkx=rdx251*(tk(i,j,k)*terra(i,j,k)+tk(ic,j,k)*terra(ic,j,k) &
               +tk(i,j,kcu)*terra(i,j,kcu)+tk(ic,j,kcu)*terra(ic,j,kcu)) 	&
            *4./(terra(i,j,k)+terra(ic,j,k)+terra(i,j,kcu)+terra(ic,j,kcu)+1.e-10)
    fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*min(terraw(ic,j,kc),terraw(i,j,kc)) &
                   +(u(ic,j,kcu)-u(ic,j,k))*min(terrau(ic,j,kcu),terrau(ic,j,k))*dxz)
   end do 
   do i=1,nx
    ib=i-1
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(ib,j,k))
   end do  

end do 

end if 

!-------------------------
rdz=1./dz
dzx=dz/dx

uwsb(:) = 0.
vwsb(:) = 0.
	
do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2=rdz*rdz *grdf_z(k)
 rdz25=0.25*rdz2
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk(i,j,k)*terra(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz*min(terraw(i,j,kc),terraw(i,j,k))
    tkz=rdz25*(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
             *4./(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(ib,j,kc))*min(terraw(i,j,kc),terraw(ib,j,kc))*dzx)*rhow(kc) 	
    fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))*iadzw*rhow(kc)
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
end do

do i=1,nx
 k=k_terra(i,j)
 tkz=rdz2*grdf_z(nzm)*tk(i,j,nzm)
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
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
  end do
end do ! k

do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
  end do
end do ! k


end subroutine diffuse_mom2D


