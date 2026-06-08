
subroutine diffuse_mom3D_DNS
	
!        momentum tendency due to molecular viscosity in DNS mode
!        Based on original SAM subroutine diffuse_mom3D()
!        MK, Sep 2023

use vars
use sgs, only: DIFF_DNS, PR_DNS
use params, only: dowallx, dowally
use terrain
implicit none

real rdx2,rdy2,rdz2,rdz
real dxy,dxz,dyx,dyz,dzx,dzy

integer i,j,k,ic,ib,jb,jc,kc,kcu
real rhoi, iadzw, iadz
real coef1, coef2
real fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

rdx2=DIFF_DNS*PR_DNS/(dx*dx)
rdy2=DIFF_DNS*PR_DNS/(dy*dy)

dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz

rdz=1./dz
rdz2 = DIFF_DNS*PR_DNS*rdz*rdz 
dzx=dz/dx
dzy=dz/dy

!-----------------------------------------
! Horizontal diffusion

do k=1,nzm
 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 dyz=dy/(dz*adzw(kc))
 do j=1,ny
   jb=j-1
   do i=0,nx
    ic=i+1
    fu(i,j,k)=-2.*rdx2*(u(ic,j,k)-u(i,j,k))
    coef1 = 2.-terrav(ic,j,k)*terrav(i,j,k)
    fv(i,j,k)=-rdx2*coef1*(v(ic,j,k)-v(i,j,k)+(u(ic,j,k)-u(ic,jb,k))*dxy)
    coef2 = 2.-terraw(ic,j,kc)*terraw(i,j,kc)
    fw(i,j,k)=-rdx2*coef2*(w(ic,j,kc)-w(i,j,kc)+(u(ic,j,kcu)-u(ic,j,k))*dxz)
   end do 
   do i=1,nx
    ib=i-1
    dudtd(i,j,k)=-(fu(i,j,k)-fu(ib,j,k))
    dvdtd(i,j,k)=-(fv(i,j,k)-fv(ib,j,k))
    dwdtd(i,j,kc)=-(fw(i,j,k)-fw(ib,j,k))
   end do  
 end do 

  do j=0,ny
   jc=j+1
   do i=1,nx
    ib=i-1
    fv(i,j,k)=-2.*rdy2*(v(i,jc,k)-v(i,j,k))
    coef1 = 2.-terrau(i,jc,k)*terrau(i,j,k)
    fu(i,j,k)=-rdy2*coef1*(u(i,jc,k)-u(i,j,k)+(v(i,jc,k)-v(ib,jc,k))*dyx)
    coef2 = 2.-terraw(i,jc,kc)*terraw(i,j,kc)
    fw(i,j,k)=-rdy2*coef2*(w(i,jc,kc)-w(i,j,kc)+(v(i,jc,kcu)-v(i,jc,k))*dyz)
   end do 
  end do 
  do j=1,ny
    jb=j-1
    do i=1,nx	    
     dudtd(i,j,k)=dudtd(i,j,k)-(fu(i,j,k)-fu(i,jb,k))
     dvdtd(i,j,k)=dvdtd(i,j,k)-(fv(i,j,k)-fv(i,jb,k))
     dwdtd(i,j,kc)=dwdtd(i,j,kc)-(fw(i,j,k)-fw(i,jb,k))
   end do 
  end do 

end do 
 
!-------------------------
! Verical diffusion

do k=1,nzm-1
 kc=k+1
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
  do j=1,ny
   jb=j-1
   do i=1,nx
    ib=i-1
    fw(i,j,kc)=-2.*rdz2*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz
    fu(i,j,kc)=-rdz2*( (u(i,j,kc)-u(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(ib,j,kc))*dzx)*rhow(kc) 	
    fv(i,j,kc)=-rdz2*( (v(i,j,kc)-v(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(i,jb,kc))*dzy)*rhow(kc)
  end do 
 end do
end do

fw(:,:,1) = 0.
fu(:,:,1) = 0.
fv(:,:,1) = 0.

do j=1,ny
 do i=1,nx
   fw(i,j,k_terra(i,j)+1) = 0.
   fu(i,j,k_terrau(i,j))=fluxbu(i,j) * rdz * rhow(k_terrau(i,j))
   fv(i,j,k_terrav(i,j))=fluxbv(i,j) * rdz * rhow(k_terrav(i,j))
   fw(i,j,nz)= 0.
   fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
   fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
  end do
 end do

if (dostatis) then
 do k=1,nzm-1
   kc=k+1
   do j=1,ny
    jb=j-1
    do i=1,nx
     uwsb(kc)=uwsb(kc)+fu(i,j,kc)*wgtu(j,kc)*terrau(i,j,kc)
     vwsb(kc)=vwsb(kc)+fv(i,j,kc)*wgtv(j,kc)*terrav(i,j,kc)
   end do
  end do
 end do
 do j=1,ny
  do i=1,nx
    uwsb(k_terrau(i,j)) = uwsb(k_terrau(i,j)) + fu(i,j,k_terrau(i,j)) &
                             *wgtu(j,k_terrau(i,j))*terrau(i,j,k_terrau(i,j))
    vwsb(k_terrav(i,j)) = vwsb(k_terrav(i,j)) + fv(i,j,k_terrav(i,j)) &
                             *wgtv(j,k_terrav(i,j))*terrav(i,j,k_terrav(i,j))
  end do
 end do
end if

 do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do j=1,ny	  
   do i=1,nx
    dudtd(i,j,k)=dudtd(i,j,k)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdtd(i,j,k)=dvdtd(i,j,k)-(fv(i,j,kc)-fv(i,j,k))*rhoi
   end do
  end do
 end do ! k

 do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do j=1,ny
   do i=1,nx	 
    dwdtd(i,j,k)=dwdtd(i,j,k)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
   end do
  end do
 end do ! k


end subroutine diffuse_mom3D_DNS
