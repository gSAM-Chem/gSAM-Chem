subroutine diffuse_mom3D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs
use params, only: docolumn, dowallx, dowally, doimplicitdiff
use check_energy
use stat_coars
use microphysics, only: nmicro_fields
use terrain, only: terra,terrau, terrav, terraw ! for coarse diagnostic compatibility
implicit none

integer i,j,k,ic,ib,jb,jc,kc,kcu
real tkx,tky,tkz,rhoi,rdx2,rdy2,rdz2,rdz,rdz25
real fux(0:nx),fvx(0:nx),fwx(0:nx)
real fuy(nx,0:ny),fvy(nx,0:ny),fwy(nx,0:ny)
real fuz(nx,ny,nz),fvz(nx,ny,nz),fwz(nx,ny,nz)
real rdx2u,rdx2v,rdx2w,rdy2u,rdy2v,rdy2w,rat

call sumenergy(1,'diff_mom3D')

!-----------------------------------------
if(docolumn) goto 111  ! skip horizontal diffusion for single-column model

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)

!---------------------------------------------------

do k=1,nzm
  kc=k+1
  kcu=min(kc,nzm)
  do j=1,ny
   jb=j-1
   rdx2u=rdx2*imu(j)**2
   rdx2v=rdx2*imuv(j)**2
   rdx2w=rdx2*imu(j)**2
   do i=0,nx
    ic=i+1
    tkx=tk(i,j,k)
    tkx=max(tkminu(j,k),min(tkmaxu(j,k),tkx))
    fux(i)=-rdx2u*tkx*(u(ic,j,k)-u(i,j,k))
    tkx=0.25*(tk(i,j,k)+tk(i,jb,k)+tk(ic,j,k)+tk(ic,jb,k)) 
    tkx=max(tkminv(j,k),min(tkmaxv(j,k),tkx))
    fvx(i)=-rdx2v*tkx*(v(ic,j,k)-v(i,j,k))
    tkx=0.25*(tk(i,j,k)+tk(ic,j,k)+tk(i,j,kcu)+tk(ic,j,kcu)) 
    tkx=max(tkminw(j,kcu),min(tkmaxw(j,kcu),tkx))
    fwx(i)=-rdx2w*tkx*(w(ic,j,kc)-w(i,j,kc))
   end do
   do i=1,nx
    ib=i-1
    dudtd(i,j,k) = -(fux(i)-fux(ib))
    dvdtd(i,j,k) = -(fvx(i)-fvx(ib))
    dwdtd(i,j,kc) = -(fwx(i)-fwx(ib))
   end do  
  end do 
end do

do k=1,nzm
  kc=k+1
  kcu=min(kc,nzm)
  do j=0,ny
   jc=j+1
   rdy2u=rdy2/adyv(jc)*muv(jc)
   rdy2v=rdy2/ady(j)*mu(j)
   rdy2w=rdy2/adyv(jc)*muv(jc)
   do i=1,nx
    ib=i-1
    tky=0.25*(tk(i,j,k)+tk(ib,j,k)+tk(i,jc,k)+tk(ib,jc,k)) 
    tky=max(tkminu(j,k),min(tkmaxu(j,k),tky))
    fuy(i,j)=-rdy2u*tky*(u(i,jc,k)-u(i,j,k))
    tky=tk(i,j,k)
    tky=max(tkminv(j,k),min(tkmaxv(j,k),tky))
    fvy(i,j)=-rdy2v*tky*(v(i,jc,k)-v(i,j,k))
    tky=0.25*(tk(i,j,k)+tk(i,jc,k)+tk(i,j,kcu)+tk(i,jc,kcu)) 
    tky=max(tkminw(j,kcu),min(tkmaxw(j,kcu),tky))
    fwy(i,j)=-rdy2w*tky*(w(i,jc,kc)-w(i,j,kc))
   end do 
  end do 
  if(dowally) then
   if(rank.lt.nsubdomains_x) then
     fvy(1:nx,0) = 0.
   end if
   if(rank.gt.nsubdomains-nsubdomains_x-1) then
     fvy(1:nx,ny) = 0.
   end if
  end if
  do j=1,ny
    jb=j-1
    rdy2u=1./(ady(j)*mu(j))
    rdy2v=1./(adyv(j)*muv(j))
    rdy2w=1./(ady(j)*mu(j))
    do i=1,nx    
     dudtd(i,j,k) = dudtd(i,j,k) - rdy2u*(fuy(i,j)-fuy(i,jb))
     dvdtd(i,j,k) = dvdtd(i,j,k) - rdy2v*(fvy(i,j)-fvy(i,jb))
     dwdtd(i,j,kc) = dwdtd(i,j,kc) - rdy2w*(fwy(i,j)-fwy(i,jb))
   end do 
  end do 
end do 

!-------------------------

111 continue

if(doimplicitdiff) goto 333

rdz=1./dz

do k=1,nzm-1
 kc=k+1
 rdz2 = rdz*rdz*rho(k)/adz(k) 
 rdz25 = rdz*rdz*rhow(kc)/adzw(kc)
  do j=1,ny
   jb=j-1
   do i=1,nx
    ib=i-1
    tkz=0.25*(tk(i,j,k)+tk(ib,j,k)+tk(i,j,kc)+tk(ib,j,kc)) 
    tkz=max(tkminu(j,k),min(tkmaxu(j,k),tkz))
    fuz(i,j,kc)=-rdz25*tkz*(u(i,j,kc)-u(i,j,k))
    tkz=0.25*(tk(i,j,k)+tk(i,jb,k)+tk(i,j,kc)+tk(i,jb,kc)) 
    tkz=max(tkminv(j,k),min(tkmaxv(j,k),tkz))
    fvz(i,j,kc)=-rdz25*tkz*(v(i,j,kc)-v(i,j,k))
    tkz=tk(i,j,k)
    tkz=max(tkminw(j,k),min(tkmaxw(j,k),tkz))
    fwz(i,j,kc)=-rdz2*tkz*(w(i,j,kc)-w(i,j,k)) 
  end do 
 end do
end do

fwz(:,:,1) = 0.
fuz(:,:,1) = 0.
fvz(:,:,1) = 0.

do j=1,ny
 do i=1,nx
   fwz(i,j,2) = 0.
   fuz(i,j,1)=fluxbu(i,j) * rdz * rhow(1)
   fvz(i,j,1)=fluxbv(i,j) * rdz * rhow(1)
!   tkz=rdz2*tk(i,j,nzm)
!   fwz(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)
   fwz(i,j,nz)= 0.
   fuz(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
   fvz(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
  end do
 end do

if (dostatis) then
 do k=1,nzm-1
   kc=k+1
   do j=1,ny
    jb=j-1
    do i=1,nx
     uwsb(kc)=uwsb(kc)+fuz(i,j,kc)*wgtu(j,kc)
     vwsb(kc)=vwsb(kc)+fvz(i,j,kc)*wgtv(j,kc)
   end do
  end do
 end do
 do j=1,ny
  do i=1,nx
    uwsb(1) = uwsb(1) + fuz(i,j,1)*wgtu(j,1)
    vwsb(1) = vwsb(1) + fvz(i,j,1)*wgtv(j,1)
  end do
 end do
end if
 
 do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do j=1,ny
   do i=1,nx
    dudtd(i,j,k) = dudtd(i,j,k)-(fuz(i,j,kc)-fuz(i,j,k))*rhoi
    dvdtd(i,j,k) = dvdtd(i,j,k)-(fvz(i,j,kc)-fvz(i,j,k))*rhoi
   end do
  end do
 end do ! k

 do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do j=1,ny
   do i=1,nx 
    dwdtd(i,j,k) = dwdtd(i,j,k)-(fwz(i,j,k+1)-fwz(i,j,k))*rhoi
   end do
  end do
 end do ! k

if(collect_coars) then

  call coars_fld(fuz(1:nx,1:ny,1:nzm)*dz,mu(1:ny),ady(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(0:nx-1,1:ny,1:nzm)),fld_flux_sgs(:,:,:,1))
  call coars_fld(fvz(1:nx,1:ny,1:nzm)*dz,muv(1:ny),adyv(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(1:nx,0:ny-1,1:nzm)),fld_flux_sgs(:,:,:,2))
  call coars_fld(fwz(1:nx,1:ny,1:nzm)*dz,mu(1:ny),ady(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(1:nx,1:ny,2:nz)),fld_flux_sgs(:,:,:,3))
  rat = cfl_diff_max/cfl_diffsc_max
  do k=1,nzm
   kc=k+1
   do j=1,ny
    do i=1,nx
      misc1(i,j,kc)=min(rat*tkmax(j,k),tk(i,j,k))
    end do
   end do
  end do
  call coars_fld(misc1(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny), &
        terra(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+nmicro_fields+1))

end if

333 continue

call sumenergy(-1,'diff_mom3D')


end subroutine diffuse_mom3D
