! diffusion of scalar specifically for DNS case.
! Specifically there are molecular diffusion fluxes at all solid boundaries
! including walls of the terrain using solid boundaries temperatures.
! MK 2023

subroutine diffuse_scalar3D_DNS (field,fluxb,fluxt,tkh,rho,rhow,flux,dosubtr)

use grid
use params, only: docolumn,dosgs
use sgs, only: tkmax, DIFF_DNS
use terrain
use vars, only: t, gamaz, misc, misc1
implicit none
! input	
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkh(0:nxp1,1-YES3D:nyp1,nzm)	! eddy conductivity
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real rho(nzm)
real rhow(nz)
real flux(nz)
logical dosubtr
! local        
real flx(0:nx,0:ny,0:nzm)
real dfdt(nx,ny,nz)
real(8) rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5
real rhoi,ttt,fmin,tmp
integer i,j,k,ib,ic,jb,jc,kc,kb

if(.not.dosgs) return

rdx2=DIFF_DNS/(dx*dx)
rdy2=DIFF_DNS/(dy*dy)
rdz2=DIFF_DNS/(dz*dz)
rdz = 1./dz

dfdt(:,:,:)=0.

if(docolumn) goto 111 ! no horizontal diffusion for single-column model

if(dosubtr) then
  fmin=minval(field(0:nx+1,1-YES3D:ny+YES3D,1:nzm))
  field=field-fmin
else
  fmin=0.
end if
!-----------------------------------------

!  Horizontal diffusion:

do k=1,nzm

 do j=1,ny
  do i=0,nx
    ic=i+1
    tmp = 2.-terra(ic,j,k)*terra(i,j,k) ! adjust molecular flux at walls (half gridstep)
    flx(i,j,k)=-rdx2*tmp*(field(ic,j,k)-field(i,j,k))
  end do 
  do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
  end do 
 end do 

 do j=0,ny
  jc=j+1
  rdy5 = rdy2/adyv(jc)         
  do i=1,nx
   tmp = 2.-terra(i,jc,k)*terra(i,j,k) ! adjust molecular flux at walls (half gridstep)
   flx(i,j,k)=-rdy5*tmp*(field(i,jc,k)-field(i,j,k))
  end do 
 end do
 do j=1,ny
  jb=j-1
  rdy5 = 1./ady(j)
  do i=1,nx	    
    dfdt(i,j,k)=dfdt(i,j,k)-rdy5*(flx(i,j,k)-flx(i,jb,k))
  end do 
 end do 
 
end do ! k

111 continue

!  Vertical explicit diffusion:

do k=1,nzm-1
 kc=k+1	
 rdz5 = rdz2*rhow(kc)/adzw(kc)
 do j=1,ny
  do i=1,nx
    flx(i,j,k)=-rdz5*(field(i,j,kc)-field(i,j,k))  
    flux(kc) = flux(kc) + flx(i,j,k)*wgtw(j,kc)
  end do 
 end do
end do

tmp=1./adzw(nz)
do j=1,ny
 do i=1,nx	
   k=k_terra(i,j)
   flx(i,j,k-1)=fluxb(i,j)*rdz*rhow(k)
   flux(k) = flux(k) + flx(i,j,k-1)*wgtw(j,k)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
   flx(i,j,0:k-2) = 0.
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do j=1,ny
  do i=1,nx		 
   field(i,j,k)=field(i,j,k)+ &
          dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)*terra(i,j,k) + fmin
  end do 
 end do 	 
end do 

if(collect_coars) then
   misc(1:nx,1:ny,1:nzm) = flx(1:nx,1:ny,0:nzm-1)*dz
   do k=1,nzm
   kc=k+1
   do j=1,ny
    do i=1,nx
      misc1(i,j,kc)=min(tkmax(j,k),tkh(i,j,k))
    end do
   end do
  end do

end if


end subroutine diffuse_scalar3D_DNS
