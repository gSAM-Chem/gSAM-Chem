subroutine diffuse_scalar3D (field,fluxb,fluxt,tkh,rho,rhow,flux)

use grid
use params, only: docolumn,dowallx,dowally,dosgs
use sgs, only: grdf_x,grdf_y,grdf_z
use terrain
implicit none
! input	
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkh(0:nxp1,1-YES3D:nyp1,nzm)	! eddy conductivity
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real rho(nzm)
real rhow(nz)
real flux(nz)
! local        
real flx(0:nx,0:ny,0:nzm)
real dfdt(nx,ny,nz)
real rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdy5m,rdz5,rhoi,tmp
real tkx,tky,tkz
integer i,j,k,ib,ic,jb,jc,kc,kb


if(.not.dosgs) return

dfdt(:,:,:)=0.

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=1,ny
         field(0,j,k) = field(1,j,k)
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=1,ny
         field(nx+1,j,k) = field(nx,j,k)
     end do
    end do
  end if

end if

if(dowally) then

  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         field(i,1-YES3D,k) = field(i,1,k)
       end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,nx
         field(i,ny+YES3D,k) = field(i,ny,k)
       end do
    end do
  end if

end if


!-----------------------------------------


!  Horizontal diffusion:

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)
rdz2=1./(dz*dz)
rdz=1./dz


do k=1,nzm
	
 rdy5=0.5*rdy2*grdf_y(k)
 rdx5=0.5*rdx2*grdf_x(k)

 do j=1,ny
  do i=0,nx
    ic=i+1
    tkx=rdx5*(tkh(i,j,k)+tkh(ic,j,k)) 	
    flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))*terrau(ic,j,k)
  end do 
  do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))*terra(i,j,k)
  end do 
 end do 

 do j=0,ny
  jc=j+1
  rdy5m = rdy5*muv(jc)
  do i=1,nx
   tky=rdy5m*(tkh(i,j,k)+tkh(i,jc,k)) 	
   flx(i,j,k)=-tky*(field(i,jc,k)-field(i,j,k))*terrav(i,jc,k)
  end do 
 end do
 do j=1,ny
  jb=j-1
  do i=1,nx	    
    dfdt(i,j,k)=dfdt(i,j,k)-imu(j)*(flx(i,j,k)-flx(i,jb,k))*terra(i,j,k)
  end do 
 end do 
 
end do ! k


!  Vertical diffusion:

flux(:) = 0.
tmp=1./adzw(nz)
do j=1,ny
 do i=1,nx	
   k=k_terra(i,j)
   flx(i,j,k-1)=fluxb(i,j)*rdz*rhow(k)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
   flux(k) = flux(k) + flx(i,j,k-1)
 end do
end do


do k=1,nzm-1
 kc=k+1	
 flux(kc)=0. 
 rhoi = rhow(kc)/adzw(kc)
 rdz5=0.5*rdz2 * grdf_z(k)
 do j=1,ny
  do i=1,nx
    tkz=rdz5*(tkh(i,j,k)+tkh(i,j,kc))
    flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi*terraw(i,j,k)
    flux(kc) = flux(kc) + flx(i,j,k)
  end do 
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do j=1,ny
  do i=1,nx		 
   dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
   field(i,j,k)=field(i,j,k)+dfdt(i,j,k)*terra(i,j,k)
  end do 
 end do 	 
end do 

end subroutine diffuse_scalar3D
