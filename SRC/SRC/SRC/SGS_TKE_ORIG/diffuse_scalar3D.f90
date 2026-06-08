subroutine diffuse_scalar3D (field,fluxb,fluxt,tkh,rho,rhow,flux)

use grid
use params, only: docolumn,dosgs,doimplicitdiff,tk_factor,cfl_diff_max,cfl_diffsc_max
use sgs, only: grdf_x,grdf_y,tkmax
use terrain
use vars, only: misc, misc1
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
real rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp,rat
real tkx,tky,tkz,rhoi,tk0
integer i,j,k,ib,ic,jb,jc,kc,kb

if(.not.dosgs) return

dfdt(:,:,:)=0.

!-----------------------------------------

rat = cfl_diffsc_max/cfl_diff_max

!  Horizontal diffusion:

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)
rdz2=1./(dz*dz)
rdz=1./dz

do k=1,nzm

! damping in stratosphere
 if(z(k).gt.20000.) then
  tk0 = 0.001*(dz*adz(k))**2/dtn
 else
  tk0 = 0.
 end if

 do j=1,ny
  rdx5=rdx2*imu(j)**2 *grdf_x(j,k)
  do i=0,nx
    ic=i+1
    tkx=max(tk0,0.5*(tkh(i,j,k)+tkh(ic,j,k))) 
    tkx=min(rat*tkmax(j,k),tkx)
    flx(i,j,k)=-rdx5*tkx*(field(ic,j,k)-field(i,j,k))*terrau(ic,j,k)
  end do 
  do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))*terra(i,j,k)
  end do 
 end do 

 do j=0,ny
  jc=j+1
  rdy5 = rdy2/adyv(jc)*muv(jc) *grdf_y(j,k)         
  do i=1,nx
   tky=max(tk0,0.5*(tkh(i,j,k)+tkh(i,jc,k)))
   tky=min(rat*tkmax(j,k),tky)
   flx(i,j,k)=-rdy5*tky*(field(i,jc,k)-field(i,j,k))*terrav(i,jc,k)
  end do 
 end do
 do j=1,ny
  jb=j-1
  rdy5 = 1./(ady(j)*mu(j))
  do i=1,nx	    
    dfdt(i,j,k)=dfdt(i,j,k)-rdy5*(flx(i,j,k)-flx(i,jb,k))*terra(i,j,k)
  end do 
 end do 
 
end do ! k

! skip vertical diffusion if implicit scheme in vertical:
if(doimplicitdiff) then
  do k=1,nzm
    do j=1,ny
     do i=1,nx
       field(i,j,k)=field(i,j,k)+dtn*dfdt(i,j,k)*terra(i,j,k)
     end do
    end do
  end do
  return
end if

!  Vertical diffusion:

do k=1,nzm-1
 kc=k+1	
 rdz5 = rdz2*rhow(kc)/adzw(kc)
 do j=1,ny
  do i=1,nx
    tkz=min(rat*tkmax(j,k),0.5*(tkh(i,j,k)+tkh(i,j,kc)))  
    flx(i,j,k)=-rdz5*tkz*(field(i,j,kc)-field(i,j,k))*terraw(i,j,kc)  
    flux(kc) = flux(kc) + flx(i,j,k)*wgtw(j,kc)*terraw(i,j,kc)
  end do 
 end do
end do

tmp=1./adzw(nz)
do j=1,ny
 do i=1,nx	
   k=k_terra(i,j)
   flx(i,j,k-1)=fluxb(i,j)*rdz*rhow(k)
   flux(k) = flux(k) + flx(i,j,k-1)*wgtw(j,k)*terraw(i,j,k)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do j=1,ny
  do i=1,nx		 
   field(i,j,k)=field(i,j,k)+ &
          dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)*terra(i,j,k)
  end do 
 end do 	 
end do 

if(collect_coars) then
   misc(1:nx,1:ny,1:nzm) = flx(1:nx,1:ny,0:nzm-1)*dz
   do k=1,nzm
   kc=k+1
   do j=1,ny
    do i=1,nx
      misc1(i,j,kc)=min(rat*tkmax(j,k),tkh(i,j,k))
    end do
   end do
  end do

end if


end subroutine diffuse_scalar3D
