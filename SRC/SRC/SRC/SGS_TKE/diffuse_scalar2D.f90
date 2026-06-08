subroutine diffuse_scalar2D (field,fluxb,fluxt,tkh,rho,rhow,flux)

use grid
use params, only: docolumn, dosgs, doimplicitdiff
use terrain
implicit none
	
! input
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkh(0:nxp1, 1-YES3D:nyp1, nzm)	! eddy conductivity
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real rho(nzm)
real rhow(nz)
real flux(nz)
    
! local        
real flx(0:nx,1,0:nzm)
real dfdt(nx,ny,nzm) 
real rdx2,rdz2,rdz,rdx5,rdz5,tmp
real dxz,dzx,tkx,tkz,rhoi
integer i,j,k,ib,ic,kc,kb

if(.not.dosgs) return

rdx2=1./(dx*dx)
rdz2=1./(dz*dz)
rdz=1./dz
dxz=dx/dz
dzx=dz/dx

j=1

dfdt(:,:,:)=0.

if(docolumn) goto 111

do k=1,nzm
	
  rdx5=0.5*rdx2     

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

111 continue

! skip vertical diffusion if implicit diffusion  scheme in vertical is used:

if(doimplicitdiff) then

  do k=1,nzm
     do i=1,nx
       field(i,j,k)=field(i,j,k)+dtn*dfdt(i,j,k)*terra(i,j,k)
     end do
  end do

  return

end if


tmp=1./adzw(nz)
do i=1,nx	
   k=k_terra(i,j)
   flx(i,j,k-1)=fluxb(i,j)*rdz*rhow(k)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
   flux(k) = flux(k) + flx(i,j,k-1)
   flx(i,j,0:k-2) = 0.
end do


do k=1,nzm-1
 kc=k+1	
 rhoi = rhow(kc)/adzw(kc)
 rdz5=0.5*rdz2
 do i=1,nx
    tkz=rdz5*(tkh(i,j,k)+tkh(i,j,kc))
    flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi*terraw(i,j,kc)
    flux(kc) = flux(kc) + flx(i,j,k)
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do i=1,nx		 
  field(i,j,k)=field(i,j,k) &
         + dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)*terra(i,j,k)
 end do 
end do 

end subroutine diffuse_scalar2D
