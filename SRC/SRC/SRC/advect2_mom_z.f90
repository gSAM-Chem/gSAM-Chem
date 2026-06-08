
subroutine advect2_mom_z

!       momentum tendency due to the 2nd-order-central vertical advection

use vars
use terrain
use params, only: gamma_RAVE
use stat_coars

implicit none


real fuz(nx,ny,nz),fvz(nx,ny,nz),fwz(nx,ny,nzm)
integer i, j, jb, k, kc, kb
real www
real(8) coef

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
 kb = k-1
 uwle(k) = 0.
 vwle(k) = 0.
 do j=1,ny
  jb = j-1
  do i=1,nx
   fuz(i,j,k) = 0.25*(w1(i,j,k)+w1(i-1,j,k))*(u(i,j,k)+u(i,j,kb)) 
   fvz(i,j,k) = 0.25*(w1(i,j,k)+w1(i,jb,k))*(v(i,j,k)+v(i,j,kb)) 
   uwle(k) = uwle(k)+fuz(i,j,k)*wgtu(j,k)*terrau(i,j,k)
   vwle(k) = vwle(k)+fvz(i,j,k)*wgtv(j,k)*terrav(i,j,k)
  end do
 end do
end do

else

do k=2,nzm
 kb = k-1
 uwle(k) = 0.
 vwle(k) = 0.
 do j=1,ny
  do i=1,nx
    www = 0.25*(w1(i,j,k)+w1(i-1,j,k))
    fuz(i,j,k) = www*(u(i,j,k)+u(i,j,kb)) 
    fvz(i,j,k) = www*(v(i,j,k)+v(i,j,kb)) 
    uwle(k) = uwle(k)+fuz(i,j,k)
    vwle(k) = vwle(k)+fvz(i,j,k)
  end do
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
 kc = k+1
 do j=1,ny
  do i=1,nx
   fwz(i,j,k)=0.25*(w1(i,j,kc)+w1(i,j,k))*(w(i,j,kc)+w(i,j,k))
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

if(collect_coars) then

 do k=2,nzm
   kb = k-1
   do j=1,ny
    jb = j-1
    do i=1,nx
     fuz(i,j,k) = 0.25*(w(i,j,k)+w(i-1,j,k))*(u(i,j,k)+u(i,j,kb))*rhow(k)
     fvz(i,j,k) = 0.25*(w(i,j,k)+w(i,jb,k))*(v(i,j,k)+v(i,j,kb))*rhow(k)
    end do
   end do
  end do

  do k=1,nzm
   kc = k+1
   do j=1,ny
    do i=1,nx
     fwz(i,j,k)=0.25*(w(i,j,kc)+w(i,j,k))*(w(i,j,kc)+w(i,j,k))*rho(k)
    end do
   end do
  end do

  call coars_fld(fuz(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(0:nx-1,1:ny,1:nzm)),fld_flux(:,:,:,1))
  call coars_fld(fvz(1:nx,1:ny,1:nzm),muv(1:ny),adyv(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(1:nx,0:ny-1,1:nzm)),fld_flux(:,:,:,2))
  call coars_fld(fwz(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny), &
        0.5*(terraw(1:nx,1:ny,1:nzm)+terraw(1:nx,1:ny,2:nz)),fld_flux(:,:,:,3))
  
end if


end subroutine advect2_mom_z

