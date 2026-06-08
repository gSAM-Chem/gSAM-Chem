
subroutine diffuse_mom3D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_y, grdf_z
use params, only: docolumn, dowallx, dowally
use terrain
implicit none

real rdx2,rdy2,rdz2,rdz,rdx25,rdy25
real rdx21,rdy21,rdx251,rdy251,rdz25
real dxy,dxz,dyx,dyz,dzx,dzy

integer i,j,k,ic,ib,jb,jc,kc,kcu
real tkx, tky, tkz, rhoi, iadzw, iadz
real fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=1,ny
         v(0,j,k) = v(1,j,k)
         w(0,j,k) = w(1,j,k)
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=1,ny
         v(nx+1,j,k) = v(nx,j,k)
         w(nx+1,j,k) = w(nx,j,k)
     end do
    end do
  end if

end if

if(dowally) then

  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         u(i,1-YES3D,k) = u(i,1,k)
         w(i,1-YES3D,k) = w(i,1,k)
       end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,nx
         u(i,ny+YES3D,k) = u(i,ny,k)
         w(i,ny+YES3D,k) = w(i,ny,k)
       end do
    end do
  end if

end if

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)

rdx25=0.25*rdx2
rdy25=0.25*rdy2

dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz


do k=1,nzm
 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 dyz=dy/(dz*adzw(kc))
  rdy21=rdy2    * grdf_y(k)
  rdy251=rdy25  * grdf_y(k)
  do j=1,ny
   jb=j-1
   rdx21=rdx2*imu(j)**2*grdf_x(k)
   rdx251=rdx25*imu(j)**2*grdf_x(k)
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk(i,j,k)*imu(j)*terra(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))*min(terrau(ic,j,k),terrau(i,j,k))
    tkx=rdx251*imuv(j)*(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
               +tk(ic,j,k)*terra(ic,j,k)+tk(ic,jb,k)*terra(ic,jb,k)) &
              *4./(terra(i,j,k)+terra(i,jb,k)+terra(ic,j,k)+terra(ic,jb,k)+1.e-10)
    fv(i,j,k)=-tkx*((v(ic,j,k)-v(i,j,k))*min(terrav(ic,j,k),terrav(i,j,k)) &
                   +(u(ic,j,k)-u(ic,jb,k))*min(terrau(ic,j,k),terrau(ic,jb,k))*dxy)
    tkx=rdx251*imu(j)*(tk(i,j,k)*terra(i,j,k)+tk(ic,j,k)*terra(ic,j,k) &
               +tk(i,j,kcu)*terra(i,j,kcu)+tk(ic,j,kcu)*terra(ic,j,kcu)) &
              *4./(terra(i,j,k)+terra(ic,j,k)+terra(i,j,kcu)+terra(ic,j,kcu)+1.e-10) 
    fw(i,j,k)=-tkx*((w(ic,j,kc)-w(i,j,kc))*min(terraw(ic,j,kc),terraw(i,j,kc)) &
                   +(u(ic,j,kcu)-u(ic,j,k))*min(terrau(ic,j,kcu),terrau(ic,j,k))*dxz)
   end do 
   do i=1,nx
    ib=i-1
    dudt(i,j,k,na)=dudt(i,j,k,na)-imu(j)*(fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-imuv(j)*(fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-imu(j)*(fw(i,j,k)-fw(ib,j,k))
   end do  
  end do 

  do j=0,ny
   jc=j+1
   do i=1,nx
    ib=i-1
    tky=rdy21*tk(i,j,k)*mu(j)*terra(i,j,k)
    fv(i,j,k)=-2.*tky*(v(i,jc,k)-v(i,j,k))*min(terrav(i,jc,k),terrav(i,j,k))
    tky=rdy251*muv(jc)*(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
               +tk(i,jc,k)*terra(i,jc,k)+tk(ib,jc,k)*terra(ib,jc,k)) &
             *4./(terra(i,j,k)+terra(ib,j,k)+terra(i,jc,k)+terra(ib,jc,k)+1.e-10)
    fu(i,j,k)=-tky*((u(i,jc,k)-u(i,j,k))*min(terrau(i,jc,k),terrau(i,j,k)) &
                   +(v(i,jc,k)-v(ib,jc,k))*min(terrav(i,jc,k),terrav(ib,jc,k))*dyx)
    tky=rdy251*muv(jc)*(tk(i,j,k)*terra(i,j,k)+tk(i,jc,k)*terra(i,jc,k) &
               +tk(i,j,kcu)*terra(i,j,kcu)+tk(i,jc,kcu)*terra(i,jc,kcu)) &
             *4./(terra(i,j,k)+terra(i,jc,k)+terra(i,j,kcu)+terra(i,jc,kcu)+1.e-10)
    fw(i,j,k)=-tky*((w(i,jc,kc)-w(i,j,kc))*min(terraw(i,jc,kc),terraw(i,j,kc)) &
                   +(v(i,jc,kcu)-v(i,jc,k))*min(terrav(i,jc,kcu),terrav(i,jc,k))*dyz)
   end do 
  end do 
  do j=1,ny
    jb=j-1
    do i=1,nx	    
     dudt(i,j,k,na)=dudt(i,j,k,na)-imu(j)*(fu(i,j,k)-fu(i,jb,k))
     dvdt(i,j,k,na)=dvdt(i,j,k,na)-imuv(j)*(fv(i,j,k)-fv(i,jb,k))
     dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-imu(j)*(fw(i,j,k)-fw(i,jb,k))
   end do 
  end do 

end do 
 
!-------------------------
rdz=1./dz
dzx=dz/dx
dzy=dz/dy

uwsb(:) = 0.
vwsb(:) = 0.
	
do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2 = rdz*rdz * grdf_z(k)
 rdz25 = 0.25*rdz2
  do j=1,ny
   jb=j-1
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk(i,j,k)*terra(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*min(terraw(i,j,kc),terraw(i,j,k))*rho(k)*iadz
    tkz=rdz25*(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
             *4./(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(ib,j,kc))*min(terraw(i,j,kc),terraw(ib,j,kc))*dzx)*rhow(kc) 	
    tkz=rdz25*(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
             *4./(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
    fv(i,j,kc)=-tkz*( (v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(i,jb,kc))*min(terraw(i,j,kc),terraw(i,jb,kc))*dzy)*rhow(kc)
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
 end do
end do

do j=1,ny
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
 end do
 
 do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do j=1,ny	  
   do i=1,nx
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
   end do
  end do
 end do ! k

 do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do j=1,ny
   do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
   end do
  end do
 end do ! k


end subroutine diffuse_mom3D
