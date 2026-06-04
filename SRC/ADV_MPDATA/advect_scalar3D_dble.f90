
subroutine advect_scalar3D_dble (f4, u, v, w, rho, rhow, flux)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
use terrain, only: terra
use vars, only: misc
implicit none


real f4(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real v(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
real w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
real rho(nzm)
real rhow(nz)
real flux(nz)
	
real(8) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real(8) mx (0:nxp1,0:nyp1,nzm)
real(8) mn (0:nxp1,0:nyp1,nzm)
real(8) uuu(-1:nxp3,-1:nyp2,nzm)
real(8) vvv(-1:nxp2,-1:nyp3,nzm)
real(8) www(-1:nxp2,-1:nyp2,nz)

real(8), parameter :: eps = 1.e-9
real(8) g(-1:nyp2,nzm),ig(-1:nyp2,nzm)
integer i,j,k,ic,ib,jc,jb,kc,kb
logical nonos
real(8) rhoz,dtdz


real(8) x1, x2, a, b, y
real aa, a1, a2
real(8) andiff,across,pp,pn
andiff(x1,x2,aa,b)=(abs(aa)-2._8*aa*aa/b)*0.5_8*(x2-x1)
across(x1,a,a1,a2)=0.0625_8*x1*a1*a2/a
pp(y)= max(0._8,y)
pn(y)=-min(0._8,y)
	
nonos = .true.

www(:,:,nz)=0.

do j=-1,nyp2
 do k=1,nzm
  g(j,k) = rho(k)*ady(j)*adz(k)*mu(j)
  ig(j,k) = 1./g(j,k)
 end do
end do

f(:,:,:) = f4(:,:,:)

!-----------------------------------------
	 	 
if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do j=0,nyp1
   jb=j-1
   jc=j+1
   do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),f(i,jb,k),f(i,jc,k))
    mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),f(i,jb,k),f(i,jc,k))
   end do
  end do
 end do
	 
end if  ! nonos

 do k=1,nzm
  do j=-1,nyp2
   do i=-1,nxp3
    uuu(i,j,k)=max(0.,u(i,j,k))*f(i-1,j,k)+min(0.,u(i,j,k))*f(i,j,k)
   end do
  end do
 end do
 do k=1,nzm
  do j=-1,nyp3
   do i=-1,nxp2
    vvv(i,j,k)=max(0.,v(i,j,k))*f(i,j-1,k)+min(0.,v(i,j,k))*f(i,j,k)
   end do
  end do
 end do
 do k=1,nzm
  kb=max(1,k-1)
  do j=-1,nyp2
   do i=-1,nxp2
    www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
   end do
  end do
  flux(k) = 0.
  do j=1,ny
   do i=1,nx
    flux(k) = flux(k) + www(i,j,k)
   end do
  end do
 end do
 if(collect_coars) misc(1:nx,1:ny,1:nzm) = www(1:nx,1:ny,1:nzm)


 do k=1,nzm
  do j=-1,nyp2
   do i=-1,nxp2
      f(i,j,k)=max(0._8,f(i,j,k) &
              -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k)+ &
                www(i,j,k+1)-www(i,j,k))*ig(j,k))
   end do
  end do
 end do 

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do j=0,nyp1
   do i=0,nxp2
    uuu(i,j,k)=andiff(f(i-1,j,k),f(i,j,k),u(i,j,k),g(j,k)+g(j,k)) &
              -across(f(i-1,j+1,k)+f(i,j+1,k)-f(i-1,j-1,k)-f(i,j-1,k), g(j,k)+g(j,k), & 
                      u(i,j,k), v(i-1,j,k)+v(i-1,j+1,k)+v(i,j+1,k)+v(i,j,k)) &
              -across(f(i-1,j,kc)+f(i,j,kc)-f(i-1,j,kb)-f(i,j,kb), g(j,k)+g(j,k), & 
                      u(i,j,k), w(i-1,j,k)+w(i-1,j,k+1)+w(i,j,k)+w(i,j,k+1))
   end do
  end do
 end do


 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do j=0,nyp2
   do i=0,nxp1
    vvv(i,j,k)=andiff(f(i,j-1,k),f(i,j,k),v(i,j,k),g(j,k)+g(j-1,k))
    vvv(i,j,k)=vvv(i,j,k) &
        -across(f(i+1,j-1,k)+f(i+1,j,k)-f(i-1,j-1,k)-f(i-1,j,k), g(j,k)+g(j-1,k), &
         v(i,j,k), u(i,j-1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k)) 
   end do
  end do
 end do
 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do j=0,nyp2
   do i=0,nxp1
    vvv(i,j,k)=vvv(i,j,k) &
        -across(f(i,j-1,kc)+f(i,j,kc)-f(i,j-1,kb)-f(i,j,kb), g(j,k)+g(j-1,k), & 
                      v(i,j,k), w(i,j-1,k)+w(i,j,k)+w(i,j,k+1)+w(i,j-1,k+1)) 
   end do
  end do
 end do

 do k=1,nzm
  kb=max(1,k-1)
  do j=0,nyp1
   jb=j-1
   jc=j+1
   do i=0,nxp1
    www(i,j,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),g(j,k)+g(j,kb)) &
             -across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), g(j,k)+g(j,kb), &
                     w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) &
             -across(f(i,jc,k)+f(i,jc,kb)-f(i,jb,k)-f(i,jb,kb), g(j,k)+g(j,kb), &
                     w(i,j,k), v(i,j,kb)+v(i,jc,kb)+v(i,jc,k)+v(i,j,k))
   end do
  end do
 end do

www(:,:,1) = 0.

!---------- non-osscilatory option ---------------

if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
  do j=0,nyp1
   jc=j+1
   do i=0,nxp1
    ic=i+1
     mx(i,j,k)=(mx(i,j,k)-f(i,j,k))*g(j,k)/ &
                       (pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+ &
                       (pn(vvv(i,jc,k)) + pp(vvv(i,j,k)))+ &
                       (pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)
     mn(i,j,k)=(f(i,j,k)-mn(i,j,k))*g(j,k)/ &
                       (pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+ &
                       (pp(vvv(i,jc,k)) + pn(vvv(i,j,k)))+ &
                       (pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)	
   end do
  end do
 end do


 do k=1,nzm
  do j=1,ny
   do i=1,nxp1
    ib=i-1
    uuu(i,j,k)=pp(uuu(i,j,k))*min(1._8,mx(i,j,k), mn(ib,j,k)) &
             - pn(uuu(i,j,k))*min(1._8,mx(ib,j,k),mn(i,j,k))
   end do
  end do
 end do

 do k=1,nzm
  do j=1,nyp1
   jb=j-1
   do i=1,nx
    vvv(i,j,k)=pp(vvv(i,j,k))*min(1._8,mx(i,j,k), mn(i,jb,k)) &
             - pn(vvv(i,j,k))*min(1._8,mx(i,jb,k),mn(i,j,k))
   end do
  end do
 end do

 do k=1,nzm
  kb=max(1,k-1)
  do j=1,ny
   do i=1,nx
    www(i,j,k)=pp(www(i,j,k))*min(1._8,mx(i,j,k), mn(i,j,kb)) &
             - pn(www(i,j,k))*min(1._8,mx(i,j,kb),mn(i,j,k))
    flux(k) = flux(k) + www(i,j,k)
   end do
  end do
 end do


if(collect_coars) then
  dtdz = dtn/dz
  do k=1,nzm
    do j=1,ny
     rhoz = 1./(dtdz*ady(j)*mu(j))
     do i=1,nx
       misc(i,j,k) = (misc(i,j,k)+www(i,j,k))*rhoz
     end do
    end do
  end do
end if

endif ! nonos

do k=1,nzm
 kc=k+1
 do j=1,ny
  do i=1,nx
 ! MK: added fix for very small negative values (relative to positive values) 
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values is
 !     most likely truncation error.

   f4(i,j,k)=max(0._8,f(i,j,k) &
                  -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k)+ &
                    www(i,j,k+1)-www(i,j,k))*ig(j,k)) !*terra(i,j,k))
  end do
 end do
end do 

end subroutine advect_scalar3D_dble


