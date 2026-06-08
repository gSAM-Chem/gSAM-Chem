
subroutine advect_scalar2D (f, u, w, rho, rhow, flux)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
implicit none


real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
real rho(nzm)
real rhow(nz)
real flux(nz)
	
real mx (0:nxp1,nzm)
real mn (0:nxp1,nzm)
real mx1 (0:nxp1,nzm)
real mn1 (0:nxp1,nzm)
real uuu(-1:nxp3,nzm)
real www(-1:nxp2,nz)

real, parameter :: eps = 1.e-9
real g(nzm),ig(nzm)
integer i,j,k,ic,ib,kc,kb
logical nonos

real x1, x2, a, b, a1, a2, y
real andiff,across,pp,pn,fmin
andiff(x1,x2,a,b)=(abs(a)-2.*a*a/b)*0.5*(x2-x1)
across(x1,a,a1,a2)=0.0625*x1*a1*a2/a
pp(y)= max(0.,y)
pn(y)=-min(0.,y)
	
nonos = .true.

www(:,nz)=0.

j=1

do k=1,nzm
  g(k) = rho(k)*adz(k)
  ig(k) = 1./g(k)
end do

!-----------------------------------------
	 	 
if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
   do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
    mn(i,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
    mx1(i,k) = mx(i,k)
    mn1(i,k) = mn(i,k)
   end do
 end do
	 
end if  ! nonos

 do k=1,nzm
   do i=-1,nxp3
    uuu(i,k)=max(0.,u(i,j,k))*f(i-1,j,k)+min(0.,u(i,j,k))*f(i,j,k)
   end do
 end do
 do k=1,nzm
  kb=max(1,k-1)
   do i=-1,nxp2
    www(i,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
   end do
  flux(k) = 0.
  do i=1,nx
    flux(k) = flux(k) + www(i,k)
  end do
 end do


 do k=1,nzm
   do i=-1,nxp2
      f(i,j,k)=max(0.,f(i,j,k) &
              -(uuu(i+1,k)-uuu(i,k)+www(i,k+1)-www(i,k))*ig(k))
   end do
 end do 

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
   do i=0,nxp2
    ib=i-1
    uuu(i,k)=andiff(f(ib,j,k),f(i,j,k),u(i,j,k),g(k)+g(k)) &
              -across(f(ib,j,kc)+f(i,j,kc)-f(ib,j,kb)-f(i,j,kb), g(k)+g(k), & 
                      u(i,j,k), w(ib,j,k)+w(ib,j,k+1)+w(i,j,k)+w(i,j,k+1))
   end do
 end do


 do k=1,nzm
  kb=max(1,k-1)
   do i=0,nxp1
    ib=i-1
    ic=i+1
    www(i,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),g(k)+g(kb)) &
             -across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), g(k)+g(kb), &
                     w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) 
   end do
 end do

www(:,1) = 0.

!---------- non-osscilatory option ---------------

if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
   do i=0,nxp1
    ic=i+1
     mx(i,k)=(mx(i,k)-f(i,j,k))*g(k)/ &
                       (pn(uuu(ic,k)) + pp(uuu(i,k))+ &
                       (pn(www(i,kc)) + pp(www(i,k)))+eps)
     mn(i,k)=(f(i,j,k)-mn(i,k))*g(k)/ &
                       (pp(uuu(ic,k)) + pn(uuu(i,k))+ &
                       (pp(www(i,kc)) + pn(www(i,k)))+eps)	
   end do
 end do


 do k=1,nzm
   do i=1,nxp1
    ib=i-1
    uuu(i,k)=pp(uuu(i,k))*min(1.,mx(i,k), mn(ib,k)) &
             - pn(uuu(i,k))*min(1.,mx(ib,k),mn(i,k))
   end do
 end do

 do k=1,nzm
  kb=max(1,k-1)
   do i=1,nx
    www(i,k)=pp(www(i,k))*min(1.,mx(i,k), mn(i,kb)) &
             - pn(www(i,k))*min(1.,mx(i,kb),mn(i,k))
    flux(k) = flux(k) + www(i,k)
   end do
 end do


endif ! nonos


do k=1,nzm
 kc=k+1
  do i=1,nx
 ! MK: added fix for very small negative values (relative to positive values) 
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values is
 !     most likely truncation error.

   f(i,j,k)=max(0.,f(i,j,k) &
                  -(uuu(i+1,k)-uuu(i,k)+(www(i,k+1)-www(i,k)))*ig(k))
  end do
end do 

end subroutine advect_scalar2D


