
subroutine advect2_mom_xy		
	
!        momentum tendency due to 2nd-order-central horizontal advection

use vars
use params, only: gamma_RAVE, dowally

implicit none
	
real fu(0:nx,1-YES3D:ny) 
real fv(0:nx,1-YES3D:ny)
real fw(0:nx,1-YES3D:ny)

integer i, j, k, kc, kcu, ic, jb, ib, jc

if(RUN3D) then

do k = 1,nzm	
 kc= k+1
 kcu =min(kc, nzm)
	 
 do j = 1, ny	
  jb = j-1
  do i = 0, nx 
   ic = i+1			
   fu(i,j)=0.25*(u1(ic,j,k)+u1(i,j,k))*(u(i,j,k)+u(ic,j,k)) 
   fv(i,j)=0.25*(u1(ic,j,k)+u1(ic,jb,k))*(v(i,j,k)+v(ic,j,k))   
   fw(i,j)=0.25*(u1(ic,j,k)+u1(ic,j,kcu))*(w(i,j,kc)+w(ic,j,kc)) 
  end do 
  do i = 1, nx	  
   ib = i-1
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - igu(j,k)*(fu(i,j)-fu(ib,j))
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - igv(j,k)*(fv(i,j)-fv(ib,j))
    dwdt(i,j,kc,na) = dwdt(i,j,kc,na) - igw(j,kcu)*(fw(i,j)-fw(ib,j))*gamma_RAVE**2
  end do 
 end do 

 do j = 0, ny 
  jc = j+1	
  do i = 1, nx
   ib = i-1
   fu(i,j)=0.25*(v1(i,jc,k)+v1(ib,jc,k))*(u(i,j,k)+u(i,jc,k)) 
   fv(i,j)=0.25*(v1(i,jc,k)+v1(i,j,k))*(v(i,j,k)+v(i,jc,k)) 
   fw(i,j)=0.25*(v1(i,jc,k)+v1(i,jc,kcu))*(w(i,j,kc)+w(i,jc,kc)) 
  end do
 end do 

 do j = 1,ny	
  jb = j-1
  do i = 1, nx
   dudt(i,j,k,na) = dudt(i,j,k,na) - igu(j,k)*(fu(i,j) - fu(i,jb))
   dvdt(i,j,k,na) = dvdt(i,j,k,na) - igv(j,k)*(fv(i,j) - fv(i,jb))
   dwdt(i,j,kc,na)= dwdt(i,j,kc,na)- igw(j,kcu)*(fw(i,j)-fw(i,jb))*gamma_RAVE**2
  end do
 end do 

end do ! k

else

j=1

do k = 1,nzm
 kc= k+1
 kcu =min(kc, nzm)
  do i = 0, nx
   ic = i+1
   fu(i,j)=0.25*(u1(ic,j,k)+u1(i,j,k))*(u(i,j,k)+u(ic,j,k)) 
   fv(i,j)=0.25*(u1(ic,j,k)+u1(i,j,k))*(v(i,j,k)+v(ic,j,k)) 
   fw(i,j)=0.25*(u1(ic,j,k)+u1(ic,j,kcu))*(w(i,j,kc)+w(ic,j,kc)) 
  end do
  do i = 1, nx   
   ib = i-1
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - igu(j,k)*(fu(i,j)-fu(ib,j))
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - igv(j,k)*(fv(i,j)-fv(ib,j))
    dwdt(i,j,kc,na) = dwdt(i,j,kc,na) - igw(j,kcu)*(fw(i,j)-fw(ib,j))
  end do

end do ! k

endif

end subroutine advect2_mom_xy

