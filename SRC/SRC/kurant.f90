
subroutine kurant

use vars
use sgs, only: kurant_sgs
use params, only: ncycle_max, cfl_max, dooceanonly, dosgs

implicit none

integer i, j, k
real cfl,cfll,cflz1,cflh1,cflz,cflh,coef,idx,idy,idz
real buf1(4),buf2(4)

call t_startf ('kurant')

cfl = 0.
cflz = 0.
cflh = 0.
do k = 1,nzm
 idz = dtn/(dz*adzw(k))
 do j=1,ny
  idx = imu(j)*dtn/dx 
  idy = YES3D*dtn/(dy*ady(j))
  do i=1,nx
   cflz1 = abs(w(i,j,k))*idz
   if(cflz1.gt.cflz_max_xy(i,j)) then
      cflz_max_xy(i,j)=cflz1
      zcflz_max_xy(i,j)=zi(k)
   end if
!   cflh1 = abs(u(i,j,k))*idx+abs(v(i,j,k))*idy
   cflh1 = sqrt((u(i,j,k)*idx)**2+(v(i,j,k)*idy)**2)
   if(cflh1.gt.cflh_max_xy(i,j)) then
      cflh_max_xy(i,j)=cflh1
      zcflh_max_xy(i,j)=z(k)
   end if
!   cfll = cflh1+cflz1
   cfll = sqrt(cflh1**2+cflz1**2)
   if(cfll.ge.cfl_max.and.cfll.gt.cfl_max_xy(i,j)) then
      cfl_max_xy(i,j)=cfll
      zcfl_max_xy(i,j)=z(k)
   end if
   cfl = max(cfl,cfll)
   cflz = max(cflz,cflz1)
   cflh = max(cflh,cflh1)
  end do
 end do
end do
w_max=max(w_max,maxval(w(1:nx,1:ny,1:nz)))
u_max=max(u_max,sqrt(maxval(u(1:nx,1:ny,1:nzm)**2+YES3D*v(1:nx,1:ny,1:nzm)**2)))

if(dosgs) then
 call kurant_sgs(cfl_sgs)
else
 cfl_sgs = 0.
end if

cfl_adv=cfl
cfl_advh=cflh
cfl_advz=cflz
if(dompi) then
  buf1(1)=cfl_adv
  buf1(2)=cfl_advh
  buf1(3)=cfl_advz
  buf1(4)= cfl_sgs
  call task_max_real(buf1,buf2,4)
!  if(buf2(1).eq.buf1(1)) print*,'rank of cflmax-',rank,buf1
  cfl_adv=buf2(1)
  cfl_advh=buf2(2)
  cfl_advz=buf2(3)
  cfl_sgs=buf2(4)
end if
dtn = min(dt,dtn*cfl_max/(1.e-10+cfl_adv))
if(ncycle.gt.ncycle_max) then
   if(masterproc) print *,'the number of cycles exceeded ', ncycle_max
   call stepout(-1)
   call write_fields2D()
   call write_fields2DM()
   call write_fields3D()
   call task_abort()
end if

call t_stopf ('kurant')

end subroutine kurant	
