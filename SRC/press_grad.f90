
subroutine press_grad
	
!       pressure term of the momentum equations

use vars
use params, only: gamma_RAVE
implicit none
	
real(8) rdx,rdy,rdz
integer i,j,k,kb,jb,ib
real(8) igam2
real, allocatable ::  du(:,:,:,:)

igam2 = 1._8/gamma_RAVE**2

do k=1,nzm
 kb=max(1,k-1)
 rdz = 1./(dz*adzw(k))
 do j=1,ny
  jb=j-YES3D
  rdx=imu(j)/dx
  rdy=1./(dy*adyv(j))
  do i=1,nx
   ib=i-1 
   u(i,j,k)=u(i,j,k)-dt3(na)*at*(p(i,j,k,na)-p(ib,j,k,na))*rdx
   v(i,j,k)=v(i,j,k)-dt3(na)*at*(p(i,j,k,na)-p(i,jb,k,na))*rdy
   w(i,j,k)=w(i,j,k)-dt3(na)*igam2*at*(p(i,j,k,na)-p(i,j,kb,na))*rdz
  end do ! i
 end do ! j	
end do ! k

call task_barrier()
	
if(dostatis) then

  allocate(du(nx,ny,nz,3))
  do k=1,nzm
   kb=max(1,k-1)
   rdz = 1./(dz*adzw(k))
   do j=1,ny
    jb=j-YES3D
    rdx=imu(j)/dx
    rdy=1./(dy*adyv(j))
    do i=1,nx
     ib=i-1
     du(i,j,k,1) = (p(i,j,k,na)-p(ib,j,k,na))*rdx
     du(i,j,k,2) = (p(i,j,k,na)-p(i,jb,k,na))*rdy
     du(i,j,k,3) = igam2*at*(p(i,j,k,na)-p(i,j,kb,na))*rdz
    end do ! i
   end do ! j
  end do ! k
  du(:,:,nz,3) = 0.

  call stat_tke(du,tkelepress)
  call stat_mom(du,momlepress)
  call setvalue(twlepres,nzm,0.)
  call setvalue(qwlepres,nzm,0.)
  call stat_sw1(du,twlepres,qwlepres)
 
  deallocate (du)

endif


end subroutine press_grad



	



