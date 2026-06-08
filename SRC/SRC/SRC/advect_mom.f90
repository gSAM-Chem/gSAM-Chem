subroutine advect_mom

use vars
use params, only: docolumn, nadv_mom, doterrain
use check_energy
use check_mom

implicit none
integer i,j,k
real, allocatable ::  du(:,:,:,:)

if(docolumn) then
  tkeleadv = 0.
  momleadv = 0.
  twleadv = 0.
  qwleadv = 0.
  twleadv = 0.
  qwleadv = 0.
  uwle = 0.
  vwle = 0.
  return
end if

call t_startf ('advect_mom')

if(dostatis) then
	
  allocate(du(nx,ny,nz,3))
 
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na) 
     du(i,j,k,2)=dvdt(i,j,k,na) 
     du(i,j,k,3)=dwdt(i,j,k,na) 
    end do
   end do
  end do

endif

call sumenergy(1,'energy:adv_mom')
call summom(1,'mom:adv_mom')

! compute grid Jacobians:

do k=1,nzm
 do j=1,ny
   gu(j,k) = mu(j)*rho(k)*ady(j)*adz(k)
   gv(j,k) = max(1.e-5_8,muv(j))*rho(k)*adyv(j)*adz(k)
   gw(j,k) = mu(j)*rhow(k)*ady(j)*adzw(k)
   igu(j,k) = 1._8/gu(j,k)
   igv(j,k) = 1._8/gv(j,k)
   igw(j,k) = 1._8/gw(j,k)
 end do
end do

if(nadv_mom.eq.3) then
  call advect23_mom_xy()
  call advect23_mom_z()
else if(nadv_mom.eq.23) then
  call advect23_mom_xy()
  call advect23_mom_z()
else 
  call advect2_mom_xy()
  call advect2_mom_z()
end if

call sumenergy(-1,'energy:adv_mom')
call summom(-1,'mom:adv_mom')

if(dostatis) then
	
  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)-du(i,j,k,1)
     du(i,j,k,2)=dvdt(i,j,k,na)-du(i,j,k,2)
     du(i,j,k,3)=dwdt(i,j,k,na)-du(i,j,k,3)
    end do
   end do
  end do
  du(:,:,nz,3)=0.

  call stat_tke(du,tkeleadv)
  call stat_mom(du,momleadv)
  call setvalue(twleadv,nzm,0.)
  call setvalue(qwleadv,nzm,0.)
  call stat_sw1(du,twleadv,qwleadv)

  deallocate (du)

endif

call t_stopf ('advect_mom')

end subroutine advect_mom

