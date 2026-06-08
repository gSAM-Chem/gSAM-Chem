subroutine diffuse_mom

!  Interface to the diffusion routines

use vars
use sgs, only: dodns
use params, only: doterrain, dodebug
use check_energy

implicit none
integer i,j,k
real, allocatable ::  du(:,:,:,:)

call t_startf ('diffuse_mom')

if(dostatis) then
	
  allocate (du(nx,ny,nz,3))

  do k=1,nzm
   do j=1,ny
    do i=1,nx
     du(i,j,k,1)=dudt(i,j,k,na)
     du(i,j,k,2)=dvdt(i,j,k,na)
     du(i,j,k,3)=dwdt(i,j,k,na)
    end do
   end do
  end do
  du(:,:,nz,3)=0.

endif

if(dodebug) then
 dudtd = 0. 
 dvdtd = 0. 
 dwdtd = 0. 
end if
call sumenergy(1,'energy:diff_mom')

if(RUN3D) then
 if(dodns) then
   call diffuse_mom3D_DNS()
 else
   if(doterrain) then
      call diffuse_mom3D_TERR()
   else
      call diffuse_mom3D()
   end if   
 end if
else
 call diffuse_mom2D()
endif


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

  call stat_tke(du,tkelediff)
  call stat_mom(du,momlediff)
  call setvalue(twlediff,nzm,0.)
  call setvalue(qwlediff,nzm,0.)
  call stat_sw1(du,twlediff,qwlediff)

  deallocate(du)

endif

call sumenergy(-1,'energy:diff_mom')

call t_stopf ('diffuse_mom')

end subroutine diffuse_mom

