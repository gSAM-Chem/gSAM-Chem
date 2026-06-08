subroutine debug(title)
use vars
implicit none
character*(*) title
integer i,j,k

if(masterproc) print*,title
call fminmax_print('dudt:',dudt(:,:,:,na),1,nx,1,ny,nzm)
call fminmax_print('dvdt:',dvdt(:,:,:,na),1,nx,1,ny,nzm)
call fminmax_print('dwdt:',dwdt(:,:,:,na),1,nx,1,ny,nzm)
!call fminmax_print('u1:',u1(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!call fminmax_print('v1:',v1(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!call fminmax_print('w1:',w1(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
!call fminmax_print('rho:',rho(1:nzm),1,1,1,1,nzm)
!call fminmax_print('rhow:',rhow(1:nz),1,1,1,1,nz)
!call fminmax_print('adz:',adz(1:nzm),1,1,1,1,nzm)
!call fminmax_print('adzw:',adzw(1:nz),1,1,1,1,nz)
!call fminmax_print('ady:',real(ady(1:ny)),1,1,1,ny,1)
!call fminmax_print('adyv:',real(adyv(1:ny+1)),1,1,1,ny+1,1)
!call fminmax_print('mu:',real(mu(1:ny)),1,1,1,ny,1)
!call fminmax_print('muv:',real(muv(1:ny+1)),1,1,1,ny+1,1)
if(masterproc) print*


end

