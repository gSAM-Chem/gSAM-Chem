
subroutine AdamsB
	
!       Adams-Bashorth scheme Part B
!       apply explicit pressure gradient terms after implicit step

use grid
use vars, only: u, v, w, p
use params, only: gamma_RAVE, dodebug, donodynamics, dofixdynamics
use terrain, only: terrau, terrav, terraw
implicit none
	
real(8) rdx,rdy,rdz
integer i,j,k,kb,jb,ib
real(8) igam2

if(donodynamics.or.dofixdynamics) return

call t_startf('adamsB')

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
   u(i,j,k)=u(i,j,k)-dt3(na)*(bt*(p(i,j,k,nb)-p(ib,j,k,nb))+ct*(p(i,j,k,nc)-p(ib,j,k,nc)))*rdx
   v(i,j,k)=v(i,j,k)-dt3(na)*(bt*(p(i,j,k,nb)-p(i,jb,k,nb))+ct*(p(i,j,k,nc)-p(i,jb,k,nc)))*rdy
   w(i,j,k)=w(i,j,k)-dt3(na)*igam2*(bt*(p(i,j,k,nb)-p(i,j,kb,nb))+ct*(p(i,j,k,nc)-p(i,j,kb,nc)))*rdz
   !   if(isnan(u(i,j,k)).or.isnan(v(i,j,k)).or.isnan(w(i,j,k))) then
   !     print*,'adansB NaN!',rank,i,j,k,u(i,j,k),v(i,j,k),w(i,j,k)
   !     call task_abort()
   !   end if
  end do ! i
 end do ! j	
end do ! k

if(dodebug) then
  if(masterproc) print*,'after AdamsB:'
  call fminmax_print('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz)
end if
call t_stopf('adamsB')


end subroutine adamsB



