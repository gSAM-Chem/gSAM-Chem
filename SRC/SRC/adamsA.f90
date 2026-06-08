
subroutine adamsA

!       Adams-Bashforth scheme
!       compute povisional velocities due to explit terms

use grid, only: na, nb, nc, at, bt, ct, dt3, nx, ny, nzm
use vars, only: u, v, w, dudt, dvdt, dwdt, dudtd, dvdtd, dwdtd
use terrain, only: terrau, terrav, terraw
use params, only: dofixdynamics, donodynamics, gamma_RAVE 

implicit none

integer i,j,k,ib,jb,kb
real(8) igam2

call t_startf('adamsA')

if(.not.(donodynamics.or.dofixdynamics)) then

 igam2 = 1._8/gamma_RAVE**2

 do k=1,nzm
   do j=1,ny
    do i=1,nx
      u(i,j,k) = terrau(i,j,k)*(u(i,j,k)+dt3(na) &
               *(at*dudt(i,j,k,na)+bt*dudt(i,j,k,nb)+ct*dudt(i,j,k,nc)+dudtd(i,j,k)))
      v(i,j,k) = terrav(i,j,k)*(v(i,j,k)+dt3(na) &
               *(at*dvdt(i,j,k,na)+bt*dvdt(i,j,k,nb)+ct*dvdt(i,j,k,nc)+dvdtd(i,j,k))) 
      w(i,j,k) = terraw(i,j,k)*(w(i,j,k)+dt3(na)*igam2 &
               *(at*dwdt(i,j,k,na)+bt*dwdt(i,j,k,nb)+ct*dwdt(i,j,k,nc)+dwdtd(i,j,k)))
    end do
   end do
 end do



else
 
  if(dofixdynamics) call setfixdynamics()

end if

call t_stopf ('adamsA')

end subroutine adamsA

	
