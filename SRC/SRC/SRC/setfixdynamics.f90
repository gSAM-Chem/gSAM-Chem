subroutine setfixdynamics()

! overwrite dynamics with prescribed fields if dosetdynamics = .true.

use grid
use vars, only: u,v,w,rho,rhow
use params, only: dofixdynamics, fixdynamics_type
implicit none

real top, x, x0, xc, zc, pi
real psi(nx+1,nzm+1), buf(1),buf1(1),w_max
integer i,j,k,it,jt

select case (fixdynamics_type)

 case(1)

  ! single prescribed 2D cumulus cloud

  top = 15000. ! height of cloud top
  w_max = 10.
  pi = acos(-1.)

  call task_rank_to_index(rank,it,jt)

  x0 = dx*nx_gl/2.
  psi = 0.
  do k=1,nzm+1
    if(zi(k).le.top) then
     do i=1,nx+1
       x = dx*(i+it-nx_gl/2-1)/x0
       if(x.gt.0.) then
        psi(i,k) = (2.5*abs(x))**0.2*exp(-(2.5*x)**2)*sin(pi*zi(k)/top) 
       else
        psi(i,k) = -(2.5*abs(x))**0.2*exp(-(2.5*x)**2)*sin(pi*zi(k)/top) 
       end if
     end do
    end if
  end do
  if(1+it.eq.1) psi(1,:) = 0.
  if(nx+1+it.eq.nx_gl+1) psi(nx+1,:) = 0.
  do k=1,nzm
   do i=1,nx
    u(i,:,k) = -(psi(i,k+1)-psi(i,k))/(rho(k)*adz(k)*dz)
    w(i,:,k) = (psi(i+1,k)-psi(i,k))/(rhow(k)*dx)
   end do
  end do
  w(:,:,nz) = 0.
  v(:,:,:) = 0.
  buf(1) = maxval(w(1:nx,1:ny,1:nzm)) 
  if(dompi) then
   buf1(1) = buf(1)
   call task_max_real(buf1,buf,1)
  end if
  w(:,:,:) = w(:,:,:)*w_max/buf(1)
  u(:,:,:) = u(:,:,:)*w_max/buf(1)

 case(2)

  ! rectangular cells (rotation along horizontal axes

  top = 15000. ! height of cloud top
  w_max = 10.
  pi = acos(-1.)

  call task_rank_to_index(rank,it,jt)

  xc = dx*(nx_gl/2+1) ! center's x
  zc = top/2. ! center's z
  psi = 0.
  do k=1,nzm+1
    if(zi(k).le.top) then
     do i=1,nx+1
       x = dx*(i+it-1)
       psi(i,k) = sin(pi*(x-xc)/xc)*cos(pi/2.*(zi(k)-zc)/zc)
     end do
    end if
  end do
  if(1+it.eq.1) psi(1,:) = 0.
  if(nx+1+it.eq.nx_gl+1) psi(nx+1,:) = 0.
  do k=1,nzm
   do i=1,nx
    u(i,:,k) = -(psi(i,k+1)-psi(i,k))/(rho(k)*adz(k)*dz)
    w(i,:,k) = (psi(i+1,k)-psi(i,k))/(rhow(k)*dx)
   end do
  end do
  w(:,:,nz) = 0.
  v(:,:,:) = 0.
  buf(1) = maxval(w(1:nx,1:ny,1:nzm))
  if(dompi) then
   buf1(1) = buf(1)
   call task_max_real(buf1,buf,1)
  end if
  w(:,:,:) = w(:,:,:)*w_max/buf(1)
  u(:,:,:) = u(:,:,:)*w_max/buf(1)


end select

end
