
subroutine press_rhs(ppp)

!       right-hand-side of the Poisson equation for pressure

use vars
use terrain, only: terrau, terrav, terraw
use params, only: dowallx, dowally, doterrain, doregion

implicit none
         
real, intent(out) :: ppp(nx,ny,nzm)
	
real *8 dta,rdx,rdy,rdz,rup,rdn
integer i,j,k,ic,jc,kc,it,jt
real tmp


if(.not.doregion.and.dowallx.and.mod(rank,nsubdomains_x).eq.0) then

    do k=1,nzm
     do j=1,ny
      u(1,j,k) = 0.
     end do
    end do	

end if

if(.not.doregion.and.dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

    do k=1,nzm
     do i=1,nx
      v(i,1,k) = 0.
     end do
    end do	

end if

if(dompi) then
   call task_bound_uv()
else
   call bound_uv()
endif

if(doregion) then
  call task_rank_to_index(rank,it,jt)
  if(it.eq.0) u(1,1:ny,1:nzm) = u_w(1:ny,1:nzm)
  if(it+nx.eq.nx_gl) u(nx+1,1:ny,1:nzm) = u_e(1:ny,1:nzm)
  if(jt.eq.0) v(1:nx,1,1:nzm) = v_s(1:nx,1:nzm)
  if(jt+ny.eq.ny_gl) v(1:nx,ny+YES3D,1:nzm) = v_n(1:nx,1:nzm)
end if

dta=1._8/(dt3(na)*at)

if(RUN3D) then

do k=1,nzm
 kc=k+1 
 rdz=1./(adz(k)*dz)
 rup = rhow(kc)/rho(k)*rdz
 rdn = rhow(k)/rho(k)*rdz
 do j=1,ny
  jc=j+1 
  rdx=imu(j)/dx
  rdy=imu(j)/(dy*ady(j))
  do i=1,nx
   ic=i+1
   ppp(i,j,k)=( rdx*(u(ic,j,k)-u(i,j,k))+ &
                rdy*(muv(jc)*v(i,jc,k)-muv(j)*v(i,j,k))+ &
                (w(i,j,kc)*rup-w(i,j,k)*rdn) )*dta
!   if(isnan(ppp(i,j,k))) then
!      print*,'ppp is NaN! (rank,i,j,k):',rank,i,j,k, &
!             u(ic,j,k),u(i,j,k),v(i,jc,k),v(i,j,k),w(i,j,kc),w(i,j,k)
!      stop
!   end if
  end do
 end do
end do

else

j=1

rdx = 1./dx
do k=1,nzm
 kc=k+1 
 rdz=1./(adz(k)*dz)
 rup = rhow(kc)/rho(k)*rdz
 rdn = rhow(k)/rho(k)*rdz
 do i=1,nx
  ic=i+1
  ppp(i,j,k)=(rdx*(u(ic,j,k)-u(i,j,k))+(w(i,j,kc)*rup-w(i,j,k)*rdn))*dta 
 end do
end do


endif

call task_barrier()

end subroutine press_rhs
