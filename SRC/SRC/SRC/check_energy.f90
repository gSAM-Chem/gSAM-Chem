module check_energy

implicit none
real(8) sumtke_u,sumtke_v,sumtke_w

contains

subroutine sumenergy(flag,title)
use vars
use params, only: dodebug,docheckenergy
integer flag
character*(*) title
real(8) buf(3), buf1(3)
integer k,j

if(.not.(docheckenergy.or.dodebug)) return

if(title.eq.'energy:diff_mom') then

 if(flag.eq.1) then
  sumtke_u = 0.
  sumtke_v = 0.
  sumtke_w = 0.
  do k=1,nzm
   do j=1,ny
     sumtke_u = sumtke_u -sum(dble(dudtd(1:nx,j,k)*u(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     sumtke_v = sumtke_v -sum(dble(dvdtd(1:nx,j,k)*v(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     sumtke_w = sumtke_w -sum(dble(dwdtd(1:nx,j,k)*w(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
 end if
 if(flag.eq.-1) then
  do k=1,nzm
   do j=1,ny
     sumtke_u = sumtke_u +sum(dble(dudtd(1:nx,j,k)*u(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     sumtke_v = sumtke_v +sum(dble(dvdtd(1:nx,j,k)*v(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     sumtke_w = sumtke_w +sum(dble(dwdtd(1:nx,j,k)*w(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
  if(dompi) then
    buf(1) = sumtke_u
    buf(2) = sumtke_v
    buf(3) = sumtke_w
    call task_sum_real8(buf,buf1,3)
    sumtke_u = buf1(1)
    sumtke_v = buf1(2)
    sumtke_w = buf1(3)
  end if
  if(masterproc) print*,'rate of change of KE in ',title,' (u,v,w,tot):',sumtke_u,sumtke_v,sumtke_w, &
                 sumtke_u+sumtke_v+sumtke_w
  if(dodebug) call debug(title)
 end if

else

 if(flag.eq.1) then
  sumtke_u = 0.
  sumtke_v = 0.
  sumtke_w = 0.
  do k=1,nzm
   do j=1,ny
     sumtke_u = sumtke_u -sum(dble(dudt(1:nx,j,k,na)*u(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     sumtke_v = sumtke_v -sum(dble(dvdt(1:nx,j,k,na)*v(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     sumtke_w = sumtke_w -sum(dble(dwdt(1:nx,j,k,na)*w(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
 end if
 if(flag.eq.-1) then
  do k=1,nzm
   do j=1,ny
     sumtke_u = sumtke_u +sum(dble(dudt(1:nx,j,k,na)*u(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     sumtke_v = sumtke_v +sum(dble(dvdt(1:nx,j,k,na)*v(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     sumtke_w = sumtke_w +sum(dble(dwdt(1:nx,j,k,na)*w(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
  if(dompi) then
    buf(1) = sumtke_u
    buf(2) = sumtke_v
    buf(3) = sumtke_w
    call task_sum_real8(buf,buf1,3)    
    sumtke_u = buf1(1)
    sumtke_v = buf1(2)
    sumtke_w = buf1(3)
  end if
  if(masterproc) print*,'rate of change of KE in ',title,' (u,v,w,tot):',sumtke_u,sumtke_v,sumtke_w, &
                 sumtke_u+sumtke_v+sumtke_w
  if(dodebug) call debug(title)
 end if

end if 

end subroutine sumenergy

end module check_energy
