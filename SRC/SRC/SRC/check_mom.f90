module check_mom

implicit none
real(8) summom_u,summom_v,summom_w
real(8) summomt_u,summomt_v,summomt_w

contains

subroutine summom(flag,title)
use vars
use params, only: docheckmom, dodebug
use terrain, only: terrau, terrav, terraw
integer flag
character*(*) title
real(8) buf(15), buf1(15)
real(8) totmom_u,totmom_v,totmom_w
real(8) totmomt_u,totmomt_v,totmomt_w
real(8) totter_u,totter_v,totter_w
integer k,j
if(.not.(docheckmom.or.dodebug)) return
if(flag.eq.1) then
  summom_u = 0.
  summom_v = 0.
  summom_w = 0.
  summomt_u = 0.
  summomt_v = 0.
  summomt_w = 0.
  do k=1,nzm
   do j=1,ny
     summom_u = summom_u - sum(dble(dudt(1:nx,j,k,na)*terrau(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     summom_v = summom_v - sum(dble(dvdt(1:nx,j,k,na)*terrav(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     summom_w = summom_w - sum(dble(dwdt(1:nx,j,k,na)*terraw(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
     summomt_u = summomt_u - sum(dble(dudt(1:nx,j,k,na)*(1.-terrau(1:nx,j,k))))*mu(j)*rho(k)*adz(k)*ady(j)
     summomt_v = summomt_v - sum(dble(dvdt(1:nx,j,k,na)*(1.-terrav(1:nx,j,k))))*muv(j)*rho(k)*adz(k)*adyv(j)
     summomt_w = summomt_w - sum(dble(dwdt(1:nx,j,k,na)*(1.-terraw(1:nx,j,k))))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
end if
if(flag.eq.-1) then
  totmom_u = 0.
  totmom_v = 0.
  totmom_w = 0.
  totmomt_u = 0.
  totmomt_v = 0.
  totmomt_w = 0.
  do k=1,nzm
   do j=1,ny
     summom_u = summom_u + sum(dble(dudt(1:nx,j,k,na)*terrau(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     summom_v = summom_v + sum(dble(dvdt(1:nx,j,k,na)*terrav(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     summom_w = summom_w + sum(dble(dwdt(1:nx,j,k,na)*terraw(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
     summomt_u = summomt_u + sum(dble(dudt(1:nx,j,k,na)*(1.-terrau(1:nx,j,k))))*mu(j)*rho(k)*adz(k)*ady(j)
     summomt_v = summomt_v + sum(dble(dvdt(1:nx,j,k,na)*(1.-terrav(1:nx,j,k))))*muv(j)*rho(k)*adz(k)*adyv(j)
     summomt_w = summomt_w + sum(dble(dwdt(1:nx,j,k,na)*(1.-terraw(1:nx,j,k))))*mu(j)*rhow(k)*adzw(k)*ady(j)
     totmom_u = totmom_u + sum(dble(u(1:nx,j,k)*terrau(1:nx,j,k)))*mu(j)*rho(k)*adz(k)*ady(j)
     totmom_v = totmom_v + sum(dble(v(1:nx,j,k)*terrav(1:nx,j,k)))*muv(j)*rho(k)*adz(k)*adyv(j)
     totmom_w = totmom_w + sum(dble(w(1:nx,j,k)*terraw(1:nx,j,k)))*mu(j)*rhow(k)*adzw(k)*ady(j)
     totmomt_u = totmomt_u + sum(dble(u(1:nx,j,k)*(1.-terrau(1:nx,j,k))))*mu(j)*rho(k)*adz(k)*ady(j)
     totmomt_v = totmomt_v + sum(dble(v(1:nx,j,k)*(1.-terrav(1:nx,j,k))))*muv(j)*rho(k)*adz(k)*adyv(j)
     totmomt_w = totmomt_w + sum(dble(w(1:nx,j,k)*(1.-terraw(1:nx,j,k))))*mu(j)*rhow(k)*adzw(k)*ady(j)
   end do
  end do
  totmom_u = max(1.e-5_8,totmom_u) ! avoid division by 0
  totmom_v = max(1.e-5_8,totmom_v)
  totmom_w = max(1.e-5_8,totmom_w)
  totmomt_u = max(1.e-5_8,totmomt_u) ! avoid division by 0
  totmomt_v = max(1.e-5_8,totmomt_v)
  totmomt_w = max(1.e-5_8,totmomt_w)
  if(dompi) then
    buf(1) = summom_u
    buf(2) = summom_v
    buf(3) = summom_w
    buf(4) = totmom_u
    buf(5) = totmom_v
    buf(6) = totmom_w
    buf(7) = summomt_u
    buf(8) = summomt_v
    buf(9) = summomt_w
    buf(10) = totmomt_u
    buf(11) = totmomt_v
    buf(12) = totmomt_w
    buf(13) = sum(dble(terrau(1:nx,1:ny,1:nzm)))
    buf(14) = sum(dble(terrav(1:nx,1:ny,1:nzm)))
    buf(15) = sum(dble(terraw(1:nx,1:ny,1:nzm)))
    call task_sum_real8(buf,buf1,15)    
    summom_u = buf1(1)
    summom_v = buf1(2)
    summom_w = buf1(3)
    totmom_u = buf1(4)
    totmom_v = buf1(5)
    totmom_w = buf1(6)
    summomt_u = buf1(7)
    summomt_v = buf1(8)
    summomt_w = buf1(9)
    totmomt_u = buf1(10)
    totmomt_v = buf1(11)
    totmomt_w = buf1(12)
    totter_u = buf1(13)
    totter_v = buf1(14)
    totter_w = buf1(15)
  end if
  if(masterproc) then
   print*,'total fraction of terrain cells (u,v,w):',(1.-totter_u/dble(nx_gl*ny_gl*nzm)), &
                 (1.-totter_v/dble(nx_gl*ny_gl*nzm)), (1.-totter_w/dble(nx_gl*ny_gl*nzm))
   print*,'absolute change of MOM in ',title,' (u,v,w,tot):',summom_u,summom_v, &
                 summom_w, summom_u+summom_v+summom_w
   print*,'momentum  (u,v,w,tot):',totmom_u,totmom_v, &
                 totmom_w, totmom_u+totmom_v+totmom_w
   print*,'relative change of MOM in ',title,' (u,v,w,tot):',summom_u/totmom_u,summom_v/totmom_v, &
                 summom_w/totmom_w, (summom_u+summom_v+summom_w)/(totmom_u+totmom_v+totmom_w)
   print*,'TERR:absolute change of MOM in ',title,' (u,v,w,tot):',summomt_u,summomt_v, &
                 summomt_w, summomt_u+summomt_v+summomt_w
   print*,'TERR: momentum  (u,v,w,tot):',totmomt_u,totmomt_v, &
                 totmomt_w, totmomt_u+totmomt_v+totmomt_w
   print*,'TERR:relative change of MOM in ',title,' (u,v,w,tot):',summomt_u/totmomt_u,summomt_v/totmomt_v, &
                 summomt_w/totmomt_w, (summomt_u+summomt_v+summomt_w)/(totmomt_u+totmomt_v+totmomt_w)
  end if
end if

end subroutine summom

end module check_mom
