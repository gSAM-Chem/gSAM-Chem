
subroutine pressz

!
! Compute the reference pressure at height levels from temperature and
! moisture sounding. Mostly the effect of surface pressure change is
! taken into account here, which can be important for CRM's long runs.

use vars
use params
implicit none
integer k

real presr(nz)
	
presr(1)=(pres0/1000.)**(rgas/cp)
presi(1)=pres0

do k=1,nzm
 prespot(k)=(1000./pres(k))**(rgas/cp)
 tv0(k)=tabs0(k)*prespot(k)*(1.+epsv*q0(k))
 presr(k+1)=presr(k)-ggr/cp/tv0(k)*(zi(k+1)-zi(k))
 presi(k+1)=1000.*presr(k+1)**(cp/rgas)
 pres(k) = exp(log(presi(k))+log(presi(k+1)/presi(k))* &
              (z(k)-zi(k))/(zi(k+1)-zi(k)))
 prespot(k)=(1000./pres(k))**(rgas/cp)
end do

do k=1,nzm
  rho(k) = (presi(k)-presi(k+1))/(zi(k+1)-zi(k))/ggr*100.
end do

do k=2,nzm
  rhow(k) = (rho(k-1)*adz(k)+rho(k)*adz(k-1))/(adz(k)+adz(k-1))
end do
rhow(1) = 2*rhow(2) - rhow(3)
rhow(nz)= 2*rhow(nzm) - rhow(nzm-1)

end subroutine pressz
