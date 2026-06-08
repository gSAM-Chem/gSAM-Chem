
subroutine abcoefs

!      coefficients for the Adams-Bashforth scheme

use grid, only: at, bt, ct, na, nb, nc, dt3, nadams, nrestart 

implicit none

real alpha, beta
real coef
	
call t_startf("abcoefs")

if((nadams.eq.3.or.nadams.eq.23).and.(at.ne.0..and.bt.ne.0..or.nrestart.eq.2)) then
  alpha = dt3(nb) / dt3(na)
  beta = dt3(nc) / dt3(na)
  ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta)
  bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha)
  at = 1. - bt - ct
else if(nadams.eq.2.and.at.ne.0..and.bt.ne.0.or.at.ne.0..and.bt.eq.0) then
  alpha = dt3(nb) / dt3(na)
  at = (1.+2*alpha)/(2.*alpha)
  bt = -1./(2.*alpha)
  ct = 0.
else
  at = 1.
  bt = 0.
  ct = 0.
end if

call t_stopf("abcoefs")
end subroutine abcoefs


