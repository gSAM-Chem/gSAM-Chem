
subroutine cloud

!  Condensation of cloud water

use vars
use microphysics
use micro_params
use params

implicit none

integer i,j,k, kb, kc
real dtabs, tabs1, an, bn, ap, bp, om, omp
real fac1,fac2  
real fff,dfff,qsatt,dqsat,powi1,powi2,dq
real lstarn,dlstarn,lstarp,dlstarp
real, external :: esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati
integer niter

ap = 1./(tprmax-tprmin)	
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
! below homogeneous freezing temperature, saturation of water is not for pure water
! but decreases as a function of temperature
an = 1./100.	
bn = (thfreez-100.) * an


do k = 1, nzm
 powi1 = 1./b_ice
 powi2 = (3.+bv_ice(k))/(2.*b_ice)
 do j = 1, ny
  do i = 1, nx

    q(i,j,k)=max(0.,q(i,j,k))


! Initail guess for temperature assuming no cloud water/ice:


    tabs(i,j,k) = t(i,j,k)-gamaz(k)
    tabs1=(tabs(i,j,k)+fac1*qp(i,j,k)+fac_sub*qi(i,j,k))/(1.+fac2*qp(i,j,k))

    qsatt = min(1.,max(0.,an*tabs1-bn))*qsatw(tabs1,pres(k))

!  Test if condensation is possible:

    if(q(i,j,k).gt.qsatt) then

      niter=0
      dtabs = 100.
      do while(abs(dtabs).gt.0.01.and.niter.lt.10)
	   lstarn=fac_cond
	   dlstarn=0.
           om = min(1.,max(0.,an*tabs1-bn))
	   qsatt=om*qsatw(tabs1,pres(k))
	   dqsat=om*dtqsatw(tabs1,pres(k))
	if(tabs1.ge.tprmax) then
	   omp=1.
	   lstarp=fac_cond
	   dlstarp=0.
        else if(tabs1.le.tprmin) then
	   omp=0.
	   lstarp=fac_sub
	   dlstarp=0.
	else
	   omp=ap*tabs1-bp
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   dlstarp=ap*fac_fus
	endif
	fff = tabs(i,j,k)-tabs1+lstarn*(q(i,j,k)-qsatt)+lstarp*qp(i,j,k)+fac_sub*qi(i,j,k)
	dfff=-dlstarp*qp(i,j,k)-lstarn*dqsat-1.
	dtabs=-fff/dfff
	niter=niter+1
	tabs1=tabs1+dtabs
      end do   

      qsatt = qsatt + dqsat * dtabs
      qc(i,j,k) = max(0.,q(i,j,k)-qsatt)

    else

      qc(i,j,k) = 0.

    endif

    tabs(i,j,k) = tabs1
    qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case
    qv(i,j,k) = q(i,j,k) - qc(i,j,k)

    if(tabs(i,j,k).le.thfreez) then
 ! freeze all liquid water below homogeneous-freezing temperature
       qi(i,j,k) = qi(i,j,k) + qc(i,j,k)
       q(i,j,k) = max(0.,q(i,j,k) - qc(i,j,k))
       qc(i,j,k) = 0.
       qpsrc(k) = qpsrc(k) + qc(i,j,k)
    end if
    if(tabs(i,j,k).le.273.15) then
  ! first make sure that the mean size of ice is not less then minimal d_ice_min (or qi_min)
      if(qi(i,j,k).lt.qi_min(k)) then
        q(i,j,k) = max(0.,q(i,j,k)+qi(i,j,k))
        qi(i,j,k) = 0.
        qpsrc(k) = qpsrc(k) - q(i,j,k)
      end if
 ! Deposition/sublimation including Bergeron-Findeisen process
      dq = dtn*(evapi1(k)*qi(i,j,k)**powi1+evapi2(k)*qi(i,j,k)**powi2) &
                        *(qv(i,j,k)/qsati(tabs(i,j,k),pres(k))-1.)
      dq = max(-qi(i,j,k),dq)
      dq = min(q(i,j,k),dq)
      qi(i,j,k) = qi(i,j,k)+dq
      q(i,j,k) = q(i,j,k) - dq
      qpevp(k) = qpevp(k) + dq
    else
      if(doprecip) then
        qp(i,j,k)=qp(i,j,k)+qi(i,j,k)
      else
        qc(i,j,k) = qc(i,j,k)+qi(i,j,k)
        q(i,j,k) = q(i,j,k)+qi(i,j,k)
        qpsrc(k) = qpsrc(k) - qi(i,j,k)
      end if
      qi(i,j,k) = 0.
    end if
  end do
 end do
end do


end subroutine cloud

