
subroutine cloud

!  Condensation of cloud water/cloud ice.

use vars
use microphysics
use micro_params
use sgs, only: tk
use terrain, only: k_terra
use params
use rad, only: qrad

implicit none

integer i,j,k
real dtabs, tabs1, tabs2, an, bn, ap, bp, om, omp
real fac1,fac2  
real fff,dfff,dqsat
real lstarn,dlstarn,lstarp,dlstarp,rh_homo
integer niter
real, external :: qsatw,qsati,dtqsatw,dtqsati

if(donomicro) return

call t_startf('cloud')

an = 1./(tbgmax-tbgmin)	
bn = tbgmin * an
ap = 1./(tprmax-tprmin)	
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap

do k = 1, nzm
 do j = 1, ny
  do i = 1, nx

    if(k.lt.k_terra(i,j)) cycle

    q(i,j,k)=max(0.,q(i,j,k))

 ! Make new ice at very cold temperatures (below homog. freezing) only
 ! when rel. hum over ice is above rh_homo. Modeled after IFS model
 ! if ice already exists - do as usual, that is no supersaturation over ice
 ! Marat 2023

    if(tabs(i,j,k).lt.235.and.qci(i,j,k).lt.1.e-8) then
     rh_homo = 2.583 - tabs(i,j,k)/207.8  ! from IFS
    else
     rh_homo = 1.
    end if

! Initail guess for temperature assuming no cloud water/ice:

    tabs2 = tabs(i,j,k)
    tabs(i,j,k) = t(i,j,k)-gamaz(k)
    tabs1=(tabs(i,j,k)+fac1*qp(i,j,k))/(1.+fac2*qp(i,j,k))

! Warm cloud:

    if(tabs1.ge.tbgmax) then

      tabs1=tabs(i,j,k)+fac_cond*qp(i,j,k)
      qsatt(i,j,k) = qsatw(tabs1,pp(i,j,k))
!!      qsatt(i,j,k) = qsatw(tabs1,pres(k))

! Ice cloud:

    elseif(tabs1.le.tbgmin) then

      tabs1=tabs(i,j,k)+fac_sub*qp(i,j,k)
      qsatt(i,j,k) = qsati(tabs1,pp(i,j,k))*rh_homo
!!      qsatt(i,j,k) = qsati(tabs1,pres(k))*rh_homo

! Mixed-phase cloud:

    else

      om = an*tabs1-bn
      qsatt(i,j,k) = om*qsatw(tabs1,pp(i,j,k))+(1.-om)*qsati(tabs1,pp(i,j,k))*rh_homo
!!      qsatt(i,j,k) = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))*rh_homo

    endif

!  Test if condensation is possible:

!    if(pres(k).gt.70..and.q(i,j,k).gt.qsatt(i,j,k)) then
    if(q(i,j,k).gt.qsatt(i,j,k)) then

      niter=0
      tabs1 = tabs2 ! better initial guess - use previous temperature -MK Mar 2024
      dtabs = 100.
      do while(abs(dtabs).gt.0.001.and.niter.lt.100)
        if(tabs1.ge.tbgmax) then
           om=1.
           lstarn=fac_cond
           dlstarn=0.
           qsatt(i,j,k)=qsatw(tabs1,pp(i,j,k))
           dqsat=dtqsatw(tabs1,pp(i,j,k))
!!           qsatt(i,j,k)=qsatw(tabs1,pres(k))
!!           dqsat=dtqsatw(tabs1,pres(k))
        else if(tabs1.le.tbgmin) then
           om=0.
           lstarn=fac_sub
           dlstarn=0.
           qsatt(i,j,k)=qsati(tabs1,pp(i,j,k))*rh_homo
           dqsat=dtqsati(tabs1,pp(i,j,k))*rh_homo
!!           qsatt(i,j,k)=qsati(tabs1,pres(k))*rh_homo
!!           dqsat=dtqsati(tabs1,pres(k))*rh_homo
        else
           om=an*tabs1-bn
           lstarn=fac_cond+(1.-om)*fac_fus
           dlstarn=an*fac_fus
           qsatt(i,j,k)=om*qsatw(tabs1,pp(i,j,k))+(1.-om)*qsati(tabs1,pp(i,j,k))
           dqsat=om*dtqsatw(tabs1,pp(i,j,k))+(1.-om)*dtqsati(tabs1,pp(i,j,k))
!!           qsatt(i,j,k)=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))
!!           dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k))
        endif
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
        fff = tabs(i,j,k)-tabs1+lstarn*(q(i,j,k)-qsatt(i,j,k))+lstarp*qp(i,j,k)
        dfff=dlstarn*(q(i,j,k)-qsatt(i,j,k))+dlstarp*qp(i,j,k)-lstarn*dqsat-1.
        dtabs=-fff/dfff
        niter=niter+1
        tabs1=tabs1+dtabs
      end do

      qn(i,j,k) = max(0.,q(i,j,k)-qsatt(i,j,k))

    else

      qn(i,j,k) = 0.

    endif

    tabs(i,j,k) = tabs1
    if(q(i,j,k)-qn(i,j,k).lt.0.) then
      print*,'cloud.f90: q(i,j,k)-qn(i,j,k).lt.0.! rank,i,j,k,landmask,pres,gamaz,q,qn,qp,qsat, &
         t,dtabs,tabs2,tabs,sst,qrad,niter:',&
           rank,i,j,k,landmask(i,j),pp(i,j,k),gamaz(k),q(i,j,k),qn(i,j,k),qp(i,j,k),qsatt(i,j,k),t(i,j,k),dtabs,tabs2,tabs(i,j,k),sst_xy(i,j),qrad(i,j,k),niter 
      call task_abort()
    end if
    qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case
  end do
 end do
end do

call t_stopf('cloud')

end subroutine cloud

