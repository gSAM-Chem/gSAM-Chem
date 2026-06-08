! processes involving precipitation in the form of graupel and rain  

subroutine precip_proc

use vars
use microphysics
use micro_params
use params

implicit none

integer i,j,k
real autor, accrr, accrci, accrig, accrcg
real dq, dqp, dqi, omp
real powi, powg1, powg2, powr1, powr2, tmp
real qcc, qrr, qss, qgg
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
real, external :: qsatw,qsati


powr1 = (3 + b_rain) / 4.
powr2 = (5 + b_rain) / 8.
powg1 = (3 + b_grau) / 4.
powg2 = (5 + b_grau) / 8.
      
     
if(dostatis) then
        
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = q(i,j,k)
     end do
    end do
  end do
         
endif


 
do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 powi = (2.+bv_ice(k))/b_ice
 do j=1,ny
  do i=1,nx	  
	  
!-------     Autoconversion/accretion/riming 

   if(qc(i,j,k)+qp(i,j,k)+qi(i,j,k).gt.0.) then


         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))

	 if(qc(i,j,k).gt.0.) then
     
           qcc = qc(i,j,k) 

      ! Autoconversion of cloud water to rain
!           autor = 1350.*qcc**1.47/Nc0**1.79   ! Linearized drizzle autoconversion
                                                !(Khairoutdinov and Kogan 2000)
          if(qcc .gt. qcw0) then
            autor = alphaelq
          else
           autor = 0.
          endif 
      ! Accretion of cloud water by rain
           accrr = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             accrr = accrrc(k) * qrr ** powr1
           end if
      ! Riming of cloud water by ice/snow
           accrci = accric(k)*qi(i,j,k)**powi 
           accrcg = 0.
      ! Riming of cloud water by graupel
           if(omp.lt.0.999) then
             qgg = qp(i,j,k) * (1.-omp)
             tmp = qgg ** powg1
             accrcg = accrgc(k) * tmp
           endif

    !       qcc = qcc/(1.+dtn*(accrr+accrci+accrcg+autor))
    !       dqp = dtn * (accrr*qcc + autor*qcc + accrcg*qcc)
           qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrci+accrcg+autor))
           dqp = dtn * (accrr*qcc + autor*(qcc-qcw0)+ accrcg*qcc)
           dqi = dtn * accrci*qcc 
           dqi= min(dqi,qc(i,j,k))
           dq= min(dqp+dqi,qc(i,j,k))
           dqp = dq - dqi
           qp(i,j,k) = qp(i,j,k) + dqp
           q(i,j,k) = q(i,j,k) - dq
           qc(i,j,k) = qc(i,j,k) - dq
           qi(i,j,k) = qi(i,j,k) + dqi
	   qpsrc(k) = qpsrc(k) + dq

         end if
         if(qp(i,j,k).gt.qp_threshold) then

           dq = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             dq = dq + dtn*(evapr1(k)*sqrt(qrr)+evapr2(k)*qrr**powr2) &
                       *min(0.,qv(i,j,k)/qsatw(tabs(i,j,k),pres(k))-1.)
           end if
           if(omp.lt.0.999) then
             qgg = qp(i,j,k) * (1.-omp)
             dq = dq + dtn*(evapg1(k)*sqrt(qgg) + evapg2(k)*qgg**powg2) &
                       *min(0.,qv(i,j,k)/qsati(tabs(i,j,k),pres(k))-1.)
           end if
           dq = max(-qp(i,j,k),dq) 
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
	   qpevp(k) = qpevp(k) + dq

	 else
	
           q(i,j,k) = q(i,j,k) + qp(i,j,k)
	   qpevp(k) = qpevp(k) - qp(i,j,k)
           qp(i,j,k) = 0.

         endif

    endif

    dq = qp(i,j,k)
    qp(i,j,k)=max(0.,qp(i,j,k))
    q(i,j,k) = q(i,j,k) + (dq-qp(i,j,k))
    dq = qi(i,j,k)
    qi(i,j,k)=max(0.,qi(i,j,k))
    q(i,j,k) = q(i,j,k) + (dq-qi(i,j,k))

  end do
 enddo
enddo
    


if(dostatis) then
                  
  call stat_varscalar(q,df,f0,df0,q2leprec)
  call setvalue(qwleprec,nzm,0.)
  call stat_sw2(q,df,qwleprec)

endif


end subroutine precip_proc

