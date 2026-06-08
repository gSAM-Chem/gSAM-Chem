  
subroutine precip_proc

use vars
use microphysics
use micro_params
use params
use terrain, only : k_terra, terra
use stat_coars

implicit none

integer i,j,k
real autor, autos, accrr, accris, accrcs, accrig, accrcg
real dq, omn, omp, omg
real pows1, pows2, powg1, powg2, powr1, powr2, tmp
real coefice
real qii, qcc, qrr, qss, qgg, q1, qp1
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
real Ncc
real, external :: qsatw,qsati
integer iflag

if(donomicro) return

call t_startf('precip_proc')

powr1 = (3 + b_rain) / 4.
powr2 = (5 + b_rain) / 8.
pows1 = (3 + b_snow) / 4.
pows2 = (5 + b_snow) / 8.
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

if(collect_coars) call coars_fld(qp(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm), &
                                         fld_tend(:,:,:,1))
 
autos = 0.
autor = 0.
accrr = 0.
accrcs = 0.
accris = 0. 
accrcg = 0.
accrig = 0. 

do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 do j=1,ny
  do i=1,nx	  
	  
        if(q(i,j,k).lt.0..or.qn(i,j,k).lt.0..or.qn(i,j,k).gt.q(i,j,k).or.qp(i,j,k).lt.0.) then
           print*,'precip_proc.f90: >>>>! rank,i,j,k,q,qn,qp,tabs,qsat,pres:',&
           rank,i,j,k,q(i,j,k),qn(i,j,k),qp(i,j,k),tabs(i,j,k),qsatt(i,j,k),pres(k)
           call task_abort()
         end if

   iflag = 0 

   if(k.lt.k_terra(i,j)) then
  ! don't just throw away precip inside terrain for conservation sake.
    q(i,j,k) = q(i,j,k) + qp(i,j,k)
    qp(i,j,k) = 0.
    cycle
   end if
!-------     Autoconversion/accretion 

   if(qn(i,j,k)+qp(i,j,k).gt.0.) then

         iflag = 1

         omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
         omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))

	 if(qn(i,j,k).gt.0.) then
     
           qcc = qn(i,j,k) * omn
           qii = qn(i,j,k) * (1.-omn)

         if(doKKauto) then

! Khairoutdinov and Kogan (2000) parameterization for autoconversion:
! Note that the power 1.47 instead of 2.47 is because of the implicit scheme 
! (autor is multiplied later by qcc)

           qcw0 = 0.
           Ncc = landmask(i,j)*Nc_land+(1.-landmask(i,j))*NC_ocn ! from IFS model
           autor = 1350.*qcc**1.47/Ncc**1.79
           if(do_scale_dependence_of_autoconv) then
              autor = autor*min(1.,10000./(dx*mu(j)*dy*ady(j)))
           else
            ! fudge factor similar to one used in IFS Cy47r3 doc)
              autor = autor * auto_fudge
           end if
         else

! Standard SAM Kessler parameterization for autovonversion:
           if(qcc .gt. qcw0) then
            autor = alphaelq
            if(do_scale_dependence_of_autoconv) then
              autor = autor*min(1.,10000./(dx*mu(j)*dy*ady(j)))
            end if
           else
            autor = 0.
           endif 

         end if ! doKKauto

         if(doKKaccr) then

! Khairoutdinov and Kogan (2000) accretion of cloud water by rain

           accrr = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             accrr = 67. * qrr**1.15 * qcc**0.15
           end if

          else

! Standard SAM Accretion of cloud water by rain:

           accrr = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             accrr = accrrc(k) * qrr ** powr1
           end if

          end if ! doKKaccr

!-----------------------------------
          if(.not.dowarmcloud) then

! Parameterization of ice-to-snow conversion (based on Lin et al (1983))

           if(tabs(i,j,k).lt.273.15) then
             coefice = exp(0.025*(tabs(i,j,k) - 273.15))
           else
             coefice = 0.
           end if

           if(qii .gt. qci0) then
            autos = betaelq*coefice
           else
            autos = 0.
           endif 

! Accretion of cloud water/ice by snow/graupel:

           accrcs = 0.
           accris = 0. 
           if(omp.lt.0.999.and.omg.lt.0.999) then
             qss = qp(i,j,k) * (1.-omp)*(1.-omg)
             tmp = qss ** pows1
             accrcs = accrsc(k) * tmp
             accris = accrsi(k) * coefice * tmp 
           end if
           accrcg = 0.
           accrig = 0. 
           if(omp.lt.0.999.and.omg.gt.0.001) then
             qgg = qp(i,j,k) * (1.-omp)*omg
             tmp = qgg ** powg1
             accrcg = accrgc(k) * tmp
             accrig = accrgi(k) * coefice * tmp 
           endif

          end if ! .not.dowarmcloud
!--------------------------------------------
! Implementation of the implicit time-scheme for all processes:

           qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
           qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
           dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
             (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
           dq = min(dq,qn(i,j,k))
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
           qn(i,j,k) = qn(i,j,k) - dq
	   qpsrc(k) = qpsrc(k) + dq

         elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then

           iflag = 2
! Evaporation of hydrometeors outside clouds:

           dq = 0.
           q1 = q(i,j,k)
           qp1 = qp(i,j,k)
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             dq = dq + (evapr1(i,j,k)*sqrt(qrr)+evapr2(i,j,k)*qrr**powr2) &
                                       *(q(i,j,k)/qsatt(i,j,k)-1.) 
           end if
           if(omp.lt.0.999) then
             if(omg.lt.0.999) then
               qss = qp(i,j,k) * (1.-omp)*(1.-omg)
               dq = dq + (evaps1(i,j,k)*sqrt(qss)+evaps2(i,j,k)*qss**pows2) &
                                       *(q(i,j,k)/qsatt(i,j,k)-1.) 
             end if
             if(omg.gt.0.001) then
               qgg = qp(i,j,k) * (1.-omp)*omg
               dq = dq + (evapg1(i,j,k)*sqrt(qgg)+evapg2(i,j,k)*qgg**powg2) &
                                       *(q(i,j,k)/qsatt(i,j,k)-1.)
             end if
           end if
           dq = dq * dtn  
           dq = max(-0.5*qp(i,j,k),dq) 
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

    if(q(i,j,k)-qn(i,j,k).lt.0.) then
      print*,'precip_proc.f90: q(i,j,k)-qn(i,j,k).lt.0.! iflag,rank,i,j,k,landmask,q,q1,qn,qsat,dq,tabs,qp,qp1,pres:',&
           iflag,rank,i,j,k,landmask(i,j),q(i,j,k),q1,qn(i,j,k),qsatt(i,j,k),dq,tabs(i,j,k),qp(i,j,k),qp1,pres(k)
      call task_abort()
    end if

  end do
 enddo
enddo
    


if(dostatis) then
                  
  call stat_varscalar(q,df,f0,df0,q2leprec)
  call setvalue(qwleprec,nzm,0.)
  call stat_sw2(q,df,qwleprec)

endif

call t_stopf('precip_proc')

end subroutine precip_proc

