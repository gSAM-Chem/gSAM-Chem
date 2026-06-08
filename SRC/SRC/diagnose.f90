subroutine diagnose
	
! Diagnose some useful stuff

use vars
use params
use sgs, only: sgs_diagnose
use terrain
use microphysics, only: nmicro_fields,micro_field,index_water_vapor,q 
use tracers
implicit none
	
integer i,j,k,kb,kc
real(8) coefd,buffer(nzm,7), buffer1(nzm,7)
real coef,coef1,coef2,coef3
real(8) u0d(nzm),v0d(nzm),tabs0d(nzm),q0d(nzm),qn0d(nzm),qp0d(nzm),t0d(nzm)
real tmp_lwp(nx,ny), coefsum, tmp, fld
real mae700, mae950, maesfc
real, external :: qsatw 
integer, external :: index_pressure

call t_startf ('diagnose')

t01(1:nzm) = tabs0(1:nzm)
q01(1:nzm) = q0(1:nzm)

do k=1,nzm
  u0d(k)=0.
  v0d(k)=0.
  t0d(k)=0.
  tabs0d(k)=0.
  q0d(k)=0.
  qn0d(k)=0.
  qp0d(k)=0.
  do j=1,ny
    do i=1,nx
     tabs(i,j,k) = t(i,j,k)-gamaz(k)+ fac_cond * (qcl(i,j,k)+qpl(i,j,k)) + &
                                      fac_sub * (qci(i,j,k) + qpi(i,j,k))
     u0d(k) = u0d(k)+u(i,j,k)*wgtu(j,k)*terrau(i,j,k)
     v0d(k) = v0d(k)+v(i,j,k)*wgtv(j,k)*terrav(i,j,k)
     coefd = wgt(j,k)*terra(i,j,k)
     t0d(k)= t0d(k)+t(i,j,k)*coefd
     tabs0d(k) = tabs0d(k)+tabs(i,j,k)*coefd
     q0d(k) = q0d(k)+(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))*coefd
     qn0d(k) = qn0d(k) + (qcl(i,j,k) + qci(i,j,k))*coefd
     qp0d(k) = qp0d(k) + (qpl(i,j,k) + qpi(i,j,k))*coefd
    end do
  end do
end do

if(dompi) then
  coefd = 1.d0/dble(nx*ny)
  do k=1,nzm
    buffer(k,1) = u0d(k)*coefd
    buffer(k,2) = v0d(k)*coefd
    buffer(k,3) = t0d(k)*coefd
    buffer(k,4) = q0d(k)*coefd
    buffer(k,5) = tabs0d(k)*coefd
    buffer(k,6) = qn0d(k)*coefd
    buffer(k,7) = qp0d(k)*coefd
  end do
  call task_sum_real8(buffer,buffer1,nzm*7)
  coefd = 1.d0/dble(nsubdomains)
  do k=1,nzm
    u0(k)=buffer1(k,1)*coefd
    v0(k)=buffer1(k,2)*coefd
    t0(k)=buffer1(k,3)*coefd
    q0(k)=buffer1(k,4)*coefd
    tabs0(k)=buffer1(k,5)*coefd
    qn0(k)=buffer1(k,6)*coefd
    qp0(k)=buffer1(k,7)*coefd
  end do
else
  coefd = 1.d0/dble(nx*ny)
  do k=1,nzm
    u0(k)=u0d(k)*coefd
    v0(k)=v0d(k)*coefd
    t0(k)=t0d(k)*coefd
    q0(k)=q0d(k)*coefd
    tabs0(k)=tabs0d(k)*coefd
    qn0(k)=qn0d(k)*coefd
    qp0(k)=qp0d(k)*coefd
  end do
end if

qv0(:) = q0(:) - qn0(:)

if(.not.docloud) qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,1)

k200 = index_pressure(200.)
k500 = index_pressure(500.)
k700 = index_pressure(700.)
k850 = index_pressure(850.)
k950 = index_pressure(950.)

if(nstep.eq.0) goto 111

if(dtfactor.gt.0.) then

 tmp_lwp = 0.
 coefsum = 0.
 do k=1,nzm
   coef1 = rho(k)*dz*adz(k)
   coefsum = coefsum + coef1
   coef3 = coef1*dtfactor
   do j=1,ny
     do i=1,nx
      coef2 = coef3*terra(i,j,k)
      pw_xy(i,j) = pw_xy(i,j)+qv(i,j,k)*coef2
      ta_xy(i,j) = ta_xy(i,j)+tabs(i,j,k)*coef2
      cw_xy(i,j) = cw_xy(i,j)+qcl(i,j,k)*coef2
      iw_xy(i,j) = iw_xy(i,j)+qci(i,j,k)*coef2
      fse_xy(i,j) = fse_xy(i,j)+ &
          (tabs(i,j,k)+gamaz(k)+fac_cond*qv(i,j,k)-fac_fus*qci(i,j,k))*coef2
      tmp_lwp(i,j) = tmp_lwp(i,j) + (qcl(i,j,k)+qci(i,j,k))*coef2
      pws_xy(i,j) = pws_xy(i,j)+qsatw(tabs(i,j,k),pres(k))*coef2
     end do
   end do
 end do ! k
 if(nstep.eq.0) coefsum=1. ! let first call from setdata to compute mean profile
 coefsum=1./coefsum
 ta_xy(:,:) = ta_xy(:,:)*coefsum
 tmp_lwp(:,:) = tmp_lwp(:,:)/dtfactor
 if(LES_S) then ! make consistent with the way statistics is gathered in *.stat
  where(tmp_lwp(:,:).gt.0.) cld_xy(:,:) = cld_xy(:,:) + dtfactor
 else
  where(tmp_lwp(:,:).gt.cwp_threshold) cld_xy(:,:) = cld_xy(:,:) + dtfactor
 end if


 do k=1,nzm
   coef1 = rho(k)*dz*adz(k)*dtfactor
   do j=1,ny
     do i=1,nx
      coef2 = coef1*terra(i,j,k)
      rw_xy(i,j) = rw_xy(i,j)+qpl(i,j,k)*coef2
      gw_xy(i,j) = gw_xy(i,j)+qpg(i,j,k)*coef2
      sw_xy(i,j) = sw_xy(i,j)+(qpi(i,j,k)-qpg(i,j,k))*coef2
     end do
   end do
 end do ! k

 coef = dtfactor
 coef1 = 1./dx
 do j=1,ny
  coef2 = 1./(dy*adyv(j))
  do i=1,nx
   k = k_terra(i,j)
   usfc_xy(i,j) = usfc_xy(i,j) + u(i,j,k_terrau(i,j))*coef
   vsfc_xy(i,j) = vsfc_xy(i,j) + v(i,j,k_terrav(i,j))*coef
   tsfc_xy(i,j) = tsfc_xy(i,j) + tabs(i,j,k)*coef
   qsfc_xy(i,j) = qsfc_xy(i,j) + qv(i,j,k)*coef
   u200_xy(i,j) = u200_xy(i,j) + u(i,j,k200)*coef *terrau(i,j,k200) 
   v200_xy(i,j) = v200_xy(i,j) + v(i,j,k200)*coef *terrav(i,j,k200) 
   w500_xy(i,j) = w500_xy(i,j) + w(i,j,k500)*coef *terraw(i,j,k500)
   phi500_xy(i,j) = phi500_xy(i,j) + p(i,j,k500,nb)*coef *terra(i,j,k500)
   vor200_xy(i,j) = vor200_xy(i,j) + coef*imuv(j)*((v(i,j,k200)-v(i-1,j,k200))*coef1 &
                            -(mu(j)*u(i,j,k200)-mu(j-YES3D)*u(i,j-YES3D,k200))*coef2)
  end do
 end do

 if(dotracers) then
  do k=1,nzm
    coef1 = rho(k)*dz*adz(k)*dtfactor
    do j=1,ny
      do i=1,nx
       tr_xy(i,j,1:ntracers) = tr_xy(i,j,1:ntracers)+ &
             tracer(i,j,k,1:ntracers)*terra(i,j,k)*coef1
      end do
    end do
  end do
  do j=1,ny
   do i=1,nx
    k = k_terra(i,j)
    tr_acs_xy(i,j,1:ntracers) = tr_acs_xy(i,j,1:ntracers) + &
                        tracer(i,j,k_terrau(i,j),1:ntracers)*rho(k_terrau(i,j))*dtn
   end do
  end do
 end if

 ! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
 coef = dtfactor
 do j=1,ny
  do i=1,nx
   ! 850 mbar horizontal winds
   u850_xy(i,j) = u850_xy(i,j) + u(i,j,k850)*coef*terrau(i,j,k850)  
   v850_xy(i,j) = v850_xy(i,j) + v(i,j,k850)*coef*terrav(i,j,k850)  
   omega200_xy(i,j) = omega200_xy(i,j) - w(i,j,k200)*ggr*rhow(k200)*coef *terraw(i,j,k200)
   omega500_xy(i,j) = omega500_xy(i,j) - w(i,j,k500)*ggr*rhow(k500)*coef *terraw(i,j,k500)
   omega700_xy(i,j) = omega700_xy(i,j) - w(i,j,k700)*ggr*rhow(k700)*coef *terraw(i,j,k700)
   omega850_xy(i,j) = omega850_xy(i,j) - w(i,j,k850)*ggr*rhow(k850)*coef *terraw(i,j,k850)
   rh200_xy(i,j) = rh200_xy(i,j) + qv(i,j,k200)/qsatw(tabs(i,j,k200),pres(k200))*coef*terra(i,j,k200)
   rh500_xy(i,j) = rh500_xy(i,j) + qv(i,j,k500)/qsatw(tabs(i,j,k500),pres(k500))*coef*terra(i,j,k500)
   rh700_xy(i,j) = rh700_xy(i,j) + qv(i,j,k700)/qsatw(tabs(i,j,k700),pres(k700))*coef*terra(i,j,k700)
   rh850_xy(i,j) = rh850_xy(i,j) + qv(i,j,k850)/qsatw(tabs(i,j,k850),pres(k850))*coef*terra(i,j,k850)
   vor850_xy(i,j) = vor850_xy(i,j) + coef*imuv(j)*((v(i,j,k850)-v(i-1,j,k850))/dx &
                                     -(mu(j)*u(i,j,k850)-mu(j-YES3D)*u(i,j-YES3D,k850))/(dy*adyv(j)))
  end do
 end do

 do k=1,nzm
   coef1 = rho(k)*dz*adz(k)*dtfactor
   do j=1,ny
     do i=1,nx
      ! Saturated water vapor path with respect to water. Can be used
      ! with water vapor path (= pw) to compute column-average
      ! relative humidity.   
      swvp_xy(i,j) = swvp_xy(i,j)+qsatw(tabs(i,j,k),pres(k))*coef1*terra(i,j,k)
     end do
   end do
 end do ! k

 if(donudge3D.or.doregion) then
   if(allocated(uobs)) then
    if(dayfld3Dobs(1).gt.0.) then
      coef=(day-dayfld3Dobs(1))/(dayfld3Dobs(2)-dayfld3Dobs(1))
    else
      coef = 0.
    end if
    do k=1,nzm
     coef1 = rho(k)*dz*adz(k)*dtfactor
     do j=1,ny
       do i=1,nx
        coef3 = coef1*terra(i,j,k)
        fld = qobs(i,j,k,1)+(qobs(i,j,k,2)-qobs(i,j,k,1))*coef
        pwobs_xy(i,j) = pwobs_xy(i,j)+fld*coef3
       end do
     end do
    end do
    do j=1,ny
     do i=1,nx
      fld = uobs(i,j,k850,1)+(uobs(i,j,k850,2)-uobs(i,j,k850,1))*coef
      uobs850_xy(i,j) = uobs850_xy(i,j) + fld*dtfactor*terrau(i,j,k850)
      fld = vobs(i,j,k850,1)+(vobs(i,j,k850,2)-vobs(i,j,k850,1))*coef
      vobs850_xy(i,j) = vobs850_xy(i,j) + fld*dtfactor*terrav(i,j,k850)
      fld = uobs(i,j,k200,1)+(uobs(i,j,k200,2)-uobs(i,j,k200,1))*coef
      uobs200_xy(i,j) = uobs200_xy(i,j) + fld*dtfactor *terrau(i,j,k200) 
      fld = vobs(i,j,k200,1)+(vobs(i,j,k200,2)-vobs(i,j,k200,1))*coef
      vobs200_xy(i,j) = vobs200_xy(i,j) + fld*dtfactor *terrav(i,j,k200) 
     end do
    end do
   end if
 end if

! Compute the Estimated Inversion Strenth (EIS). Based on IFS doc for Cy47r3
! Here it is used only over ocean.

 do j=1,ny
  do i=1,nx
   k = k_terra(i,j)
   if(landmask(i,j).eq.0.and.k.lt.k950) then
     mae700 = tabs(i,j,k700)*(1.+5.87*(qv(i,j,k700)+qcl(i,j,k700)+qci(i,j,k700))) &
         - fac_cond*qcl(i,j,k700) - fac_sub*qci(i,j,k700) + gamaz(k700)
     mae950 = tabs(i,j,k950)*(1.+5.87*(qv(i,j,k950)+qcl(i,j,k950)+qci(i,j,k950))) &
         - fac_cond*qcl(i,j,k950) - fac_sub*qci(i,j,k950) + gamaz(k950)
     maesfc = tabs(i,j,k)*(1.+5.87*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))) &
         - fac_cond*qcl(i,j,k) - fac_sub*qci(i,j,k) + gamaz(k)
     eis(i,j) = max(mae700-mae950,mae950-maesfc)
   else 
    eis(i,j) = 0. 
   end if
   eis_xy(i,j) = eis_xy(i,j) + eis(i,j)*dtfactor
  end do
 end do

 if(mod(nstep,nsave2D).eq.0.and.nstep.ge.nsave2Dstart.and.nstep.le.nsave2Dend) then

! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.

! initially, zero out heights and set cloudtoptemp to SST
  cloudtopheight = 0.
  cloudtoptemp = sstxy(1:nx,1:ny) + t00
  tmp_lwp = 0.
  cloudcover = 0.
  do k = nzm,1,-1
    coef = rho(k)*dz*adz(k)
    do j = 1,ny
      do i = 1,nx
        if(k.lt.k_terra(i,j)) cycle
        ! FIND CLOUD TOP HEIGHT
        if(tmp_lwp(i,j).lt.0.01) then
             tmp_lwp(i,j) = tmp_lwp(i,j) + (qcl(i,j,k)+qci(i,j,k))*coef
              if (tmp_lwp(i,j).ge.0.01) then
              cloudtopheight(i,j) = z(k)
              cloudtoptemp(i,j) = tabs(i,j,k)*terra(i,j,k)
              tmp_lwp(i,j) = 1. ! done
              cloudcover(i,j) = 1.
              cycle
           end if
        end if
      end do
    end do
  end do
  echotopheight = 0.
  do k = nzm,1,-1
    do j = 1,ny
      do i = 1,nx
        ! FIND ECHO TOP HEIGHT
        if(k.lt.k_terra(i,j)) cycle
        if ((qpl(i,j,k)+qpi(i,j,k)).gt.1.e-6) then
             echotopheight(i,j) = z(k)
             cycle
        end if
      end do
    end do
  end do
 end if

end if ! dtfactor.gt.0.

!-----------------
! compute some sgs diagnostics:

call sgs_diagnose()

111 continue

call t_stopf ('diagnose')


end subroutine diagnose


integer function index_pressure(p)
 use vars, only: nzm, pres
 implicit none
 real, intent(in) ::  p
 integer kk,k,kc
 kk = nzm
 do k = 1,nzm
    kc=min(nzm,k+1)
    if((pres(kc).le.p).and.(pres(k).gt.p)) then
       if ((p-pres(kc)).lt.(pres(k)-p))then
          kk=kc
       else
          kk=k
       end if
    end if
 end do
 index_pressure = kk
end

