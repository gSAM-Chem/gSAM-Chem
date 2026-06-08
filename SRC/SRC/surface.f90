subroutine surface()

use vars
use mpi_stuff
use params
use microphysics, only: micro_field, index_water_vapor
use terrain
use simple_ocean, only: sst_evolve
use rad, only: swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy
use slm_vars, only: seaicemask
use sgs, only: tke
use cup, only: wd_cu
implicit none

real ts, qs, ws, z0s, th, qh, tah, uh, vh, uh0
real taux, tauy, xlmo
real diag_ustar, trf,qrf,urf,vrf,ra,rs,raq,wd,prate,coef
integer i,j,k,iqv
real tmp, tmp1, tmp2, tmpu(0:nx,ny), tmpv(nx,0:ny)
real zr(nx,ny),presr(nx,ny),ur(nx,ny),vr(nx,ny),tr(nx,ny),qvr(nx,ny),t_sfc(nx,ny)
real rhosf(nx,ny),pressf(nx,ny)
real qsat

integer tag(2),req(2),irank,count
real(8) buffer(2), buffer1(2)
real buffers(max(nx,ny),2)
logical flag(2)
real, external :: qsatw,qsati
integer niter, m


call t_startf ('surface')

shf_ocean = 0.
lhf_ocean = 0.
shf_land = 0.
lhf_land = 0.
evp_all = 0.
shf_all = 0.
lhf_all = 0.

iqv = index_water_vapor

if(.not.SFC_FLX_FXD) then

  if(OCEAN.or.ISLAND) then

     if(dossthomozonal) then
      niter = 3
     else
      niter = 1
     end if 
     do m=1,niter

       do j=1,ny
         do i=1,nx

          if(landmask(i,j).eq.0.and.seaicemask(i,j).eq.0) then

           k=k_terra(i,j)
           uh = (u(i+1,j,k)*terrau(i+1,j,k)+u(i,j,k)*terrau(i,j,k)) &
                /(terrau(i+1,j,k)+terrau(i,j,k)+1.e-10)+ug 
           vh = (v(i,j+YES3D,k)*terrav(i,j+YES3D,k)+v(i,j,k)*terrav(i,j,k)) &
                /(terrav(i,j+YES3D,k)+terrav(i,j,k)+1.e-10)+vg
           gustsfc_xy(i,j) = max(gustsfc_xy(i,j), uh**2+vh**2)
           tah = tabs(i,j,k)
           qh = qv(i,j,k)
           th = tabs(i,j,k)*prespot(k)
           ts = (sstxy(i,j)+t00)*prespoti(k)
           if(doseaice.and.sstxy(i,j)+t00.lt.271._8) then
              qsat = qsati(real(sstxy(i,j))+t00,presi(k))
              qs = qsat
           else
              qsat = qsatw(real(sstxy(i,j))+t00,presi(k))
              qs = salt_factor*qsat
           end if
           if(docup) then
             wd = wd_cu(i,j)
           else
             ! adding w* to surface wind 
             ! as well as gustiness due to precipitation following parameterization by
             ! Redelsperger at al (2000, J. Clim 13, 402-421)
          !   prate = precinst(i,j)*8640. ! convert precip rate from mm/s to cm/d
          !   wd = sqrt((9800./tah*max(0.,(fluxbt(i,j)+0.61*tah*fluxbq(i,j))))**0.333 &
          !       + min(3.2,log(1.+6.69*prate-0.474*prate**2))**2)
          !  From IFS documentation:
          !   wd = (1000.*ggr/tah*max(0.,fluxbt(i,j)+epsv*tah*fluxbq(i,j)))**0.3333
             wd = 0. 
           end if

           call oceflx(rhow(k), uh, vh, tah, qh, th, z(k)-zi(k), ts, qs, wd, &
                        fluxt0, fluxq0, taux, tauy, ra, raq, trf, qrf, urf, vrf)
           fluxbu(i,j) = taux/rhow(k) 
           fluxbv(i,j) = tauy/rhow(k) 
           fluxbt(i,j) = fluxt0 
           fluxbq(i,j) = fluxq0 
           trf = trf/prespoti(k)
           if(m.eq.niter) then
            t2m_xy(i,j) = t2m_xy(i,j) + trf * dtfactor
            q2m_xy(i,j) = q2m_xy(i,j) + qrf * dtfactor
            u10m_xy(i,j) = u10m_xy(i,j) + urf * dtfactor
            v10m_xy(i,j) = v10m_xy(i,j) + vrf * dtfactor
            u10ma_xy(i,j) = u10ma_xy(i,j) + urf * dtfactor
            v10ma_xy(i,j) = v10ma_xy(i,j) + vrf * dtfactor
            ra_xy(i,j) = ra_xy(i,j) + ra * dtfactor
            gust10m_xy(i,j) = max(gust10m_xy(i,j),urf**2+vrf**2)
           end if
           t2m(i,j) = trf
           q2m(i,j) = qrf
           u10m(i,j) = urf
           v10m(i,j) = vrf
           raf(i,j) = ra
           rafq(i,j) = raq
           shf_ocean(i,j) = cp*rhow(k)*fluxbt(i,j)
           lhf_ocean(i,j) = lcond*rhow(k)*fluxbq(i,j)
           shf_all(i,j) = shf_ocean(i,j)
           lhf_all(i,j) = lhf_ocean(i,j)
           evp_all(i,j) = lhf_ocean(i,j)/lcond

!           if(isnan(fluxbu(i,j)).or.isnan(fluxbv(i,j)).or.isnan(fluxbt(i,j)).or.isnan(fluxbq(i,j))) then
!             print*,'oceflx:',rank,i,j,k,'fluxes:',fluxbu(i,j),fluxbv(i,j),fluxbt(i,j),fluxbq(i,j), &
!             uh,vh,tah,qh,th,ts,qs,sstxy(i,j),'pres:',pres(k),presi(k),prespot(k),prespoti(k), &
!            'z:',z(k),zi(k),rho(k),rhow(k)
!             call task_abort()
!           end if
          end if ! landmask.eq.0

         end do
       end do

       call sst_evolve()

     end do ! m

     do j=1,ny
        do i=1,nx
          if(landmask(i,j).eq.0.and.seaicemask(i,j).eq.0) then
            sst_xy(i,j)=sst_xy(i,j)+sstxy(i,j)*dtfactor   ! 2D output
          end if
       end do
     end do

  end if ! OCEAN .or. ISLAND

  if((LAND.or.ISLAND).and.SLM) then

   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      rhosf(i,j) = rhow(k)
      presr(i,j) = pres(k)
      pressf(i,j) = presi(k)
! make the reference level to be the same as the first grid level
! interpolate all variables to that level
      zr(i,j) = min(z(1)-zi(1),z(k)-zi(k))
      coef = (zr(i,j) - (z(k)-zi(k))) / (z(k+1)-z(k))
      ur(i,j) = (u(i+1,j,k)*terrau(i+1,j,k)+u(i,j,k)*terrau(i,j,k)) &
                          /(terrau(i+1,j,k)+terrau(i,j,k)+1.e-10)+ug
      tmp = (u(i+1,j,k+1)*terrau(i+1,j,k+1)+u(i,j,k+1)*terrau(i,j,k+1)) &
                          /(terrau(i+1,j,k+1)+terrau(i,j,k+1)+1.e-10)+ug
      ur(i,j) = ur(i,j) + (tmp - ur(i,j))*coef 
      vr(i,j) = (v(i,j+YES3D,k)*terrav(i,j+YES3D,k)+v(i,j,k)*terrav(i,j,k)) &
                          /(terrav(i,j+YES3D,k)+terrav(i,j,k)+1.e-10)+vg
      tmp = (v(i,j+YES3D,k+1)*terrav(i,j+YES3D,k+1)+v(i,j,k)*terrav(i,j,k+1)) &
                          /(terrav(i,j+YES3D,k+1)+terrav(i,j,k+1)+1.e-10)+vg
      vr(i,j) = vr(i,j) + (tmp - vr(i,j))*coef 
      tr(i,j) = tabs(i,j,k) + (tabs(i,j,k+1) - tabs(i,j,k))*coef 
      qvr(i,j) = qv(i,j,k) + (qv(i,j,k+1) - qv(i,j,k))*coef 
      t_sfc(i,j) = sstxy(i,j)+t00
      gustsfc_xy(i,j) = max(gustsfc_xy(i,j), ur(i,j)**2+vr(i,j)**2)
 

!      if(abs(i-nx/2.-0.5).lt.1..and.abs(j-ny/2.-0.5).lt.1.and.nstep.lt.8640) then
!           precinst(i,j) = 200./3600.
!      else
!           precinst(i,j) = 0.
!      end if

    end do
   end do
 


   call run_slm(ur, vr, precinst, qvr, tr, t_sfc, zr, presr, pressf, rhosf, &
              swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy, &
              precinstsoil(:,:), fluxbu, fluxbv, fluxbq, fluxbt, lhf_land, raf)
   sstxy = t_sfc-t00
   do j=1,ny
       do i=1,nx
         if(landmask(i,j).eq.1.or.seaicemask(i,j).eq.1) then
           sst_xy(i,j)=sst_xy(i,j)+sstxy(i,j)*dtfactor   ! 2D output
           ra_xy(i,j) = ra_xy(i,j) + raf(i,j)*dtfactor
           shf_land(i,j) = cp*rhosf(i,j)*fluxbt(i,j)
           shf_all(i,j) = shf_land(i,j)
           lhf_all(i,j) = lhf_land(i,j)
           evp_all(i,j) = fluxbq(i,j)*rhosf(i,j)
         end if
       end do
   end do

  end if ! SLM

! 
! compute stresses at velocity positions:

  tmpu(1:nx,1:ny) = fluxbu(1:nx,1:ny)
  tmpv(1:nx,1:ny) = fluxbv(1:nx,1:ny)
  if(dompi) then
     call task_receive_float(buffers(:,1),max(nx,ny),req(1))
     call task_receive_float(buffers(:,2),max(nx,ny),req(2))
     call task_bsend_float(rankee,fluxbu(nx,1:ny),ny,133)
     call task_bsend_float(ranknn,fluxbv(1:nx,ny),nx,134)
     count = 2  
     flag(1:2) = .false.
     do while(count.gt.0)
      do k=1,2
        if(.not.flag(k)) then
          call task_test(req(k),flag(k),irank,tag(k))
          if(flag(k)) then
            if(tag(k).eq.133) then
             tmpu(0,1:ny) = buffers(1:ny,k)  
             count=count-1
            else if(tag(k).eq.134) then
             tmpv(1:nx,0) = buffers(1:nx,k)  
             count=count-1
            else
              print*,'surface: wrong tag. Should be 133 or 134, got ',tag
              call task_abort
            end if
          end if
        end if 
      end do
     end do
  else
     tmpu(0,1:ny) = fluxbu(nx,1:ny)
     tmpv(1:nx,0) = fluxbv(1:nx,ny)
  end if
  do j=1,ny
   do i=1,nx
    k = k_terrau(i,j)
    fluxbu(i,j) = (tmpu(i,j)*terra(i,j,k)+tmpu(i-1,j)*terra(i-1,j,k)) &
                  /(terra(i,j,k)+terra(i-1,j,k)+1.e-10)
   end do
  end do

  do j=1,ny
   do i=1,nx
    k = k_terrav(i,j)
    fluxbv(i,j) = (tmpv(i,j)*terra(i,j,k)+tmpv(i,j-YES3D)*terra(i,j-YES3D,k)) &
                /(terra(i,j,k)+terra(i,j-YES3D,k)+1.e-10)
   end do
  end do

  call task_barrier()

end if! .not.SFC_FLX_FXD

if(SFC_FLX_FXD) then

  uh0 = max(0.1,sqrt((u0(1)+ug)**2+(v0(1)+vg)**2))

  if(.not.SFC_TAU_FXD) then

    if(OCEAN) z0 = 0.0001  ! for LAND z0 should be set in namelist (default z0=0.035)

    tau0 = diag_ustar(z(1),  &
                bet(1)*(fluxt0+epsv*(t0(1)-gamaz(1))*fluxq0),uh0,z0)**2  

  end if ! .not.SFC_TAU_FXD

  fluxbu(:,:) = -(u(1:nx,1:ny,1)+ug)/uh0*tau0
  fluxbv(:,:) = -(v(1:nx,1:ny,1)+vg)/uh0*tau0

  fluxbt(:,:) = fluxt0
  fluxbq(:,:) = fluxq0

  shf_all(:,:) = cp*rhow(1)*fluxbt(:,:)
  lhf_all(:,:) = lcond*rhow(1)*fluxbq(:,:)

end if ! SFC_FLX_FXD


!
! Homogenize the surface scalar fluxes if needed for sensitivity studies
!
   if(dosfchomo) then

	fluxt0 = 0.
	fluxq0 = 0.
	do j=1,ny
         do i=1,nx
	   fluxt0 = fluxt0 + fluxbt(i,j)
	   fluxq0 = fluxq0 + fluxbq(i,j)
         end do
        end do
	fluxt0 = fluxt0 / float(nx*ny)
	fluxq0 = fluxq0 / float(nx*ny)
        if(dompi) then
            buffer(1) = fluxt0
            buffer(2) = fluxq0
            call task_sum_real8(buffer,buffer1,2)
	    fluxt0 = buffer1(1) /float(nsubdomains)
	    fluxq0 = buffer1(2) /float(nsubdomains)
        end if ! dompi
	fluxbt(:,:) = fluxt0
	fluxbq(:,:) = fluxq0

   end if

shf_xy(:,:) = shf_xy(:,:) + fluxbt(:,:) * dtfactor
lhf_xy(:,:) = lhf_xy(:,:) + fluxbq(:,:) * dtfactor
taux_xy(:,:) = taux_xy(:,:) + fluxbu(:,:) * dtfactor
tauy_xy(:,:) = tauy_xy(:,:) + fluxbv(:,:) * dtfactor


!-------------------------------------------------------------------------------

do j=1,ny
  do i=1,nx
    sst_min_xy(i,j)=min(sst_min_xy(i,j),real(sstxy(i,j))+t00)
    sst_max_xy(i,j)=max(sst_max_xy(i,j),real(sstxy(i,j))+t00)
  end do
end do

call t_stopf ('surface')

end




! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx .ne. 0.0) then 
        do iterate=1,4
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = z*rlmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------

