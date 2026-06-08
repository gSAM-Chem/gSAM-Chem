module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!------------------------------------------------------------

use grid
implicit none

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type

CONTAINS


SUBROUTINE set_sst()

use vars, only: sstxy,t00,landmask,latitude,longitude
use params, only: tabs_s, delta_sst, ocean_type, sst_mean, pi, &
                   readsst, earth_factor, dolatlon, shift_sst

! parameters of the sinusoidal SST destribution 
! along the X for Walker-type simulatons( ocean-type = 1):

real(8) tmpx(nx), lx, ly, yy, phi, lamda
integer i,j, it,jt,nx1,ny1

if(readsst) return

select case (ocean_type)

   case(0) ! fixed constant SST

     do j=1,ny
       do i=1,nx
         if(landmask(i,j).eq.0) sstxy(i,j) = tabs_s - t00
       end do
     end do
     sst_mean = tabs_s

   case(1) ! Sinusoidal distribution along the x-direction:

     lx = float(nx_gl)*dx
     do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
     end do
     do j=1,ny
       do i=1,nx
         if(landmask(i,j).eq.0) sstxy(i,j) = tabs_s-delta_sst*cos(2.*pi*tmpx(i)/lx) - t00
       end do
     end do
   
   case(2) ! Sinusoidal distribution along the y-direction:
     
     call task_rank_to_index(rank,it,jt)
     
     lx = y_gl(ny_gl)
     do j=1,ny
        yy = y_gl(j+jt)
       do i=1,nx
         if(landmask(i,j).eq.0) sstxy(i,j) = tabs_s+delta_sst*(2.*cos(pi*yy/lx)-1.) - t00
       end do
     end do 


   case(3) !  distribution along the y-direction:

     call task_rank_to_index(rank,it,jt)

     ly = y_gl(ny_gl)
     do j=1,ny
        yy = y_gl(j+jt)
       do i=1,nx
         if(landmask(i,j).eq.0.) sstxy(i,j) = tabs_s+delta_sst/(1+(yy/ly/0.4)**6)-t00
       end do
     end do

   case(4)

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
        yy = y_gl(j+jt)
       do i=1,nx
        if(landmask(i,j).eq.0.) then
         ! approx SST disribution in January averaged over 60E-90E
         sstxy(i,j) = tabs_s+delta_sst*cos(pi*min(1._8,abs(yy+333333.)/ly)**1.2)-t00
         tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
         sstxy(i,j) = sstxy(i,j)-0.5*exp(-(2.*yy/ly)**2)*delta_sst*cos(2.*pi*tmpx(i)/lx)
        end if
       end do
     end do

   case(5) !  Aquaplanet Experiment (APE) "Control" (Neale and Hoskins 2001)

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
       if(dolatlon) then
        phi = (latitude(1,j)-shift_sst)*pi/180.
       else
        yy = y_gl(j+jt)
        phi = pi/180.*(yy*2.5e-8*earth_factor*360.-shift_sst)
       end if
       do i=1,nx
         if(landmask(i,j).eq.0.) then
           if(abs(phi).le.pi/3.) then
               sstxy(i,j) = tabs_s+delta_sst*(1-sin(1.5*phi)**2)-t00
           else
               sstxy(i,j) = tabs_s - t00
           end if
         end if
       end do
     end do

   case(6) !  Aquaplanet Experiment (APE) "QOBS" (Neale and Hoskins 2001)

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
       if(dolatlon) then
        phi = (latitude(1,j)-shift_sst)*pi/180.
       else
        yy = y_gl(j+jt)
        phi = pi/180.*(yy*2.5e-8*earth_factor*360.-shift_sst)
       end if
       do i=1,nx
         if(landmask(i,j).eq.0.) then
           if(abs(phi).le.pi/3.) then
              sstxy(i,j) = tabs_s+delta_sst*(1-0.5*(sin(1.5*phi)**2+sin(1.5*phi)**4))-t00
           else
              sstxy(i,j) = tabs_s - t00
           end if
         end if
       end do
     end do

   case(7) !  Aquaplanet Experiment (APE) FLAT 

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
       if(dolatlon) then
        phi = (latitude(1,j)-shift_sst)*pi/180.
       else
        yy = y_gl(j+jt)
        phi = pi/180.*(yy*2.5e-8*earth_factor*360.-shift_sst)
       end if
       do i=1,nx
         if(landmask(i,j).eq.0.) then
           if(abs(phi).le.pi/3.) then
             sstxy(i,j) =  tabs_s+delta_sst*(1-sin(1.5*phi)**4)-t00
           else
             sstxy(i,j) = tabs_s - t00
           end if
         end if
       end do
     end do

   case(8,81) !  super FLAT for MJO experiments

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
       if(dolatlon) then
        phi = (latitude(1,j)-shift_sst)*pi/180.
       else
        yy = y_gl(j+jt)
        phi = pi/180.*(yy*2.5e-8*earth_factor*360.-shift_sst)
       end if
       do i=1,nx
         if(landmask(i,j).eq.0.) then
             sstxy(i,j) =  tabs_s+delta_sst*(1-sin(phi)**6)-t00
            ! add warm pool
             if(ocean_type.eq.81) then
              if(dolatlon) then
                lamda = longitude(i,j)*pi/180.
              else
                lamda = (i+it-1)*dx*pi/180.*2.5e-8*earth_factor*360.
              end if
              if(abs(lamda-pi/2.).le.(pi/6.).and.abs(phi).le.pi/12.) then
                 sstxy(i,j) = sstxy(i,j)+1.*cos(pi/2.*(lamda-pi/2)/(pi/6.))**2*cos(pi/2.*phi/(pi/12.))**2
              end if
             end if
         end if
       end do
     end do

   case(9) !  Aquaplanet Experiment (APE) "QOBS" and Warm Pool (Neale and Hoskins 2001)

     call task_rank_to_index(rank,it,jt)

     lx = float(nx_gl)*dx
     ly = 0.5*yv_gl(ny_gl+1)
     do j=1,ny
       if(dolatlon) then
        phi = (latitude(1,j)-shift_sst)*pi/180.
       else
        yy = y_gl(j+jt)
        phi = pi/180.*(yy*2.5e-8*earth_factor*360.-shift_sst)
       end if
       do i=1,nx
         if(landmask(i,j).eq.0.) then
           if(dolatlon) then
             lamda = longitude(i,j)*pi/180.
           else
             lamda = (i+it-1)*dx*pi/180.*2.5e-8*earth_factor*360.
           end if
           if(abs(phi).le.pi/3.) then
              sstxy(i,j) = tabs_s+delta_sst*(1-0.5*(sin(1.5*phi)**2+sin(1.5*phi)**4))-t00
           else
              sstxy(i,j) = tabs_s - t00
           end if
           ! add warm pool
           if(abs(lamda-pi/2.).le.(pi/6.).and.abs(phi).le.pi/12.) then
              sstxy(i,j) = sstxy(i,j)+3.*cos(pi/2.*(lamda-pi/2)/(pi/6.))**2*cos(pi/2.*phi/(pi/12.))**2
           end if
         end if
       end do
     end do

   case(10) ! Sinusoidal distribution along the y-direction for RCEMIP Phase II

     call task_rank_to_index(rank,it,jt)

     ly = pi
     do j=1,ny
        phi = (latitude(1,j)-shift_sst)*pi/180.
       do i=1,nx
         if(landmask(i,j).eq.0.) sstxy(i,j) = tabs_s+0.5*delta_sst*cos(2*pi*phi/ly) - t00
       end do
     end do

   case default

     if(masterproc) then
         print*, 'unknown ocean type in set_sst. Exitting...'
         call task_abort
     end if

end select

end subroutine set_sst



SUBROUTINE sst_evolve
 use vars, only: sstxy, t00, fluxbt, fluxbq, rhow, qocean_xy, landmask, prespoti, raf, rafq, &
                 tabs, seaice_h
 use params, only: doslabocean, cp, lcond, tabs_s, ocean_type, dossthomo, &
                   depth_slab_ocean, Szero, deltaS, timesimpleocean, dossthomozonal, &
                   sst_climo, tau_ocean, dosstclimo, sst_mean, salt_factor, doequilocean, &
                   doseaice, seaicethickness, doseaiceevol
 use rad, only: swnsxy, lwnsxy, lwdsxy, initrad
 use terrain, only: k_terra
 use slm_vars, only: seaicemask, icemask
 use slm_params, only: cond_ice
 use consts, only: lcond, sigmaSB, emis_water

 real, parameter :: rhor = 1000. ! density of water (kg/m3)
 real, parameter :: cw = 4187.   ! Liquid Water heat capacity = 4187 J/kg/K
 real factor_cp, factor_lc, qoceanxy
 real tmpx(nx), lx
 real(8) sss(2),ssss(2),sstnudge
 real, external :: dtqsatw, dtqsati
 real f(nx,ny),df(nx,ny)
 real dtqsat
 real icecoef, tmp
 real fm(ny),dfm(ny)
 integer i,j,k,m

   if(initrad) return

   if(doslabocean.or.doequilocean) then
      if(time.lt.timesimpleocean) return
      lx = float(nx_gl)*dx
      do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
      end do
      do j=1,ny
       do i=1,nx
         if(landmask(i,j).eq.0.and.seaicemask(i,j).eq.0) then
           if(dosstclimo) then
             qoceanxy = -(sst_mean-sst_climo)/tau_ocean*rhor*cw*depth_slab_ocean
           else
             qoceanxy = Szero + deltaS*abs(2.*tmpx(i)/lx - 1)
           end if
           qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy*dtfactor
         end if
       end do
      end do
      if(doslabocean) then
      ! Use forward Euler to integrate the differential equation
      ! for the ocean mixed layer temperature: dT/dt = S - E.
      ! The source: CPT?GCSS WG4 idealized Walker-circulation 
      ! RCE Intercomparison proposed by C. Bretherton.
        do j=1,ny
         do i=1,nx
           k = k_terra(i,j)
           if(landmask(i,j).eq.0.and.seaicemask(i,j).eq.0) then
             factor_cp = rhow(k)*cp
             factor_lc = rhow(k)*lcond
             sstxy(i,j) = sstxy(i,j) &
                   + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
                   - lwnsxy(i,j)               & ! LW Radiative Heating
                   - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                   - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                   + qoceanxy)            & ! Ocean Heating
                   /(rhor*cw*depth_slab_ocean)         ! Convert W/m^2 Heating to K/s
           end if
         end do
        end do
      end if ! doslabocean

      if(doequilocean) then
    ! Use Neuton's method to itarate the SST to ocean surface equilibrium
    ! MK 2020
       icecoef = cond_ice/seaicethickness
       seaice_h(:,:) = 0.
       if(dossthomozonal) then
        factor_cp = rhow(1)*cp
        factor_lc = rhow(1)*lcond
        do j=1,ny
         do i=1,nx
            f(i,j) = &
                   swnsxy(i,j)          & ! SW net solar flux
                 + lwdsxy(i,j)               & ! downwelling LW Flux
                 - emis_water*sigmaSB*(sstxy(i,j)+t00)**4  & ! upwelling LW Flux
                 - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                 - factor_lc*fluxbq(i,j)      ! Latent Heat Flux
            if(doseaice.and.sstxy(i,j)+t00.lt.271._8) then
              dtqsat = dtqsati(real(sstxy(i,j))+t00,presi(1))
              icemask(i,j) = 1.
              f(i,j) = f(i,j) + icecoef*(271._8 - sstxy(i,j)-t00)
              if(doseaiceevol) then
  ! parameterization of relationship between seaice thinkness 
  ! and temperature above sea-ice. Loosely based on parameterization by
  ! Levedev (1938) as described by Y. Yu and R. W. Lindsay (JGR 2003, doi:10.1029/2002JC001319)
  ! assuming 90 days for ice accumulation time scale. -MK
                seaice_h(i,j) = 0.0133*(90.*max(0.,271.-tabs(i,j,1)))**0.58
                tmp = cond_ice/max(0.1,seaice_h(i,j)) 
              else
                tmp = icecoef
              end if
            else
              dtqsat = salt_factor*dtqsatw(real(sstxy(i,j))+t00,presi(1))
              icemask(i,j) = 0.
              tmp = 0.
            end if
            df(i,j) = &
                     rhow(1)*(cp*prespoti(1)/raf(i,j) &
                   + lcond/rafq(i,j)*dtqsat) &
                   + 4.*sigmaSB*(sstxy(i,j)+t00)**3 + tmp
         end do
        end do
        call mean_x_2D(f(1:nx,1:ny), fm(1:ny))
        call mean_x_2D(df(1:nx,1:ny), dfm(1:ny))
        do j=1,ny
          sstxy(:,j) = sstxy(:,j) + fm(j)/dfm(j)
!          if(any(isnan(sstxy(:,j)))) then
!            print*,'NAN sstxy>>>',j,minval(f(:,j)),maxval(f(:,j)),fm(j), &
!                              minval(swnsxy(:,j)),maxval(swnsxy(:,j)), &
!                              minval(lwdsxy(:,j)),maxval(lwdsxy(:,j)), &
!                              minval(factor_cp*fluxbt(:,j)),maxval(factor_cp*fluxbt(:,j)), &
!                              minval(factor_lc*fluxbq(:,j)),maxval(factor_lc*fluxbq(:,j)), &
!                              minval(raf(:,j)),maxval(raf(:,j)), &
!                              minval(df(:,j)),maxval(df(:,j)),dfm(j)
!            call task_abort()
!          end if
        end do
       else
        do j=1,ny
         do i=1,nx
          k = k_terra(i,j)
          if(landmask(i,j).eq.0) then
            factor_cp = rhow(k)*cp
            factor_lc = rhow(k)*lcond
            if(doseaice.and.sstxy(i,j)+t00.lt.271._8) then
              dtqsat = dtqsati(real(sstxy(i,j))+t00,presi(k))
              icemask(i,j) = 1.
            else
              dtqsat = salt_factor*dtqsatw(real(sstxy(i,j))+t00,presi(k))
              icemask(i,j) = 0.
            end if
            sstxy(i,j) = sstxy(i,j) &
                 + (swnsxy(i,j)          & ! SW net solar flux
                 + lwdsxy(i,j)               & ! downwelling LW Flux
                 - emis_water*sigmaSB*(real(sstxy(i,j))+t00)**4  & ! upwelling LW Flux
                 - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                 - factor_lc*fluxbq(i,j))     & ! Latent Heat Flux
                 /(rhow(k)/raf(i,j)*(cp*prespoti(k) &
                   + lcond*dtqsat) &
                   + 4.*sigmaSB*(real(sstxy(i,j))+t00)**3)
          end if
         end do
        end do
       end if ! dossthomozonal
      end if

      if(dossthomo.or.dosstclimo) then
        sss = 0.
        do j=1,ny
         do i=1,nx
           if(landmask(i,j).eq.0) then
            sss(1) = sss(1) + sstxy(i,j)
            sss(2) = sss(2) + 1
           end if
         end do
        end do
        if(dompi) then
            call task_sum_real8(sss,ssss,2)
            sss = ssss
        end if ! dompi
        if(sss(2).eq.0.) then
           if(masterproc) print*,'no ocean points found in dynamic ocean!'
           call task_abort
        end if
        sst_mean = sss(1) / sss(2) + t00
        if(dossthomo) then
          sstxy = sst_mean - t00
          tabs_s = sst_mean
        end if
      end if

   end if ! doslabocean.or.doequilocean

end subroutine sst_evolve


end module simple_ocean
