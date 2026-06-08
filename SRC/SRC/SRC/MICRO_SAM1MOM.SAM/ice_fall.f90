
subroutine ice_fall()


! Sedimentation of ice:

use vars
use microphysics, only: micro_field, index_cloud_ice, total_water, icefall_fudge
!use micro_params
use params
use terrain, only: k_terra,terraw
use stat_coars

implicit none

integer i,j,k, kb, kc, kmax, kmin, ici
real coef,coef1,dqi,lat_heat,vt_ice
real omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
real delta,fudge(ny),fudge1km,fudge4km 
real fz(nx,ny,nz)

if(index_cloud_ice.eq.-1) return

call t_startf ('ice_fall')

total_water_prec = total_water_prec + total_water()

kmax=0
kmin=nzm+1

do k = 1,nzm
 do j = 1, ny
  do i = 1, nx
      if(qcl(i,j,k)+qci(i,j,k).gt.0..and. tabs(i,j,k).lt.273.15) then
        kmin = min(kmin,k)
        kmax = max(kmax,k)
      end if
  end do
 end do
end do

do k = 1,nzm
   qifall(k) = 0.
   tlatqi(k) = 0.
end do


if(doglobalpresets) then
  ! compute terminal velocity fudge factor depending on local grid spacing
  fudge1km = 0.8 ! fudge for 1 km
  fudge4km = 0.3 ! fudge for 4 km
  do j=1,ny
    delta = min(8.,0.001*sqrt(dx*mu(j)*dy*ady(j))) ! compute effective local grid spacing in km)
    fudge(j) = max(0.2,min(1.,0.3333*((fudge4km-fudge1km)*delta +4.*fudge1km-fudge4km)))
  end do
else
  fudge(:) = icefall_fudge 
end if

fz = 0.

! Compute cloud ice flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002). 
do k = max(1,kmin-1),kmax
   ! Set up indices for x-y planes above and below current plane.
   kc = min(nzm,k+1)
   kb = max(1,k-1)
   ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
   coef = dtn/(0.5*(adz(kb)+adz(k))*dz)
   do j = 1,ny
      do i = 1,nx
        if(k.ge.k_terra(i,j)) then
         ! Compute cloud ice density in this cell and the ones above/below.
         ! Since cloud ice is falling, the above cell is u (upwind),
         ! this cell is c (center) and the one below is d (downwind). 

         qiu = rho(kc)*qci(i,j,kc)
         qic = rho(k) *qci(i,j,k) 
         qid = rho(kb)*qci(i,j,kb) 

         ! Ice sedimentation velocity depends on ice content. The fiting is
         ! based on the data by Heymsfield (JAS,2003). -Marat
         ! 0.1 m/s low bound was suggested by Chris Bretherton 
!          vt_ice = 0.0
!         vt_ice = max(0.1,0.5*log10(qic+1.e-12)+3.) ! based on Heymsfield's figure
         vt_ice = fudge(j)*8.66*(max(0.,qic)+1.e-10)**0.24   ! Heymsfield (JAS, 2003, p.2607)
!         vt_ice = 0.1
         vt_ice = vt_ice/gamma_RAVE  ! MK - slow down sedimentation, predominantly in deep anvils

         ! Use MC flux limiter in computation of flux correction.
         ! (MC = monotonized centered difference).
         if (abs(qic-qid).lt.1.e-6) then
            tmp_phi = 0.
         else
            tmp_theta = (qiu-qic)/(qic-qid)
            tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
         end if

         ! Compute limited flux.
         ! Since falling cloud ice is a 1D advection problem, this
         ! flux-limited advection scheme is monotonic.
         fz(i,j,k) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid))
       end if
      end do
   end do
end do
fz(:,:,nz) = 0.

ici = index_cloud_ice

coef1 = dt/dz 
do k=max(1,kmin-2),kmax
   coef=dtn/(dz*adz(k)*rho(k))
   do j=1,ny
      do i=1,nx
       if(k.ge.k_terra(i,j)) then
         ! The cloud ice increment is the difference of the fluxes.
         dqi=coef*(fz(i,j,k)-fz(i,j,k+1))
         ! Add this increment to both non-precipitating and total water.
         micro_field(i,j,k,ici)  = micro_field(i,j,k,ici)  + dqi
         ! Include this effect in the total moisture budget.
         qifall(k) = qifall(k) + dqi

         ! The latent heat flux induced by the falling cloud ice enters
         ! the liquid-ice static energy budget in the same way as the
         ! precipitation.  Note: use latent heat of sublimation. 
         lat_heat  = (fac_cond+fac_fus)*dqi
         ! Add divergence of latent heat flux to liquid-ice static energy.
         t(i,j,k)  = t(i,j,k)  - lat_heat
         ! Add divergence to liquid-ice static energy budget.
         tlatqi(k) = tlatqi(k) - lat_heat
         precflux(k) = precflux(k) - fz(i,j,k)*coef1*wgt(j,k)
       end if
      end do
   end do
end do

do j=1,ny
   do i=1,nx
     precinst(i,j) = precinst(i,j) - fz(i,j,k_terra(i,j))
     prectot = prectot - wgty(j)*fz(i,j,k_terra(i,j))*dtn ! For statistics
     if(dostatis) then
        precsfc(i,j) = precsfc(i,j) - fz(i,j,k_terra(i,j))*dtn ! For statistics
        precsfcsnow(i,j) = precsfcsnow(i,j) - fz(i,j,k_terra(i,j))*dtn ! For statistics
     end if
     prec_xy(i,j) = prec_xy(i,j) - fz(i,j,k_terra(i,j))*dtn ! For 2D output
     preca_xy(i,j) = preca_xy(i,j) - fz(i,j,k_terra(i,j))*dtn ! For 2D output
   end do
end do

if(collect_coars) then
     call coars_fld(fz,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm), &
                                     fld_flux(:,:,:,4+nmicro_fields+3))
     tmp_phi = fac_cond+fac_fus
     call coars_fld(fz*tmp_phi,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm), &
                                     fld_flux(:,:,:,4+nmicro_fields+4))
end if


total_water_prec = total_water_prec - total_water()

call t_stopf ('ice_fall')

end subroutine ice_fall

