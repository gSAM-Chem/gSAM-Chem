
subroutine cloud_fall()

! Sedimentation of cloud water:

use vars
use microphysics, only: micro_field, total_water, donomicro, sigmag, Nc_land, Nc_ocn
use params
use terrain, only: k_terra,terraw
use stat_coars

implicit none

integer i,j,k, kb, kc, kmax(nx,ny), kmin, ici, n, nsub
real dtc,coef,coef1,dqi,lat_heat,vt_cl,coef_cl,www
real omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
real fz(nz), Nc0, qcl_max(nx,ny), vt_cl_max 

if(donomicro) return

call t_startf ('cloud_fall')

total_water_prec = total_water_prec + total_water()

kmax = 0
qcl_max = 0.

do k = 1,nzm
 do j = 1, ny
  do i = 1, nx
      if(qcl(i,j,k).gt.0.) then
        kmax(i,j) = max(kmax(i,j),k)
        qcl_max(i,j) = max(qcl_max(i,j),qcl(i,j,k))
      end if
  end do
 end do
end do

do k = 1,nzm
   qcfall(k) = 0.
   tlatqc(k) = 0.
end do

! assume lognormal distribution of cloud droplets by size
! following Bretherton, C. S., P. N. Blossey, and J. Uchida (GRL,2007)
coef_cl = 1.19e8*(3./(4.*3.1415*1000.*1.e6))**(2./3.)*exp(5.*log(sigmag)**2)

ici = 1  ! cloud water index 
coef1 = dt/dz 

! Compute cloud flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002). 
do j = 1,ny
 do i = 1,nx
  fz = 0.
  Nc0 = landmask(i,j)*Nc_land + (1.-landmask(i,j))*Nc_ocn
  vt_cl_max = coef_cl*(qcl_max(i,j)/Nc0)**(2./3.)
  nsub = max(1,nint(0.9*(vt_cl_max*dtn)/(minval(adz(:))*dz)))
  dtc = dtn/nsub
  www = 1./nsub
  kmin = k_terra(i,j) 
  do n=1,nsub
    do k = kmin,kmax(i,j)
     ! Set up indices for x-y planes above and below current plane.
     kc = min(nzm,k+1)
     kb = max(kmin,k-1)
     ! CFL number based on grid spacing interpolated to interface k-1/2
     coef = dtc/(0.5*(adz(kb)+adz(k))*dz)
     ! Compute cloud water density in this cell and the ones above/below.
     ! Since cloud water is falling, the above cell is u (upwind),
     ! this cell is c (center) and the one below is d (downwind). 

     qiu = max(0.,rho(kc)*qcl(i,j,kc))
     qic = max(0.,rho(k) *qcl(i,j,k))
     qid = max(0.,rho(kb)*qcl(i,j,kb)) 

     vt_cl = min(vt_cl_max,coef_cl*(qic/Nc0)**(2./3.))

     ! Use MC flux limiter in computation of flux correction.
     ! (MC = monotonized centered difference).
     if (abs(qic-qid).lt.1.e-6) then
        tmp_phi = 0.
     else
        tmp_theta = (qiu-qic)/(qic-qid)
        tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
     end if

     ! Compute limited flux.
     ! Since falling cloud water is a 1D advection problem, this
     ! flux-limited advection scheme is monotonic.
     fz(k) = -vt_cl*(qic - 0.5*(1.-coef*vt_cl)*tmp_phi*(qic-qid))
    end do ! k
    do k = kmin,kmax(i,j)
     coef=dtc/(dz*adz(k)*rho(k))
     ! The cloud water increment is the difference of the fluxes.
     dqi=coef*(fz(k)-fz(k+1))
     ! Add this increment to both non-precipitating and total water.
     micro_field(i,j,k,ici)  = micro_field(i,j,k,ici)  + dqi
     qcl(i,j,k)  = qcl(i,j,k)  + dqi
     ! Include this effect in the total moisture budget.
     qcfall(k) = qcfall(k) + dqi

     ! The latent heat flux induced by the falling cloud water enters
     ! the liquid-ice static energy budget in the same way as the
     ! precipitation.  
     lat_heat  = fac_cond*dqi
     ! Add divergence of latent heat flux to liquid-ice static energy.
     t(i,j,k)  = t(i,j,k)  - lat_heat
     ! Add divergence to liquid-ice static energy budget.
     tlatqc(k) = tlatqc(k) - lat_heat
     precflux(k) = precflux(k) - fz(k)*coef1*wgt(j,k)
    end do
    precinst(i,j) = precinst(i,j) - fz(k_terra(i,j))*www
    prectot = prectot - wgty(j)*fz(kmin)*dtc ! For statistics
    if(dostatis) then
      precsfc(i,j) = precsfc(i,j) - fz(kmin)*dtc ! For statistics
    end if
    prec_xy(i,j) = prec_xy(i,j) - fz(kmin)*dtc ! For 2D output
    preca_xy(i,j) = preca_xy(i,j) - fz(kmin)*dtc ! For 2D output
  end do ! n
 end do
end do


total_water_prec = total_water_prec - total_water()

call t_stopf ('cloud_fall')

end subroutine cloud_fall

