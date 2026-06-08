! Initialize 2D output

subroutine stat_2Dinit()

use vars
use params, only: dotracers
use tracers, only: tr_xy, tr_ac_xy, tr_acs_xy
implicit none

     dtp = 0.
     prec_xy(:,:) = 0.
     precs_xy(:,:) = 0.
     taux_xy(:,:) = 0.
     tauy_xy(:,:) = 0.
     taux_top_xy(:,:) = 0.
     tauy_top_xy(:,:) = 0.
     shf_xy(:,:) = 0.
     lhf_xy(:,:) = 0.
     shf_top_xy(:,:) = 0.
     lhf_top_xy(:,:) = 0.
     pw_xy(:,:) = 0.
     pwobs_xy(:,:) = 0.
     pws_xy(:,:) = 0.
     fse_xy(:,:) = 0.
     ta_xy(:,:) = 0.
     cw_xy(:,:) = 0.
     iw_xy(:,:) = 0.
     rw_xy(:,:) = 0.
     sw_xy(:,:) = 0.
     gw_xy(:,:) = 0.
     cld_xy(:,:) = 0.
     eis_xy(:,:) = 0.
     u200_xy(:,:) = 0.
     v200_xy(:,:) = 0.
     uobs200_xy(:,:) = 0.
     vobs200_xy(:,:) = 0.
     usfc_xy(:,:) = 0.
     vsfc_xy(:,:) = 0.
     tsfc_xy(:,:) = 0.
     qsfc_xy(:,:) = 0.
     w500_xy = 0.
     phi500_xy = 0.
     ZdBZ_xy = 0.
     vor200_xy(:,:) = 0.
     vor850_xy(:,:) = 0.
     omega200_xy = 0.
     omega500_xy = 0.
     omega700_xy = 0.
     omega850_xy = 0.
     rh200_xy = 0.
     rh500_xy = 0.
     rh700_xy = 0.
     rh850_xy = 0.
     lwnt_xy(:,:) = 0.
     swnt_xy(:,:) = 0.
     lwntc_xy(:,:) = 0.
     swntc_xy(:,:) = 0.
     lwns_xy(:,:) = 0.
     swns_xy(:,:) = 0.
     lwnsc_xy(:,:) = 0.
     swnsc_xy(:,:) = 0.
     lwds_xy(:,:) = 0.
     swds_xy(:,:) = 0.
     lwdsc_xy(:,:) = 0.
     swdsc_xy(:,:) = 0.
     solin_xy(:,:) = 0.
     qocean_xy(:,:) = 0.
     soil_wet_xy(:,:) = 0.
     snow_depth_xy(:,:) = 0.
     snow_melt_xy(:,:) = 0.
     sst_xy(:,:) = 0.
     alb_xy(:,:) = 0.
     albc_xy(:,:) = 0.
     if(mod(nstep,nsave2D).eq.0) then
      dtpa = 0.
      preca_xy(:,:) = 0.
      precsa_xy(:,:) = 0.
      lwnta_xy(:,:) = 0.
      swnta_xy(:,:) = 0.
      lwnsa_xy(:,:) = 0.
      swnsa_xy(:,:) = 0.
      lwdsa_xy(:,:) = 0.
      swdsa_xy(:,:) = 0.
      lwntca_xy(:,:) = 0.
      swntca_xy(:,:) = 0.
      lwnsca_xy(:,:) = 0.
      swnsca_xy(:,:) = 0.
      lwdsca_xy(:,:) = 0.
      swdsca_xy(:,:) = 0.
      solina_xy(:,:) = 0.
      gust10m_xy(:,:) = 0.
      u10ma_xy(:,:) = 0.
      v10ma_xy(:,:) = 0.
     end if
     if(nstep.eq.0) then
       precac_xy = 0.
       precsac_xy = 0.
       tr_ac_xy = 0.
       tr_acs_xy = 0.
       lwntac_xy = 0.
       swntac_xy = 0. 
       lwnsac_xy = 0.
       swnsac_xy = 0.
       lwdsac_xy = 0.
       swdsac_xy = 0.
     end if
     t2m_xy(:,:) = 0.
     q2m_xy(:,:) = 0.
     tsoil_xy(:,:) = 0.
     grnd_xy(:,:) = 0.
     ra_xy(:,:) = 0.
     rc_xy(:,:) = 0.
     u10m_xy(:,:) = 0.
     v10m_xy(:,:) = 0.
     if(dotracers) then
       tr_xy(:,:,:) = 0.
     end if

!===================================
! UW ADDITIONS: MOSTLY 2D STATISTICS

    !bloss: store initial profiles for computation of storage terms in budgets
     ustor(:) = u0(1:nzm)
     vstor(:) = v0(1:nzm)
     tstor(:) = t0(1:nzm)
     qstor(:) = q0(1:nzm)

     utendcor(:) = 0.
     vtendcor(:) = 0.

     psfc_xy(:,:) = 0.

     u850_xy(:,:) = 0.
     v850_xy(:,:) = 0.
     uobs850_xy(:,:) = 0.
     vobs850_xy(:,:) = 0.

     swvp_xy(:,:) = 0.

     cloudtopheight(:,:) = 0.
     echotopheight(:,:) = 0.
     cloudtoptemp(:,:) = 0.

! END UW ADDITIONS

end

