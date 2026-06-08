! interface to ocean model

subroutine dyn_ocean()

use vars
use rad, only: swnsxy, lwnsxy
use slm_vars, only: seaicemask
use params 
implicit none

real coef

coef = 1./float(ncycle)
shf_ocn(:,:) = shf_ocn(:,:) + fluxbt(:,:)*coef
lhf_ocn(:,:) = lhf_ocn(:,:) + fluxbq(:,:)*coef
taux_ocn(:,:) = taux_ocn(:,:) + fluxbu(:,:)*coef
tauy_ocn(:,:) = tauy_ocn(:,:) + fluxbv(:,:)*coef
prec_ocn(:,:) = prec_ocn(:,:) + precinst(:,:)*coef
lw_ocn(:,:) = lw_ocn(:,:) + lwnsxy(:,:)*coef
sw_ocn(:,:) = sw_ocn(:,:) + swnsxy(:,:)*coef

if(mod(nstep,ncallocean).eq.0.and.icycle.eq.ncycle) then
  coef = 1./float(ncallocean)
  where (seaicemask.eq.0) 
    shf_ocn(:,:) = shf_ocn(:,:)*coef*rhow(1)*cp
    lhf_ocn(:,:) = lhf_ocn(:,:)*coef*rhow(1)*lcond
    taux_ocn(:,:) = taux_ocn(:,:)*rhow(1)*coef
    tauy_ocn(:,:) = tauy_ocn(:,:)*rhow(1)*coef
    prec_ocn(:,:) = prec_ocn(:,:)*coef
    lw_ocn(:,:) = lw_ocn(:,:)*coef
    sw_ocn(:,:) = sw_ocn(:,:)*coef
  elsewhere
    shf_ocn(:,:) = 0.
    lhf_ocn(:,:) = 0.
    taux_ocn(:,:) = 0.
    tauy_ocn(:,:) = 0.
    prec_ocn(:,:) = 0.
    lw_ocn(:,:) = 0.
    sw_ocn(:,:) = 0.
  end where

  call ocn(real(dt*ncallocean,4))

  shf_ocn(:,:) = 0.
  lhf_ocn(:,:) = 0.
  taux_ocn(:,:) = 0.
  tauy_ocn(:,:) = 0.
  prec_ocn(:,:) = 0.
  lw_ocn(:,:) = 0.
  sw_ocn(:,:) = 0.

end if

end
