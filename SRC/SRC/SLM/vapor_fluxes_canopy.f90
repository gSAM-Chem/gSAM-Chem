!=====================================================================================================================
! subroutine vapor_fluxes_canopy
! Note: 
! computes latent heat fluxes from canopy (for vegetated land), and combines them to return the total to SAM
!
! Input:
! qsfc: specific humidity at surface level 
!	(the surface that the atmosphere sees
! i,j : indices for the grid location
!
! History : Jungmin Lee, March, 2016
! MK Jan 2020: separated into separate subroutone for canopy and soil
!=====================================================================================================================

 SUBROUTINE vapor_fluxes_canopy(qsfc,pressf,rhosf,dtn_iter,i,j)

 use slm_vars, only : DBL,nsoil, soilw, w_s_FC, q_gr, pii, r_d, r_litter, r_b, &
     w_s_WP,evapo_s, vege_YES, wet_canop, qsat_canop, t_canop, mw, tfriz, &
     mw_mx, evapo_wet,rootF, evapo_dry, dtn_dbl, lhf_canop, lcond, lsub, r_c, r_soil, mw_inc, &
     icemask, evp_canop
 IMPLICIT NONE

 REAL (KIND=DBL), INTENT(IN) :: qsfc 
 REAL, INTENT(IN) :: pressf,rhosf 
 REAL (KIND=DBL), INTENT(IN) :: dtn_iter 
 INTEGER,INTENT(IN) :: i,j
        
 INTEGER :: k
 REAL (KIND=DBL) :: flag_ice ! flag to switch off evapotranspiration when freezing T
 REAL (KIND=DBL) :: evapo_dry0  ! actual evapotranspitation (from soil layers above wilting point)
 real, external :: qsatw,qsati

!=====================================================================================================================
! Evaporation from canopy
!=====================================================================================================================
 if(icemask(i,j).eq.0.and.t_canop(i,j).ge.tfriz) then
    qsat_canop = qsatw(real(t_canop(i,j)),pressf)
    flag_ice = 1.
 else
    qsat_canop = qsati(real(t_canop(i,j)),pressf)
    flag_ice = 0.
 end if

 ! direct evaporation from the water held on canopy
 ! evaporation/dew only possible if canopy temperature is above freezing
   evapo_wet = min(mw(i,j)/dtn_iter, &
             flag_ice*(qsat_canop-qsfc)*rhosf/(2.0_DBL*r_b))*vege_YES(i,j)

 ! increment/decrement of the water amount held on leaves following the direct evaporation/dew formation
 mw_inc = -dtn_iter*evapo_wet  ! evapo_wet [kg/m2s=mm/s] 
 mw(i,j) = mw(i,j) +  mw_inc
 wet_canop = min(1.0_DBL,mw(i,j)/mw_mx(i,j))
 evapo_wet = wet_canop * evapo_wet

 ! Transpiration - only occures when qsat_canop > qsfc 
 evapo_dry = max(0.,flag_ice*(qsat_canop-qsfc)*rhosf* &
             (1.0_DBL-wet_canop)/(2.0_DBL*r_b+r_c))*vege_YES(i,j)
 if(evapo_dry.gt.0.) then
  evapo_dry0 = evapo_dry       
  ! Check soil moisture availability for transpiration
  do k = 1, nsoil
    if(soilw(k,i,j).lt.0.05) then
       evapo_dry = evapo_dry - evapo_dry0*rootF(k,i,j)
    end if
  end do! k
 end if

!=====================================================================================================================
! Convert evaporation (kg/m2/s) to latent heat flux (W/m2)
!=====================================================================================================================
 evp_canop(i,j) = evapo_wet+evapo_dry
 lhf_canop(i,j) = (flag_ice*lcond+(1.-flag_ice)*lsub)*(evapo_wet+evapo_dry)


END SUBROUTINE vapor_fluxes_canopy
