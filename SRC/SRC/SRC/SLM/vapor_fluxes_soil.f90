!=====================================================================================================================
! subroutine vapor_fluxes_soil
! Note: 
! computes latent heat fluxes from topsoil, canopy (for vegetated land), and combines them to return the total to SAM
!
! Input:
! qsfc: specific humidity at surface level 
!	(the surface that the atmosphere sees: for instance, topsoil for baresoil case, and canopy for vegetated land)
! qr:   specific humidity at the reference level (input from the lowest level of SAM) 
! i,j : indices for the grid location
!	(Indices are passed to access other variables of the same grid point, in case landmask includes ocean)
!
! History : Jungmin Lee, March, 2016
! MK Jani 2020, separated original vapor_fluxes into two for canopy and soil.
!=====================================================================================================================

 SUBROUTINE vapor_fluxes_soil(qsfc,qr,pressf,rhosf,sdew,i,j)

 use slm_vars, only : DBL,nsoil, vegetated, soilw, w_s_FC, q_gr, pii, soilt, r_a, r_d, r_litter, &
     w_s_WP,evapo_s, vege_YES, wet_canop, qsat_canop, IMPERV, &
     mw_mx, evapo_wet,rootF, evapo_dry, dtn_dbl, lhf_canop, lhf_soil, lhf_air, lcond, lsub, r_soil, mw_inc, &
     icemask, snow_mass, evp_canop, evp_soil, evp_air, tfriz
 IMPLICIT NONE

 REAL (KIND=DBL), INTENT(IN) :: qsfc 
 REAL, INTENT(IN) :: qr,pressf,rhosf 
 REAL(KIND=DBL), INTENT(IN) :: sdew  ! 0 - soil is not saturated, 1 - saturated (dew possible) ! -Marat 03/19
 INTEGER,INTENT(IN) :: i,j
        
 INTEGER :: k
 REAL (KIND=DBL) :: totalR_soil,& ! total aerodynamic resistance over soil surface
   qref_tmp, & 
   soil_diff      ! soil diffusion factor for the topsoil evaporation from soil pore space to the overlying air
 REAL (KIND=DBL) :: ice_flag ! flag to switch off evapotranspiration when freezing T
 real, external :: qsatw,qsati

!=====================================================================================================================
! Evaporation from the Top soil layer 
!=====================================================================================================================

!=====================================================================================================================
! New reference level, which is the end point for water vapor turbulent flux, is defined over the vegetated land 
! 	for baresoil soil => reference level q stays to be qr
! 	for vegetated => vegetation layer q, qsfc, become the new reference level q
!=====================================================================================================================
 IF(vegetated(i,j)) THEN ! for baresoil
   qref_tmp = qsfc
 ELSE
   qref_tmp = DBLE(qr)
 END IF

!=========================================================================================================================
! make sure a correct formula is used when freezing temperature or ice surface (MK 11.20)
 if(icemask(i,j).eq.1.or.snow_mass(i,j).gt.0..or.soilt(1,i,j).lt.tfriz) then
  ice_flag = 0.
 else
  ice_flag = 1.
 end if
!=====================================================================================================================
! Soil diffusion factor (soil_diff) follows "Sakaguchi and Zeng, 2009", but the snow factor is discarded in SLM adoptation
! Maximum evaporation is assumed when dew formation and/or soil top layer wetness is greater than field capacity
!=====================================================================================================================
 IF((soilw(1,i,j).ge.w_s_FC(1,i,j)).or.(qref_tmp.gt.q_gr)) THEN
   soil_diff = 1._DBL
 ELSE
   soil_diff = 0.25_DBL*((1.0_DBL-cos(pii*max(0.01,soilw(1,i,j))/w_s_FC(1,i,j)))**2)
 END IF

!=====================================================================================================================
! total resistance for the topsoil evaporation = moisture diffusion factor + aerodynamic resistance
! r_soil base value is applied to prevent too active evaporation and roughly represents the ground litter resistance
! r_litter at this version is not in use, contains zero
!=====================================================================================================================
 IF(.not.vegetated(i,j)) then ! for baresoil
   if(icemask(i,j).eq.1.or.snow_mass(i,j).gt.0.) then
     r_soil = 0.
   else 
     r_soil = min(10000._DBL,max(100._DBL,r_a*(1.0_DBL/soil_diff-1.0_DBL))) 
   end if
   totalR_soil = r_soil + r_a
 ELSE
   r_soil = min(10000._DBL,max(50._DBL,r_d*(1.0_DBL/soil_diff-1.0_DBL)))
   totalR_soil = r_soil + r_d + r_litter
 END IF

!=====================================================================================================================
! Evaporation from soil top underneath the canopy
!=====================================================================================================================
 evapo_s  = vege_YES(i,j)*(q_gr-qsfc)*rhosf/totalR_soil
 if(evapo_s.lt.0.) evapo_s = evapo_s*sdew ! allow negative evaporation only when soil is saturated (dew) ! -Marat

!=====================================================================================================================
! Convert evaporation (kg/m2/s) to latent heat flux (W/m2)
!=====================================================================================================================
  
 IF(.not.vegetated(i,j)) THEN
 ! 	For baresoil case, lhf_air is equal to lhf_soil, lhf_air is updated to include soil diffusion factor 
   evp_air(i,j) = (qsfc-DBLE(qr))*rhosf/(r_soil+r_a) * (1.-IMPERV(i,j))
   if(evp_air(i,j).lt.0.) evp_air(i,j) = evp_air(i,j)*sdew
   evapo_s = evp_air(i,j)
   lhf_air(i,j) = (ice_flag*lcond+(1.-ice_flag)*lsub)*evp_air(i,j)
   evp_soil(i,j) = evp_air(i,j)
 else 
   evp_soil(i,j) = evapo_s * (1.-IMPERV(i,j))
   evp_air(i,j) = evp_canop(i,j) + evp_soil(i,j)
   lhf_air(i,j) = lhf_canop(i,j) + (ice_flag*lcond+(1.-ice_flag)*lsub)*evp_soil(i,j)
 end if
 lhf_soil(i,j) = (ice_flag*lcond+(1.-ice_flag)*lsub)*evp_soil(i,j)



!energy_bal(i,j) = energy_bal(i,j) - lhf_air(i,j) - shf_air(i,j)
!energy_bal_c(i,j) = energy_bal_c(i,j) - lhf_canop(i,j) - shf_canop(i,j)
!energy_bal_s(i,j) = energy_bal_s(i,j) - lhf_soil(i.j) - shf_soil(i.j)

END SUBROUTINE vapor_fluxes_soil
