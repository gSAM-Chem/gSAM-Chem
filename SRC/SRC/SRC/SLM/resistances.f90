!=============================================================================
! subroutine resistances
! Note:
! 	compute aerodynamic resistances (r_b, r_d), stomatal resistance (r_c)
!=============================================================================

SUBROUTINE resistances(pressf,i,j)
use slm_vars, only : DBL, r_a, r_b, r_c, r_d, r_litter, vegetated, z0_soil, ustar, LAI, ggr, soilt, &
Rgl, LAI, Rc_min, Rc_max, hs_rc, t_cas, q_cas, T_opt, &
rootF, soilw, w_s_FC, s_depth, w_s_WP, theta_FC, theta_WP, poro_soil, &
ztop, nsoil, dtfactor_dbl, snow_mass, rho_snow, phi_1, phi_2
use rad, only: swdsxy
IMPLICIT NONE

REAL, INTENT(IN) :: pressf
INTEGER, INTENT(IN) :: i,j

REAL (KIND=DBL) :: &
   Cs,      & ! Total Turbulent transfer coef between canopy and underlying soil
   Cs_bare ,& !  turbulent transfer coefficient over baresoil
   Cs_dense   !  turbulent transfer coefficient under dense canopy : Sakaguchi and Zeng (2009)
REAL (KIND=DBL) :: z_m
REAL (KIND=DBL) :: rb_correc_fac, rd_correc_fac, temp_diff
REAL (KIND=DBL) :: &
   rc_fac_rad,&  ! stomatal resistance factor from radiation
   rc_fac_vpd,&  ! stomatal resistance factor from vapor pressure deficit
   rc_fac_t,&    ! stomatal rersistance factor from temperature
   rc_fac_sw     ! stomatal resistance factor from soil moisture

REAL (KIND=DBL) :: tmp, tmp2, temp, d_root, e_cas
REAL (KIND=DBL) :: k_beer, f_shade, lai_sun, lai_shade, sw_sun, sw_shade, tmp_sun, tmp_shade, rc_fac_sun, rc_fac_shade, factor
INTEGER :: k
real, external :: esatw

!==============================================================================
! Note:
! Initialize resistances with large number to prevent "divide by zero error" 
! when surface heat and water vapor fluxes are computed
!==============================================================================

if(vegetated(i,j)) then
 !=================================================================================
 ! Aerodynamic resistance for heat and vapor transfer under canopy space : r_d 
 !=================================================================================
 ! temp_diff > 0 :: stable undercanopy
 ! temp_diff < 0 :: unstable undercanopy
 temp_diff = t_cas(i,j) - soilt(1,i,j)
 
 IF(temp_diff.lt.0.0_DBL) THEN  
   factor = 1.
 ELSE                   
   ! rd_correc_fac : undercanopy stability parameter, in effect only for stable undercanopy
   rd_correc_fac = ggr*ztop(i,j)*max(0.0_DBL,temp_diff)/soilt(1,i,j)/(ustar(i,j)**2)
   factor = 1./(1.0_DBL+0.5_DBL*MIN(10.0_DBL,rd_correc_fac))
 END IF
 
 ! turbulent transfer coefficient under dense canopy
 Cs_dense = 0.004_DBL
 ! correct for stability 
 Cs_dense = Cs_dense * factor

 ! turbulent transfer coefficient over the exposed topsoil
 ! typical value of Cs_bare ~0.2
 Cs_bare = 0.4_DBL/0.13_DBL*(z0_soil*ustar(i,j)/(1.5e-5_DBL))**(-0.45_DBL) 
 
 ! Turbulence transfer coefficient undercanopy 
 ! LAI-weighed sum of Cs_bare & Cs_dense
 factor = exp(-1.0_DBL*LAI(i,j))
 Cs = Cs_bare*factor+Cs_dense*(1.0_DBL-factor)
 
 ! Undercanopy aerodynamic resistance depends on the weighed sum of the dense canopy covered soil 
 ! and baresoil turbulent transfer coefficient and friction velocity 
 ! Reference :[Oleson et al., 2004] [Zeng et al., 2005]
 r_d = min(400._DBL,1.0_DBL/ustar(i,j)/Cs) ! cap the r_d from too large values in stable conditions
! r_d = 1.0_DBL/ustar(i,j)/Cs
 
 !==================================================================================
 ! Leaf boundary layer resistance : r_b
 !==================================================================================
 ! turbulent transfer coefficient between canopy surface and canopy air : Cv = 0.01m/s^-0.5
 ! characteristic dimension of the elaves in the direction of wind flox : d_leaf = 0.04m
 !  r_b = 1.0_DBL/0.01_DBL*((ustar(i,j)/0.04_DBL)**(-0.5_DBL))/LAI(i,j)
 ! the exression above seems to overestimate the canopy LHF
                 ! especially for wet evaporation when r_c=0.. Replacing 2.*r_b with simply r_a
                 ! follows IFS model description - MK Jan 2023
 r_b = 0.5*r_a 
 !=================================================================================
 ! Stomatal resistance : r_c
 !=================================================================================
   ! radiation factor
! original code:
!   tmp = 0.55_DBL*swdsxy(i,j)*2.0_DBL/Rgl(i,j)/LAI(i,j)
!   rc_fac_rad = (Rc_min(i,j)/Rc_max + tmp)/(1.0_DBL+tmp)
! Modification assuming that parts of leaves are in shadow of other leaves for dense canopies -MK 2025
! Parameters
   k_beer = phi_1(i,j)+phi_2(i,j)   ! extinction coefficient 
   f_shade = 0.1_DBL  ! fraction of radiation reaching shaded leaves (empirical)
! Partition LAI into sunlit and shaded components
   lai_sun = (1.0_DBL - exp(-k_beer * LAI(i,j))) / k_beer
   lai_shade = LAI(i,j) - lai_sun
   lai_sun = max(1.e-6_DBL, lai_sun)
   lai_shade = max(0._DBL, lai_shade)
! Radiation reaching each
   sw_sun = swdsxy(i,j)
   sw_shade = f_shade * swdsxy(i,j)
! Compute radiation factors separately
   tmp_sun   = 0.55_DBL * sw_sun   * 2.0_DBL / Rgl(i,j) / lai_sun
   tmp_shade = 0.55_DBL * sw_shade * 2.0_DBL / Rgl(i,j) / max(1.e-6_DBL, lai_shade)
   rc_fac_sun   = (Rc_min(i,j)/Rc_max + tmp_sun) / (1.0_DBL + tmp_sun)
   rc_fac_shade = (Rc_min(i,j)/Rc_max + tmp_shade) / (1.0_DBL + tmp_shade)
! Combine weighted by LAI
   rc_fac_rad = (lai_sun * rc_fac_sun + lai_shade * rc_fac_shade) / LAI(i,j)


   ! vapor pressure deficit factor
   e_cas = q_cas(i,j)*pressf/(0.622+0.388*q_cas(i,j)) ! vapor pressure in hPa
   rc_fac_vpd = exp(-hs_rc(i,j)*(esatw(real(t_cas(i,j)))-e_cas))
   
   ! temperature factor
   rc_fac_t = max(0.,1._DBL-0.0016_DBL*(T_opt-t_cas(i,j))**2)
   
   ! rootzone soil moisture factor
   rc_fac_sw = 0._DBL
   d_root = 0.
   DO k = 1,nsoil
     if(rootF(k,i,j).gt.0.0_DBL) then ! soil layer with the root
      if(soilw(k,i,j).gt.w_s_FC(k,i,j)) then
       ! no water stree if water level exceed the value at field capacity
       temp = 1.*s_depth(k,i,j) * rootF(k,i,j)
       d_root = d_root + s_depth(k,i,j) * rootF(k,i,j)
      else if(soilw(k,i,j).lt.w_s_WP(k,i,j)) then
      ! below wilting point
       temp = 0.0_DBL
      else   ! otherwise
       temp = s_depth(k,i,j) & 
              *(soilw(k,i,j)*poro_soil(k,i,j)-theta_WP(k,i,j)) &
              /(theta_FC(k,i,j)-theta_WP(k,i,j)) * rootF(k,i,j)
       d_root = d_root + s_depth(k,i,j) * rootF(k,i,j)
      end if
      rc_fac_sw = rc_fac_sw + temp
!   if(swdsxy(i,j).gt.900.) then
!      print*,'>>>',k,rootF(k,i,j),soilw(k,i,j),w_s_FC(k,i,j),d_root,temp,rc_fac_sw
!   end if
     end if
   END DO
   rc_fac_sw = rc_fac_sw/max(1.e-6,d_root)
  
   tmp2 = max(1.e-6_DBL,rc_fac_rad*rc_fac_vpd*rc_fac_t*rc_fac_sw)
   r_c = min(Rc_max,(Rc_min(i,j)/LAI(i,j)/tmp2)) 
!   if(swdsxy(i,j).gt.900.) then
!      print*,'>>>>>>',rc_fac_rad ,rc_fac_vpd,rc_fac_t,rc_fac_sw,tmp2,r_c 
!      stop
!   end if
 !=================================================================================
 ! r_litter : litter resistance : not in use in this version
 !=================================================================================
 ! Ref.[Sakaguchi and Zeng 2009]
 ! set litter LAI as 0.5
 !r_litter = 1.0_DBL/0.004_DBL/ustar*(1._DBL-exp(-0.5_DBL))		
  r_litter = 0.  

else ! not vegetated land surfaces

 r_d = 0. 
 r_c = 0.
 r_b = 0.
 r_litter = 0.  

end if ! vegetated

! overwrite r_d if snow is present: 
! Marat, March 2019

if(snow_mass(i,j).gt.0.) then
 r_d = max(r_d,10000.*snow_mass(i,j)/rho_snow)  
end if

END SUBROUTINE resistances

