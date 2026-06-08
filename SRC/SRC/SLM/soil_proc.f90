!=================================================================================================
! module soil_proc
! Note: module contains subroutines for soil temperature and moisture calculations
!=================================================================================================

MODULE soil_proc
use slm_vars, only : DBL,soilw, soilw_inc, soilt, precip_in, precip_sfc, ks, run_off_sfc,nsoil,&
                     sh_eff_cond, s_depth_mm,Bconst,m_pot_sat, poro_soil, sh_eff_vel, evapo_dry, &
                     dtn_dbl,sst_cond, sst_capa,st_cond, st_capa, s_depth, st_eff_cond, alpha, beta, &
                     evp_soil, rootF, lcond, lsub, grflux0, w_s_WP, dtfactor_dbl, IMPERV, &
                     icemask, seaicemask, tfriz, tfrizs, &
                     drainage_flux, mws, mws_mx, snow_mass, snowt, &
                     cond_water, cond_ice, cond_snow, capa_water,  &
                     capa_ice, rho_water, rho_ice, rho_snow, w_s_FC, landtype, dorunoff
use params, only: lfus
IMPLICIT NONE
INTEGER :: k
REAL (KIND=DBL) :: aa, bb, cc, dd
PRIVATE :: k, aa, bb, cc, dd
PUBLIC :: soil_water, soil_temperature

CONTAINS

SUBROUTINE soil_water(tr,i,j)

REAL,INTENT(IN) :: tr ! reference level temperature 
INTEGER, INTENT(IN) :: i,j
REAL(KIND=DBL) :: sw_wgt(nsoil), drain, capa_soil, capa, deltas, evap


if(icemask(i,j).eq.1) then

  soilw(:,i,j) = 0.
  soilw_inc(:) = 0.
  precip_in = 0.
  run_off_sfc = 0.
  drainage_flux = 0.
  sh_eff_cond(:) = 0.
  sh_eff_vel(:) = 0.

else

  soilw_inc(:) = soilw(:,i,j)

! Snow melting parameterization:
! Marat, Apr 2017
! see if some snow should be melted to keep the top soil level at freezing temparture:
  if(snow_mass(i,j).gt.0..and.soilt(1,i,j).gt.tfriz) then
    
!   heat_capacity of first soil layer and snow:
    capa = ((1.-poro_soil(1,i,j))*sst_capa(1,i,j) &
          + (917.0_DBL*2030.0_DBL)*soilw(1,i,j)*poro_soil(1,i,j)) * s_depth(1,i,j)
! compute mass of snow (mm) melted to bring top soil temperature  to tfriz
    deltas = capa*(soilt(1,i,j)-tfriz)/lfus
    if(deltas.lt.snow_mass(i,j)) then
! not all snow is melted:
     snow_mass(i,j) = snow_mass(i,j) - deltas 
     mws(i,j) = mws(i,j) + deltas
     soilt(1,i,j) = tfriz
    else
! all remainig snow is melted:
     soilt(1,i,j) = soilt(1,i,j) - lfus*snow_mass(i,j)/capa
     mws(i,j) = mws(i,j) + snow_mass(i,j)
     snow_mass(i,j) = 0.
    end if
  end if

! Add some snow if appropriate:

  if(tr.lt.tfriz.and. &
    (snow_mass(i,j).gt.0..or.snow_mass(i,j).eq.0.and.soilt(1,i,j).lt.tfriz)) then
        snow_mass(i,j) = snow_mass(i,j) + precip_sfc*dtn_dbl
        precip_sfc = 0.
        snow_mass(i,j) = snow_mass(i,j) + mws(i,j)
        mws(i,j) = 0.
  end if
        
! Evaporate/vapor-deposite snow:
! MK Dec, 2019
  if(snow_mass(i,j).gt.0.) then
    snow_mass(i,j) = max(0.,snow_mass(i,j) - evp_soil(i,j)*dtn_dbl)
  else
! Evaporate standing water:
! MK Dec, 2019
    mws(i,j) = max(0.,mws(i,j) - evp_soil(i,j)*dtn_dbl)
  end if

!=================================================================================================
! Calculate precipitation infiltration rate into the first soil layer
! Note:
! MK 2025: Infiltration occurs only if any soil layer is unsaturated and ground temperature
! is above freezing point (before it was only top soil layer)
! Max. infiltration rate is equal to top soil hydraulic conductivity @ sat.(ks).
! Marat: modification Feb 2017:
! before running off, the excess of water accumulates in puddle until maximum
! amount is reached.
! puddle water continues to infiltrate even  when there is no rain.
! Also use puddle parameterization for simple snow model: Marat, Apr 2017
!=================================================================================================
  IF(any(soilw(:,i,j).lt.1._DBL).and.(soilt(1,i,j).ge.tfriz)) THEN
      precip_in = (1.-IMPERV(i,j))*MIN(precip_sfc+mws(i,j)/dtn_dbl, ks(1,i,j))  ! unit : kg/m2/s
  ELSE
      precip_in = 0.0_DBL
  END IF
!=================================================================================================
! calculate diffusion coefficient and velocity for soil moisture transfer 
! at each interfacial layer(= between adjacent soil layers)
! Note:
! depth-weighed averaging method is used
! sh_eff_cond : diffusion coefficient in diffusion equation ; mm*mm/s
! sh_eff_vel : effective velocity in advection term at interface layer : mm/s
! ks = hydraulic conductivity at saturation
! Bconst  = constant
! m_pot_sat = moisture potential at saturation
! soilw = soil wetness		  
! poro_soil = porosity		
! above parameter values are calculated based on percentage of SAND and CLAY as specified in input file
!=================================================================================================
  DO k = 1, nsoil-1
    if(soilt(k,i,j).ge.tfriz.and.soilt(k+1,i,j).ge.tfriz) then
      sh_eff_cond(k) &
      = (s_depth_mm(k,i,j)*(soilw(k,i,j)**(Bconst(k,i,j)+2._DBL))&
      +s_depth_mm(k+1,i,j)*(soilw(k+1,i,j)**(Bconst(k,i,j)+2._DBL)))&
      /(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) &
        *ks(k,i,j)*Bconst(k,i,j)*ABS(m_pot_sat(k,i,j))/poro_soil(k,i,j)

      sh_eff_vel(k) = &
         (s_depth_mm(k,i,j)*(soilw(k,i,j)**(2._DBL*Bconst(k,i,j)+2._DBL))&
         +s_depth_mm(k+1,i,j)*(soilw(k+1,i,j)**(2._DBL*Bconst(k,i,j)+2._DBL)))&
        /(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) &
        *ks(k,i,j)/poro_soil(k,i,j) 
    else ! no water movment between two layers one of which is frozen
      sh_eff_cond(k) = 0.
      sh_eff_vel(k) = 0.
    end if
  END DO

!=================================================================================================
! from FDE by the implicit method, Thomas algorithm is applied
! aa: terms related with soilw(k-1)
! bb: terms related with soilw(k)
! cc: terms related with soilw(k+1)
! dd: current soil wetness - sink + source
!=================================================================================================

! make saturated soil for wetlands: -MK
  if(landtype(i,j).eq.11) soilw(:,i,j) = w_s_FC(:,i,j)

  sw_wgt = 1.
  where(soilw(:,i,j).lt.w_s_WP(:,i,j)) sw_wgt = 0.

  ! no evaporation and precipitation infiltration from the soil itself if there is snow on in
  ! MK Dec, 2019
  ! Also, no soil top layer evaporation if there is standing water on the soil surface -MK Feb, 2025
  if(mws(i,j).gt.0.or.snow_mass(i,j).gt.0..and.soilt(1,i,j).lt.tfriz) then
   evap = 0.
  else
   evap = 1.
  end if

  aa = 0._DBL
  cc = -1.*2.*sh_eff_cond(1)*dtn_dbl/(s_depth_mm(1,i,j)*(s_depth_mm(1,i,j)+s_depth_mm(2,i,j)))  
  dd = soilw(1,i,j)& ! soil wetness at time step n
   -(max(0.,evp_soil(i,j))*evap & ! mm/s ; soil top evaporation
   +rootF(1,i,j)*evapo_dry*sw_wgt(1) & ! mm/s ; transpiration
   -precip_in) & ! mm/s 		! precipitation infiltration
!   -precip_in*evap) & ! mm/s 		! bug: MK 2025 precipitation infiltration
   *dtn_dbl/poro_soil(1,i,j)/s_depth_mm(1,i,j) ! s/mm

  bb = 1._DBL-cc + sh_eff_vel(1)*dtn_dbl/(s_depth_mm(1,i,j))
  alpha(1) = cc/bb
  beta(1) = dd/bb

! for 2-(nsoil-1) layer
  DO k = 2, nsoil-1
    aa = -1.*dtn_dbl/s_depth_mm(k,i,j) & ! s/mm
       * (sh_eff_vel(k-1) & ! mm/s
       + sh_eff_cond(k-1)*2._DBL/(s_depth_mm(k-1,i,j)+s_depth_mm(k,i,j))) ! mm/s
    cc = -1.*dtn_dbl/s_depth_mm(k,i,j) & ! s/mm
        * 2._DBL*sh_eff_cond(k)/(s_depth_mm(k,i,j)+s_depth_mm(k+1,i,j)) ! mm/s
    bb = 1._DBL-cc &
      + dtn_dbl/s_depth_mm(k,i,j)*sh_eff_vel(k) &
      + 2._DBL*dtn_dbl/s_depth_mm(k,i,j)*sh_eff_cond(k-1)/(s_depth_mm(k-1,i,j)+s_depth_mm(k,i,j))
    dd = soilw(k,i,j) & ! current time step
     - rootF(k,i,j)*evapo_dry*sw_wgt(k) & ! tranpiration mm/s
     *dtn_dbl/poro_soil(k,i,j)/s_depth_mm(k,i,j)  ! s/mm

    alpha(k) = cc/(bb-aa*alpha(k-1))
    beta(k) = (dd-aa*beta(k-1))/(bb-aa*alpha(k-1))
  END DO

  ! for bottom layer
  aa = -1.*dtn_dbl/s_depth_mm(nsoil,i,j) & ! s/mm
     * (sh_eff_vel(nsoil-1) & ! mm/s
     + sh_eff_cond(nsoil-1)*2._DBL/(s_depth_mm(nsoil-1,i,j)+s_depth_mm(nsoil,i,j))) ! mm/s
  cc = 0._DBL
  bb = 1._DBL &
     + (dtn_dbl/s_depth_mm(nsoil,i,j) &  ! s/mm
     * 2._DBL*sh_eff_cond(nsoil-1)/(s_depth_mm(nsoil-1,i,j)+s_depth_mm(nsoil,i,j))) ! mm/s

! drainge when it exceed 1. mm/s
  drainage_flux = max(soilw(nsoil,i,j)-1._DBL, 0._DBL)*poro_soil(nsoil,i,j)*s_depth_mm(nsoil,i,j)/dtn_dbl  
  dd = soilw(nsoil,i,j) & !current time step
     - (rootF(nsoil,i,j)*evapo_dry*sw_wgt(nsoil) & ! tranpiration mm/s
     + drainage_flux) & ! drainge when it exceed 1. mm/s
     *dtn_dbl/poro_soil(nsoil,i,j)/s_depth_mm(nsoil,i,j)  ! s/mm
  alpha(nsoil) = 0._DBL
  beta(nsoil) = (dd-aa*beta(nsoil-1))/(bb-aa*alpha(nsoil-1))

  ! (n+1) time step soil wetness
  soilw(nsoil,i,j) = beta(nsoil)
  DO k = nsoil-1, 1, -1
    soilw(k,i,j) = max(0.,beta(k)-alpha(k)*soilw(k+1,i,j))
  END DO 

  dd = 0.
  if(any(soilw(:,i,j).gt.1._DBL)) then
  ! fix the levels where wetness exceeds 1 preserving total water:
  ! compute the excess:
    do k=1,nsoil
     if(soilw(k,i,j).gt.1._dbl) then
       dd = dd + (soilw(k,i,j)-1._DBL)*poro_soil(k,i,j)*s_depth_mm(k,i,j) 
       soilw(k,i,j) = 1._DBL 
     end if
    end do
  ! distribute  the excess among layer into deepest layers first:
  !  do k=nsoil,1,-1
    do k=1,nsoil ! change the order for rain infiltration from the top when - MK 2025
     if(soilw(k,i,j).lt.1._DBL) then
       cc = min(dd,(1.-soilw(k,i,j))*poro_soil(k,i,j)*s_depth_mm(k,i,j))
       soilw(k,i,j) = soilw(k,i,j) + cc/(poro_soil(k,i,j)*s_depth_mm(k,i,j))
       dd = dd  - cc
       if(dd.lt.0._DBL) exit
     end if
    end do
  end if

! fixing the issue with rain infiltration even if soil is completely saturated
! move all the access to mws and modify the precip_in accordinally

  if(dd.gt.0.) then ! still some water left after saturating all the soil layers
                    ! modify the precip_in so that the access of water is moved
                    ! to the surface water
    precip_in = precip_in-dd/dtn_dbl
  end if

  mws(i,j) = max(0.,mws(i,j) + (precip_sfc - precip_in)*dtn_dbl)
  
  if(.not.dorunoff.and.mws(i,j).gt.mws_mx(i,j)) then
     drain = (mws(i,j)-mws_mx(i,j))/dtn_dbl
     mws(i,j) = mws_mx(i,j)
  else
     drain = 0.
  end if

  run_off_sfc = drain

  soilw_inc(:) = soilw(:,i,j)-soilw_inc(:)


end if ! icemask

END SUBROUTINE soil_water

!=============================================================

SUBROUTINE soil_temperature(i,j)

use vars, only: landmask
use grid, only: rank
INTEGER, INTENT(IN) :: i,j
! List of Local variables
REAL (KIND=DBL) :: k_dry, k_sat, Ke
REAL (KIND=DBL) :: temp
REAL(KIND=DBL) :: snow_depth, st_eff_snow, snow_capa
real, parameter :: snow_min_depth = 0.01
integer it,jt
real deltas

if(icemask(i,j).eq.1) then

  DO k = 1,nsoil
     st_cond(k,i,j) = cond_ice ! thermal conductivity of ice
     st_capa(k,i,j) = rho_ice*capa_ice ! ice heat capacity
  END DO


else


  DO k = 1,nsoil
  ! dry soil density [kg/m3] ! 2700kg/m3 = soilds unit weight 
  ! Empirical method from FIFE site observation
  temp = (1.0_DBL-poro_soil(k,i,j))*2700.0_DBL 
 
  ! dry thermal conductivity [W/mK]
  k_dry = (0.135_DBL*temp+64.7_DBL)/(2700._DBL-0.947_DBL*temp) 

  ! Saturated soil thermal conductivity : unsing percent of sand and clay content
  if(soilt(k,i,j).ge.tfriz) then
      k_sat = ((sst_cond(k,i,j))**(1.0_DBL-poro_soil(k,i,j))) &
             *(cond_water**poro_soil(k,i,j)) 
  else
      k_sat = ((sst_cond(k,i,j))**(1.0_DBL-poro_soil(k,i,j))) &
             *(cond_ice**poro_soil(k,i,j)) 
  end if

  ! Weighing factor between dry and saturated soil thermal conductivity : Kersten number
  Ke = log10(max(0.1_DBL,soilw(k,i,j)))+1.0_DBL

  ! Total soil thermal conductivity at each node_z
  ! Johansen (1975) 
  st_cond(k,i,j) = Ke*(k_sat - k_dry) + k_dry

  !soil volumetric heat capacity at each node_z depth
  if(soilt(k,i,j).ge.tfriz) then
     st_capa(k,i,j) = (1.-poro_soil(k,i,j))*sst_capa(k,i,j) &  
          + (rho_water*capa_water)*soilw(k,i,j)*poro_soil(k,i,j) ! water heat capacity
  else
     st_capa(k,i,j) = (1.-poro_soil(k,i,j))*sst_capa(k,i,j) &  
          + (rho_ice*capa_ice)*soilw(k,i,j)*poro_soil(k,i,j) ! ice heat capacity
  end if
  END DO

end if

DO k = 1,nsoil-1
! Calculate effective conductiviy at the adjacent soil layer interface [at interface_z]
! depth weighed average
   st_eff_cond(k) = (st_cond(k+1,i,j)*s_depth(k+1,i,j) + st_cond(k,i,j)*s_depth(k,i,j))/ &
                    (s_depth(k+1,i,j)+s_depth(k,i,j))
END DO! k

! from FDE by the implicit method, Thomas algorithm is applied
! aa: terms related with T(k-1)
! bb: terms related with T(k)
! cc: terms related with T(k+1)
! dd: current soil temperature + (additional source/sink on soil top)
! for first layer
if(snow_mass(i,j).gt.rho_snow*snow_min_depth) then
  snow_capa = rho_snow*capa_ice  ! snow heat capacity
  snow_depth = snow_mass(i,j)/rho_snow
!  st_eff_snow = (st_cond(1,i,j)*s_depth(1,i,j) + snow_cond*snow_depth)/(s_depth(1,i,j)+snow_depth)
  st_eff_snow = cond_snow   ! conductivity over snow/soil interface is for snow - MK
  aa = 0._DBL
  cc = -1.*2.*st_eff_snow*dtn_dbl/snow_capa/(snow_depth*(snow_depth+s_depth(1,i,j)))
  bb = 1.-cc
  dd = snowt(i,j)&  ! snow temperature at time step n
      -grflux0*dtn_dbl/snow_depth/snow_capa ! soil top boundary layer contidion
  alpha(0) = cc/bb
  beta(0) = dd/bb
  aa = -2._DBL*st_eff_snow*dtn_dbl/s_depth(1,i,j)/st_capa(1,i,j)/(snow_depth+s_depth(1,i,j))
  cc = -2._DBL*st_eff_cond(1)*dtn_dbl/s_depth(1,i,j)/st_capa(1,i,j)/(s_depth(1,i,j)+s_depth(2,i,j))
  bb = 1._DBL-cc-aa
  dd = soilt(1,i,j) ! current time step
  alpha(1) = cc/(bb-aa*alpha(0))
  beta(1) = (dd-aa*beta(0))/(bb-aa*alpha(0))
else
  aa = 0._DBL
  cc = -1.*2.*st_eff_cond(1)*dtn_dbl/st_capa(1,i,j)/(s_depth(1,i,j)*(s_depth(1,i,j)+s_depth(2,i,j)))
  bb = 1.-cc 
  dd = soilt(1,i,j)&  ! soil temperature at time step n
      -grflux0*dtn_dbl/s_depth(1,i,j)/st_capa(1,i,j) ! soil top boundary layer contidion
  alpha(1) = cc/bb
  beta(1) = dd/bb
end if


! for 2-(nsoil-1) layer
DO k = 2, nsoil-1
  aa = -2._DBL*st_eff_cond(k-1)*dtn_dbl/s_depth(k,i,j)/st_capa(k,i,j)/(s_depth(k-1,i,j)+s_depth(k,i,j))
  cc = -2._DBL*st_eff_cond(k)*dtn_dbl/s_depth(k,i,j)/st_capa(k,i,j)/(s_depth(k,i,j)+s_depth(k+1,i,j))
  bb = 1._DBL-cc-aa 
  dd = soilt(k,i,j) ! current time step
  alpha(k) = cc/(bb-aa*alpha(k-1))
  beta(k) = (dd-aa*beta(k-1))/(bb-aa*alpha(k-1))
END DO
! for bottom layer
aa = -2._DBL*st_eff_cond(nsoil-1)*dtn_dbl/s_depth(nsoil,i,j)/st_capa(nsoil,i,j)/ &
     (s_depth(nsoil-1,i,j)+s_depth(nsoil,i,j))
cc = 0._DBL
bb = 1._DBL-aa 
dd = soilt(nsoil,i,j)

alpha(nsoil) = 0._DBL
if(seaicemask(i,j).eq.0) then
  beta(nsoil) = (dd-aa*beta(nsoil-1))/(bb-aa*alpha(nsoil-1))
else
  beta(nsoil) = tfrizs ! lowest ice layer is in contact with seawater
end if

!call task_rank_to_index(rank,it,jt)
!if(i+it.eq.7068.and.j+jt.eq.3630) then
!  print*,'>>>',rank,i,j,i+it,j+jt
!  print*,'seaicemask=',seaicemask(i,j),icemask(i,j),landmask(i,j),grflux0
!  do k=1,nsoil
!    print*,k,soilt(k,i,j)
!  end do
!end if
! (n+1) time step soil temperature

soilt(nsoil,i,j) = beta(nsoil)
DO k = nsoil-1, 1, -1
  soilt(k,i,j) = beta(k)-alpha(k)*soilt(k+1,i,j)
END DO 

if(icemask(i,j).eq.1) then
  DO k=1,nsoil
    soilt(k,i,j) = min(tfriz,soilt(k,i,j))
  end do
end if

if(snow_mass(i,j).gt.rho_snow*snow_min_depth) then
 snowt(i,j) = beta(0)-alpha(0)*soilt(1,i,j)
else
 snowt(i,j) = min(tfriz,soilt(1,i,j))
end if

!melt some snow if it's temperature is above freezing: (MK 2020)

if(snowt(i,j) .gt.tfriz) then

    deltas = snow_capa*snow_depth*(snowt(i,j)-tfriz)/lfus
    if(deltas.lt.snow_mass(i,j)) then
! not all snow is melted:
     snow_mass(i,j) = snow_mass(i,j) - deltas
     mws(i,j) = mws(i,j) + deltas
     snowt(i,j) = tfriz
    else
! all remainig snow is melted:
     snowt(i,j) = tfriz
     mws(i,j) = mws(i,j) + snow_mass(i,j)
     snow_mass(i,j) = 0.
    end if
end if

if(landtype(i,j).eq.11) mws(i,j) = 0. ! for wetland, surface water is already implied. -MK

end subroutine soil_temperature
 
END MODULE soil_proc
