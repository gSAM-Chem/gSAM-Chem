!=============================================================================================
! subroutine run_slm 
! Note:
!	main() for SLM
!	called from SAM/SRC/surface()
!
! History: March, 2016
!=============================================================================================

SUBROUTINE run_slm(ur, vr, precip_ref, qr, tr, ts, zr, pref, pressf, rhosf, &
                   swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, coszrsxy, &
                   prsfc, flbu, flbv, flbq, flbt, flbqe, ra)

use slm_vars
use soil_proc, only : soil_water, soil_temperature
use vars, only: u10m_xy, v10m_xy, t2m_xy, q2m_xy, u10m, v10m, t2m, q2m, & 
                u10ma_xy, v10ma_xy, gust10m_xy, soilflux, dtfactor
use params, only: donorestart
use runoff, only: move_water


IMPLICIT NONE
!=============================================================================================
! Input : variables at reference level
!=============================================================================================
REAL , INTENT(IN) , DIMENSION(nx,ny) :: &
        ur, & ! u-wind
        vr, & ! v-wind
        qr, & ! qv at reference level
        tr, & ! t  at reference level
        precip_ref, & ! precipitation rate from reference level, mm/s
        swdsvisxy, & ! downward shortwave direct visible flux
        swdsnirxy, & ! downward shortwave direct near-IR flux
        swdsvisdxy, & ! downward shortwave diffuse visible flux
        swdsnirdxy, & ! downward shortwave idiffuse near-IR flux
        lwdsxy, & ! downward longwave flux
        coszrsxy, &  ! zenith angle
        zr, &  ! height of ref. level
        pref, &   ! pressure at ref. level
        pressf, &   ! surface pressure 
                                        rhosf   ! surface air density
!==============================================================================================
! Output : 
!==============================================================================================
REAL , INTENT(OUT), DIMENSION(nx,ny) :: &
        flbu, & ! surface u-momentum
        flbv, & ! surface v-momentum
        flbq, & ! surface vapor flux (m/s)
        flbqe, & ! surface latent heat flux (W/m2)
        flbt, & ! surface heat flux (Km/s)
        ra,   & ! surface scalar flux resistance
        prsfc   ! precipitation rate reaching to soil surface

!====================================================================================
! INPUT & OUTPUT : surface temperature 
!====================================================================================
REAL , INTENT(INOUT), DIMENSION(nx,ny) :: ts  ! surface temperature

!====================================================================================
! Local Variables
!====================================================================================
REAL (KIND=DBL) :: t_sfc, q_sfc,cp_vege_tot, coef
REAL (KIND=DBL) :: sdew
INTEGER :: i,j,k
integer iter,niter
REAL (KIND=DBL) :: dtn_iter, evapo_dry0, evapo_wet0, drain0, shf0, lhf0, evp0, tr_pot,tsfc_pot

real, external :: qsatw,qsati
integer it,jt
!====================================================================================
! note:
! Downward longwave (lwdsxy_slm), shortwave radiation (swdsxy_slm), 
! and cosine of zenith angle (coszrs_slm) are read from the equivalent SAM variables 
! at each time step.
! When the run restarts, restarts swdsxy_slm, lwdsxy_slm, coszrs_slm are overwritten 
! from the values stored in the LAND restart file. 
! Otherwise values will be zero as radiation is calculated after surface() in SAM.
!====================================================================================

if(flag_vege_init.and. (nrestart.eq.0)) then  
  ! initialize vegetation and soil parameters 
    if(masterproc) write(*,*) 'initializing SLM canopy if available ......'
    call init_slm_vars(tr, ts, qr, pressf)  
    if(masterproc) write(*,*) 'canopy initialized'
    flag_vege_init = .false.
end if

dtn_dbl = dtn

call task_rank_to_index(rank,it,jt)

!====================================================================================
! local variable initialization
cnp_mw_drip = 0._DBL ! dripping from canopy water storage when the storage exceeds mw_mx : mm/s
drain = 0._DBL ! drainige rate from canopy

DO j = 1,ny
DO i = 1,nx

IF(landmask(i,j).eq.1.or.seaicemask(i,j).eq.1) THEN ! not ocean    

!====================================================================================
! Note:
! Nudging soil wetness and soil temperature profile to initial over time period, tau_soil.
! soil_relax_hgt (in the range of 0 - 1) variable determines the effectiveness of nudging as set in the 
! initial soil sounding. 
!====================================================================================
  if(dosoilwnudging) then
        DO k = 1, nsoil
        soilw(k,i,j) = soilw(k,i,j)-(soilw(k,i,j)-soilw_obs(k,i,j))*dtn_dbl/tau_soil(i,j)*soil_relax_hgt(k,i,j)
        soilw_nudge(k,i,j) = (soilw(k,i,j)-soilw_obs(k,i,j))*soil_relax_hgt(k,i,j)/tau_soil(i,j)
        END DO ! k
  end if

  if(dosoiltnudging) then
        DO k = 1, nsoil
        soilt(k,i,j) = soilt(k,i,j)-(soilt(k,i,j)-soilt_obs(k,i,j))*dtn_dbl/tau_soil(i,j)*soil_relax_hgt(k,i,j)
        soilt_nudge(k,i,j) = (soilt(k,i,j)-soilt_obs(k,i,j))*soil_relax_hgt(k,i,j)/tau_soil(i,j)
        END DO ! k
  end if

!====================================================================================
! Note: 
! Calculate net radiation absorbed by canopy and soil surface  
!====================================================================================
   call radiative_fluxes(i,j)

!====================================================================================
! Note:
! precip 	: precipitation interception rate at canoppy
! precip_ref 	: precipiptation rate at the reference height, zr
! precip_extinc : extinction coefficient for rainfall penetration through 
!		  vegetation layer
! LAI 		: leaf area index
! For baresoil, precip = 0, as LAI = 0
!====================================================================================
   precip = precip_ref(i,j)*(1.-exp(-1.*precip_extinc(i,j)*LAI(i,j)))
!====================================================================================
! Note:
! drain 	: drainage rate from canopy 
! cnp_mw_drip 	: dripping rate from canopy
! precip_sfc 	: precipitation reaching to soil surface
! mw 		: water storage of canopy [kg/m^2]
! mw_mx 	: maximum water storage of canopy
! mw_inc 	: (prognostic variable) increment of moisture storage on canopy
!====================================================================================
   if(mw(i,j).lt.mw_mx(i,j)) then  
       precip = min((mw_mx(i,j)-mw(i,j))/dtn_DBL,precip)
       drain =  0.
   else
       ! When water holding storage exceeds its maximum, 
       ! no precipiataion is intercepted
       drain = precip + (mw(i,j)-mw_mx(i,j))/dtn_DBL
   end if
                
   precip_sfc = precip_ref(i,j)-precip+drain 

 ! canapy heat capacity: assume 0.001 m leaf thickness, basal area in sq.feet/acre=43560 m2/m2,
 ! 900 kg/m3 density of leaves and wood, 2800 J/kg/K specific heat capacity.
!  cp_vege(i,j) = LAI(i,j)*0.001*900.*2800.
  cp_vege(i,j) = (LAI(i,j)*leaf_thickness*0.001 + ztop(i,j)*BAI(i,j)/43560.)*900.*2800.
  cp_vege_tot = cp_vege(i,j) + mw(i,j)*1.e-3_DBL*cp_water

  
  mw_inc = dtn_dbl*(precip-drain)
  mw(i,j) = mw(i,j) + mw_inc

! intercepted precip cools the canopy:
! temperature of intercepted rain is the same as reference level
! note it cools even when water storage is full as old water on leaves is replaced by new rain water
! MK 2020
  if(vegetated(i,j)) then
   t_canop(i,j) = (cp_vege_tot*t_canop(i,j)+tr(i,j)*precip*dtn_dbl*1.e-3_DBL*cp_water)/ &
                  (cp_vege_tot+precip*dtn_dbl*1.e-3_DBL*cp_water)
  end if
!====================================================================================
! Note:
! Assign appropriate "surface level" values for each land type, for the calculation  
! of surface turbulent fluxes
! 
! For baresoil,   surface level = soil surface
! For vegetation, surface level = canopy level
!====================================================================================
! specific humidity at top soil
  if(icemask(i,j).eq.0.and.soilt(1,i,j).gt.tfriz) then
   if(mws(i,j).gt.0.) then
     q_gr = qsatw(real(soilt(1,i,j)),pressf(i,j))
     sdew = 1.
   else
     fh = fh_calc(soilt(1,i,j),m_pot_sat(1,i,j), soilw(1,i,j), Bconst(1,i,j)) 
     if(fh.gt.0.99) then
       sdew = 1.
     else
       sdew = 0.
     end if
     q_gr = fh*qsatw(real(soilt(1,i,j)),pressf(i,j))
   end if
  else
    if(snow_mass(i,j).gt.0.) then
     q_gr = qsati(real(snowt(i,j)),pressf(i,j))
!     q_gr = qsati(real(0.5*(snowt(i,j)+tr(i,j))),pressf(i,j))
    else
      q_gr = qsati(real(soilt(1,i,j)),pressf(i,j))
    end if
    sdew = 1.
  end if

  IF(vegetated(i,j)) THEN ! for vegetated surfaces
     t_sfc = t_cas(i,j)
     q_sfc = q_cas(i,j)
  ELSE    ! for baresoil
   ! correct for snow
   ! MK Dec 2019
    if(snow_mass(i,j).gt.0.) then
     t_sfc = snowt(i,j)
    else
     t_sfc = soilt(1,i,j)
    end if
     q_sfc = q_gr
  END IF

! when there is no dew formation (sdew=0), there is non-positive evaporation flux
! affecting other diagnostics like q at 2-m
  if (sdew.lt.0.1 .and. q_sfc.lt.qr(i,j)) q_sfc = qr(i,j) 

!====================================================================================
! Note:
! Calculate turbulent transfer coefficient between reference level and surface
! Input : pref  	reference level pressure
!	  tr    	reference level temperature
!	  qr    	reference level specific humidity
!	  ur,vr 	reference level zonal, meridional velocities
!	  zr		reference level height
!	  t_sfc 	surface level temperature
!	  q_sfc 	surface level specific humidity
!	  z0_sfc 	surface roughness length
!	  disp_hgt 	displacement height
!====================================================================================
  tsfc_pot = t_sfc*(1000./pressf(i,j))**(rgas/cp)
  tr_pot = tr(i,j)*(1000./pref(i,j))**(rgas/cp)

  call transfer_coef(tsfc_pot, tr_pot, qr(i,j), q_sfc,  &
                     ur(i,j), vr(i,j), zr(i,j), z0_sfc(i,j), disp_hgt(i,j),i,j)

  temp_2m = temp_2m*(pressf(i,j)/1000.)**(rgas/cp)  ! convert from pot temp to temp
!====================================================================================
! Note:
! Calculate surface momentum fluxes
! Surface Stress : kg/m/s2
!====================================================================================
  taux_sfc = -1._DBL*mom_trans_coef*vel_m*ur(i,j)*rhosf(i,j)
  tauy_sfc = -1._DBL*mom_trans_coef*vel_m*vr(i,j)*rhosf(i,j)

!====================================================================================
! Note:
! Calculate aerodynamic resistances + stomatal resistance
!====================================================================================
call resistances(pressf(i,j),i,j)

!====================================================================================
! Subcycle in time to avoid swings of canopy temperature
! because of small heat capacity of vegetation in curtain places and seasons
! Of course, the best would be to use implicit scheme, but for now,
! subcycling seems like a good fix.
! MK Jan 2020
!====================================================================================
if(vegetated(i,j)) then

  niter = max(1,nint(dtn_dbl/1._DBL))
  !niter = 1
  dtn_iter = dtn_dbl/niter
  shf0 = 0
  lhf0 = 0
  evp0 = 0
  evapo_dry0=0.
  evapo_wet0=0.
  drain0 = 0.

  do iter=1,niter

!====================================================================================
    shf_canop(i,j) =(t_canop(i,j)-t_sfc)*rhosf(i,j)*cp/r_b
    shf0 = shf0 + shf_canop(i,j)
!====================================================================================
! Note:
! Calculate latent heat fluxes
! Input: q_sfc          surface level specific humidity
!        qr             reference level specific humidity
!        i,j            current index of the grids
!====================================================================================
    call vapor_fluxes_canopy(q_sfc,pressf(i,j),rhosf(i,j),dtn_iter,i,j)
    lhf0 = lhf0 + lhf_canop(i,j)
    evp0 = evp0 + evp_canop(i,j)
    evapo_dry0 = evapo_dry0+evapo_dry
    evapo_wet0 = evapo_wet0+evapo_wet
!====================================================================================
! Note:
! Update vegetation moisture storage
!====================================================================================
    if(mw(i,j).gt.mw_mx(i,j)) then
     drain0 = drain0 + (mw(i,j)-mw_mx(i,j))/dtn_iter  ! dripping excess of dew
     mw(i,j) = mw_mx(i,j)
    end if
!====================================================================================
! Note:
! Update vegetatation temperature
! cp_vege_tot [J/m2/K]  : total vegetation heat capacity
!         = heat capacity from vegetation + heat capacity of water held on vegetation    
!====================================================================================
    cp_vege_tot = cp_vege(i,j) + mw(i,j)*1.e-3_DBL*cp_water
  
   ! max(1.e-3, vege_cp) is for the baresoil case where vege_cp = 0.(to prevent divide by zero)
   ! vege_YES keeps t_canop_inc zero for baresoil case
    t_canop_inc = dtn_iter/max(1.e-3_DBL,cp_vege_tot)&
                      *(net_rad(1)-shf_canop(i,j)-lhf_canop(i,j))*vege_YES(i,j)
    t_canop(i,j) = min(t_canop_max,t_canop(i,j) +  t_canop_inc)

  end do ! iter
  shf_canop(i,j) = shf0/real(niter)
  lhf_canop(i,j) = lhf0/real(niter)
  evp_canop(i,j) = evp0/real(niter)
  precip_sfc = precip_sfc + drain0/real(niter)
  drain = drain + drain0/real(niter)
  evapo_wet = evapo_wet0/real(niter)
  evapo_dry = evapo_dry0/real(niter)

!  if(lhf_canop(i,j).lt.-50000.) then
!    print*,'Landtype=',landtype(i,j),landtype(i-1,j),landtype(i,j-1)
!    print*,'LAI=',lai(i,j),lai(i-1,j),lai(i,j-1)
!    print*,'LHF_CANOP is insane! =',lhf_canop(i,j),lhf_canop(i-1,j),lhf_canop(i,j-1)
!    print*,'SHF =',shf_canop(i,j),shf_canop(i-1,j),shf_canop(i,j-1)
!    print*,'r_b=',r_b
!    print*,'r_c=',r_c
!    print*,'rank, i, j=',rank, i,j
!    print*,'t_canop=',t_canop(i,j),t_canop(i-1,j),t_canop(i,j-1)
!    print*,'tr=',tr(i,j)
!    print*,'ts=',ts(i,j)
!    print*,'qr=',qr(i,j)
!    print*,'t_cas=',t_cas(i,j)
!    print*,'q_cas=',q_cas(i,j)
!    print*,'t_sfc=',t_sfc
!    print*,'q_sfc=',q_sfc
!    print*,'t_soil=',soilt(1,i,j)
!    print*,'evapo_dry=',evapo_dry
!    print*,'evapo_wet=',evapo_wet
!    print*,'wet_canop=',wet_canop
!    print*,'precip_sfc=',precip_sfc
!    print*,'precip=',precip
!    print*,'precip_ref=',precip_ref(i,j)
!    print*,'drain=',drain
!    print*,'mw=',mw(i,j),'mw_mx=',mw_mx(i,j)
!    call task_abort()
!  end if

else
  shf_canop(i,j) = 0.
  lhf_canop(i,j) = 0.
  evp_canop(i,j) = 0.
  wet_canop = 0.
  mw(i,j) = 0.
  evapo_wet = 0.
  evapo_dry = 0.
  t_canop(i,j) = tr(i,j)
end if
!====================================================================================
! Sensible heat fluxes
!====================================================================================
! temperature of snow is average between air above and soil below
! MK Dec 2019

  IF(vegetated(i,j)) THEN 
     if(snow_mass(i,j).gt.0.) then
        shf_soil(i,j) = (snowt(i,j)-t_sfc)*rhosf(i,j)*cp/r_d
     else
        shf_soil(i,j) = (soilt(1,i,j)-t_sfc)*rhosf(i,j)*cp/r_d
     end if
     shf_air(i,j) = shf_canop(i,j)+shf_soil(i,j)
  ELSE
     shf_air(i,j) = (tsfc_pot-tr_pot)*rhosf(i,j)*cp/r_a
     shf_soil(i,j) = shf_air(i,j)
  END IF

  call vapor_fluxes_soil(q_sfc,qr(i,j),pressf(i,j),rhosf(i,j),sdew,i,j)
!====================================================================================
! Note:
! Calculate soil moisture increment
!====================================================================================
  call soil_water(tr(i,j),i,j)

!====================================================================================
! Note:
! Calculate soil temperature increment
! grflux0:	Total heat flux on top of soil surface : top boundary condition
! grflux0 < 0 => heat is added into soil top
!====================================================================================
  grflux0 = -1.0_DBL*(net_rad(2) - shf_soil(i,j) - lhf_soil(i,j))  
  call soil_temperature(i,j)

!====================================================================================
! Note:
! Compute diagnostic variables at canopy air space level
!====================================================================================

 ! Calculate heat conducatnces ; non-zero only for canopy land type
  if(vegetated(i,j)) then
    cond_heat = 1._DBL/r_a + 1._DBL/r_b + 1._DBL/r_d
    cond_href = 1._DBL/r_a/cond_heat
    cond_hcnp = 1._DBL/r_b/cond_heat
    cond_hundercnp = 1._DBL/r_d/cond_heat
  else
    cond_heat = 1._DBL/r_a
    cond_href = 1._DBL
    cond_hcnp = 0._DBL
    cond_hundercnp = 0._DBL
  end if

 ! Calculate Vapor conductances
  if(vegetated(i,j)) then
    cond_vapor = 1._DBL/r_a + wet_canop/(2.*r_b) + (1.-wet_canop)/(2.*r_b+r_c) &
                 + 1._DBL/(r_d+r_litter+r_soil)
    cond_vref  = 1._DBL/r_a/cond_vapor
    cond_vcnp  = (wet_canop/(2.*r_b) + (1.-wet_canop)/(2.*r_b+r_c))/cond_vapor
    cond_vundercnp = (1._DBL/(r_d+r_litter+r_soil))/cond_vapor
  else
    cond_vapor = 1._DBL/(r_a+r_soil)
    cond_vref  = 1._DBL
    cond_vcnp  = 0._DBL
    cond_vundercnp = 0._DBL
  end if

  if(snow_mass(i,j).eq.0.) then
    t_cas(i,j) = tr(i,j)*cond_href &
              + t_canop(i,j) *cond_hcnp &
              + soilt(1,i,j) *cond_hundercnp
  else
    t_cas(i,j) = tr(i,j)*cond_href &
              + t_canop(i,j) *cond_hcnp &
              + snowt(i,j) *cond_hundercnp
  end if
  if(icemask(i,j).eq.0.and.t_canop(i,j).ge.tfriz) then
    qsat_canop = qsatw(real(t_canop(i,j)),pressf(i,j))
  else
    qsat_canop = qsati(real(t_canop(i,j)),pressf(i,j))
  end if
  if(icemask(i,j).eq.0.) then
    if(snow_mass(i,j).eq.0.and.soilt(1,i,j).ge.tfriz) then
      fh = fh_calc(soilt(1,i,j),m_pot_sat(1,i,j),soilw(1,i,j),Bconst(1,i,j))
      q_gr = fh * qsatw(real(soilt(1,i,j)),pressf(i,j))
    else
      fh = 1.
      q_gr = qsati(real(snowt(i,j)),pressf(i,j))
    end if
  else
    fh = 1.
    q_gr = qsati(real(soilt(1,i,j)),pressf(i,j))
  end if
  q_cas(i,j) = qr(i,j)*cond_vref &
             + qsat_canop*cond_vcnp &
             + q_gr*cond_vundercnp

!  if(lhf_air(i,j).gt.8000.) then 
!    print*,'>>>','it+i jt+j=',i+it,j+jt,'q_cas=',q_cas(i,j),'qr=',qr(i,j),'cond_vref=',cond_vref, &
!   'qsat_canop=',qsat_canop,'cond_vcnp=',cond_vcnp,'q_gr=',q_gr,'cond_vundercnp=',cond_vundercnp, &
!    't_cas=',t_cas(i,j),'t_canop',t_canop(i,j),'soilt=',soilt(1,i,j),'pressf=',pressf(i,j), &
!    'r_a=',r_a,'r_b=',r_b,'r_c=',r_c,'r_d=',r_d,'LAI=',LAI(i,j),'wet_canop=',wet_canop,&
!    'cond_vapor=',cond_vapor,'mw=',mw(i,j),'shf_canop=',shf_canop(i,j),'lhf_canop=',lhf_canop(i,j), &
!    'shf_air=',shf_air(i,j),'lhf_air=',lhf_air(i,j),'shf_soil=',shf_soil(i,j),'lhf_soil=',lhf_soil(i,j), &
!    'cp_vege_tot=',cp_vege_tot,"vegetated=",vegetated(i,j),"evp_soil=",evp_soil(i,j)
!    stop
!  end if
!====================================================================================
! output variables
!====================================================================================
!  ts(i,j)   = real(t_sfc)
  ts(i,j)   = real(t_skin(i,j))
  flbu(i,j) = real(taux_sfc)/rhosf(i,j)  ! Surface Stress in x-direction 
  flbv(i,j) = real(tauy_sfc)/rhosf(i,j)  ! Surface Stress in y-direction
  flbq(i,j) = real(evp_air(i,j))/rhosf(i,j) ! latent heat flux [kg/kg m/s]
  flbqe(i,j) = real(lhf_air(i,j)) ! latent heat flux [W/m2]
  flbt(i,j) = real(shf_air(i,j))/(cp*rhosf(i,j))    ! sensible heat flux [W/m2] -> [Km/s]
  prsfc(i,j) = real(precip_sfc)
  ra(i,j) = real(r_a)
! diagnistic variables:

  u10m_xy(i,j) = u10m_xy(i,j) + u_10m * dtfactor
  v10m_xy(i,j) = v10m_xy(i,j) + v_10m * dtfactor
  u10ma_xy(i,j) = u10ma_xy(i,j) + u_10m * dtfactor
  v10ma_xy(i,j) = v10ma_xy(i,j) + v_10m * dtfactor
  gust10m_xy(i,j) = max(gust10m_xy(i,j), u_10m**2+v_10m**2)
  t2m_xy(i,j) = t2m_xy(i,j) + temp_2m * dtfactor
  q2m_xy(i,j) = q2m_xy(i,j) + q_2m * dtfactor
  u10m(i,j) = u_10m
  v10m(i,j) = v_10m
  t2m(i,j) = temp_2m
  q2m(i,j) = q_2m
  soilflux(i,j) = grflux0
!====================================================================================
! Collect statistics
!====================================================================================

 ! collect 2D stat variables
  call collect_2D_stat_vars(i,j)

!==================================================================================
else  ! ocean points

lhf_canop(i,j) = 0.
shf_canop(i,j) = 0.
lhf_soil(i,j) = 0.
shf_soil(i,j) = 0.
lhf_air(i,j) = 0.
shf_air(i,j) = 0.
soilflux(i,j) = 0.
evp_canop(i,j) = 0.
evp_soil(i,j) = 0.
evp_air(i,j) = 0.

END IF ! landmask.eq.1.ot.seaicemak.eq.1
END DO
END DO
 
if(dorunoff) call move_water()

s_precip_ref = s_precip_ref + precip_ref*dtfactor

! Save land variables to restart file
!if(mod(nstep,nstat*(1+nrestart_skip)).eq.0) then
!if(mod(nstep,nstat*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) then
!
!  if(.not.donorestart) call write_statement_slm() ! save restart file
!
!end if

END SUBROUTINE run_slm


