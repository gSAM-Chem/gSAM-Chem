!===============================================================================
!	subroutine radiative_fluxes
!	Note:
!	computes radiation flux transfer over land surace
!===============================================================================

SUBROUTINE radiative_fluxes(i, j)

use slm_vars, only : DBL,vege_YES, LAI_min, &
                 albedovis_v, albedonir_v, albedovis_s, albedonir_s, snow_mass, &
		 net_rad, net_sw, net_lw, net_swup, net_swdn, net_lwdn, net_lwup, t_skin, &
		 t_canop, soilt, IR_emis_vege, sigma, IR_emis_grnd, phi_1, phi_2, LAI, soilw, &
                 IR_emis_ice, IR_emis_snow, snowt, snow_mass, icemask, NSTEP
use params, only: dolongwave, doshortwave
use rad, only: swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy
IMPLICIT NONE

INTEGER,INTENT(IN) :: i, j ! incices
 
REAL (KIND=DBL) :: fdn1, & ! downwelling flux  on top
	           fdn2, & ! Transmitted flux (downwelling flux below)
	           fup1, & ! upwelling flux from top
	           fup2    ! upwelling flux from below
!===================================== 
! Note: 
! Definition of fdn1, fdn2, fup1, fup2
!=====================================
!	|                  ^
!      fdn1                |
!       |                 fup1   [fluxes above a certain layer, either canopy or topsoil]
!       V                  |
!====================================
!   Canopy / soil top
!====================================
!       |                  ^ 
!      fdn2                |     [fluxes below a certain layer, either canopy or topsoil]
!       |                 fup2         
!       V                  |
! Albedo is computed as alb = alb_v*(1-exp(-kLAI))+alb_s*exp(-kLAI)

REAL (KIND=DBL), DIMENSION(2) :: tir 	! Thermal IR from canopy, tir(1)
					! Thermal IR from soil top, tir(2)
REAL (KIND=DBL) :: ka ! optical depth
REAL (KIND=DBL) explai, explai0, wetfactor

!===============================================================================
! Compute shortwave radiation transfer between land surface and reference level 
! sum visible, NearIR for direct and diffuse radiation fluxes.
!===============================================================================

net_sw = 0.
net_swup = 0.
net_swdn = 0.
net_lw = 0.
net_lwup = 0.
net_lwdn = 0.
net_rad = 0.

if(snow_mass(i,j).gt.0.) then
  t_skin(i,j) = snowt(i,j) ! provisional 
else
  t_skin(i,j) = soilt(1,i,j) ! provisional 
end if

if(doshortwave.and.coszrsxy(i,j).gt.0.0_DBL) then

!  optical depth of the direct beam per unit leaf area 
  ka = phi_1(i,j)/max(0.01,coszrsxy(i,j)) + phi_2(i,j)
  explai = exp(-ka*LAI(i,j)) ! for direct radiation
!  optical depth of the diffuse  beam per unit leaf area 
  ka = phi_1(i,j) + phi_2(i,j)
  explai0 = exp(-ka*LAI(i,j)) !  for diffuse radiation
  
!  print*,'coszrsxy(i,j)=',coszrsxy(i,j)
!  print*,'explai=',explai
!  print*,'explai0=',explai0

!  net_rad(1) = net absorbed shortwave radiation by canopy

  net_rad(1) = net_rad(1)+swdsvisxy(i,j)*(1.-albedovis_v(i,j)*(1.-explai)-explai)
  net_rad(1) = net_rad(1)+swdsvisdxy(i,j)*(1.-albedovis_v(i,j)*(1.-explai0)-explai)
  net_rad(1) = net_rad(1)+swdsnirxy(i,j)*(1.-albedonir_v(i,j)*(1.-explai)-explai)
  net_rad(1) = net_rad(1)+swdsnirdxy(i,j)*(1.-albedonir_v(i,j)*(1.-explai0)-explai)

!  print*,'net_rad(1)=',net_rad(1)
  
  net_swdn(1) = swdsvisxy(i,j)+swdsvisdxy(i,j)+swdsnirxy(i,j)+swdsnirdxy(i,j)
  net_swup(1) = net_rad(1) - net_swdn(1)*(1.-explai)   

!  print*,'net_swdn=',net_swdn
!  print*,'net_swup=',net_swup
  
!  net_rad(2) = net absorbed shortwave radiation by soil

  if(snow_mass(i,j).eq.0) then
    wetfactor = 1.-0.5*soilw(1,i,j) ! soil wetnes factor: assume that wet soil is twice as dark
    net_rad(2) = net_rad(2)+swdsvisxy(i,j)*(1.-albedovis_s(i,j)*wetfactor)*explai
    net_rad(2) = net_rad(2)+swdsvisdxy(i,j)*(1.-albedovis_s(i,j)*wetfactor)*explai0
    net_rad(2) = net_rad(2)+swdsnirxy(i,j)*(1.-albedonir_s(i,j)*wetfactor)*explai
    net_rad(2) = net_rad(2)+swdsnirdxy(i,j)*(1.-albedonir_s(i,j)*wetfactor)*explai0
  else
    net_rad(2) = net_rad(2)+swdsvisxy(i,j)*(1.-0.75)*explai
    net_rad(2) = net_rad(2)+swdsvisdxy(i,j)*(1.-0.45)*explai0
    net_rad(2) = net_rad(2)+swdsnirxy(i,j)*(1.-0.75)*explai
    net_rad(2) = net_rad(2)+swdsnirdxy(i,j)*(1.-0.45)*explai0
  end if

!  print*,'net_rad(2)=',net_rad(2)

  net_swdn(2) = net_swdn(1)*explai
  net_swup(2) = net_rad(2) -net_swdn(2)

!  print*,'net_swdn(2)=',net_swdn(2)
!  print*,'net_swup(2)=',net_swup(2)

! Store net absorbed SW

  net_sw = net_rad

end if

if(dolongwave) then

!===================================================
! LongWave Radiation 
!===================================================

!===================================================
! tir: Emitted thermal infrared radiation
!===================================================
! Note:
! For no vegetation : 
! IR_trans becomes zero => tir(1) automatically becomes zero
!===================================================
tir(1) = IR_emis_vege(i,j)*sigma*(t_canop(i,j)**4)
if(snow_mass(i,j).gt.0.) then
  tir(2) = IR_emis_snow*sigma*(snowt(i,j)**4)
else
  if(icemask(i,j).eq.1.) then
    tir(2) = IR_emis_ice*sigma*(soilt(1,i,j)**4)
  else
    tir(2) = IR_emis_grnd(i,j)*sigma*(soilt(1,i,j)**4)
  end if
end if

!===================================================
! downwelling LW on canopy top: input
!===================================================
fdn1 = lwdsxy(i,j)
net_lwdn(1) = fdn1

!===================================================
! downwelling LW below canopy layer
!===================================================
! Note:
! below canopy layer : 
!  incoming LW (fdn1) that is not absorbed by canopy 
! + emitted thermal IR by canopy (tir(1)) toword soil surface)
! with no canopy : fdn2 is computed to be fdn1
! (1-IR_emis) = area of canopy gap (skyview factor)
!===================================================
fdn2 = (1.-IR_emis_vege(i,j))*fdn1 + tir(1)

net_lw(1) = fdn1-fdn2

!===================================================
! Note:
! At this stage, 
! net_rad(1,:,:) = net absorbed SW by canopy 
!		 + downwelling LW on canopy top
!		 - transmitted LW through canopy layer
!		 - emitted TIR toward soil surface
!===================================================
net_rad(1) = net_rad(1) + fdn1 - fdn2

!===================================================
! downwelling LW for soil surface
!===================================================
fdn1 = fdn2 
net_lwdn(2) = fdn1

! no fluxes below topsoil
fdn2 = 0.0_DBL
fup2 = 0.0_DBL

!===================================================
! Note:
! Emitted LW from topsoil = emitted tir from topsoil
!			  + portion of incoming LW that is reflected back toward canopy
!===================================================
fup1 = (tir(2)+(1.-IR_emis_grnd(i,j))*fdn1)
net_lwup(2) = fup1

!===================================================
! Note:
! At this stage, 
! net_rad(2,:,:) = net absorbed SW by soil surface 
!		 + net absorbed LW by soil surface
!===================================================
net_rad(2) = net_rad(2) + fdn1 - fdn2 - fup1 + fup2

! net_lw(2,:,:) = net absorbed LW by soil surface
net_lw(2) = fdn1-fdn2-fup1+fup2

!===================================================
! Incoming LW from below canopy
! Note:
! incoming LW from below canopy = upwelling flux at topsoil
!===================================================
fup2 = fup1

!===================================================
! Upwelling LW from canopy top
! Note:
! fup1 = portion of fup2 that is not absorbed + tir emitted from canopy
! for no canopy : fup1 = fup2
!===================================================
fup1 = (1.-IR_emis_vege(i,j))*fup2 + tir(1)

! total upward lw from surface (for canopy cover- from canopy top, for no canopy - from soil surface)
net_lwup(1) = fup1
t_skin(i,j) = (fup1/sigma)**0.25

net_lw(1) = net_lw(1) + fup2 - fup1

!if(i.eq.4.and.j.eq.4) then
!  print*,'swds:>>>',swdsvisxy(i,j),swdsvisdxy(i,j),swdsnirxy(i,j),swdsnirdxy(i,j),swdsxy(i,j)
!  print*,'swup/dn:>>>',net_swup,net_swdn
!  print*,'lwup/dn:>>>',net_lwup,net_lwdn
!  print*,'rad:>>>',' sw:',net_sw,'  lw:',net_lw
!end if
!===================================================
! Note:
! At this stage, 
! net_rad(1,:,:) = net absorbed SW by canopy 
!		 + downwelling LW on canopy top
!		 - transmitted LW through canopy layer (downward direction)
!		 - emitted TIR from canopy toward soil surface
!		 - emitted TIR from canopy toward atmosphere
!	         + upwelling LW from topsoil
!		 - transmitted LW through canopy layer (upward direction)	
!===================================================
net_rad(1) = net_rad(1) - fup1 + fup2

end if ! dolongwave

!energy_bal(i,j) = net_swdn(1,:,:)-net_swup(1,:,:)+net_lwdn(1,:,:)-net_lwup(1,:,:)
!energy_bal_c(i,j) = net_sw(1,:,:)+net_lw(1,:,:)
!energy_bal_s(i,j) = net_sw(2,:,:)+net_lw(2,:,:)

return
END SUBROUTINE radiative_fluxes

