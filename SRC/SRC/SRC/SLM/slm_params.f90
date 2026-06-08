! some constants and parameters used in SLM

module slm_params

INTEGER, PARAMETER :: DBL = 4 !SIZEOF(REAL)
!INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=12)
REAL (KIND=DBL), PARAMETER :: FC = -3.3_DBL ! field capacity ~-0.33bar=-3.3m
REAL (KIND=DBL), PARAMETER :: WP = -150._DBL ! wilting point
REAL (KIND=DBL), PARAMETER :: pii = 3.141592653589793 ! acos(-1)

! parameters for stomatal_resistance
REAL (KIND=DBL), PARAMETER :: Rc_max = 5000.0_DBL ! maximum stomatal resistance
REAL (KIND=DBL), PARAMETER :: T_opt = 298.0_DBL   ! optimum temperature for transpiration

! Heat capacity of water

REAL (KIND=DBL), PARAMETER :: sigma = 5.67e-8_DBL   ! Stefen- Boltzmann constant [W/m^2/K^4]

REAL (KIND=DBL), PARAMETER :: z0_soil = 0.005_DBL   ! baresoil roughness length (m) ! tuning to ERA5 fluxes
!REAL (KIND=DBL), PARAMETER :: z0_soil = 0.01_DBL   ! baresoil roughness length (m)
REAL (KIND=DBL), PARAMETER :: z0_ice = 0.001_DBL      ! ice roughness length (m)
!REAL (KIND=DBL), PARAMETER :: z0_soil = 0.0387_DBL   ! baresoil roughness length (m)
!REAL (KIND=DBL), PARAMETER :: z0_ice = 0.01_DBL      ! ice roughness length (m)
REAL (KIND=DBL), PARAMETER :: seaicedepth = 1.5_DBL  ! prescribed thickness of seaice (m)
REAL (KIND=DBL), PARAMETER :: tfriz = 273.15_DBL     ! freezing temperature
REAL (KIND=DBL), PARAMETER :: tfrizs = 271.35_DBL    ! seawater freezing temperature
REAL (KIND=DBL), PARAMETER :: mws_mx0 = 50._DBL      ! maximum puddle water storage (mm)
REAL (KIND=DBL), PARAMETER :: IR_emis_urban = 0.90_DBL    ! IR emissivity of urban area/buildings
REAL (KIND=DBL), PARAMETER :: IR_emis_leaf = 0.97_DBL    ! IR emissivity of a leaf
REAL (KIND=DBL), PARAMETER :: IR_emis_soil = 0.96_DBL    ! IR emissivity of bare soil
REAL (KIND=DBL), PARAMETER :: IR_emis_snow = 0.985_DBL    ! IR emissivity of snow
!REAL (KIND=DBL), PARAMETER :: IR_emis_snow = 0.85_DBL    ! IR emissivity of snow
REAL (KIND=DBL), PARAMETER :: IR_emis_ice = 0.97_DBL    !  IR emissivity of ice

REAL (KIND=DBL), PARAMETER :: rho_water = 998._DBL    ! density of snow (kg/m3)
REAL (KIND=DBL), PARAMETER :: rho_ice = 917._DBL    ! density of snow (kg/m3)
REAL (KIND=DBL), PARAMETER :: rho_snow = 100._DBL    ! density of snow (kg/m3)
REAL (KIND=DBL), PARAMETER :: cond_water = 0.57_DBL    ! heat conductivity of ice (W/mK)
REAL (KIND=DBL), PARAMETER :: cond_ice = 2.2_DBL       ! heat conductivity of ice (W/mK)
REAL (KIND=DBL), PARAMETER :: cond_snow = 0.01_DBL     ! heat conductivity of snow (W/mK)
REAL (KIND=DBL), PARAMETER :: capa_water = 4182._DBL    ! heat capacity of water (J/kgK)
REAL (KIND=DBL), PARAMETER :: capa_ice = 2030._DBL    ! heat capacity of ice (J/kgK)
REAL (KIND=DBL), PARAMETER :: cp_water = 4.1796e6_DBL ! [J/m3K]

REAL (KIND=DBL), PARAMETER :: LAI_min = 0.1_DBL     !  minimum LAI
REAL (KIND=DBL), PARAMETER :: t_canop_max = 343._DBL ! maximum plant temperature
REAL (KIND=DBL), PARAMETER :: leaf_thickness = 1.0 ! typical leaf thickness (mm)-used for heat capacity

! albedos of snow are for dry snow in Ebert abd Curry (1993)
REAL, parameter :: albvis_snow = 0.98 ! visible
REAL, parameter :: albnir_snow = 0.65 ! NIR

! albedo of glaciers:
REAL, parameter :: albvis_ice = 0.91 !0.96 ! visible
REAL, parameter :: albnir_ice = 0.65 !0.69 ! NIR

! albedo of bare soil:
REAL, parameter :: albvis_soil = 0.208 ! visible
REAL, parameter :: albnir_soil = 0.344 ! NIR
!REAL, parameter :: albvis_soil = 0.23 ! visible
!REAL, parameter :: albnir_soil = 0.38 ! NIR

! albedo of urban area:
REAL, parameter :: albvis_urban = 0.15 ! visible
REAL, parameter :: albnir_urban = 0.25 ! NIR


end module slm_params
