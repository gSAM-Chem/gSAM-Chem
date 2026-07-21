module emissions

    use grid, only : nx, ny, nz, nzm, z, dx, dy, dz, day, time, dt, pres, nstep, dimx1_s, dimx2_s, dimy1_s, dimy2_s               ! pres is in mbar!
    use cloudchem_Parameters, only : ind_NO, ind_ISOP, NVAR
    use chemistry_params, only : do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, do_CTG_lightning, do_IC_lightning, CTG_decaria_reflectivity, CTG_price_and_rind, IC_decaria, tropopause_index, increase_in_soil_moisture_linear, tau_decay_in_soil_moisture
    use vars, only : dtn
    
    implicit none

    ! Common constants
    real, parameter :: pi = 3.1415927
    real, parameter :: universal_gas_constant = 8.314                                  ! gas constant, in J mol-1 K-1
    real, parameter :: avogadro_number = 6.022e23       ! molecules/mol
    real, parameter :: isop_molar_mass = 68.12   ! in g/mol
    real, parameter :: no_molar_mass = 30.01     ! in g/mol
    logical :: isallocatedEmissions = .false.

    ! Conversion factors for downward shortwave -> PAR
    real, parameter :: J_to_mol = 4.6                                           ! Approximate unit conversion between W m-2 and umol photons m-2 s-1
    real, parameter :: frac_PAR = 0.5                                           ! Fraction of downward shortwave that is PAR

    ! Allocatable arrays for emissions
    real, allocatable, dimension(:,:) :: soil_NOx_activity_factor ! soil NOx array (x,y)
    real, allocatable, dimension(:,:) :: ppfd ! to calculate PAR-dependence of isoprene emissions
    real, allocatable, dimension(:,:) :: radiation_activity_factor ! to calculate PAR-dependence of isoprene emissions
    real, allocatable, dimension(:,:) :: temperature_activity_factor ! to calculate PAR-dependence of isoprene emissions
    real, allocatable, dimension(:,:) :: delta_ppv         ! isoprene to add to gridbox

    ! Allocatable arrays for lightning
    real, allocatable, dimension(:) :: vertical_function_profile            ! f(x) for CTG in DeCaria
    real, allocatable, dimension(:,:,:) :: change_in_mixing_ratio                   ! delta q_NO(z) in DeCaria
    real, allocatable, dimension(:) :: number_of_20dbz_per_altitude
    real, allocatable, dimension(:,:,:) :: refl_10cm                            ! reflectivity at 10 cm
    real, allocatable, dimension(:) :: refl_10cm_col                        ! reflectivity at 10 cm

    real, allocatable, dimension(:, :) :: cloud_top_heights                  ! for use in Price & Rind parametrization
    real, allocatable, dimension(:, :) :: cloud_top_temps                    ! for use in Price & Rind parametrization
    real, allocatable, dimension(:, :) :: price_and_rind_flash_rates         ! for use in Price & Rind parametrization
    real, allocatable, dimension(:) :: x_area_of_storm                                                 ! area the storm takes up

contains

    subroutine emissions_init()
        if ( .not. isallocatedEmissions ) then 
            if ( do_megan_isoprene ) then 
                allocate(radiation_activity_factor(nx, ny), temperature_activity_factor(nx, ny), ppfd(nx, ny), delta_ppv(nx, ny))
            endif

            if ( do_bdsnp_no ) then 
                allocate(soil_NOx_activity_factor(nx, ny))
            endif

            if ( do_CTG_lightning .or. do_IC_lightning ) then 
                allocate(vertical_function_profile(nzm), change_in_mixing_ratio(nx, ny, nzm))
            endif

            if ( do_IC_lightning ) then 
                allocate(x_area_of_storm(nzm))
            endif

            if ( CTG_price_and_rind ) then 
                allocate(cloud_top_heights(nx, ny), cloud_top_temps(nx, ny), price_and_rind_flash_rates(nx, ny))
            endif

            if ( CTG_decaria_reflectivity ) then 
                allocate(number_of_20dbz_per_altitude(nzm), refl_10cm_col(nzm), refl_10cm(nx, ny, nzm))
            endif



            isallocatedEmissions = .true.
        endif
    end subroutine 

    subroutine emissions_finalize()

        if ( do_megan_isoprene ) then 
            deallocate(radiation_activity_factor, temperature_activity_factor, ppfd, delta_ppv)
        endif

        if ( do_bdsnp_no ) then
            deallocate(soil_NOx_activity_factor)
        endif

        isallocatedEmissions = .false.
    end subroutine 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! SURFACE EMISSIONS !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! do_surface_emission_flux() provides a source of surface emissions to the bottommost model box
    ! you must send by reference the fluxbch (flux) (x,y) array to populate
    ! fluxbch is a variable in the chemistry module
    subroutine surface_emission_flux_driver(gchem_field, fluxbch, canopy_index, M_profile, isop_emission_flux, soil_NO_emission_flux)
        use vars, only : tabs
        use rad, only : swDownSurface, lwDownSurface, insolation_TOA                    ! eventually converted to PAR

        implicit none

        real, intent(inout) :: gchem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, NVAR)     ! Must specify size due to interface issue!
        real, intent(inout) :: fluxbch(:,:,:)                      ! parameter containing the surface fluxes
        integer, intent(in) :: canopy_index      
        real, intent(in)  :: M_profile(:)                      ! parameter containing the surface fluxes
        real, intent(out) :: isop_emission_flux(:), soil_NO_emission_flux(:)  ! parameter containing the surface fluxes

        ! Average fluxes
        real, parameter :: iso_avg_flux = 3.15e-10 * 1.75 ! From midpoint of Sarkar et al. (1 mg C m-2 hr-1)
        ! real :: conversion_between_isoprene_and_nox = 0.1765 * 0.8 * 0.15 !  0.0154454 !0.00077227 ! 0.0072    ! To scale average isoprene flux [ratio of kg m-2 s-1 to kg m-2 s-1]
        real, parameter :: no_avg_flux = 2.001e-12 * 1.35

        ! Other variables
        real, parameter :: LDF_i = 1                                                ! Set for isoprene
        real :: layer_height
        integer :: i, j
        real :: iso_avg_flux_ppv_s  ! From midpoint of Sarkar et al. (1 mg C m-2 hr-1)

        !! For non-MEGAN ISOP diurnal cycle
        real :: t_solar_peak = 0.167                                                ! days
        real :: frac_of_peak

        ! Default values for emissions
        fluxbch(:, :, ind_ISOP) = 0.
        fluxbch(:, :, ind_NO) = 0.
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Isoprene Flux Module !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        if ( do_megan_isoprene ) then 
            radiation_activity_factor = 0.
            temperature_activity_factor = 0.
            ppfd = 0.

            !! MEGAN Isoprene Radiation Module !!
            ! Convert radiation from W m-2 to umol m-2 s-1 (https://search.r-project.org/CRAN/refmans/bigleaf/html/Rg.to.PPFD.html)
            ! if ( ALLOCATED(swDownSurface) ) then
            ppfd(:,:) = swDownSurface(:,:) * J_to_mol * frac_PAR
            radiation_activity_factor(:,:) = calculate_megan_BVOC_radiation(ppfd(:,:), LDF_i)
            ! endif

            !! MEGAN Isoprene Temperature Module !!
            ! Calculate temperature activity factor for isoprene
            temperature_activity_factor = calculate_megan_BVOC_temperature(tabs(:,:,1), LDF_i)
            
            if ( canopy_index .eq. 1 ) then
                fluxbch(:,:,ind_ISOP) =  ( iso_avg_flux * 1000 / isop_molar_mass / M_profile(1) * avogadro_number / 100**3 ) * temperature_activity_factor * radiation_activity_factor     ! Calculate isoprene fluxes using MEGAN
                isop_emission_flux(1) = SUM(fluxbch(:,:,ind_ISOP))

            else 
                if (canopy_index < nzm) then
                    layer_height = (z(canopy_index + 1) - z(canopy_index)) * 100    ! in cm    
                else
                    layer_height = (z(canopy_index) - z(canopy_index - 1)) * 100
                endif

                iso_avg_flux_ppv_s = ( iso_avg_flux * 1000 / isop_molar_mass / M_profile(canopy_index) * avogadro_number / 100**2 / layer_height * dtn) 

                do i = 1, nx
                    do j = 1, ny
                        delta_ppv(i,j) = iso_avg_flux_ppv_s * temperature_activity_factor(i, j) * radiation_activity_factor(i, j)     ! Calculate isoprene fluxes using MEGAN
                        gchem_field(i, j, canopy_index, ind_ISOP) = gchem_field(i, j, canopy_index, ind_ISOP) + delta_ppv(i, j)
                    end do
                enddo
                      
                isop_emission_flux(canopy_index) = SUM(delta_ppv) / dtn
            endif

        elseif ( do_surface_Isoprene_diurnal ) then
            frac_of_peak = MAX(0., cos((day - t_solar_peak)*2*pi))
            fluxbch(:,:,ind_ISOP) =  frac_of_peak * iso_avg_flux * (pi/2.)
        else
            fluxbch(:,:,ind_ISOP) = iso_avg_flux
        endif
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Soil NO Flux Module !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        if ( do_bdsnp_no ) then
            soil_NOx_activity_factor = 0.                       ! Initialize to zero

            call calculate_bdsnp_NO()
            fluxbch(:, :, ind_NO) = ( no_avg_flux * 1000 / no_molar_mass / M_profile(1) * avogadro_number / 100**3 ) * soil_NOx_activity_factor(:, :)
            soil_NO_emission_flux(1) = SUM(fluxbch(:, :, ind_NO))
        endif 
    end subroutine surface_emission_flux_driver

    function calculate_megan_BVOC_temperature(surface_temperature, LDF_i) result(gamma_T)
        ! Returns the value of the MEGAN scaling factor (gamma_T) for temperature
        ! See Guenther et al. (2012) for parametrization
        implicit none

        ! Parameters
        real, intent(in) :: surface_temperature(nx, ny)                  ! temperature at bottommost layer
        real, intent(in) :: LDF_i                                                    ! light-dependent fraction, = 1 for isoprene

        real :: gamma_T_ldf(nx, ny)                                     ! temperature activity factor, light-dependent
        real :: gamma_T_lif(nx, ny)                                     ! temperature activity factor, light-independent
        real :: gamma_T(nx, ny)                                          ! final temperature activity factor, to return

        ! MEGAN - light-independent fraction
        real, parameter :: beta_i = 0.13                                            ! empirically determined coefficient for each VOC, tuned to ISOP
        real, parameter :: T_s = 297                                                ! standard temperature conditions for leaf temperature [K]

        ! MEGAN - light-dependent fraction
        real, parameter :: C_eo = 2                                                 ! Changes with species (currently set at isoprene)
        real, parameter :: C_t1 = 95                                                ! Changes with species (currently set at isoprene)
        real, parameter :: C_t2 = 230                                               ! Empirical coefficient
        real, parameter :: T_24 = 298                                               ! Average leaf temperature of past 24 hours
        real, parameter :: T_240 = 298                                              ! Average leaf temperature of past 240 hours

        ! intermediate calculations
        real :: T_opt(nx,ny)
        real :: E_opt(nx, ny)
        real :: x(nx, ny)

        gamma_T_lif = exp(beta_i * (surface_temperature(:,:) - T_s))                                ! Light-independent -- similar to monoterpene flux

        T_opt = 313 + (0.6 * (T_240 - T_s))                                                         ! Optimal temperature
        E_opt = C_eo * exp(0.05 * (T_24 - T_s)) * exp(0.05 * (T_240 - T_s))
        x = ((1.0 / T_opt) - (1.0 / surface_temperature)) / 0.00831

        gamma_T_ldf = E_opt * (C_t2 * exp(C_t1 * x) / (C_t2 - C_t1 * (1.0 - exp(C_t2 * x))))          ! light dependent emission activity factor
        gamma_T = (1.0 - LDF_i) * gamma_T_lif + (LDF_i * gamma_T_ldf)                                 ! MEGAN emission activity factor, accounts for light dependent and light independent factors    
    end function calculate_megan_BVOC_temperature


    function calculate_megan_BVOC_radiation(ppfd, LDF_i) result(gamma_p)
        ! Returns the value of the MEGAN scaling factor (gamma_P) for temperature
        ! See Guenther et al. (2012) for parametrization
        implicit none

        ! parameters
        real, intent(in) :: ppfd(nx, ny)                ! photosynthetic photon flux density, in umol m-2 s-1
        real, intent(in) :: LDF_i                    ! light-dependent fraction, = 1 for isoprene

        ! constants
        real, parameter :: p_s = 200               ! Standard conditions for PPFD, equal to 200 umol m-2 s-1 for sunlit leaves, 50 for shaded leaves
        real, parameter :: p_24 = 310              ! average PPFD of past 24 hours
        real, parameter :: p_240 = 310             ! average PPFD of past 240 hours

        ! intermediate calculations
        real :: c_p 
        real :: alpha
        real :: gamma_p_ldf(nx, ny)
        real :: gamma_p(nx, ny)

        alpha = 0.004 - ( 0.0005 * log(p_240) )
        c_p = 0.0468 * exp(0.0005 * ( p_24 - p_s )) * p_240**(0.6)

        gamma_p_ldf = c_p * ((alpha * ppfd(:,:)) / (1.0 + (alpha**2  * ppfd(:,:)**2))**0.5)       ! light dependent emission activityfactor
        gamma_p = (1.0 - LDF_i) + ( LDF_i * gamma_p_ldf )                                         ! MEGAN emission activity factor, accounts for light dependent and light independent factors    
    end function calculate_megan_BVOC_radiation


    subroutine calculate_bdsnp_NO()
        ! Calculates soil NOx
        ! See Hudman et al. (2012) and Wang et al. (2021) for parametrization
        use vars, only : prec_xy, precinstsoil, precsfc, precinst, tabs, interactive_soil_wetness ! precsfc (x,y), tabs (x,y,z) (JY), and interactive_soil_wetness (only used in chem_flux) (JY)

        implicit none

        ! BDSNP Parameters
        real, parameter :: a_bdsnp = 1.65
        real, parameter :: b_bdsnp = 3.3
        real :: temperature_in_celsius
        real :: rain_threshold = 1e-5
        
        ! Conversion factors in soil moisture calculations
        real :: conversion_between_hour_and_second = 3600       ! MERRA2 is on an hourly time grid versus seconds for SAM

        integer :: i, j                                                  ! Counter variables

        ! Loop through the horizontal axes


        do i = 1, nx
            do j = 1, ny

                ! If there is precipitation, increase soil moisture
                ! Precsfc is outputted as mm/day, but in the model it is kg m-2 s-1!
                if ( precinst(i, j) < 0.01 ) then        ! Check if precsfc is unreasonably high; if so, don't update!
                    ! if ( interactive_soil_wetness(i, j) <= SM_threshold ) then 
                    interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) + ( increase_in_soil_moisture_linear * precinst(i,j) * dt / conversion_between_hour_and_second ) - ( ( 1 / tau_decay_in_soil_moisture ) * interactive_soil_wetness(i, j) * dt / conversion_between_hour_and_second )

                    ! elseif ( interactive_soil_wetness(i, j) > SM_threshold ) then
                    !     interactive_soil_wetness(i, j) = interactive_soil_wetness(i, j) + ( increase_in_soil_moisture_saturated * precsfc(i,j) * dt / conversion_between_hour_and_second ) - ( ( 1 / tau_decay_in_soil_moisture ) * interactive_soil_wetness(i, j) * dt / conversion_between_hour_and_second )
                    ! endif
                endif

                if ( interactive_soil_wetness(i, j) .gt. 1 ) then
                    interactive_soil_wetness(i, j) = 1
                elseif ( interactive_soil_wetness(i, j) .lt. 0 ) then
                    interactive_soil_wetness(i, j) = 0
                endif

                !! Now, calculate soil NOx !!
                temperature_in_celsius = tabs(i, j, 1) - 273.15

                if ( temperature_in_celsius < 20 ) then
                    soil_NOx_activity_factor(i, j) = exp( 0.103 * temperature_in_celsius ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                
                elseif ( ( temperature_in_celsius >= 20 ) .and. ( temperature_in_celsius <= 40 ) ) then
                    soil_NOx_activity_factor(i, j) = ( -0.009 * temperature_in_celsius**3 + 0.837 * temperature_in_celsius**2 - 22.52 * temperature_in_celsius + 196.149 ) * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                
                else
                    soil_NOx_activity_factor(i, j) = 58.269 * a_bdsnp * interactive_soil_wetness(i,j) * exp(-1 * b_bdsnp * interactive_soil_wetness(i,j)**2)
                endif

            end do
        end do
    end subroutine calculate_bdsnp_NO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! LIGHTNING EMISSIONS !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lightning_decaria_ctg(gchem_field)
        ! Parametrizes (as close as possible!) to DeCaria et al. (2000), 
        ! A cloud-scale model study of lightning-generated NOx...
        
        use grid, only : z
        use vars, only : tabs, qpl, qpi                       ! tabs (x,y,z) contains absolute temperature in K
        use params, only : land, ocean
        ! use microphysics, only: dBZ_cloudradar    ! not every microphysics module has a dBZ_cloudradar variable! 

        implicit none

        real, intent(inout) :: gchem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, NVAR)     ! Must specify size due to interface issue!

        ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
        real, parameter :: time_between_flashes = 180                                      ! 3 minutes
        integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)
        real :: moist_adiabatic_lapse_rate = 5e-3                               ! an estimate of the moist adiabatic lapse rate, in C / m

        ! real, parameter:: radar_threshold = 1e-4                                           ! to match 20 dBZ, approximate qp mixing ratio (kg/kg) based on regression
        real, parameter:: radar_threshold = 20

        ! Cloud-to-ground lightning
        real, parameter :: cloud_to_ground_isotherm = -15                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
        integer :: vertical_index_for_isotherm                                  ! index of the CTG mean isotherm
        real :: mu
        real :: std_dev                                                         ! in meters
        real :: temperature_of_closest_altitude

        ! Production rates
        real :: no_production_per_flash = 460                                   ! CTG NO production per flash, in mol/flash (460 is CTG from Price & Rind)
        real :: integral_of_vertical_function                                   ! Used to get LC from N_tot

        real, parameter :: cloud_top_height_threshold = 5000.                               ! in meters
        integer :: i, j, k                                                          ! Counter variables
        real, parameter :: lightning_scaling = 250.                            ! Scaling on lightning emissions (part of quotient)
        real, parameter :: storm_extent_y_m = 20000.                            ! Spatial extent of storm in y-direction (based on Decaria)

        ! Only for Price & Rind, eventually remove
        real :: x_area_of_storm

        lightning_time_step = time_between_flashes / dt

        ! If the desired time between flashes has passed, then do lightning
        if ( mod(nstep, lightning_time_step) == 0 .and. nstep .gt. 0 ) then
            change_in_mixing_ratio = 0.

            do i = 1, nx
                do j = 1, ny 
                    ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                    vertical_index_for_isotherm = minloc(abs((tabs(i,j,:) - 273.15) - cloud_to_ground_isotherm), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                    mu = z(vertical_index_for_isotherm)
                    temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm) - 273.15 

                    ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                    if ( abs( temperature_of_closest_altitude - cloud_to_ground_isotherm ) > 1.5 ) then
                        mu = mu - ( cloud_to_ground_isotherm - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                    endif
                    
                    std_dev = mu / 3.
                    vertical_function_profile = 1 / (sqrt(2 * pi) * std_dev) * exp(-1 * (z - mu)**2 / (2 * std_dev**2))
                    
                    ! Integration for denominator
                    integral_of_vertical_function = 0.
                    do k = 1, nzm
                        integral_of_vertical_function = integral_of_vertical_function + ( vertical_function_profile(k) * pres(k) * 100 ) * dz           ! 100 is to convert from mbar to Pa
                    enddo
                    
                    if ( integral_of_vertical_function > 0 ) then
                            change_in_mixing_ratio(i,j,:) = ( no_production_per_flash / integral_of_vertical_function ) * universal_gas_constant * vertical_function_profile * tabs(i, j, :)
                    endif

                enddo
            enddo

            if ( CTG_price_and_rind ) then
                cloud_top_heights = 0.
                cloud_top_temps = 0.

                call calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)

                x_area_of_storm = 0.
                price_and_rind_flash_rates = 0.

                do i = 1, nx
                    do j = 1, ny 
                        if ( cloud_top_heights(i, j) <= cloud_top_height_threshold ) then
                            cloud_top_heights(i, j) = 0.
                        else
                            x_area_of_storm = x_area_of_storm + 1

                            if ( land ) then
                                    price_and_rind_flash_rates(i,j) = 3.44e-5 * ( cloud_top_heights(i,j) / 1000 )**4.9
                            elseif ( ocean ) then
                                    price_and_rind_flash_rates(i,j) = 6.4e-4 * ( cloud_top_heights(i,j) / 1000 )**1.73
                            endif
                        endif
                    end do
                end do

                do i = 1, nx
                    do j = 1, ny 
                        if ( cloud_top_heights(i, j) < 1 ) then 
                            change_in_mixing_ratio(i,j,:)= 0.
                        else
                            do k = 1, nzm
                                change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( x_area_of_storm * dx * dx ) ! eventually replace dx with dy
                                change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) * price_and_rind_flash_rates(i,j) / sum( price_and_rind_flash_rates(:,:) )
                            enddo
                        endif
                    enddo
                enddo
                gchem_field(:,:,:, ind_NO) = gchem_field(:,:,:, ind_NO) + change_in_mixing_ratio(:,:,:)
            endif

            if ( CTG_decaria_reflectivity ) then
                number_of_20dbz_per_altitude = 0.
                refl_10cm = 0.
                call calculate_radar_reflectivity()

                do k = 1,nzm 
                    do i = 1,nx
                        do j = 1,ny
                            if ( refl_10cm(i, j, k) >= radar_threshold ) then
                                number_of_20dbz_per_altitude(k) = number_of_20dbz_per_altitude(k) + 1
                            else
                                change_in_mixing_ratio(i,j,k) = 0.
                            endif
                        enddo
                    enddo
                enddo

                do i = 1, nx
                    do j = 1, ny 
                        do k = 1, nzm
                            if ( number_of_20dbz_per_altitude(k) > 0. ) then 
                                change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( number_of_20dbz_per_altitude(k) * dx * storm_extent_y_m * lightning_scaling ) ! eventually replace dx with dy; ADDED *2 11/19
                            endif
                        enddo
                    enddo
                enddo

                ! , "*************** Maximum CTG lightning = ", MAXVAL(change_in_mixing_ratio)
                ! print*, "*************** Minimum CTG lightning = ", MINVAL(change_in_mixing_ratio)

                ! print*, "*************** Maximum NO  = ", MAXVAL(gchem_field(:, :, :, ind_NO))
                ! print*, "*************** Minimum NO  = ", MINVAL(gchem_field(:, :, :, ind_NO))

                do i = 1, nx
                    do j = 1, ny 
                        do k = 1, nzm
                            gchem_field(i, j, k, ind_NO) = gchem_field(i, j, k, ind_NO) + change_in_mixing_ratio(i, j, k)
                        enddo
                    enddo
                enddo                
            endif
        endif   
    end subroutine lightning_decaria_ctg

    subroutine lightning_decaria_ic(gchem_field)
        ! Parametrizes (as close as possible!) to DeCaria et al. (2000), A cloud-scale model study of lightning-generated NOx...
        use grid, only : z
        use vars, only : tabs, qcl, qpl, qci, qpi                                                   ! tabs (x,y,z) contains absolute temperature in K
        use params, only : land, ocean
        ! use microphysics, only: dBZ_cloudradar

        implicit none

        ! real :: radar_threshold = 1e-4                                           ! to match 20 dBZ, approximate qp mixing ratio (kg/kg) based on regression
        real, parameter :: radar_threshold = 20

        ! General parameters for both cloud-to-ground (CTG) and intracloud (IC) lightning
        real :: time_between_flashes = 180                                      ! 3 minutes
        integer :: lightning_time_step                                          ! How many time steps to skip before calculating lightning (function of dt)

        real, intent(inout) :: gchem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, NVAR)     ! Must specify size due to interface issue!
        real :: moist_adiabatic_lapse_rate = 5e-3                               ! an estimate of the moist adiabatic lapse rate, in C / m

        ! Cloud-to-ground lightning
        real, parameter :: ic_isotherm_bottom = -15                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
        integer :: vertical_index_for_isotherm_bottom                                  ! index of the CTG mean isotherm

        real, parameter :: ic_isotherm_top = -45                                  ! isoterm of the Gaussian mean for CTG lightning, in Celsius
        integer :: vertical_index_for_isotherm_top                                  ! index of the CTG mean isotherm

        real :: mu_bottom
        real :: mu_top

        real :: std_dev_bottom                                                         ! in meters
        real :: std_dev_top                                                         ! in meters

        real :: temperature_of_closest_altitude

        ! Production rates
        real :: no_production_per_flash = 460                                   ! CTG NO production per flash, in mol/flash (460 is CTG from Price & Rind)
        real :: integral_of_vertical_function                                   ! Used to get LC from N_tot
        
        real :: ic_mixing_ratio_threshold = 0.01                                ! in g/kg
        integer :: i, j, k                                                          ! Counter variables
        real, parameter :: lightning_scaling = 250.
        real, parameter :: storm_extent_y_m = 20000.                            ! Spatial extent of storm in y-direction (based on Decaria)

        lightning_time_step = time_between_flashes / dt

        ! If the desired time between flashes has passed, then do lightning
        if ( mod(nstep, lightning_time_step) == 0 ) then
            change_in_mixing_ratio = 0.

            do i = 1, nx
                do j = 1, ny 

                    ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                    vertical_index_for_isotherm_bottom = minloc(abs((tabs(i,j,:) - 273.15) - ic_isotherm_bottom), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                    mu_bottom = z(vertical_index_for_isotherm_bottom)
                    temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm_bottom) - 273.15 

                    ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                    if ( abs( temperature_of_closest_altitude - ic_isotherm_bottom ) > 1.5 ) then
                        mu_bottom = mu_bottom - ( ic_isotherm_bottom - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                    endif

                    std_dev_bottom = mu_bottom / 3.

                    ! Gives the closest to the isotherm, but there is no guarantee that this is close enough
                    vertical_index_for_isotherm_top = minloc(abs((tabs(i,j,:) - 273.15) - ic_isotherm_top), dim = 1) ! dim=1 is required to return an integer, not an 1-element array

                    mu_top = z(vertical_index_for_isotherm_top)
                    temperature_of_closest_altitude = tabs(i, j, vertical_index_for_isotherm_top) - 273.15 

                    ! If the closest temperature is farther than 1.5 deg C away from the isotherm, extrapolate using the MALR
                    if ( abs( temperature_of_closest_altitude - ic_isotherm_top ) > 1.5 ) then
                        mu_top = mu_top - ( ic_isotherm_top - temperature_of_closest_altitude ) / moist_adiabatic_lapse_rate
                    endif
                    std_dev_top = std_dev_bottom / 3.

                    vertical_function_profile = ( (0.8 / (sqrt(2 * pi) * std_dev_bottom) * exp(-1 * (z - mu_bottom)**2 / (2 * std_dev_bottom**2))) + (1 / (sqrt(2 * pi) * std_dev_top) * exp(-1 * (z - mu_top)**2 / (2 * std_dev_top**2))) )
                    
                    ! Integration for denominator
                    integral_of_vertical_function = 0.

                    do k = 1, nzm
                        integral_of_vertical_function = integral_of_vertical_function + ( vertical_function_profile(k) * pres(k) * 100 ) * dz           ! 100 is to convert from mbar to Pa
                    enddo

                    ! if ( integral_of_vertical_function > 0 ) then
                    !         change_in_mixing_ratio(i,j,:) = ( no_production_per_flash / integral_of_vertical_function ) * universal_gas_constant * vertical_function_profile * tabs(i, j, :)
                    !         print*, "change in mixing ratio @ 256 = ", change_in_mixing_ratio(i, j, 256)
                    ! endif
                enddo
            enddo

            if ( IC_decaria ) then
                x_area_of_storm = 0.
                refl_10cm = 0.
                call calculate_radar_reflectivity()

                do k = 1,nzm 
                    do i = 1,nx
                        do j = 1,ny
                            if ( ( (qcl(i, j, k) + qpl(i, j, k) + qci(i, j, k) + qpi(i, j, k)) * 1000. > ic_mixing_ratio_threshold ) .and. ( maxval( refl_10cm ) >= radar_threshold ) ) then
                                x_area_of_storm(k) = x_area_of_storm(k) + 1
                            else
                                change_in_mixing_ratio(i,j,k) = 0.
                            endif
                        enddo
                    enddo
                enddo

                do i = 1, nx
                    do j = 1, ny 
                        do k = 1, nzm
                            if ( x_area_of_storm(k) > 0. ) then 
                            change_in_mixing_ratio(i,j,k) = change_in_mixing_ratio(i,j,k) / ( x_area_of_storm(k) * dx * storm_extent_y_m * lightning_scaling ) ! eventually replace dx with dy, ADDED *10 11/19
                            endif
                        enddo
                    enddo
                enddo

                do i = 1, nx
                    do j = 1, ny 
                        do k = 1, nzm
                            gchem_field(i, j, k, ind_NO) = gchem_field(i, j, k, ind_NO) + change_in_mixing_ratio(i, j, k)
                        enddo
                    enddo
                enddo

            endif
        endif
    end subroutine lightning_decaria_ic

    subroutine calculate_cloud_top_height(cloud_top_heights, cloud_top_temps)
        ! Copy of the function that calculates cloud top height in diagnose
        ! This function is required because diagnose is run after emissions/chemistry
        
        use grid, only : adz, dz, z
        use vars, only : qcl, qci, rho, tabs
        
        implicit none

        integer :: i, j, k                                                          ! Counter variables
        real :: tmp_lwp                                                             ! Temporary variable
        real, allocatable, dimension(:,:) :: cloud_top_heights
        real, allocatable, dimension(:,:) :: cloud_top_temps

        cloud_top_heights = 0.
        cloud_top_temps = 0.

        do j = 1,ny
            do i = 1,nx
                tmp_lwp = 0.
                do k = nzm,1,-1

                    tmp_lwp = tmp_lwp + (qcl(i,j,k) + qci(i,j,k)) * rho(k) * dz * adz(k)
                
                    if (tmp_lwp.gt.0.01) then
                        cloud_top_heights(i,j) = z(k)
                        cloud_top_temps(i,j) = tabs(i,j,k)
                        exit
                    end if
                
                end do
            end do
        end do
    end subroutine calculate_cloud_top_height

    subroutine calculate_radar_reflectivity()
        use module_mp_GRAUPEL, only : calc_refl10cm
        use microphysics, only :  micro_field, iqv, iqr, iqs, iqg, inr, ins, ing
        use vars, only : tabs
        use grid, only : pres

        implicit none 

        integer i, j, k

        do i = 1, nx
            do j = 1, ny
                call calc_refl10cm( &
                    micro_field(i,j,:,iqv), &
                    micro_field(i,j,:,iqr), &
                    micro_field(i,j,:,iqs), &
                    micro_field(i,j,:,iqg), &
                    tabs(i,j,:), &
                    pres(:), &
                    refl_10cm_col, &
                    1, nzm, i, j, &
                    micro_field(i,j,:,inr), &
                    micro_field(i,j,:,ins), &
                    micro_field(i,j,:,ing))

                do k = 1, nzm
                    refl_10cm(i,j,k) = refl_10cm_col(k)
                enddo
            end do
        end do

    end subroutine calculate_radar_reflectivity

end module emissions
