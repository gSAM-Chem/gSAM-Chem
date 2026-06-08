module rad
  use shr_kind_mod, only: r4 => shr_kind_r4
  use grid

implicit none

!--------------------------------------------------------------------
!
! Variables accumulated between two calls of radiation routines


        real(r4), allocatable ::  tabs_rad(:,:,:)     ! accumulated temperature
        real(r4), allocatable ::  qc_rad(:,:,:)       ! accumulated cloud water (g/g)
        real(r4), allocatable ::  qi_rad(:,:,:)       ! accumulated cloud ice (g/g)
        real(r4), allocatable ::  qv_rad(:,:,:)       ! accumulated water vapor (g/g)
        real(r4), allocatable ::  cld_rad(:,:,:)      ! accumulated cloud fraction 
        real(r4), allocatable ::  rel_rad(:,:,:)      ! accumulated effective radius for liquid water (mkm)
        real(r4), allocatable ::  rei_rad(:,:,:)      ! accumulated effective radius for ice water (mkm)
        ! Fields for radiatively-active snow
        real(r4), allocatable ::  qs_rad  (:,:,:)       ! accumulated snow mass mixing ratio (g/g)
        real(r4), allocatable ::  res_rad  (:,:,:)      ! accumulated effective radius for snow (mkm)
        real  ozone  (nx, ny, nzm)        ! ozone 3D field (g/g) if read from 3D initfile
	real:: qrad    (nx, ny, nzm)=0 ! radiative heating(K/s)
	real:: lwnsxy  (nx, ny)=0
	real:: swnsxy  (nx, ny)=0
	real:: lwntxy  (nx, ny)=0
	real:: swntxy  (nx, ny)=0
	real:: lwntmxy (nx, ny)=0
	real:: swntmxy (nx, ny)=0
	real:: lwnscxy (nx, ny)=0
	real:: swnscxy (nx, ny)=0
	real:: lwntcxy (nx, ny)=0
	real:: swntcxy (nx, ny)=0
	real:: lwdsxy  (nx, ny)=0
	real:: swdsxy  (nx, ny)=0
	real:: lwdscxy (nx, ny)=0
	real:: swdscxy (nx, ny)=0
	real:: swdsvisxy  (nx, ny)=0
	real:: swdsnirxy  (nx, ny)=0
	real:: swdsvisdxy  (nx, ny)=0
	real:: swdsnirdxy  (nx, ny)=0
	real:: swusvisdxy  (nx, ny)=0
	real:: swusnirdxy  (nx, ny)=0
	real:: solinxy (nx, ny)=0
        real:: coszrsxy(nx,ny)=0
        real:: albvisxy(nx,ny)=0
        real:: albnirxy(nx,ny)=0


    ! Instrument simulator fields
    real tau_067   (nx, ny, nzm)
    real emis_105  (nx, ny, nzm)
    real rad_reffc (nx, ny, nzm)
    real rad_reffi (nx, ny, nzm)

    ! separate optical depths for liquid and ice for MODIS simulator
    real tau_067_cldliq   (nx, ny, nzm) 
    real tau_067_cldice   (nx, ny, nzm)
    real tau_067_snow   (nx, ny, nzm)


	logical:: initrad =.true. ! flag to initialize profiles of traces
	integer nradsteps       ! number of steps done before calling radiation routines
        character(3), parameter :: RAD_NAME="CAM"
        logical, parameter :: do_output_clearsky_heating_profiles = .true.

!       Gas traces (mass mixing ratios):  
        
        real(r4) o3(nzm)            ! Ozone
        real(r4) n2o(nzm)           ! N2O
        real(r4) ch4(nzm)           ! CH4
        real(r4) cfc11(nzm)         ! CFC11
        real(r4) cfc12(nzm)         ! CFC12

	real(r4) p_factor(nx, ny) ! perpetual-sun factor

end module rad
