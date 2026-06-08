module rad

use grid

implicit none

!--------------------------------------------------------------------
!
! Variables accumulated between two calls of radiation routines

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

        logical:: initrad =.true. ! flag to initialize profiles of traces



 ! variables that need to exist in order for simulators to compile.
 real, allocatable, dimension(:,:,:) :: tau_067, emis_105, &
                    tau_067_cldliq, tau_067_cldice, tau_067_snow, &
                    rad_reffc, rad_reffi 

 ! variables needed so that clearsky heating rates can be output from other radiation schemes.
 logical, parameter :: do_output_clearsky_heating_profiles = .false.
 real, dimension(nz) :: radqrclw, radqrcsw

 character(3), parameter :: RAD_NAME="DUM"

end module rad
