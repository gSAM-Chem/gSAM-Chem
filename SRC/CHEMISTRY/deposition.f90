module deposition
    ! This module performs dry deposition of species onto the surface.
    ! Depositing species and their dry deposition velocities should be inputted in the prm file.

    use cloudchem_Parameters, only: NVAR
    use chemistry_params, only: dry_deposition_species, dry_deposition_velocities
    use grid, only : z, nx, ny, nzm
    use cloudchem_Monitor, only: SPC_NAMES

    CONTAINS
    
    ! Note (JY): Could also use fluxbch instead of subtracting directly from gchem_field
    ! I kept it as the latter, since in the future we might want to have deposition in gridboxes not at the surface (e.g. leaves)
    subroutine dry_deposition_driver(gchem_field, g_depos_horiz_mean_tend_ISOPOOH, g_depos_horiz_mean_tend_IEPOX)
        real, allocatable, dimension(:,:,:,:) :: gchem_field
        real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_ISOPOOH      ! Old variables that were used to store ISOPOOH deposition rates
        real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_IEPOX        ! Old variables that were used to store IEPOX deposition rates

        integer :: i, j, v, v_selected
        real :: bottommost_layer_height

        bottommost_layer_height = ( z(2) - z(1) ) * 100                         ! Calculates height of lowest box; converts from m to cm b/c velocities in cm/s
        
        do v_selected = 1, NVAR                                                 ! v_selected is the index of the deposition array
            do v = 1, NVAR                                                      ! v is the index of the species list (e.g. gchem_fields)
                if ( dry_deposition_species(v_selected) == trim(SPC_NAMES(v)) ) then
                    ! Performs homogeneous deposition across the surface; no consideration for stomata/land proprties
                    do i = 1, nx
                        do j = 1, ny
                            gchem_field(i, j, 1, v) = ( 1. - ( dry_deposition_velocities(v_selected) * dtn ) / bottommost_layer_height ) * gchem_field(i, j, 1, v)
                            ! JY: units are (1 - ( cm/s * s ) / (cm))
                            ! This assumes that the lowest box is completely well-mixed
                            gchem_field(i, j, 1, v) = max(gchem_field(i, j, 1, v), 0.0)         ! Prevents negative mixing ratios
                        enddo
                    enddo

                    exit                                                         ! Species has been found, move onto the next species
                end if
            end do  
        end do
    end subroutine dry_deposition_driver

end module deposition