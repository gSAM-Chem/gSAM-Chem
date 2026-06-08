!-----------------------------------------------------------------------------
!
! This module handles the multigrid kernels, including the procedures for
! level computation, coarsening, prolongation, and V-cycles.
!
!-----------------------------------------------------------------------------
! Author: Oliver Yang and Xiangmin Jiao
! Adapted from code written by Cao Lu.
!
! Initial version: November 3, 2017
! Last updated: September 27, 2024 (Refactored to encapsulate multigrid and
! options into a single type for thread safety, updated type names to
! follow capitalization conventions, and added structure for managing
! multigrid solvers.)
!-----------------------------------------------------------------------------

MODULE gmg_kernel
   USE gmg_data, ONLY: MultigridSolver
   USE gmg_utils, ONLY: debugging
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: compute_levels, coarsen, apply_cycle

CONTAINS

   !-----------------------------------------------------------------------
   ! FUNCTION: compute_levels
   ! PURPOSE:  Computes the number of grids required for the multigrid method
   !           based on the number of cells in the fine-grid y and z directions,
   !           the desired coarse-grid size, and whether semi-coarsening is used.
   !
   ! INPUTS:
   !   y_fine_cells  - INTEGER, number of fine cells in y-direction
   !   z_fine_cells  - INTEGER, number of fine cells in z-direction
   !   coarse_size   - INTEGER, desired grid size of the coarsest grid
   !   semi_coarsen  - LOGICAL, flag for semi-coarsening
   !
   ! OUTPUT:
   !   levels        - INTEGER, number of multigrid levels
   !-----------------------------------------------------------------------
   FUNCTION compute_levels(y_fine_cells, z_fine_cells, coarse_size, semi_coarsen) RESULT(levels)
      IMPLICIT NONE

      ! Input parameters
      INTEGER, INTENT(IN) :: y_fine_cells, z_fine_cells, coarse_size
      LOGICAL, INTENT(IN) :: semi_coarsen

      ! Output parameter
      INTEGER :: levels

      ! Local variables
      INTEGER :: y_levels, z_levels
      INTEGER :: reduced_cells

      y_levels = 1
      reduced_cells = y_fine_cells

      ! Calculate the number of levels in the y-direction.
      DO WHILE (reduced_cells > coarse_size)
         reduced_cells = (reduced_cells + 1)/2
         y_levels = y_levels + 1
      END DO

      ! Calculate the number of levels in the z-direction.
      z_levels = 1
      reduced_cells = z_fine_cells
      DO WHILE (reduced_cells > coarse_size)
         reduced_cells = (reduced_cells + 1)/2
         z_levels = z_levels + 1
      END DO

      ! Determine the total levels based on semi-coarsening.
      IF (semi_coarsen) THEN
         levels = MAX(y_levels, z_levels)
      ELSE
         levels = MIN(y_levels, z_levels)
      END IF

   END FUNCTION compute_levels

   !-----------------------------------------------------------------------
   ! SUBROUTINE: coarsen
   ! PURPOSE:    Computes the data associated with one coarser grid. Handles
   !             coarsening of cells based on finer grid information. For odd
   !             number of cells, alternates the padding side for robustness.
   !
   ! INPUTS/OUTPUTS:
   !   solver                    - MultigridSolver type, input/output
   !   finer_lvl                 - INTEGER, finer multigrid level
   !   y_finer_cell_cent_coords  - REAL(8) array, finer y-center coordinates
   !   z_finer_cell_cent_coords  - REAL(8) array, finer z-center coordinates
   !   y_padside                 - INTEGER, padding side in y-direction
   !   z_padside                 - INTEGER, padding side in z-direction
   !   y_cell_edge_coords        - REAL(8) array, y-edge coordinates, input/output
   !   z_cell_edge_coords        - REAL(8) array, z-edge coordinates, input/output
   !   y_cell_widths             - REAL(8) array, cell widths in y, input/output
   !   z_cell_widths             - REAL(8) array, cell widths in z, input/output
   !   p_cell_edge               - REAL(8) array, pressure cell edges in z-direction, input/output
   !
   ! OUTPUTS:
   !   y_coarser_cell_cent_coords - REAL(8) array, coarser y-center coordinates
   !   z_coarser_cell_cent_coords - REAL(8) array, coarser z-center coordinates
   !   y_cell_widths_work         - REAL(8) array, working array for cell widths in y
   !   z_cell_widths_work         - REAL(8) array, working array for cell widths in z
   !-----------------------------------------------------------------------
   SUBROUTINE coarsen(solver, finer_lvl, y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
                      y_padside, z_padside, y_cell_edge_coords, z_cell_edge_coords, &
                      y_cell_widths, z_cell_widths, p_cell_edge, &
                      y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, &
                      y_cell_widths_work, z_cell_widths_work)
      IMPLICIT NONE

      ! Input/output parameters
      TYPE(MultigridSolver), INTENT(INOUT) :: solver
      INTEGER, INTENT(IN) :: finer_lvl
      REAL(8), INTENT(IN) :: y_finer_cell_cent_coords(:)
      REAL(8), INTENT(IN) :: z_finer_cell_cent_coords(:)
      INTEGER, INTENT(INOUT) :: y_padside, z_padside
      REAL(8), INTENT(INOUT) :: y_cell_edge_coords(:)
      REAL(8), INTENT(INOUT) :: z_cell_edge_coords(:)
      REAL(8), INTENT(INOUT) :: y_cell_widths(:)
      REAL(8), INTENT(INOUT) :: z_cell_widths(:)
      REAL(8), INTENT(INOUT) :: p_cell_edge(:)
      REAL(8), INTENT(OUT) :: y_coarser_cell_cent_coords(:)
      REAL(8), INTENT(OUT) :: z_coarser_cell_cent_coords(:)
      REAL(8), INTENT(OUT) :: y_cell_widths_work(:)
      REAL(8), INTENT(OUT) :: z_cell_widths_work(:)

      ! Local variables
      INTEGER :: i, y_remainder, z_remainder
      INTEGER :: coarser_lvl
      INTEGER :: num_y_finer_cells, num_z_finer_cells
      INTEGER :: num_y_coarser_cells, num_z_coarser_cells

      coarser_lvl = finer_lvl + 1
      ASSOCIATE (grd_finer => solver%multigrid(finer_lvl), grd_coarser => solver%multigrid(coarser_lvl), opts => solver%options)

         num_y_finer_cells = grd_finer%cells(1)
         num_z_finer_cells = grd_finer%cells(2)

         ! Copy the finer grid cell widths to the working arrays.
         y_cell_widths_work(1:num_y_finer_cells) = y_cell_widths(1:num_y_finer_cells)
         z_cell_widths_work(1:num_z_finer_cells) = z_cell_widths(1:num_z_finer_cells)

         ! Coarsen cells in the y-direction.
         IF (num_y_finer_cells > opts%coarse_size) THEN
            y_remainder = MOD(num_y_finer_cells, 2)
            num_y_coarser_cells = (num_y_finer_cells + y_remainder)/2

            ! Update padding side for y-direction
            IF (y_remainder == 1) THEN
               IF (y_padside == 0) THEN
                  y_padside = 1
               ELSE
                  y_padside = -y_padside
               END IF
            END IF

            IF (y_remainder == 0 .OR. y_padside == 1) THEN
               ! Pad the upper end
               DO i = 1, num_y_coarser_cells - y_remainder
                  y_cell_widths(i) = y_cell_widths(2*i - 1) + y_cell_widths(2*i)
                  y_cell_edge_coords(i + 1) = y_cell_edge_coords(i) + y_cell_widths(i)
                  y_coarser_cell_cent_coords(i) = 0.5D0*(y_cell_edge_coords(i) + y_cell_edge_coords(i + 1))
               END DO
               IF (y_remainder == 1) THEN
                  i = num_y_coarser_cells
                  ! Retain the width of the last finer cell for the last coarser cell.
                  y_cell_widths(i) = y_cell_widths(num_y_finer_cells)
                  y_cell_edge_coords(i + 1) = y_cell_edge_coords(num_y_finer_cells + 1)
                  y_coarser_cell_cent_coords(i) = 0.5D0*(y_cell_edge_coords(i) + y_cell_edge_coords(i + 1))
               END IF
            ELSE
               ! Pad the lower end
               ! Retain the width of the first finer cell for the first coarser cell.
               ! y_cell_widths(1) = y_cell_widths(1)
               ! y_cell_edge_coords(2) = y_cell_edge_coords(2)
               y_coarser_cell_cent_coords(1) = 0.5D0*(y_cell_edge_coords(1) + y_cell_edge_coords(2))
               DO i = 2, num_y_coarser_cells
                  y_cell_widths(i) = y_cell_widths(2*i - 2) + y_cell_widths(2*i - 1)
                  y_cell_edge_coords(i + 1) = y_cell_edge_coords(i) + y_cell_widths(i)
                  y_coarser_cell_cent_coords(i) = 0.5D0*(y_cell_edge_coords(i) + y_cell_edge_coords(i + 1))
               END DO
            END IF
         ELSE
            num_y_coarser_cells = num_y_finer_cells
            y_coarser_cell_cent_coords(1:num_y_coarser_cells) = y_finer_cell_cent_coords(1:num_y_coarser_cells)
         END IF

         ! Coarsen cells in the z-direction.
         IF (num_z_finer_cells > opts%coarse_size) THEN
            z_remainder = MOD(num_z_finer_cells, 2)
            num_z_coarser_cells = (num_z_finer_cells + z_remainder)/2

            ! Update padding side for z-direction
            IF (z_remainder == 1) THEN
               IF (z_padside == 0) THEN
                  z_padside = 1
               ELSE
                  z_padside = -z_padside
               END IF
            END IF

            IF (z_remainder == 0 .OR. z_padside == 1) THEN
               ! Pad the upper end
               DO i = 1, num_z_coarser_cells - z_remainder
                  z_cell_widths(i) = z_cell_widths(2*i - 1) + z_cell_widths(2*i)
                  z_cell_edge_coords(i + 1) = z_cell_edge_coords(i) + z_cell_widths(i)
                  z_coarser_cell_cent_coords(i) = 0.5D0*(z_cell_edge_coords(i) + z_cell_edge_coords(i + 1))
                  ! Update background pressures at cell edges in the z-direction.
                  p_cell_edge(i) = p_cell_edge(2*i - 1)
               END DO
               IF (z_remainder == 1) THEN
                  i = num_z_coarser_cells
                  ! Retain the width of the last finer cell for the last coarser cell.
                  z_cell_widths(i) = z_cell_widths(num_z_finer_cells)
                  z_cell_edge_coords(i + 1) = z_cell_edge_coords(i) + z_cell_widths(i)
                  z_coarser_cell_cent_coords(i) = 0.5D0*(z_cell_edge_coords(i) + z_cell_edge_coords(i + 1))
                  p_cell_edge(i) = p_cell_edge(num_z_finer_cells)
               END IF
               p_cell_edge(num_z_coarser_cells + 1) = p_cell_edge(num_z_finer_cells + 1)
            ELSE
               ! Pad the lower end
               ! Retain the width of the first finer cell for the first coarser cell.
               ! z_cell_widths(1) = z_cell_widths(1)
               ! z_cell_edge_coords(2) = z_cell_edge_coords(2)
               z_coarser_cell_cent_coords(1) = 0.5D0*(z_cell_edge_coords(1) + z_cell_edge_coords(2))
               ! Update background pressures at cell edges in the z-direction.
               p_cell_edge(1) = p_cell_edge(1)
               DO i = 2, num_z_coarser_cells
                  z_cell_widths(i) = z_cell_widths(2*i - 2) + z_cell_widths(2*i - 1)
                  z_cell_edge_coords(i + 1) = z_cell_edge_coords(i) + z_cell_widths(i)
                  z_coarser_cell_cent_coords(i) = 0.5D0*(z_cell_edge_coords(i) + z_cell_edge_coords(i + 1))
                  ! Update background pressures at cell edges in the z-direction.
                  p_cell_edge(i) = p_cell_edge(2*i - 2)
               END DO
               p_cell_edge(num_z_coarser_cells + 1) = p_cell_edge(num_z_finer_cells + 1)
            END IF
         ELSE
            num_z_coarser_cells = num_z_finer_cells
            z_coarser_cell_cent_coords(1:num_z_coarser_cells) = z_finer_cell_cent_coords(1:num_z_coarser_cells)
         END IF

         ! Update the number of cells for the coarser grid.
         grd_coarser%cells(1) = num_y_coarser_cells
         grd_coarser%cells(2) = num_z_coarser_cells
         grd_coarser%cells(3) = num_y_coarser_cells*num_z_coarser_cells

         ! Compute the prolongation operator from coarser grid to finer grid.
         CALL build_prolongator(solver, coarser_lvl, y_padside, z_padside, &
                                y_cell_widths, z_cell_widths, &
                                y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, &
                                y_finer_cell_cent_coords, z_finer_cell_cent_coords)

         ! Compute the restriction operator from finer grid to coarser grid.
         CALL build_restrictor_fvm(solver, finer_lvl, y_padside, z_padside, &
                                   y_cell_widths_work, z_cell_widths_work, y_cell_widths, z_cell_widths, &
                                   y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
                                   y_coarser_cell_cent_coords, z_coarser_cell_cent_coords)
      END ASSOCIATE

   END SUBROUTINE coarsen

   !-----------------------------------------------------------------------
   ! FUNCTION: determine_GS_type
   ! PURPOSE:
   !   Determines the type of Gauss-Seidel (GS) smoothing to apply based
   !   on the grid configuration and user options.
   !
   ! INPUTS:
   !   gs_option         - CHARACTER string specifying the GS type option
   !                       (e.g., 'adapt', 'yline', 'zline').
   !   zero_alpha        - LOGICAL indicating whether the reaction coefficient
   !   multigrid_ycells - INTEGER representing the number of cells in the
   !                      first dimension of the finest grid.
   !   multigrid_zcells - INTEGER representing the number of cells in the
   !                      second dimension of the finest grid.
   !   finer_ycells     - INTEGER representing the number of cells in the
   !                      first dimension of the current finer grid.
   !   finer_zcells     - INTEGER representing the number of cells in the
   !                      second dimension of the current finer grid.
   !
   ! RETURNS:
   !   gs_type         - CHARACTER string indicating the GS type to use
   !                     ('yline', or 'zline').
   !-----------------------------------------------------------------------
   FUNCTION determine_GS_type(gs_option, zero_alpha, multigrid_ycells, multigrid_zcells, &
                              finer_ycells, finer_zcells) RESULT(gs_type)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: gs_option
      LOGICAL, INTENT(IN)    :: zero_alpha
      INTEGER, INTENT(IN)    :: multigrid_ycells
      INTEGER, INTENT(IN)    :: multigrid_zcells
      INTEGER, INTENT(IN)    :: finer_ycells
      INTEGER, INTENT(IN)    :: finer_zcells
      CHARACTER(LEN=5)        :: gs_type

      ! Check if the GS type option is set to 'adapt'
      IF (TRIM(gs_option) == 'adapt') THEN
         IF (multigrid_ycells >= multigrid_zcells ) THEN
            gs_type = 'zline'
         ELSE
            gs_type = 'yline'
         END IF
      ELSE
         ! Use the GS type specified in options, trimming any trailing spaces
         gs_type = TRIM(gs_option)
      END IF

   END FUNCTION determine_GS_type

   !-----------------------------------------------------------------------
   ! SUBROUTINE: apply_cycle
   ! PURPOSE:
   !   Performs a single V-cycle or Full Multigrid (FMG) cycle.
   !
   ! INPUTS/OUTPUTS:
   !   solver     - TYPE(MultigridSolver), input/output
   !   sol_dble   - REAL(8) array, solution vector, input/output
   !   alpha      - REAL(8), reaction coefficient
   !   fmg        - LOGICAL. If present and true, perform FMG cycle.
   !-----------------------------------------------------------------------
   SUBROUTINE apply_cycle(solver, sol_dble, alpha, fmg, mixed_precision)
      USE gmg_data, ONLY: MultigridSolver, Grid
      USE gmg_utils, ONLY: KIND_dble, KIND_sngl
      IMPLICIT NONE

      ! Input/output parameters
      TYPE(MultigridSolver), INTENT(INOUT) :: solver
      REAL(8), INTENT(IN) :: alpha
      REAL(8), INTENT(INOUT) :: sol_dble(:)
      LOGICAL, INTENT(IN) :: fmg
      LOGICAL, INTENT(IN) :: mixed_precision

      ! Local variable
      INTEGER :: i

      IF (mixed_precision) THEN
         CALL apply_cycle_sngl(solver, sol_dble, REAL(alpha), fmg)
      ELSE
         CALL apply_cycle_dble(solver, sol_dble, alpha, fmg)
      END IF
   END SUBROUTINE apply_cycle

!-----------------------------------------------------------------------
! SUBROUTINE: build_prolongator
! PURPOSE:    Computes the prolongation operator on a pair of grids. Assumes
!             cells are numbered in z-major order.
!
! INPUTS:
!   solver                     - MultigridSolver type, input/output
!   coarser_lvl                - INTEGER, coarser multigrid level
!   y_padside                  - INTEGER, padding side in y-direction
!   z_padside                  - INTEGER, padding side in z-direction
!   y_coarser_cell_widths      - REAL(8) array, cell widths of coarser grid in y
!   z_coarser_cell_widths      - REAL(8) array, cell widths of coarser grid in z
!   y_coarser_cell_cent_coords - REAL(8) array, coarser y-center coordinates
!   z_coarser_cell_cent_coords - REAL(8) array, coarser z-center coordinates
!   y_finer_cell_cent_coords   - REAL(8) array, finer y-center coordinates
!   z_finer_cell_cent_coords   - REAL(8) array, finer z-center coordinates
!
! OUTPUTS:
!   solver%multigrid(coarser_lvl)%P_row_ptr - INTEGER array, prolongation row pointers
!   solver%multigrid(coarser_lvl)%P_col_ind - INTEGER array, prolongation column indices
!   solver%multigrid(coarser_lvl)%P_val_dble - REAL(8) array, prolongation values
!   solver%multigrid(coarser_lvl)%P_val_sngl - SINGLE PRECISION array, prolongation values
!-----------------------------------------------------------------------
   SUBROUTINE build_prolongator(solver, coarser_lvl, y_padside, z_padside, &
                                y_coarser_cell_widths, z_coarser_cell_widths, y_coarser_cell_cent_coords, &
                                z_coarser_cell_cent_coords, y_finer_cell_cent_coords, z_finer_cell_cent_coords)
      IMPLICIT NONE

      ! Input/output parameters
      TYPE(MultigridSolver), INTENT(INOUT) :: solver
      INTEGER, INTENT(IN) :: coarser_lvl, y_padside, z_padside
      REAL(8), INTENT(IN) :: y_coarser_cell_widths(:), z_coarser_cell_widths(:)
      REAL(8), INTENT(IN) :: y_finer_cell_cent_coords(:), z_finer_cell_cent_coords(:)
      REAL(8), INTENT(IN) :: y_coarser_cell_cent_coords(:), z_coarser_cell_cent_coords(:)

      ! Local variables
      INTEGER :: i, y_remainder, z_remainder, finer_lvl, y_multiple, z_multiple
      INTEGER :: cell_ind, y_finer_ind, z_finer_ind
      INTEGER :: y_coarser_ind, z_coarser_ind
      INTEGER :: entry_ind, nnz, y_side, z_side
      REAL(8) :: factor, y_factor, z_factor
      INTEGER :: y_finer_cells, z_finer_cells
      INTEGER :: num_finer_cells, num_coarser_cells

      ASSOCIATE (grd_coarser => solver%multigrid(coarser_lvl), &
                 grd_finer => solver%multigrid(coarser_lvl - 1), &
                 opts => solver%options)
         finer_lvl = coarser_lvl - 1

         y_finer_cells = grd_finer%cells(1)
         z_finer_cells = grd_finer%cells(2)
         num_finer_cells = y_finer_cells*z_finer_cells
         num_coarser_cells = grd_coarser%cells(3)

         y_remainder = MOD(y_finer_cells, 2)
         z_remainder = MOD(z_finer_cells, 2)

         y_multiple = (y_finer_cells + 1)/grd_coarser%cells(1)
         z_multiple = (z_finer_cells + 1)/grd_coarser%cells(2)

         ! Allocate row pointer array for the prolongation matrix.
         ALLOCATE (grd_coarser%P_row_ptr(num_finer_cells + 1))
         grd_coarser%P_row_ptr(1) = 1

         ! Compute P_row_ptr considering multiple factors for coarsening.
         DO y_finer_ind = 1, y_finer_cells
            DO z_finer_ind = 1, z_finer_cells
               cell_ind = (y_finer_ind - 1)*z_finer_cells + z_finer_ind
               nnz = 4  ! Default number of non-zero entries

               ! Adjust nnz for boundaries and corners
               IF (y_finer_cell_cent_coords(y_finer_ind) <= y_coarser_cell_cent_coords(1) .OR. &
                   y_finer_cell_cent_coords(y_finer_ind) >= y_coarser_cell_cent_coords(grd_coarser%cells(1)) .OR. &
                   y_multiple == 1) THEN
                  nnz = nnz/2
               END IF
               IF (z_finer_cell_cent_coords(z_finer_ind) <= z_coarser_cell_cent_coords(1) .OR. &
                   z_finer_cell_cent_coords(z_finer_ind) >= z_coarser_cell_cent_coords(grd_coarser%cells(2)) .OR. &
                   z_multiple == 1) THEN
                  nnz = nnz/2
               END IF

               grd_coarser%P_row_ptr(cell_ind + 1) = grd_coarser%P_row_ptr(cell_ind) + nnz
            END DO
         END DO

         ! Allocate column indices and values arrays for prolongation matrix.
         ALLOCATE (grd_coarser%P_col_ind(grd_coarser%P_row_ptr(num_finer_cells + 1) - 1))
         ALLOCATE (grd_coarser%P_val_dble(grd_coarser%P_row_ptr(num_finer_cells + 1) - 1))

         ! Fill P_col_ind and P_val_dble for the prolongation matrix.
         DO y_finer_ind = 1, y_finer_cells
            IF (y_remainder == 0 .OR. y_padside == 1 .OR. y_multiple == 1) THEN
               y_coarser_ind = (y_finer_ind - 1)/y_multiple + 1
            ELSE
               ! Pad the lower end
               IF (y_finer_ind == 1) THEN
                  y_coarser_ind = 1
               ELSE
                  y_coarser_ind = y_finer_ind/y_multiple + 1
               END IF
            END IF
            IF (y_finer_cell_cent_coords(y_finer_ind) >= y_coarser_cell_cent_coords(grd_coarser%cells(1))) THEN
               y_coarser_ind = grd_coarser%cells(1)
            END IF

            DO z_finer_ind = 1, z_finer_cells
               IF (z_remainder == 0 .OR. z_padside == 1 .OR. z_multiple == 1) THEN
                  z_coarser_ind = (z_finer_ind - 1)/z_multiple + 1
               ELSE
                  ! Pad the lower end
                  IF (z_finer_ind == 1) THEN
                     z_coarser_ind = 1
                  ELSE
                     z_coarser_ind = z_finer_ind/z_multiple + 1
                  END IF
               END IF
               IF (z_finer_cell_cent_coords(z_finer_ind) >= z_coarser_cell_cent_coords(grd_coarser%cells(2))) THEN
                  z_coarser_ind = grd_coarser%cells(2)
               END IF

               cell_ind = (y_finer_ind - 1)*z_finer_cells + z_finer_ind
               entry_ind = grd_coarser%P_row_ptr(cell_ind)
               nnz = grd_coarser%P_row_ptr(cell_ind + 1) - grd_coarser%P_row_ptr(cell_ind)

               SELECT CASE (nnz)
               CASE (1)
                  ! Injection for a single neighboring cell (corner).
                  grd_coarser%P_col_ind(entry_ind) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                  grd_coarser%P_val_dble(entry_ind) = 1.0D0

               CASE (2)
                  ! Linear interpolation for edge cells.
                  IF (y_finer_cell_cent_coords(y_finer_ind) <= y_coarser_cell_cent_coords(1) .OR. &
                      y_finer_cell_cent_coords(y_finer_ind) >= y_coarser_cell_cent_coords(grd_coarser%cells(1)) .OR. &
                      y_multiple == 1) THEN
                     ! Interpolate along z-direction
                     IF (z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(z_coarser_ind)) THEN
                        z_side = 1
                        factor = 2.0D0
                     ELSE
                        z_side = -1
                        factor = -2.0D0
                     END IF
                     z_factor = factor*(z_finer_cell_cent_coords(z_finer_ind) - z_coarser_cell_cent_coords(z_coarser_ind))/ &
                                (z_coarser_cell_widths(z_coarser_ind) + z_coarser_cell_widths(z_coarser_ind + z_side))
                     IF (debugging .AND. (z_factor > 1 .OR. z_factor < 0)) THEN
                        WRITE (*, *) 'Error: z_factor = ', z_factor, 'is out of range'
                     END IF

                     grd_coarser%P_col_ind(entry_ind) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                     grd_coarser%P_val_dble(entry_ind) = 1.0D0 - z_factor

                     grd_coarser%P_col_ind(entry_ind + 1) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                     grd_coarser%P_val_dble(entry_ind + 1) = z_factor
                  ELSE
                     IF (debugging .AND. z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(1) .AND. &
                         z_finer_cell_cent_coords(z_finer_ind) < z_coarser_cell_cent_coords(grd_coarser%cells(2)) .AND. &
                         z_multiple > 1) THEN
                        WRITE (*, *) 'Error: Invalid edge case at z_finer_ind = ', z_finer_ind
                     END IF

                     IF (y_finer_cell_cent_coords(y_finer_ind) > y_coarser_cell_cent_coords(y_coarser_ind)) THEN
                        y_side = 1
                        factor = 2.0D0
                     ELSE
                        y_side = -1
                        factor = -2.0D0
                     END IF

                     ! Interpolate along y-direction
                     y_factor = factor*(y_finer_cell_cent_coords(y_finer_ind) - y_coarser_cell_cent_coords(y_coarser_ind))/ &
                                (y_coarser_cell_widths(y_coarser_ind) + y_coarser_cell_widths(y_coarser_ind + y_side))

                     IF (debugging .AND. (y_factor > 1 .OR. y_factor < 0)) THEN
                        WRITE (*, *) 'Error: y_factor = ', y_factor, 'is out of range'
                     END IF
                     grd_coarser%P_col_ind(entry_ind) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                     grd_coarser%P_val_dble(entry_ind) = 1.0D0 - y_factor

                     grd_coarser%P_col_ind(entry_ind + 1) = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind
                     grd_coarser%P_val_dble(entry_ind + 1) = y_factor
                  END IF
               CASE DEFAULT
                  IF (debugging .AND. nnz /= 4) THEN
                     WRITE (*, *) 'Error: Invalid number of non-zero entries = ', nnz
                     STOP
                  END IF
                  ! Bilinear interpolation for four neighboring cells.
                  IF (z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(z_coarser_ind)) THEN
                     z_side = 1
                     factor = 2.0D0
                  ELSE
                     z_side = -1
                     factor = -2.0D0
                  END IF
                  z_factor = factor*(z_finer_cell_cent_coords(z_finer_ind) - z_coarser_cell_cent_coords(z_coarser_ind))/ &
                             (z_coarser_cell_widths(z_coarser_ind) + z_coarser_cell_widths(z_coarser_ind + z_side))
                  IF (debugging .AND. (z_factor > 1 .OR. z_factor < 0)) THEN
                     WRITE (*, *) 'Error: z_factor = ', z_factor, 'is out of range'
                  END IF

                  IF (y_finer_cell_cent_coords(y_finer_ind) > y_coarser_cell_cent_coords(y_coarser_ind)) THEN
                     y_side = 1
                     factor = 2.0D0
                  ELSE
                     y_side = -1
                     factor = -2.0D0
                  END IF

                  y_factor = factor*(y_finer_cell_cent_coords(y_finer_ind) - y_coarser_cell_cent_coords(y_coarser_ind))/ &
                             (y_coarser_cell_widths(y_coarser_ind) + y_coarser_cell_widths(y_coarser_ind + y_side))
                  IF (debugging .AND. (y_factor > 1 .OR. y_factor < 0)) THEN
                     WRITE (*, *) 'Error: y_factor = ', y_factor, 'is out of range'
                  END IF

                  ! Base corner.
                  grd_coarser%P_col_ind(entry_ind) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                  grd_coarser%P_val_dble(entry_ind) = (1 - y_factor)*(1 - z_factor)

                  ! Y-shifted corner.
                  grd_coarser%P_col_ind(entry_ind + 1) = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind
                  grd_coarser%P_val_dble(entry_ind + 1) = y_factor*(1 - z_factor)

                  ! Z-shifted corner.
                  grd_coarser%P_col_ind(entry_ind + 2) = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                  grd_coarser%P_val_dble(entry_ind + 2) = (1 - y_factor)*z_factor

                  ! Y- and Z-shifted corner.
                  grd_coarser%P_col_ind(entry_ind + 3) = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                  grd_coarser%P_val_dble(entry_ind + 3) = y_factor*z_factor
               END SELECT
            END DO
         END DO

         ! Convert double precision values to single precision if needed.
         ALLOCATE (grd_coarser%P_val_sngl(grd_coarser%P_row_ptr(num_finer_cells + 1) - 1))
         grd_coarser%P_val_sngl = REAL(grd_coarser%P_val_dble)

      END ASSOCIATE

   END SUBROUTINE build_prolongator

!-----------------------------------------------------------------------
! SUBROUTINE: build_restrictor_fvm
! PURPOSE:    Computes the restriction operator on a pair of grids in FVM. Assumes
!             cells are numbered in z-major order. The restriction operator
!             has a non-zero pattern equal to the transpose of the prolongation
!             operator but uses different scaling factors appropriate for
!             linear and bilinear interpolation. The scaling factors are computed
!             as the area ratio of the finer cell to the coarser cell since we
!             are using FVM and multiplying each row by the cell area recovers
!             partition of unity of test functions.
!
! INPUTS:
!   solver                     - MultigridSolver type, input/output
!   finer_lvl                  - INTEGER, finer multigrid level
!   y_padside                  - INTEGER, padding side in y-direction (-1, 0, or 1)
!   z_padside                  - INTEGER, padding side in z-direction (-1, 0, or 1)
!   y_finer_cell_cent_coords   - REAL(8) array, finer y-center coordinates
!   z_finer_cell_cent_coords   - REAL(8) array, finer z-center coordinates
!   y_coarser_cell_cent_coords - REAL(8) array, coarser y-center coordinates
!   z_coarser_cell_cent_coords - REAL(8) array, coarser z-center coordinates
!
! OUTPUTS:
!   solver%multigrid(coarser_lvl)%R_row_ptr   - INTEGER array, restriction row pointers
!   solver%multigrid(coarser_lvl)%R_col_ind   - INTEGER array, restriction column indices
!   solver%multigrid(coarser_lvl)%R_val_dble  - REAL(8) array, restriction values
!   solver%multigrid(coarser_lvl)%R_val_sngl  - SINGLE PRECISION array, restriction values (if mixed precision)
!-----------------------------------------------------------------------
   SUBROUTINE build_restrictor_fvm(solver, finer_lvl, y_padside, z_padside, &
                                   y_finer_cell_widths, z_finer_cell_widths, &
                                   y_coarser_cell_widths, z_coarser_cell_widths, &
                                   y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
                                   y_coarser_cell_cent_coords, z_coarser_cell_cent_coords)
      IMPLICIT NONE

      ! Input/output parameters
      TYPE(MultigridSolver), INTENT(INOUT) :: solver
      INTEGER, INTENT(IN) :: finer_lvl, y_padside, z_padside
      REAL(8), INTENT(IN) :: y_finer_cell_widths(:), z_finer_cell_widths(:)
      REAL(8), INTENT(IN) :: y_coarser_cell_widths(:), z_coarser_cell_widths(:)
      REAL(8), INTENT(IN) :: y_finer_cell_cent_coords(:), z_finer_cell_cent_coords(:)
      REAL(8), INTENT(IN) :: y_coarser_cell_cent_coords(:), z_coarser_cell_cent_coords(:)

      ! Local variables
      INTEGER :: i, j, k
      INTEGER :: coarser_lvl
      INTEGER :: y_finer_cells, z_finer_cells, num_finer_cells
      INTEGER :: y_coarser_cells, z_coarser_cells, num_coarser_cells
      INTEGER :: y_multiple, z_multiple, y_remainder, z_remainder
      INTEGER :: y_finer_ind, z_finer_ind
      INTEGER :: y_coarser_ind, z_coarser_ind
      INTEGER :: cell_ind_coarse, cell_ind_fine
      INTEGER :: entry_ind, entry_start, entry_end
      REAL(8) :: y_factor, z_factor, y_overlap, z_overlap
      REAL(8) :: area_total, area_coarse, area_fine
      INTEGER :: nnz, num_entries
      INTEGER :: y_side, z_side
      INTEGER :: idx, m, n
      INTEGER :: y_fine, z_fine
      INTEGER, DIMENSION(2) :: y_f_inds, z_f_inds
      REAL(8) :: y_diff, z_diff, y_range, z_range
      REAL(8) :: weight, total_weight

      ! Access grids at coarser and finer levels
      coarser_lvl = finer_lvl + 1

      ASSOCIATE (grd_coarser => solver%multigrid(coarser_lvl), &
                 grd_finer => solver%multigrid(finer_lvl), &
                 opts => solver%options)

         y_finer_cells = grd_finer%cells(1)
         z_finer_cells = grd_finer%cells(2)
         num_finer_cells = y_finer_cells*z_finer_cells

         y_coarser_cells = grd_coarser%cells(1)
         z_coarser_cells = grd_coarser%cells(2)
         num_coarser_cells = y_coarser_cells*z_coarser_cells

         y_remainder = MOD(y_finer_cells, 2)
         z_remainder = MOD(z_finer_cells, 2)

         y_multiple = (y_finer_cells + 1)/(y_coarser_cells)
         z_multiple = (z_finer_cells + 1)/(z_coarser_cells)

         num_entries = SIZE(grd_coarser%P_col_ind)

         ! Allocate row pointer array for the restriction operator.
         ALLOCATE (grd_finer%R_row_ptr(num_coarser_cells + 1))
         ALLOCATE (grd_finer%R_col_ind(num_entries))
         ALLOCATE (grd_finer%R_val_dble(num_entries))

         grd_finer%R_row_ptr = 0

         ! Compute number of nonzeros per row in the restriction operator.
         DO j = 1, num_finer_cells
            DO k = grd_coarser%P_row_ptr(j), grd_coarser%P_row_ptr(j + 1) - 1
               i = grd_coarser%P_col_ind(k)
               grd_finer%R_row_ptr(i + 1) = grd_finer%R_row_ptr(i + 1) + 1
            END DO
         END DO

         ! Compute R_row_ptr
         grd_finer%R_row_ptr(1) = 1
         DO i = 1, num_coarser_cells
            grd_finer%R_row_ptr(i + 1) = grd_finer%R_row_ptr(i) + grd_finer%R_row_ptr(i + 1)
         END DO

         ! Fill R_col_ind and R_val_dble for the restriction oprator using R_row_ptr as index.
         DO y_finer_ind = 1, y_finer_cells
            IF (y_remainder == 0 .OR. y_padside == 1 .OR. y_multiple == 1) THEN
               y_coarser_ind = (y_finer_ind - 1)/y_multiple + 1
            ELSE
               ! Pad the lower end
               IF (y_finer_ind == 1) THEN
                  y_coarser_ind = 1
               ELSE
                  y_coarser_ind = y_finer_ind/y_multiple + 1
               END IF
            END IF
            IF (y_finer_cell_cent_coords(y_finer_ind) >= y_coarser_cell_cent_coords(grd_coarser%cells(1))) THEN
               y_coarser_ind = grd_coarser%cells(1)
            END IF

            DO z_finer_ind = 1, z_finer_cells
               IF (z_remainder == 0 .OR. z_padside == 1 .OR. z_multiple == 1) THEN
                  z_coarser_ind = (z_finer_ind - 1)/z_multiple + 1
               ELSE
                  ! Pad the lower end
                  IF (z_finer_ind == 1) THEN
                     z_coarser_ind = 1
                  ELSE
                     z_coarser_ind = z_finer_ind/z_multiple + 1
                  END IF
               END IF
               IF (z_finer_cell_cent_coords(z_finer_ind) >= z_coarser_cell_cent_coords(grd_coarser%cells(2))) THEN
                  z_coarser_ind = grd_coarser%cells(2)
               END IF

               cell_ind_fine = (y_finer_ind - 1)*z_finer_cells + z_finer_ind
               nnz = grd_coarser%P_row_ptr(cell_ind_fine + 1) - grd_coarser%P_row_ptr(cell_ind_fine)
               area_fine = (y_finer_cell_widths(y_finer_ind)*z_finer_cell_widths(z_finer_ind))

               SELECT CASE (nnz)
               CASE (1)
                  cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                  entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)

                  ! Injection for a single neighboring cell (corner).
                  grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                  area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind)
                  grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse
                  grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1
               CASE (2)
                  ! Linear interpolation for edge cells.
                  IF (y_finer_cell_cent_coords(y_finer_ind) <= y_coarser_cell_cent_coords(1) .OR. &
                      y_finer_cell_cent_coords(y_finer_ind) >= y_coarser_cell_cent_coords(grd_coarser%cells(1)) .OR. &
                      y_multiple == 1) THEN
                     ! Interpolate along z-direction
                     IF (z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(z_coarser_ind)) THEN
                        z_side = 1
                     ELSE
                        z_side = -1
                     END IF
                     z_factor = (z_finer_cell_cent_coords(z_finer_ind) - z_coarser_cell_cent_coords(z_coarser_ind))/ &
                                (z_coarser_cell_cent_coords(z_coarser_ind + z_side) - z_coarser_cell_cent_coords(z_coarser_ind))
                     IF (debugging .AND. (z_factor > 1 .OR. z_factor < 0)) THEN
                        WRITE (*, *) 'Error: z_factor = ', z_factor, 'is out of range'
                     END IF

                     cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                     entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                     grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                     area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind)
                     grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*(1.0D0 - z_factor)
                     grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1

                     cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                     entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                     grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                     area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind + z_side)
                     grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*z_factor
                     grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1
                  ELSE
                     IF (debugging .AND. z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(1) .AND. &
                         z_finer_cell_cent_coords(z_finer_ind) < z_coarser_cell_cent_coords(grd_coarser%cells(2)) .AND. &
                         z_multiple > 1) THEN
                        WRITE (*, *) 'Error: Invalid edge case at z_finer_ind = ', z_finer_ind
                     END IF

                     IF (y_finer_cell_cent_coords(y_finer_ind) > y_coarser_cell_cent_coords(y_coarser_ind)) THEN
                        y_side = 1
                     ELSE
                        y_side = -1
                     END IF

                     ! Interpolate along y-direction
                     y_factor = (y_finer_cell_cent_coords(y_finer_ind) - y_coarser_cell_cent_coords(y_coarser_ind))/ &
                                (y_coarser_cell_cent_coords(y_coarser_ind + y_side) - y_coarser_cell_cent_coords(y_coarser_ind))

                     IF (debugging .AND. (y_factor > 1 .OR. y_factor < 0)) THEN
                        WRITE (*, *) 'Error: y_factor = ', y_factor, 'is out of range'
                     END IF

                     cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                     entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                     grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                     area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind)
                     grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*(1.0D0 - y_factor)
                     grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1

                     cell_ind_coarse = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind
                     entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                     grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                     area_coarse = y_coarser_cell_widths(y_coarser_ind + y_side)*z_coarser_cell_widths(z_coarser_ind)
                     grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*y_factor
                     grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1
                  END IF
               CASE DEFAULT
                  IF (debugging .AND. nnz /= 4) THEN
                     WRITE (*, *) 'Error: Invalid number of non-zero entries = ', nnz
                     STOP
                  END IF
                  ! Bilinear interpolation for four neighboring cells.
                  IF (z_finer_cell_cent_coords(z_finer_ind) > z_coarser_cell_cent_coords(z_coarser_ind)) THEN
                     z_side = 1
                  ELSE
                     z_side = -1
                  END IF
                  z_factor = (z_finer_cell_cent_coords(z_finer_ind) - z_coarser_cell_cent_coords(z_coarser_ind))/ &
                             (z_coarser_cell_cent_coords(z_coarser_ind + z_side) - z_coarser_cell_cent_coords(z_coarser_ind))
                  IF (debugging .AND. (z_factor > 1 .OR. z_factor < 0)) THEN
                     WRITE (*, *) 'Error: z_factor = ', z_factor, 'is out of range'
                  END IF

                  IF (y_finer_cell_cent_coords(y_finer_ind) > y_coarser_cell_cent_coords(y_coarser_ind)) THEN
                     y_side = 1
                  ELSE
                     y_side = -1
                  END IF

                  y_factor = (y_finer_cell_cent_coords(y_finer_ind) - y_coarser_cell_cent_coords(y_coarser_ind))/ &
                             (y_coarser_cell_cent_coords(y_coarser_ind + y_side) - y_coarser_cell_cent_coords(y_coarser_ind))
                  IF (debugging .AND. (y_factor > 1 .OR. y_factor < 0)) THEN
                     WRITE (*, *) 'Error: y_factor = ', y_factor, 'is out of range'
                  END IF

                  ! Base corner.
                  cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind
                  entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                  grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                  area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind)
                  grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*(1 - y_factor)*(1 - z_factor)
                  grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1

                  ! Y-shifted corner.
                  cell_ind_coarse = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind
                  entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                  grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                  area_coarse = y_coarser_cell_widths(y_coarser_ind + y_side)*z_coarser_cell_widths(z_coarser_ind)
                  grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*y_factor*(1 - z_factor)
                  grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1

                  ! Z-shifted corner.
                  cell_ind_coarse = (y_coarser_ind - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                  entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                  grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                  area_coarse = y_coarser_cell_widths(y_coarser_ind)*z_coarser_cell_widths(z_coarser_ind + z_side)
                  grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*(1 - y_factor)*z_factor
                  grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1

                  ! Y- and Z-shifted corner.
                  cell_ind_coarse = (y_coarser_ind + y_side - 1)*grd_coarser%cells(2) + z_coarser_ind + z_side
                  entry_ind = grd_finer%R_row_ptr(cell_ind_coarse)
                  grd_finer%R_col_ind(entry_ind) = cell_ind_fine
                  area_coarse = y_coarser_cell_widths(y_coarser_ind + y_side)*z_coarser_cell_widths(z_coarser_ind + z_side)
                  grd_finer%R_val_dble(entry_ind) = area_fine/area_coarse*y_factor*z_factor
                  grd_finer%R_row_ptr(cell_ind_coarse) = grd_finer%R_row_ptr(cell_ind_coarse) + 1
               END SELECT
            END DO
         END DO

         ! Recover the correct row pointers
         DO i = num_coarser_cells, 1, -1
            grd_finer%R_row_ptr(i + 1) = grd_finer%R_row_ptr(i)
         END DO
         grd_finer%R_row_ptr(1) = 1

         ! Convert double precision values to single precision if needed.
         ALLOCATE (grd_finer%R_val_sngl(num_entries))
         grd_finer%R_val_sngl = REAL(grd_finer%R_val_dble)
      END ASSOCIATE

   END SUBROUTINE build_restrictor_fvm

#define USE_DOUBLE_PRECISION
#include "gmg_kernel.inc"
#undef USE_DOUBLE_PRECISION

#define USE_SINGLE_PRECISION
#include "gmg_kernel.inc"
#undef USE_SINGLE_PRECISION

END MODULE gmg_kernel
