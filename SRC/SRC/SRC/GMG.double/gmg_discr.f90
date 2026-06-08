!-----------------------------------------------------------------------------
!
! This module handles the discretization on each level.
!
!-----------------------------------------------------------------------------

MODULE gmg_discr
   USE gmg_data, ONLY: Grid, Options
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: discretize

CONTAINS

   !-----------------------------------------------------------------------
   ! SUBROUTINE: discretize
   ! PURPOSE:    Computes the PDE discretization matrix on one grid. The y and
   !             z directions correspond to latitude and altitude, respectively.
   !             The cells are numbered in z-major order.
   !
   ! INPUTS:
   !   grd                  - GRID, current multigrid level (input/output)
   !   opts                 - OPTIONS, solver options
   !   y_cell_edge_coords   - DOUBLE PRECISION array, y-edge coordinates
   !   p_cell_edge          - DOUBLE PRECISION array, pressure cell edges
   !   y_cell_widths        - DOUBLE PRECISION array, cell widths in y
   !   z_cell_widths        - DOUBLE PRECISION array, cell widths in z
   !   y_cell_cent_coords   - DOUBLE PRECISION array, y-center coordinates
   !   mu_discr             - CHARACTER, discretization type for mu
   !   rho_discr            - CHARACTER, discretization type for rho
   !   mu_is_const          - LOGICAL, flag for constant mu
   !   rad_reduc_factor     - DOUBLE PRECISION, radial reduction factor
   !
   ! OUTPUTS:
   !   grd%A_row_ptr        - INTEGER array, CRS row pointers
   !   grd%A_col_ind        - INTEGER array, CRS column indices
   !   grd%A_val_dble       - DOUBLE PRECISION array, CRS values
   !   grd%A_diag_val_dble  - DOUBLE PRECISION array, diagonal values (if mixed_precision)
   !-----------------------------------------------------------------------
   SUBROUTINE discretize(grd, opts, y_cell_edge_coords, p_cell_edge, &
                         y_cell_widths, z_cell_widths, y_cell_cent_coords, &
                         mu_discr, rho_discr, mu_is_const, rad_reduc_factor)
      IMPLICIT NONE

      ! Input/output parameters
      TYPE(Grid), INTENT(INOUT) :: grd
      TYPE(Options), INTENT(IN) :: opts
      DOUBLE PRECISION, INTENT(IN) :: y_cell_edge_coords(grd%cells(1) + 1)
      DOUBLE PRECISION, INTENT(IN) :: p_cell_edge(grd%cells(2) + 1)
      DOUBLE PRECISION, INTENT(IN) :: y_cell_widths(grd%cells(1)), z_cell_widths(grd%cells(2))
      DOUBLE PRECISION, INTENT(IN) :: y_cell_cent_coords(grd%cells(1))
      CHARACTER(len=1), INTENT(IN) :: mu_discr
      CHARACTER(len=1), INTENT(IN) :: rho_discr
      LOGICAL, INTENT(IN) :: mu_is_const
      DOUBLE PRECISION, INTENT(IN) :: rad_reduc_factor

      ! Local variables
      INTEGER :: i, j
      DOUBLE PRECISION, ALLOCATABLE :: mu_cent(:), mu_edge(:)
      DOUBLE PRECISION, ALLOCATABLE :: rho_cent(:)
      INTEGER :: cent_ind
      INTEGER :: ind
      INTEGER :: total_nnz
      DOUBLE PRECISION :: cent_val
      INTEGER :: val_cent_ind

      ! Constants
      INTEGER :: n_cells_y, n_cells_z, n_total_cells

      ! Assign sizes
      n_cells_y = grd%cells(1)
      n_cells_z = grd%cells(2)
      n_total_cells = n_cells_y*n_cells_z

      ! Allocate and compute diffusion coefficients at cell centers.
      ALLOCATE (mu_cent(n_cells_y))
      DO i = 1, n_cells_y
         mu_cent(i) = y_cell_cent_coords(i)
      END DO
      CALL diffusion_coeff(n_cells_y, mu_is_const, opts%rad_earth, rad_reduc_factor, mu_cent)

      ! Allocate and compute air densities at cell centers.
      ALLOCATE (rho_cent(n_cells_z))
      DO i = 1, n_cells_z
         rho_cent(i) = -(p_cell_edge(i + 1) - p_cell_edge(i))/z_cell_widths(i)
      END DO

      ! Allocate and compute diffusion coefficients at cell edges.
      ALLOCATE (mu_edge(n_cells_y + 1))
      DO i = 1, n_cells_y + 1
         mu_edge(i) = y_cell_edge_coords(i)
      END DO
      CALL diffusion_coeff(n_cells_y + 1, mu_is_const, opts%rad_earth, rad_reduc_factor, mu_edge)

      ! Allocate row pointer array for CRS format.
      ALLOCATE (grd%A_row_ptr(grd%cells(3) + 1))
      DO i = 1, grd%cells(3) + 1
         grd%A_row_ptr(i) = 1
      END DO
      cent_ind = 1

      ! Compute the number of non-zero entries in each row.
      DO i = 1, n_cells_y
         DO j = 1, n_cells_z
            ! Check for neighboring cells and increment row_nnz accordingly.
            IF (i > 1) THEN
               grd%A_row_ptr(cent_ind + 1) = grd%A_row_ptr(cent_ind + 1) + 1
            END IF
            IF (i < n_cells_y) THEN
               grd%A_row_ptr(cent_ind + 1) = grd%A_row_ptr(cent_ind + 1) + 1
            END IF
            IF (j > 1) THEN
               grd%A_row_ptr(cent_ind + 1) = grd%A_row_ptr(cent_ind + 1) + 1
            END IF
            IF (j < n_cells_z) THEN
               grd%A_row_ptr(cent_ind + 1) = grd%A_row_ptr(cent_ind + 1) + 1
            END IF
            cent_ind = cent_ind + 1
         END DO
      END DO

      ! Adjust row pointers to get actual indices.
      DO i = 1, grd%cells(3)
         grd%A_row_ptr(i + 1) = grd%A_row_ptr(i + 1) + grd%A_row_ptr(i)
      END DO
      total_nnz = grd%A_row_ptr(n_total_cells + 1) - 1

      ! Allocate CRS matrix arrays.
      ALLOCATE (grd%A_col_ind(total_nnz))
      ALLOCATE (grd%A_val_dble(total_nnz))
      ALLOCATE (grd%A_diag_val_dble(n_total_cells))

      ! Initialize cent_ind for cell indexing.
      ind = 1
      cent_ind = 1

      ! Fill column indices and values for the CRS matrix.
      DO i = 1, n_cells_y
         DO j = 1, n_cells_z
            grd%A_col_ind(ind) = cent_ind
            cent_val = 0.0D+0
            val_cent_ind = ind

            SELECT CASE (rho_discr)
            CASE ('B')
               ! Set the column index and entry value for the cell below.
               IF (j > 1) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind - 1
                  grd%A_val_dble(ind) = 2.0D+0*mu_cent(i)**2* &
                                        (rho_cent(j - 1)*z_cell_widths(j) + rho_cent(j)*z_cell_widths(j - 1))/ &
                                        ((z_cell_widths(j - 1) + z_cell_widths(j))**2*z_cell_widths(j)*rho_cent(j))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
               ! Set the column index and entry value for the cell above.
               IF (j < n_cells_z) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind + 1
                  grd%A_val_dble(ind) = 2.0D+0*mu_cent(i)**2* &
                                        (rho_cent(j)*z_cell_widths(j + 1) + rho_cent(j + 1)*z_cell_widths(j))/ &
                                        ((z_cell_widths(j) + z_cell_widths(j + 1))**2*z_cell_widths(j)*rho_cent(j))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
            CASE DEFAULT
               ! Set the entry values based on the defined discretization.
               IF (j > 1) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind - 1
                  grd%A_val_dble(ind) = mu_cent(i)**2* &
                                        (rho_cent(j - 1) + rho_cent(j))/ &
                                        ((z_cell_widths(j - 1) + z_cell_widths(j))*z_cell_widths(j)*rho_cent(j))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
               IF (j < n_cells_z) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind + 1
                  grd%A_val_dble(ind) = mu_cent(i)**2* &
                                        (rho_cent(j) + rho_cent(j + 1))/ &
                                        ((z_cell_widths(j) + z_cell_widths(j + 1))*z_cell_widths(j)*rho_cent(j))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
            END SELECT

            ! Set the entry values based on the defined diffusion coefficient.
            SELECT CASE (mu_discr)
            CASE ('B')
               ! Set the column index and entry value for the cell to the left.
               IF (i > 1) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind - n_cells_z
                  grd%A_val_dble(ind) = 2.0D+0*(mu_cent(i - 1)*y_cell_widths(i) + &
                                                mu_cent(i)*y_cell_widths(i - 1))*mu_cent(i)/ &
                                        ((y_cell_widths(i - 1) + y_cell_widths(i))**2*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
               ! Set the column index and entry value for the cell to the right.
               IF (i < n_cells_y) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind + n_cells_z
                  grd%A_val_dble(ind) = 2.0D+0*(mu_cent(i)*y_cell_widths(i + 1) + &
                                                mu_cent(i + 1)*y_cell_widths(i))*mu_cent(i)/ &
                                        ((y_cell_widths(i) + y_cell_widths(i + 1))**2*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
            CASE ('C')
               IF (i > 1) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind - n_cells_z
                  grd%A_val_dble(ind) = mu_edge(i)*mu_cent(i)*2.0D+0/ &
                                        ((y_cell_widths(i - 1) + y_cell_widths(i))*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
               IF (i < n_cells_y) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind + n_cells_z
                  grd%A_val_dble(ind) = mu_edge(i + 1)*mu_cent(i)*2.0D+0/ &
                                        ((y_cell_widths(i) + y_cell_widths(i + 1))*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
            CASE DEFAULT
               IF (i > 1) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind - n_cells_z
                  grd%A_val_dble(ind) = (mu_cent(i - 1) + mu_cent(i))*mu_cent(i)/ &
                                        ((y_cell_widths(i - 1) + y_cell_widths(i))*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
               IF (i < n_cells_y) THEN
                  ind = ind + 1
                  grd%A_col_ind(ind) = cent_ind + n_cells_z
                  grd%A_val_dble(ind) = (mu_cent(i) + mu_cent(i + 1))*mu_cent(i)/ &
                                        ((y_cell_widths(i) + y_cell_widths(i + 1))*y_cell_widths(i))
                  cent_val = cent_val - grd%A_val_dble(ind)
               END IF
            END SELECT

            ! Set the diagonal value.
            grd%A_val_dble(val_cent_ind) = cent_val
            grd%A_diag_val_dble(cent_ind) = cent_val
            ind = ind + 1
            cent_ind = cent_ind + 1
         END DO
      END DO

      ! Deallocate temporary arrays.
      DEALLOCATE (mu_cent, rho_cent, mu_edge)

   END SUBROUTINE discretize

!-----------------------------------------------------------------------
! Subroutine: diffusion_coeff
! Purpose:    Computes the diffusion coefficients, either setting them to a constant
!             or applying a radial reduction based on the Earth's radius and reduction factor.
! Inputs:
!   length           - INTEGER, Length of the grid (number of grid points).
!   mu_is_const      - LOGICAL, If TRUE, all coefficients are set to 1.
!   rad_earth        - DOUBLE PRECISION, Earth's radius.
!   rad_reduc_factor - DOUBLE PRECISION, Scaling factor for radial reduction.
! Outputs:
!   mu               - DOUBLE PRECISION array, Updated diffusion coefficients (size length).
!-----------------------------------------------------------------------
   SUBROUTINE diffusion_coeff(length, mu_is_const, rad_earth, rad_reduc_factor, mu)
      IMPLICIT NONE

      ! Input/output parameters
      INTEGER, INTENT(IN) :: length
      LOGICAL, INTENT(IN) :: mu_is_const
      DOUBLE PRECISION, INTENT(IN) :: rad_earth
      DOUBLE PRECISION, INTENT(IN) :: rad_reduc_factor
      DOUBLE PRECISION, INTENT(INOUT) :: mu(length)

      ! Constants
      DOUBLE PRECISION :: irad

      ! Local variables
      INTEGER :: i

      irad = 1.0D0/rad_earth

      ! Apply constant or variable diffusion coefficients.
      IF (mu_is_const) THEN
         mu = 1.0D+0
      ELSE
         DO i = 1, length
            mu(i) = COS(mu(i)*irad*rad_reduc_factor)
         END DO
      END IF

   END SUBROUTINE diffusion_coeff
END MODULE gmg_discr
