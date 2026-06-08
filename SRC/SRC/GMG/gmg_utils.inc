#ifdef USE_DOUBLE_PRECISION

#define KIND_generic KIND_dble
#define GS_generic GS_dble
#define tridiag_solv_zline_generic tridiag_solv_zline_dble
#define tridiag_solv_yline_generic tridiag_solv_yline_dble
#define resid_generic resid_dble
#define mat_vec_mult_crs_generic mat_vec_mult_crs_dble
#define trans_mat_vec_mult_crs_generic trans_mat_vec_mult_crs_dble
#define mat_vec_mult_add_crs_generic mat_vec_mult_add_crs_dble
#define apply_scaled_prol_generic apply_scaled_prol_dble

! Private functions
#define GS_zline_generic GS_zline_dble
#define GS_yline_generic GS_yline_dble
#define GS_point_generic GS_point_dble

#elif defined(USE_SINGLE_PRECISION)

#define KIND_generic KIND_sngl
#define GS_generic GS_sngl
#define tridiag_solv_zline_generic tridiag_solv_zline_sngl
#define tridiag_solv_yline_generic tridiag_solv_yline_sngl
#define resid_generic resid_sngl
#define mat_vec_mult_crs_generic mat_vec_mult_crs_sngl
#define trans_mat_vec_mult_crs_generic trans_mat_vec_mult_crs_sngl
#define mat_vec_mult_add_crs_generic mat_vec_mult_add_crs_sngl
#define apply_scaled_prol_generic apply_scaled_prol_sngl

! Private functions
#define GS_zline_generic GS_zline_sngl
#define GS_yline_generic GS_yline_sngl
#define GS_point_generic GS_point_sngl

#endif

!-----------------------------------------------------------------------
! subroutine: GS_generic
! Purpose:    Applies Gauss-Seidel iteration on a Helmholtz discretization matrix in given
!             precision. Supports 'zline', 'yline', and 'point' types of GS relaxation.
! Inputs:
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size y_cells*z_cells + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   rhs            - REAL array, Right-hand-side vector (size y_cells*z_cells).
!   GS_type        - CHARACTER, Type of Gauss-Seidel relaxation ('zline', 'yline', or 'point').
!   iter           - INTEGER, Number of iterations.
! Input/Outputs:
!   sol            - REAL array, Updated solution vector (size y_cells*z_cells), input/output.
! Outputs:
!   GS_rhs         - REAL array, Workspace for the RHS (size MAX(y_cells, z_cells)).
!   GS_up_diag     - REAL array, Workspace for modified upper diagonals (size MAX(y_cells, z_cells)).
!-----------------------------------------------------------------------
SUBROUTINE GS_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, rhs, &
                      GS_type, iter, sol, GS_rhs, GS_up_diag)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: A_row_ptr(y_cells*z_cells + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: rhs(y_cells*z_cells)
   CHARACTER(LEN=*), INTENT(IN) :: GS_type
   INTEGER, INTENT(IN) :: iter

   ! Output variables
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(y_cells*z_cells)
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_rhs(MAX(z_cells, y_cells))
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_up_diag(MAX(z_cells, y_cells))

   ! Local variables
   INTEGER :: i

   ! Debugging mode: check CRS consistency
   IF (debugging) THEN
      CALL check_crs_consistency(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind)
   END IF

   ! Gauss-Seidel relaxation based on GS_type
   SELECT CASE (TRIM(GS_type))

   CASE ('zline')
      DO i = 1, iter
         CALL GS_zline_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, rhs, &
                               sol, GS_rhs, GS_up_diag)
      END DO

   CASE ('yline')
      DO i = 1, iter
         CALL GS_yline_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, rhs, &
                               sol, GS_rhs, GS_up_diag)
      END DO

   CASE DEFAULT
      ! Performs point-wise Gauss-Seidel.
      IF (debugging .AND. TRIM(GS_type) /= 'point') THEN
         PRINT *, 'Warning: Unrecognized GS_type "', TRIM(GS_type), '". Defaulting to point-wise Gauss-Seidel.'
      END IF

      ! Perform point-wise Gauss-Seidel relaxation
      DO i = 1, iter
         CALL GS_point_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, rhs, sol)
      END DO

   END SELECT

END SUBROUTINE GS_generic

!-----------------------------------------------------------------------
! subroutine: GS_zline_generic
! Purpose:    Handles Gauss-Seidel relaxation along the z-direction (z-line).
!             The routine solves a tridiagonal system along each z-line for
!             each y-block. It applies line relaxation for the multigrid
!             method used in solving partial differential equations.
! Inputs:
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size y_cells*z_cells + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   rhs            - REAL array, Right-hand side vector (size y_cells*z_cells).
! Input/Outputs:
!   sol            - REAL array, Updated solution vector (size y_cells*z_cells), input/output.
! Outputs:
!   GS_rhs         - REAL array, Workspace for RHS (size MAX(z_cells, y_cells)).
!   GS_up_diag     - REAL array, Workspace for upper diagonal elements and
!                    modification to the diagonal entries (size MAX(z_cells, y_cells)).
!-----------------------------------------------------------------------
SUBROUTINE GS_zline_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, &
                            rhs, sol, GS_rhs, GS_up_diag)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: A_row_ptr(y_cells*z_cells + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: rhs(y_cells*z_cells)

   ! Output variables
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(y_cells*z_cells)
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_rhs(MAX(z_cells, y_cells))
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_up_diag(MAX(z_cells, y_cells))

   ! Local variables
   INTEGER :: j, k, offset, start_ind, cell_index
   REAL(KIND=KIND_generic) :: sum_offdiag, rho

   rho = 0.1  ! Over-relaxation factor

   ! Loop over each y-block (from 1 to y_cells)
   DO j = 1, y_cells
      start_ind = (j - 1)*z_cells  ! Starting index for y-line j

      ! Setup GS_rhs for the current z-line by removing the effects in the y-direction
      IF (j == 1 .OR. j == y_cells) THEN
         IF (j == 1) THEN
            ! First y-block: only upper neighbor (y = 2)
            offset = z_cells
         ELSE
            ! Last y-block: only lower neighbor (y = y_cells - 1)
            offset = -z_cells
         END IF
         DO k = 1, z_cells
            cell_index = start_ind + k  ! (y=j, z=k)

            ! Compute the updating to the diagonal entry
            sum_offdiag = A_val(A_row_ptr(cell_index) + 1) + &
                          A_val(A_row_ptr(cell_index + 1) - 1)
            IF (A_row_ptr(cell_index + 1) - A_row_ptr(cell_index) == 4) THEN
               sum_offdiag = sum_offdiag + A_val(A_row_ptr(cell_index) + 2)
            END IF
            GS_up_diag(k) = rho*A_val(A_row_ptr(cell_index))* &
                            A_val(A_row_ptr(cell_index + 1) - 1)/sum_offdiag

            GS_rhs(k) = rhs(cell_index) - GS_up_diag(k)*sol(cell_index) - &
                        sol(cell_index + offset)*A_val(A_row_ptr(cell_index + 1) - 1)
         END DO
      ELSE
         ! Intermediate y-blocks: both lower (y-1) and upper (y+1) neighbors
         DO k = 1, z_cells
            cell_index = start_ind + k  ! (y=j, z=k)

            ! Compute the updating to the diagonal entry
            sum_offdiag = A_val(A_row_ptr(cell_index) + 1) + &
                          A_val(A_row_ptr(cell_index + 1) - 2) + &
                          A_val(A_row_ptr(cell_index + 1) - 1)
            IF (A_row_ptr(cell_index + 1) - A_row_ptr(cell_index) == 5) THEN
               sum_offdiag = sum_offdiag + A_val(A_row_ptr(cell_index) + 2)
            END IF
            GS_up_diag(k) = rho*A_val(A_row_ptr(cell_index))* &
                            (A_val(A_row_ptr(cell_index + 1) - 1) + &
                             A_val(A_row_ptr(cell_index + 1) - 2))/sum_offdiag

            GS_rhs(k) = rhs(cell_index) - GS_up_diag(k)*sol(cell_index) - &
                        sol(cell_index - z_cells)*A_val(A_row_ptr(cell_index + 1) - 2) - &
                        sol(cell_index + z_cells)*A_val(A_row_ptr(cell_index + 1) - 1)
         END DO
      END IF

      ! Solve the tridiagonal system for the current z-line
      CALL tridiag_solv_zline_generic(y_cells*z_cells, A_row_ptr, A_entries, &
                                      A_col_ind, A_val, y_cells, z_cells, &
                                      j, GS_rhs, sol, GS_up_diag)
   END DO

END SUBROUTINE GS_zline_generic

!-----------------------------------------------------------------------
! subroutine: GS_yline_generic
! Purpose:    Handles Gauss-Seidel relaxation along the y-direction (y-line).
!             The routine solves a tridiagonal system along each y-line for
!             each z-block. It applies line relaxation for the multigrid
!             method used in solving partial differential equations.
! Inputs:
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size y_cells*z_cells + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   rhs            - REAL array, Right-hand side vector (size y_cells*z_cells).
! Input/Outputs:
!   sol            - REAL array, Updated solution vector (size y_cells*z_cells), input/output.
! Outputs:
!   GS_rhs         - REAL array, Workspace for RHS (size MAX(z_cells, y_cells)).
!   GS_up_diag     - REAL array, Workspace for upper diagonal elements and
!                    modification to the diagonal entries (size MAX(z_cells, y_cells)).
!-----------------------------------------------------------------------
SUBROUTINE GS_yline_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, &
                            rhs, sol, GS_rhs, GS_up_diag)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: A_row_ptr(y_cells*z_cells + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: rhs(y_cells*z_cells)

   ! Output variables
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(y_cells*z_cells)
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_rhs(MAX(z_cells, y_cells))
   REAL(KIND=KIND_generic), INTENT(OUT) :: GS_up_diag(MAX(z_cells, y_cells))

   ! Local variables
   INTEGER :: j, k, offset, cell_index
   REAL(KIND=KIND_generic) :: sum_offdiag, rho

   rho = 0.1 ! Over-relaxation factor

   ! Loop over each z-block (from 1 to z_cells)
   DO j = 1, z_cells
      ! Setup GS_rhs for the current y-line by removing the effects in the z-direction
      IF (j == 1 .OR. j == z_cells) THEN
         IF (j == 1) THEN
            ! First z-block: only upper neighbor (z = 2)
            offset = 1
         ELSE
            ! Last z-block: only lower neighbor (z = z_cells - 1)
            offset = -1
         END IF
         DO k = 1, y_cells
            cell_index = (k - 1)*z_cells + j  ! (y=k, z=j)

            ! Compute the updating to the diagonal entry
            sum_offdiag = A_val(A_row_ptr(cell_index) + 1) + &
                          A_val(A_row_ptr(cell_index + 1) - 1)
            IF (A_row_ptr(cell_index + 1) - A_row_ptr(cell_index) == 4) THEN
               sum_offdiag = sum_offdiag + A_val(A_row_ptr(cell_index) + 2)
            END IF
            GS_up_diag(k) = rho*A_val(A_row_ptr(cell_index))* &
                            A_val(A_row_ptr(cell_index) + 1)/sum_offdiag

            GS_rhs(k) = rhs(cell_index) - GS_up_diag(k)*sol(cell_index) - &
                        sol(cell_index + offset)*A_val(A_row_ptr(cell_index) + 1)
         END DO
      ELSE
         ! Intermediate z-blocks: both lower (z-1) and upper (z+1) neighbors
         DO k = 1, y_cells
            cell_index = (k - 1)*z_cells + j  ! (y=k, z=j)

            ! Compute the updating to the diagonal entry
            sum_offdiag = A_val(A_row_ptr(cell_index) + 1) + &
                          A_val(A_row_ptr(cell_index) + 2) + &
                          A_val(A_row_ptr(cell_index + 1) - 1)
            IF (A_row_ptr(cell_index + 1) - A_row_ptr(cell_index) == 5) THEN
               sum_offdiag = sum_offdiag + A_val(A_row_ptr(cell_index + 1) - 2)
            END IF
            GS_up_diag(k) = rho*A_val(A_row_ptr(cell_index))* &
                            (A_val(A_row_ptr(cell_index) + 1) + &
                             A_val(A_row_ptr(cell_index) + 2))/sum_offdiag

            GS_rhs(k) = rhs(cell_index) - GS_up_diag(k)*sol(cell_index) - &
                        sol(cell_index - 1)*A_val(A_row_ptr(cell_index) + 1) - &
                        sol(cell_index + 1)*A_val(A_row_ptr(cell_index) + 2)
         END DO
      END IF

      ! Solve the tridiagonal system for the current y-line
      CALL tridiag_solv_yline_generic(y_cells*z_cells, A_row_ptr, A_entries, &
                                      A_col_ind, A_val, y_cells, z_cells, &
                                      j, GS_rhs, sol, GS_up_diag)
   END DO

END SUBROUTINE GS_yline_generic

!-----------------------------------------------------------------------
! subroutine: GS_point_generic
! Purpose:    Performs point-wise Gauss-Seidel relaxation in given precision.
! Inputs:
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size y_cells*z_cells + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   rhs            - REAL array, Right-hand-side vector (size y_cells*z_cells).
! Input/Outputs:
!   sol            - REAL array, Updated solution vector (size y_cells*z_cells), input/output.
!-----------------------------------------------------------------------
SUBROUTINE GS_point_generic(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, A_val, rhs, sol)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: A_row_ptr(y_cells*z_cells + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: rhs(y_cells*z_cells)

   ! Output variable
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(y_cells*z_cells)

   ! Local variables
   INTEGER :: i, k
   REAL(KIND=KIND_generic) :: temp

   ! Point-wise Gauss-Seidel
   DO i = 1, y_cells*z_cells
      temp = rhs(i)

      ! Subtract contributions from neighboring solutions
      DO k = A_row_ptr(i) + 1, A_row_ptr(i + 1) - 1
         temp = temp - A_val(k)*sol(A_col_ind(k))
      END DO

      ! Check for zero diagonal
      IF (debugging .AND. A_val(A_row_ptr(i)) == 0.0) THEN
         PRINT *, 'Error: Zero diagonal element encountered at row ', i
         STOP
      END IF

      ! Update solution
      sol(i) = temp/A_val(A_row_ptr(i))
   END DO

END SUBROUTINE GS_point_generic

!-----------------------------------------------------------------------
! subroutine: tridiag_solv_zline_generic
! Purpose:    Solves a tridiagonal system along the z-line using the Thomas
!             algorithm for a system stored in CRS format (REAL(KIND=KIND_generic)).
!             The system is of the form A * sol = rhs for each y-block.
! Inputs:
!   A_size         - INTEGER, Total number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_size + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   bloc           - INTEGER, Block index corresponding to the y-block currently being processed.
!   block_rhs      - REAL array, Right-hand-side vector for the current block (size MAX(y_cells, z_cells)).
! Inputs/Outputs:
!   sol            - REAL array, Solution vector (size A_size).
!   GS_up_diag     - REAL array, Workspace for upper diagonal elements and modification
!                    to the diagonal entries (size MAX(z_cells, y_cells)).
!-----------------------------------------------------------------------
SUBROUTINE tridiag_solv_zline_generic(A_size, A_row_ptr, A_entries, A_col_ind, A_val, &
                                      y_cells, z_cells, bloc, block_rhs, sol, block_up_diag)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_size
   INTEGER, INTENT(IN) :: A_row_ptr(A_size + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: bloc  ! Block index corresponding to the current y-block.
   REAL(KIND=KIND_generic), INTENT(IN) :: block_rhs(MAX(y_cells, z_cells))

   ! Input/Output variables
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(A_size)
   REAL(KIND=KIND_generic), INTENT(INOUT) :: block_up_diag(MAX(y_cells, z_cells))

   ! Local variables
   INTEGER :: i
   INTEGER :: row
   REAL(KIND=KIND_generic) :: a_i, b_i, c_i, factor

   IF (debugging) THEN
      CALL check_crs_consistency(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind)
   END IF

   ! Process the first row in 'zline'
   row = (bloc - 1)*z_cells + 1

   ! Extract coefficients for the first row
   b_i = A_val(A_row_ptr(row)) - block_up_diag(1)   ! Diagonal element with splitting
   c_i = A_val(A_row_ptr(row) + 1)                  ! Upper diagonal element

   ! Check for zero diagonal
   IF (debugging .AND. b_i == 0.0) THEN
      PRINT *, 'Error: Zero diagonal element encountered in tridiag_solv_zline at row ', row
      STOP
   END IF

   ! Forward Elimination Phase
   ! Initialize workspace and solution vector
   block_up_diag(1) = c_i/b_i
   sol(row) = block_rhs(1)/b_i

   ! Forward sweep for remaining rows
   DO i = 2, z_cells - 1
      row = row + 1

      ! Extract coefficients
      a_i = A_val(A_row_ptr(row) + 1)                 ! Lower diagonal element
      b_i = A_val(A_row_ptr(row)) - block_up_diag(i)  ! Diagonal element with splitting
      c_i = A_val(A_row_ptr(row) + 2)                 ! Upper diagonal element

      ! Compute elimination factor
      factor = b_i - a_i*block_up_diag(i - 1)

      ! Check for zero diagonal after elimination
      IF (debugging .AND. factor == 0.0) THEN
         PRINT *, 'Error: Zero diagonal element encountered during elimination at row ', row
         STOP
      END IF

      ! Update workspace and solution vector
      block_up_diag(i) = c_i/factor
      sol(row) = (block_rhs(i) - a_i*sol(row - 1))/factor
   END DO

   ! Process the last row in 'zline'
   row = row + 1
   a_i = A_val(A_row_ptr(row) + 1)                        ! Lower diagonal element
   b_i = A_val(A_row_ptr(row)) - block_up_diag(z_cells)   ! Diagonal element with splitting

   factor = b_i - a_i*block_up_diag(z_cells - 1)

   IF (debugging .AND. factor == 0.0) THEN
      PRINT *, 'Error: Zero diagonal element encountered during elimination at row ', row
      STOP
   END IF

   sol(row) = (block_rhs(z_cells) - a_i*sol(row - 1))/factor

   ! Back Substitution Phase
   ! row = bloc*z_cells
   DO i = z_cells - 1, 1, -1
      row = row - 1
      sol(row) = sol(row) - block_up_diag(i)*sol(row + 1)
   END DO

END SUBROUTINE tridiag_solv_zline_generic

!-----------------------------------------------------------------------
! subroutine: tridiag_solv_yline_generic
! Purpose:    Solves a tridiagonal system along the y-line using the Thomas
!             algorithm for a system stored in CRS format (REAL(KIND=KIND_generic)).
!             The system is of the form A * sol = rhs for each z-block.
! Inputs:
!   A_size         - INTEGER, Total number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_size + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   y_cells        - INTEGER, Number of cells in the y-direction.
!   z_cells        - INTEGER, Number of cells in the z-direction.
!   bloc           - INTEGER, Block index corresponding to the z-block currently being processed.
!   block_rhs      - REAL array, Right-hand-side vector for the current block (size MAX(y_cells, z_cells)).
! Inputs/Outputs:
!   sol            - REAL array, Solution vector (size A_size).
!   GS_up_diag     - REAL array, Workspace for upper diagonal elements and modification
!                    to the diagonal entries (size MAX(z_cells, y_cells)).
!-----------------------------------------------------------------------
SUBROUTINE tridiag_solv_yline_generic(A_size, A_row_ptr, A_entries, A_col_ind, A_val, &
                                      y_cells, z_cells, bloc, block_rhs, sol, block_up_diag)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_size
   INTEGER, INTENT(IN) :: A_row_ptr(A_size + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: y_cells, z_cells
   INTEGER, INTENT(IN) :: bloc  ! Block index corresponding to the current z-block.
   REAL(KIND=KIND_generic), INTENT(IN) :: block_rhs(MAX(y_cells, z_cells))

   ! Input/Output variables
   REAL(KIND=KIND_generic), INTENT(INOUT) :: sol(A_size)
   REAL(KIND=KIND_generic), INTENT(INOUT) :: block_up_diag(MAX(y_cells, z_cells))

   ! Local variables
   INTEGER :: i
   INTEGER :: row
   REAL(KIND=KIND_generic) :: a_i, b_i, c_i, factor

   IF (debugging) THEN
      CALL check_crs_consistency(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind)
   END IF

   ! Process the first row in 'yline'
   row = bloc

   ! Extract coefficients for the first row
   b_i = A_val(A_row_ptr(row)) - block_up_diag(1)  ! Diagonal element with splitting
   c_i = A_val(A_row_ptr(row + 1) - 1)             ! Upper diagonal element

   IF (debugging .AND. b_i == 0.0) THEN
      PRINT *, 'Error: Zero diagonal element encountered in tridiag_solv_yline at row ', row
      STOP
   END IF

   ! Forward Elimination Phase
   ! Initialize workspace and solution vector
   block_up_diag(1) = c_i/b_i
   sol(row) = block_rhs(1)/b_i

   ! Forward Elimination Phase
   DO i = 2, y_cells - 1
      row = row + z_cells

      a_i = A_val(A_row_ptr(row + 1) - 2)            ! Lower diagonal element
      b_i = A_val(A_row_ptr(row)) - block_up_diag(i) ! Diagonal element with splitting
      c_i = A_val(A_row_ptr(row + 1) - 1)            ! Upper diagonal element

      factor = b_i - a_i*block_up_diag(i - 1)

      IF (debugging .AND. factor == 0.0) THEN
         PRINT *, 'Error: Zero diagonal element encountered during elimination at row ', row
         STOP
      END IF

      block_up_diag(i) = c_i/factor
      sol(row) = (block_rhs(i) - a_i*sol(row - z_cells))/factor
   END DO

   ! Process the last row in 'yline'
   row = row + z_cells
   a_i = A_val(A_row_ptr(row + 1) - 1)                   ! Lower diagonal element
   b_i = A_val(A_row_ptr(row)) - block_up_diag(y_cells)  ! Diagonal element with splitting

   factor = b_i - a_i*block_up_diag(y_cells - 1)

   IF (debugging .AND. factor == 0.0) THEN
      PRINT *, 'Error: Zero diagonal element encountered during elimination at row ', row
      STOP
   END IF

   sol(row) = (block_rhs(y_cells) - a_i*sol(row - z_cells))/factor

   ! Back Substitution Phase
   ! row = (y_cells - 1)*z_cells + bloc
   DO i = y_cells - 1, 1, -1
      row = row - z_cells
      sol(row) = sol(row) - block_up_diag(i)*sol(row + z_cells)
   END DO

END SUBROUTINE tridiag_solv_yline_generic

!-----------------------------------------------------------------------
! subroutine: resid_generic
! Purpose:    Computes the residual for a linear system using double precision.
!             The residual is the difference between the right-hand-side
!             and the matrix-vector product A * sol.
! Inputs:
!   A_rows         - INTEGER, Number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_rows + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   A_cols         - INTEGER, Number of columns in matrix A (size of solution vector sol).
!   sol            - REAL array, Solution vector (size A_cols).
!   rhs            - REAL array, Right-hand-side vector (size A_rows).
! Outputs:
!   res            - REAL array, Residual vector (size A_rows).
!-----------------------------------------------------------------------
SUBROUTINE resid_generic(A_rows, A_row_ptr, A_entries, A_col_ind, A_val, A_cols, sol, rhs, res)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_rows
   INTEGER, INTENT(IN) :: A_row_ptr(A_rows + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: A_cols
   REAL(KIND=KIND_generic), INTENT(IN) :: sol(A_cols)
   REAL(KIND=KIND_generic), INTENT(IN) :: rhs(A_rows)

   ! Output variable
   REAL(KIND=KIND_generic), INTENT(OUT) :: res(A_rows)

   ! Local variables
   INTEGER :: i

   ! Compute matrix-vector product A * sol and store it in the residual.
   CALL mat_vec_mult_crs_generic(A_rows, A_row_ptr, A_entries, A_col_ind, A_val, A_cols, sol, res)

   ! Calculate the residual: res = rhs - A * sol.
   DO i = 1, A_rows
      res(i) = rhs(i) - res(i)
   END DO

END SUBROUTINE resid_generic

!-----------------------------------------------------------------------
! subroutine: mat_vec_mult_crs_generic
! Purpose:    Performs matrix-vector multiplication in CRS format using REAL(KIND=KIND_generic).
!             This computes the product A * x = b where A is stored in CRS format.
! Inputs:
!   A_rows         - INTEGER, Number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_rows + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   A_cols         - INTEGER, Number of columns in matrix A.
!   x              - REAL array, Input vector (size A_cols).
! Outputs:
!   b              - REAL array, Output vector (size A_rows).
!-----------------------------------------------------------------------
SUBROUTINE mat_vec_mult_crs_generic(A_rows, A_row_ptr, A_entries, A_col_ind, A_val, A_cols, x, b)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_rows
   INTEGER, INTENT(IN) :: A_row_ptr(A_rows + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: A_cols
   REAL(KIND=KIND_generic), INTENT(IN) :: x(A_cols)

   ! Output variable
   REAL(KIND=KIND_generic), INTENT(OUT) :: b(A_rows)

   ! Local variables
   INTEGER :: i, j
   INTEGER :: ind1, ind2
   REAL(KIND=KIND_generic) :: temp

   ind1 = A_row_ptr(1)
   DO i = 1, A_rows
      temp = 0.0
      ind2 = A_row_ptr(i + 1) - 1
      ! Compute the dot product for the current row.
      DO j = ind1, ind2
         temp = temp + A_val(j)*x(A_col_ind(j))
      END DO
      b(i) = temp
      ind1 = ind2 + 1
   END DO

END SUBROUTINE mat_vec_mult_crs_generic

!-----------------------------------------------------------------------
! subroutine: trans_mat_vec_mult_crs_generic
! Purpose:    Performs transposed matrix-vector multiplication in CRS format (REAL(KIND=KIND_generic)).
!             Computes the product b = A^T * x where A is stored in CRS format.
! Inputs:
!   A_rows         - INTEGER, Number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_rows + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   A_cols         - INTEGER, Number of columns in matrix A.
!   x              - REAL array, Input vector (size A_rows).
! Outputs:
!   b              - REAL array, Output vector (size A_cols).
!-----------------------------------------------------------------------
SUBROUTINE trans_mat_vec_mult_crs_generic(A_rows, A_row_ptr, A_entries, A_col_ind, A_val, A_cols, x, b)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_rows
   INTEGER, INTENT(IN) :: A_row_ptr(A_rows + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: A_cols
   REAL(KIND=KIND_generic), INTENT(IN) :: x(A_rows)

   ! Output variable
   REAL(KIND=KIND_generic), INTENT(OUT) :: b(A_cols)

   ! Local variables
   INTEGER :: i, j
   INTEGER :: ind1, ind2
   REAL(KIND=KIND_generic) :: temp

   ! Initialize the result vector to zero.
   b(:) = 0.0

   ind1 = A_row_ptr(1)
   DO i = 1, A_rows
      ind2 = A_row_ptr(i + 1) - 1
      temp = x(i)
      ! Loop over the non-zero entries in the i-th row of A.
      DO j = ind1, ind2
         b(A_col_ind(j)) = b(A_col_ind(j)) + A_val(j)*temp
      END DO
      ind1 = ind2 + 1
   END DO

END SUBROUTINE trans_mat_vec_mult_crs_generic

!-----------------------------------------------------------------------
! subroutine: mat_vec_mult_add_crs_generic
! Purpose:    Performs matrix-vector multiplication followed by vector addition
!             in CRS format (REAL(KIND=KIND_generic)). Computes the operation y = A * x + y.
! Inputs:
!   A_rows         - INTEGER, Number of rows in matrix A.
!   A_row_ptr      - INTEGER array, Row pointer array in CRS format (size A_rows + 1).
!   A_entries      - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind      - INTEGER array, Column index array in CRS format (size A_entries).
!   A_val          - REAL array, Non-zero values in CRS format (size A_entries).
!   A_cols         - INTEGER, Number of columns in matrix A.
!   x              - REAL array, Input vector (size A_cols).
! Outputs:
!   y              - REAL array, Result vector (size A_rows).
!-----------------------------------------------------------------------
SUBROUTINE mat_vec_mult_add_crs_generic(A_rows, A_row_ptr, A_entries, A_col_ind, A_val, A_cols, x, y)
   IMPLICIT NONE

   ! Input variables
   INTEGER, INTENT(IN) :: A_rows
   INTEGER, INTENT(IN) :: A_row_ptr(A_rows + 1)
   INTEGER, INTENT(IN) :: A_entries
   INTEGER, INTENT(IN) :: A_col_ind(A_entries)
   REAL(KIND=KIND_generic), INTENT(IN) :: A_val(A_entries)
   INTEGER, INTENT(IN) :: A_cols
   REAL(KIND=KIND_generic), INTENT(IN) :: x(A_cols)

   ! Output variable
   REAL(KIND=KIND_generic), INTENT(INOUT) :: y(A_rows)

   ! Local variables
   INTEGER :: i, j
   INTEGER :: ind1
   INTEGER :: ind2
   REAL(KIND=KIND_generic) :: temp

   ind1 = A_row_ptr(1)
   DO i = 1, A_rows
      temp = y(i)
      ind2 = A_row_ptr(i + 1) - 1
      ! Compute the dot product for the current row and add it to y.
      DO j = ind1, ind2
         temp = temp + A_val(j)*x(A_col_ind(j))
      END DO
      y(i) = temp
      ind1 = ind2 + 1
   END DO

END SUBROUTINE mat_vec_mult_add_crs_generic

#undef KIND_generic
#undef GS_generic
#undef tridiag_solv_zline_generic
#undef tridiag_solv_yline_generic
#undef resid_generic
#undef mat_vec_mult_crs_generic
#undef trans_mat_vec_mult_crs_generic
#undef mat_vec_mult_add_crs_generic
#undef apply_scaled_prol_generic

#undef GS_zline_generic
#undef GS_yline_generic
#undef GS_point_generic
