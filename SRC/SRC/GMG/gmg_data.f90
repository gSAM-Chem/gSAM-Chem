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
!
! This module defines the data structures and options used by the multigrid
! solver. It includes types for grid information, solver options, and the
! multigrid solver itself.
!
!-----------------------------------------------------------------------------

MODULE gmg_data
   IMPLICIT NONE

   !-----------------------------------------------------------------------------
   ! TYPE: Grid
   !
   ! Represents a single grid level in the multigrid hierarchy, containing
   ! information about the grid cells, the discretized PDE matrix in CRS format,
   ! solution vectors, right-hand side vectors, smoother buffers, residuals,
   ! and prolongation matrices.
   !-----------------------------------------------------------------------------
   TYPE :: Grid
      ! Number of cells in y and z directions, and total number of cells
      INTEGER, DIMENSION(3) :: cells = 0

      ! CRS format matrix for the discretized PDE
      INTEGER, ALLOCATABLE :: A_row_ptr(:), A_col_ind(:)
      REAL(8), ALLOCATABLE :: A_val_dble(:), A_diag_val_dble(:)
      REAL, ALLOCATABLE :: A_val_sngl(:), A_diag_val_sngl(:)

      ! Solution arrays
      REAL(8), ALLOCATABLE :: sol_dble(:)
      REAL, ALLOCATABLE :: sol_sngl(:)

      ! Right-hand side
      REAL(8), ALLOCATABLE :: rhs_dble(:)
      REAL, ALLOCATABLE :: rhs_sngl(:)

      ! Smoother buffer
      REAL(8), ALLOCATABLE :: GS_rhs_dble(:)
      REAL, ALLOCATABLE :: GS_rhs_sngl(:)

      ! Tridiagonal solver buffer
      REAL(8), ALLOCATABLE :: GS_up_diag_dble(:)
      REAL, ALLOCATABLE :: GS_up_diag_sngl(:)

      ! Residual
      REAL(8), ALLOCATABLE :: res_dble(:)
      REAL, ALLOCATABLE :: res_sngl(:)

      ! Prolongation operator
      INTEGER, ALLOCATABLE :: P_row_ptr(:), P_col_ind(:)
      REAL(8), ALLOCATABLE :: P_val_dble(:)
      REAL, ALLOCATABLE :: P_val_sngl(:)

      ! Restriction operator
      INTEGER, ALLOCATABLE :: R_row_ptr(:), R_col_ind(:)
      REAL(8), ALLOCATABLE :: R_val_dble(:)
      REAL, ALLOCATABLE :: R_val_sngl(:)
   END TYPE Grid

   !-----------------------------------------------------------------------------
   ! TYPE: Options
   !
   ! Contains all the configurable parameters for the multigrid solver, such
   ! as the number of smoothing iterations, types of smoothing, precision
   ! settings, coarse grid size, and semi-coarsening options.
   !-----------------------------------------------------------------------------
   TYPE :: Options
      ! The radius of the earth in meters. Default is 6.37122e6.
      REAL(8) :: rad_earth = 6.371229D+6

      ! Number of pre-smoothing iterations. Default is 2.
      INTEGER :: pre_smooth_iter = 2

      ! Number of post-smoothing iterations. Default is 2.
      INTEGER :: post_smooth_iter = 2

      ! Type of Gauss-Seidel smoothing. Default is 'adapt'.
      CHARACTER(len=5) :: GS_type = 'adapt'

      ! Null space removal method. Default is 'cycl'.
      CHARACTER(len=4) :: remov_null = 'cycl'

      ! Whether to use mixed single-double precision when alpha/=0. Default is .TRUE.
      LOGICAL :: mixed_precision = .TRUE.

      ! Maximum size of the coarse grid. Default is 8.
      INTEGER :: coarse_size = 8

      ! Maximum number of levels. Default is 20.
      INTEGER :: max_levels = 20

      ! Semi-coarsening option. Default is .TRUE.
      LOGICAL :: semi_coarsen = .TRUE.

      ! Number of iterations for coarse-grid solver if LU/QR is not used. Default is 2.
      INTEGER :: coarse_solve_iter = 2

      ! Enables the use of full multigrid (FMG) cycles. Default is .TRUE., which
      ! enables FMG for the first cycle if alpha /= 0 and every cycle if alpha = 0.
      LOGICAL :: full_mg = .TRUE.

#if defined(DEBUG)
      ! Whether to enable debugging output. Default is .TRUE. if DEBUG is defined.
      LOGICAL :: debug = .TRUE.
#else
      LOGICAL :: debug = .FALSE.
#endif
   END TYPE Options

   !-----------------------------------------------------------------------------
   ! TYPE: MultigridSolver
   !
   ! Encapsulates both the multigrid solver's grid hierarchy and its options.
   ! It includes residual storage for termination conditions.
   !-----------------------------------------------------------------------------
   TYPE :: MultigridSolver
      TYPE(Grid), ALLOCATABLE :: multigrid(:)   ! Array of grid structures
      TYPE(Options) :: options                  ! Multigrid solver options

      ! Declares a copy of the residual on the fine-grid from the previous
      ! v_cycle in order to compute the first termination condition.
      REAL(8), ALLOCATABLE :: res_old(:)

      ! Declares the transpose of the discretization matrix on the fine-grid
      ! applied to the right-hand-side and the residual in order to compute
      ! the second termination condition.
      REAL(8), ALLOCATABLE :: rhs_A_trans(:)
      REAL(8), ALLOCATABLE :: res_A_trans(:)

   END TYPE MultigridSolver

CONTAINS
   !-----------------------------------------------------------------------
   ! SUBROUTINE: dump_linsys
   ! PURPOSE:
   !   Dumps the sparse matrix (A) and RHS vector from the finest level
   !   in full double-precision, suitable for loading into MATLAB.
   !
   ! INPUTS:
   !   grd      - TYPE(Grid), the finest grid
   !   alppha   - REAL(8), the reaction coefficient
   !   rhs      - REAL(8) array, the right-hand side vector
   !   filename - CHARACTER(len=*), OPTIONAL
   !              Base name for the output files. Defaults to 'linsys'.
   !
   ! OUTPUT:
   !   Two files containing the matrix and RHS vector:
   !     - <filename>_matrix.txt
   !     - <filename>_rhs.txt
   !-----------------------------------------------------------------------
   SUBROUTINE dump_linsys(grd, alpha, rhs, filename)
      IMPLICIT NONE

      TYPE(Grid), INTENT(IN) :: grd
      REAL(8), INTENT(IN) :: alpha
      REAL(8), INTENT(IN) :: rhs(grd%cells(3))
      CHARACTER(len=*), OPTIONAL, INTENT(IN) :: filename

      INTEGER :: i, j, k
      INTEGER :: row_start, row_end
      INTEGER :: unit_matrix, unit_rhs, ios
      CHARACTER(LEN=100) :: matrix_filename, rhs_filename
      INTEGER :: nrows

      ! Determine filenames
      IF (PRESENT(filename)) THEN
         matrix_filename = TRIM(filename)//'_matrix.txt'
         rhs_filename = TRIM(filename)//'_rhs.txt'
      ELSE
         matrix_filename = 'linsys_matrix.txt'
         rhs_filename = 'linsys_rhs.txt'
      END IF

      ! Check if A_row_ptr and A_col_ind are allocated
      IF (.NOT. ALLOCATED(grd%A_row_ptr) .OR. .NOT. ALLOCATED(grd%A_col_ind)) THEN
         PRINT *, 'Error: Matrix A_row_ptr or A_col_ind not allocated.'
         RETURN
      END IF

      ! Check if A_val_dble or A_val_sngl is allocated
      IF (.NOT. ALLOCATED(grd%A_val_dble) .AND. .NOT. ALLOCATED(grd%A_val_sngl)) THEN
         PRINT *, 'Error: Matrix values not allocated.'
         RETURN
      END IF

      ! Check if rhs_dble or rhs_sngl is allocated
      IF (.NOT. ALLOCATED(grd%rhs_dble) .AND. .NOT. ALLOCATED(grd%rhs_sngl)) THEN
         PRINT *, 'Error: RHS vector not allocated.'
         RETURN
      END IF

      ! Open files for writing
      OPEN (NEWUNIT=unit_matrix, FILE=matrix_filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
         PRINT *, 'Error opening file ', matrix_filename
         RETURN
      END IF

      OPEN (NEWUNIT=unit_rhs, FILE=rhs_filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
         PRINT *, 'Error opening file ', rhs_filename
         CLOSE (unit_matrix)
         RETURN
      END IF

      ! Number of rows
      nrows = grd%cells(3)

      ! Loop over rows to write matrix data
      DO i = 1, nrows
         row_start = grd%A_row_ptr(i)
         row_end = grd%A_row_ptr(i + 1) - 1
         DO k = row_start, row_end
            j = grd%A_col_ind(k)
            ! Get value
            IF (i == j) THEN
               ! Diagonal element
               WRITE (unit_matrix, '(I10,1X,I10,1X,E25.16)') i, j, grd%A_diag_val_dble(i) - alpha
            ELSE
               ! Off-diagonal element
               WRITE (unit_matrix, '(I10,1X,I10,1X,E25.16)') i, j, grd%A_val_dble(k)
            END IF
         END DO
      END DO

      CLOSE (unit_matrix)

      ! Write RHS vector
      DO i = 1, nrows
         WRITE (unit_rhs, '(I10,1X,E25.16)') i, rhs(i)
      END DO

      CLOSE (unit_rhs)

      PRINT *, 'Linear system dumped to ', TRIM(matrix_filename), ' and ', TRIM(rhs_filename)

   END SUBROUTINE dump_linsys

   !-----------------------------------------------------------------------
   ! SUBROUTINE: dump_griddata
   ! PURPOSE:
   !   Writes the multigrid data structure fields to an ASCII file for debugging.
   !
   ! INPUTS:
   !   mg       - TYPE(Grid), array of grid structures to be dumped
   !   filename - CHARACTER(len=*), OPTIONAL
   !              Name of the output ASCII file. If not provided, defaults to 'dump_griddata.txt'.
   !
   ! OUTPUT:
   !   ASCII file containing the specified fields for each grid
   !
   !-----------------------------------------------------------------------
   SUBROUTINE dump_griddata(mg, filename)
      IMPLICIT NONE

      ! Interface declarations
      TYPE(Grid), INTENT(IN) :: mg(:)
      CHARACTER(len=*), OPTIONAL, INTENT(IN) :: filename

      INTEGER :: lvl
      CHARACTER(LEN=100) :: out_filename
      INTEGER :: dump_unit, ios
      INTEGER :: i
      INTEGER :: n_A_row_ptr, n_A_col_ind, n_A_val
      INTEGER :: n_P_row_ptr, n_P_col_ind, n_P_val

      ! Determine the output filename
      IF (PRESENT(filename)) THEN
         out_filename = TRIM(filename)
      ELSE
         out_filename = 'dump_griddata.txt'
      END IF

      ! Find an available unit number
      dump_unit = find_free_unit()

      IF (dump_unit == -1) THEN
         PRINT *, 'Error: No available unit number found for dumping multigrid data.'
         RETURN
      END IF

      ! Open the file for writing
      OPEN (unit=dump_unit, file=out_filename, status='replace', action='write', iostat=ios)

      ! Check for successful file opening
      IF (ios /= 0) THEN
         PRINT *, 'Error: Unable to open file ', out_filename
         RETURN
      END IF

      ! Loop over each grid level and write the specified fields
      DO lvl = 1, SIZE(mg)
         WRITE (dump_unit, '(A,I0)') 'Grid Level: ', lvl
         WRITE (dump_unit, '(A,I0,A,I0,A,I0)') 'Sizes: ', mg(lvl)%cells(1), ',', mg(lvl)%cells(2), ',', mg(lvl)%cells(3)

         ! Write A_row_ptr
         WRITE (dump_unit, '(A)') 'A_row_ptr:'
         n_A_row_ptr = SIZE(mg(lvl)%A_row_ptr)
         IF (n_A_row_ptr > 0) THEN
            DO i = 1, n_A_row_ptr
               WRITE (dump_unit, '(I0, A)', ADVANCE='NO') mg(lvl)%A_row_ptr(i), ' '
               IF (MOD(i, 10) == 0 .OR. i == n_A_row_ptr) THEN
                  WRITE (dump_unit, *)  ! New line after every 10 elements or at the end
               END IF
            END DO
         ELSE
            WRITE (dump_unit, '(A)') 'None allocated.'
         END IF

         ! Write A_col_ind
         WRITE (dump_unit, '(A)') 'A_col_ind:'
         n_A_col_ind = SIZE(mg(lvl)%A_col_ind)
         IF (n_A_col_ind > 0) THEN
            DO i = 1, n_A_col_ind
               WRITE (dump_unit, '(I0, A)', ADVANCE='NO') mg(lvl)%A_col_ind(i), ' '
               IF (MOD(i, 10) == 0 .OR. i == n_A_col_ind) THEN
                  WRITE (dump_unit, *)  ! New line after every 10 elements or at the end
               END IF
            END DO
         ELSE
            WRITE (dump_unit, '(A)') 'None allocated.'
         END IF

         IF (ALLOCATED(mg(lvl)%A_val_sngl)) THEN
            ! Write A_val_sngl
            WRITE (dump_unit, '(A)') 'A_val_sngl:'
            n_A_val = SIZE(mg(lvl)%A_val_sngl)
            IF (n_A_val > 0) THEN
               DO i = 1, n_A_val
                  WRITE (dump_unit, '(E15.8, A)', ADVANCE='NO') mg(lvl)%A_val_sngl(i), ' '
                  IF (MOD(i, 5) == 0 .OR. i == n_A_val) THEN
                     WRITE (dump_unit, *)  ! New line after every 5 elements or at the end
                  END IF
               END DO
            ELSE
               WRITE (dump_unit, '(A)') 'None allocated.'
            END IF
         ELSE
            ! Write A_val_dble
            WRITE (dump_unit, '(A)') 'A_val_dble:'
            n_A_val = SIZE(mg(lvl)%A_val_dble)
            IF (n_A_val > 0) THEN
               DO i = 1, n_A_val
                  WRITE (dump_unit, '(E15.8, A)', ADVANCE='NO') mg(lvl)%A_val_dble(i), ' '
                  IF (MOD(i, 5) == 0 .OR. i == n_A_val) THEN
                     WRITE (dump_unit, *)  ! New line after every 5 elements or at the end
                  END IF
               END DO
            ELSE
               WRITE (dump_unit, '(A)') 'None allocated.'
            END IF
         END IF

         ! Check if prolongation operator fields are allocated
         IF (ALLOCATED(mg(lvl)%P_row_ptr) .AND. ALLOCATED(mg(lvl)%P_col_ind) .AND. &
             (ALLOCATED(mg(lvl)%P_val_dble) .OR. ALLOCATED(mg(lvl)%P_val_sngl))) THEN

            ! Write P_row_ptr
            WRITE (dump_unit, '(A)') 'P_row_ptr:'
            n_P_row_ptr = SIZE(mg(lvl)%P_row_ptr)
            IF (n_P_row_ptr > 0) THEN
               DO i = 1, n_P_row_ptr
                  WRITE (dump_unit, '(I0, A)', ADVANCE='NO') mg(lvl)%P_row_ptr(i), ' '
                  IF (MOD(i, 10) == 0 .OR. i == n_P_row_ptr) THEN
                     WRITE (dump_unit, *)  ! New line after every 10 elements or at the end
                  END IF
               END DO
            ELSE
               WRITE (dump_unit, '(A)') 'None allocated.'
            END IF

            ! Write P_col_ind
            WRITE (dump_unit, '(A)') 'P_col_ind:'
            n_P_col_ind = SIZE(mg(lvl)%P_col_ind)
            IF (n_P_col_ind > 0) THEN
               DO i = 1, n_P_col_ind
                  WRITE (dump_unit, '(I0, A)', ADVANCE='NO') mg(lvl)%P_col_ind(i), ' '
                  IF (MOD(i, 10) == 0 .OR. i == n_P_col_ind) THEN
                     WRITE (dump_unit, *)  ! New line after every 10 elements or at the end
                  END IF
               END DO
            ELSE
               WRITE (dump_unit, '(A)') 'None allocated.'
            END IF

            IF (ALLOCATED(mg(lvl)%P_val_sngl)) THEN
               ! Write P_val_sngl
               WRITE (dump_unit, '(A)') 'P_val_sngl:'
               n_P_val = SIZE(mg(lvl)%P_val_sngl)
               IF (n_P_val > 0) THEN
                  DO i = 1, n_P_val
                     WRITE (dump_unit, '(E15.8, A)', ADVANCE='NO') mg(lvl)%P_val_sngl(i), ' '
                     IF (MOD(i, 5) == 0 .OR. i == n_P_val) THEN
                        WRITE (dump_unit, *)  ! New line after every 5 elements or at the end
                     END IF
                  END DO
               ELSE
                  WRITE (dump_unit, '(A)') 'None allocated.'
               END IF
            ELSE
               ! Write P_val_dble
               WRITE (dump_unit, '(A)') 'P_val_dble:'
               n_P_val = SIZE(mg(lvl)%P_val_dble)
               IF (n_P_val > 0) THEN
                  DO i = 1, n_P_val
                     WRITE (dump_unit, '(E15.8, A)', ADVANCE='NO') mg(lvl)%P_val_dble(i), ' '
                     IF (MOD(i, 5) == 0 .OR. i == n_P_val) THEN
                        WRITE (dump_unit, *)  ! New line after every 5 elements or at the end
                     END IF
                  END DO
               ELSE
                  WRITE (dump_unit, '(A)') 'None allocated.'
               END IF
            END IF
         ELSE
            WRITE (dump_unit, '(A)') 'P_row_ptr, P_col_ind, and P_val: Not allocated for this grid level.'
         END IF

         ! Separator for readability
         WRITE (dump_unit, '(A)') '--------------------------'
      END DO

      ! Close the file
      CLOSE (dump_unit)

      PRINT *, 'Multigrid data has been successfully dumped to ', out_filename

   END SUBROUTINE dump_griddata

   !-----------------------------------------------------------------------
   ! FUNCTION: find_free_unit
   ! PURPOSE:
   !   Finds an available unit number for file operations.
   !
   ! INPUTS:
   !   None
   !
   ! OUTPUT:
   !   free_unit - INTEGER, available unit number
   !
   ! RETURNS:
   !   free_unit = -1 if no available unit number is found
   !-----------------------------------------------------------------------
   INTEGER FUNCTION find_free_unit()
      IMPLICIT NONE
      INTEGER :: unit_try
      LOGICAL :: is_available

      ! Define the range of unit numbers to search
      DO unit_try = 10, 99
         ! Check if the unit is currently open
         INQUIRE (UNIT=unit_try, OPENED=is_available)
         IF (.NOT. is_available) THEN
            find_free_unit = unit_try
            RETURN
         END IF
      END DO

      ! If no free unit is found, return -1
      find_free_unit = -1
   END FUNCTION find_free_unit

END MODULE gmg_data
