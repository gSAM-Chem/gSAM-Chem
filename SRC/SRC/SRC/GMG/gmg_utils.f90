!-----------------------------------------------------------------------------
!
! This module defines some utility routines used by the multigrid solver.
!
!-----------------------------------------------------------------------------
! Author: Oliver Yang and Xiangmin Jiao
! Adapted from code written by Cao Lu.
!
! Initial version: November 3, 2017
! Last updated:
!  - Sept. 27, 2024 (Refactored to encapsulate multigrid and
!    options into a single type for thread safety, updated type names to
!    follow capitalization conventions, and added structure for managing
!    multigrid solvers.).
!  - Oct. 14, 2024 (Changed to use generic programming).
!-----------------------------------------------------------------------------

MODULE gmg_utils
   IMPLICIT NONE

   ! Define kind parameters with descriptive names
   INTEGER, PARAMETER :: KIND_sngl = KIND(1.0E0)  ! Single precision
!   INTEGER, PARAMETER :: KIND_dble = KIND(1.0D0)  ! Double precision
! for gfortran copiler, promotion of real to real(8) makes double precision real(16)
! while intel compiler's -r8  does not. To avoid specify explicit size for double precision
! regardless of the size of real. -MK 2024
   INTEGER, PARAMETER :: KIND_dble = KIND(1._8)  ! Double precision

   !-----------------------------------------------------------------------
   ! Optional Debugging Parameter
   !-----------------------------------------------------------------------
#ifdef DEBUG
   LOGICAL, PARAMETER :: debugging = .TRUE.
#else
   LOGICAL, PARAMETER :: debugging = .FALSE.
#endif

CONTAINS

! Include subroutine templates for Double Precision
#define USE_DOUBLE_PRECISION
#include "gmg_utils.inc"
#undef USE_DOUBLE_PRECISION

! Include subroutine templates for Single Precision
#define USE_SINGLE_PRECISION
#include "gmg_utils.inc"
#undef USE_SINGLE_PRECISION

!-----------------------------------------------------------------------
! Subroutine: check_crs_consistency
! Purpose:    Verifies the consistency of a CRS-formatted Helmholtz discretization matrix.
!             Ensures diagonal entries point to the correct cell, and neighbors are properly indexed.
! Inputs:
!   y_cells   - INTEGER, Number of cells in the y-direction.
!   z_cells   - INTEGER, Number of cells in the z-direction.
!   A_row_ptr - INTEGER array, Row pointer array in CRS format (size y_cells*z_cells + 1).
!   A_entries - INTEGER, Total number of non-zero entries in matrix A.
!   A_col_ind - INTEGER array, Column index array in CRS format (size A_entries).
!-----------------------------------------------------------------------
   SUBROUTINE check_crs_consistency(y_cells, z_cells, A_row_ptr, A_entries, A_col_ind)
      IMPLICIT NONE

      ! Input variables
      INTEGER, INTENT(IN) :: y_cells, z_cells
      INTEGER, INTENT(IN) :: A_row_ptr(y_cells*z_cells + 1)
      INTEGER, INTENT(IN) :: A_entries
      INTEGER, INTENT(IN) :: A_col_ind(A_entries)

      ! Local variables
      INTEGER :: j, k, cell_index, start_ind, offset

      ! Validate grid sizes
      IF (y_cells < 1 .OR. z_cells < 1) THEN
         PRINT *, 'Error: Invalid grid sizes. y_cells and z_cells must be >= 1.'
         STOP
      END IF

      ! Check z-line consistency for boundary and interior cells
      DO j = 1, y_cells
         start_ind = (j - 1)*z_cells

         IF (j == 1 .OR. j == y_cells) THEN
            IF (j == 1) THEN
               ! First y-block: only upper neighbor (y = 2)
               offset = z_cells
            ELSE
               ! Last y-block: only lower neighbor (y = y_cells - 1)
               offset = -z_cells
            END IF

            ! Loop through each z-cell in the y-line
            DO k = 1, z_cells
               cell_index = start_ind + k

               ! Check diagonal entry points to itself
               IF (A_col_ind(A_row_ptr(cell_index)) /= cell_index) THEN
                  PRINT *, 'Error: Diagonal mismatch at row ', cell_index
                  STOP
               END IF

               ! Check neighbor in the z direction
               IF (A_col_ind(A_row_ptr(cell_index + 1) - 1) /= cell_index + offset) THEN
                  IF (j == 1) THEN
                     PRINT *, 'Error: Upper diagonal index mismatch at row ', cell_index
                  ELSE
                     PRINT *, 'Error: Lower diagonal index mismatch at row ', cell_index
                  END IF
                  STOP
               END IF
            END DO
         ELSE
            ! Interior cells: check both neighbors
            DO k = 1, z_cells
               cell_index = start_ind + k

               ! Check diagonal
               IF (A_col_ind(A_row_ptr(cell_index)) /= cell_index) THEN
                  PRINT *, 'Error: Diagonal mismatch at row ', cell_index
                  STOP
               END IF

               ! Check lower and upper neighbors
               IF (A_col_ind(A_row_ptr(cell_index + 1) - 2) /= cell_index - z_cells) THEN
                  PRINT *, 'Error: Lower neighbor mismatch at row ', cell_index
                  STOP
               END IF

               IF (A_col_ind(A_row_ptr(cell_index + 1) - 1) /= cell_index + z_cells) THEN
                  PRINT *, 'Error: Upper neighbor mismatch at row ', cell_index
                  STOP
               END IF
            END DO
         END IF
      END DO

      ! Check y-line consistency for boundary and interior cells
      DO j = 1, z_cells
         IF (j == 1 .OR. j == z_cells) THEN
            IF (j == 1) THEN
               ! First z-block: only upper neighbor (z = 2)
               offset = 1
            ELSE
               ! Last z-block: only lower neighbor (z = z_cells - 1)
               offset = -1
            END IF

            ! Loop through each cell in the current z-block
            DO k = 1, y_cells
               cell_index = (k - 1)*z_cells + j  ! Cell index (y=k, z=j)

               ! Check if neighbor index is correct based on boundary condition
               IF (A_col_ind(A_row_ptr(cell_index) + 1) /= cell_index + offset) THEN
                  IF (j == 1) THEN
                     PRINT *, 'Error: Upper diagonal index mismatch at row ', cell_index
                  ELSE
                     PRINT *, 'Error: Lower diagonal index mismatch at row ', cell_index
                  END IF
                  STOP
               END IF
            END DO

         ELSE
            ! Intermediate z-blocks: both lower (z-1) and upper (z+1) neighbors
            DO k = 1, y_cells
               cell_index = (k - 1)*z_cells + j  ! Cell index (y=k, z=j)

               ! Check lower neighbor index (z-1)
               IF (A_col_ind(A_row_ptr(cell_index) + 1) /= cell_index - 1) THEN
                  PRINT *, 'Error: Lower diagonal index mismatch at row ', cell_index
                  STOP
               END IF

               ! Check upper neighbor index (z+1)
               IF (A_col_ind(A_row_ptr(cell_index) + 2) /= cell_index + 1) THEN
                  PRINT *, 'Error: Upper diagonal index mismatch at row ', cell_index
                  STOP
               END IF
            END DO
         END IF
      END DO

   END SUBROUTINE check_crs_consistency

END MODULE gmg_utils
