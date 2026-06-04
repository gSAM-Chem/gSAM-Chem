!-----------------------------------------------------------------------------
! MODULE: gmg
!
! This module implements the public interface of the pressure GMG (Geometric
! Multigrid) solver. It provides routines for initializing the solver,
! performing the multigrid solve, and finalizing the solver.
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

MODULE gmg
   USE gmg_data, ONLY: MultigridSolver
   IMPLICIT NONE

   TYPE(MultigridSolver) :: solver
CONTAINS

   !-----------------------------------------------------------------------
   ! SUBROUTINE: gmg_init
   ! PURPOSE:
   !   Sets up the multigrid data structure for the specified grid parameters
   !   and equation coefficients. Initializes the solver with default options
   !   if they are uninitialized.
   !
   ! INPUTS:
   !   solver                - TYPE(MultigridSolver), input/output
   !   y_fine_cells          - INTEGER, number of fine cells in y-direction
   !   z_fine_cells          - INTEGER, number of fine cells in z-direction
   !   y_cell_edge_coords    - DOUBLE PRECISION array, y-edge coordinates (size y_fine_cells + 1)
   !   z_cell_edge_coords    - DOUBLE PRECISION array, z-edge coordinates (size z_fine_cells + 1)
   !   p_cell_edge           - DOUBLE PRECISION array, pressure at cell edges (size z_fine_cells + 1)
   !   mu_discr              - CHARACTER(len=1), discretization type for mu
   !   rho_discr             - CHARACTER(len=1), discretization type for rho
   !   mu_is_const           - LOGICAL, flag indicating if mu is constant
   !   rad_reduc_factor      - DOUBLE PRECISION, radial reduction factor
   !
   ! OUTPUTS:
   !   solver                - TYPE(MultigridSolver), initialized solver
   !-----------------------------------------------------------------------
   SUBROUTINE gmg_init(y_fine_cells, z_fine_cells, &
                       y_cell_edge_coords, z_cell_edge_coords, p_cell_edge, mu_discr, &
                       rho_discr, mu_is_const, rad_reduc_factor)
      USE gmg_kernel, ONLY: compute_levels, coarsen
      USE gmg_discr, ONLY: discretize
      IMPLICIT NONE

      ! Declare the input variables
      INTEGER, INTENT(IN) :: y_fine_cells, z_fine_cells
      DOUBLE PRECISION, INTENT(INOUT) :: y_cell_edge_coords(y_fine_cells + 1), &
         z_cell_edge_coords(z_fine_cells + 1)
      DOUBLE PRECISION, INTENT(INOUT) :: p_cell_edge(z_fine_cells + 1)
      CHARACTER(len=1), INTENT(IN) :: mu_discr, rho_discr
      LOGICAL, INTENT(IN) :: mu_is_const
      DOUBLE PRECISION, INTENT(IN) :: rad_reduc_factor

      ! Declare local variables
      INTEGER :: i, lvls, lvl, y_padside, z_padside
      DOUBLE PRECISION, ALLOCATABLE :: y_cell_widths(:), z_cell_widths(:)
      DOUBLE PRECISION, ALLOCATABLE :: y_cell_widths_work(:), z_cell_widths_work(:)
      DOUBLE PRECISION, ALLOCATABLE :: y_finer_cell_cent_coords(:), &
         z_finer_cell_cent_coords(:), y_coarser_cell_cent_coords(:), z_coarser_cell_cent_coords(:)
      DOUBLE PRECISION :: min_diag_val, max_diag_val

      ! Error checking
      IF (y_fine_cells <= 0 .OR. z_fine_cells <= 0) THEN
         PRINT *, 'Error in gmg_init: Number of fine cells must be positive.'
         STOP
      END IF
      IF (SIZE(y_cell_edge_coords) /= y_fine_cells + 1) THEN
         PRINT *, 'Error in gmg_init: Size mismatch in y_cell_edge_coords.'
         STOP
      END IF
      IF (SIZE(z_cell_edge_coords) /= z_fine_cells + 1) THEN
         PRINT *, 'Error in gmg_init: Size mismatch in z_cell_edge_coords.'
         STOP
      END IF
      IF (SIZE(p_cell_edge) /= z_fine_cells + 1) THEN
         PRINT *, 'Error in gmg_init: Size mismatch in p_cell_edge.'
         STOP
      END IF

      y_padside = 0
      z_padside = 0
      ! Compute the number of levels based on grid size
      lvls = MIN(solver%options%max_levels, &
                 compute_levels(y_fine_cells, z_fine_cells, solver%options%coarse_size, &
                                solver%options%semi_coarsen))

      ! Allocate the array of grid variables
      ALLOCATE (solver%multigrid(lvls))

      ASSOCIATE (mg => solver%multigrid, opts => solver%options)
         ! Set GS_type based on the number of cells in each direction
         IF (TRIM(opts%GS_type) == '') THEN
            opts%GS_type = 'adapt'
         END IF

         ! Set the number of cells on the fine-grid and compute them
         mg(1)%cells(1) = y_fine_cells
         mg(1)%cells(2) = z_fine_cells
         mg(1)%cells(3) = mg(1)%cells(1)*mg(1)%cells(2)

         ! Allocate arrays for the widths of the cells in the y and z-directions
         ALLOCATE (y_cell_widths(mg(1)%cells(1)), z_cell_widths(mg(1)%cells(2)))
         ALLOCATE (y_cell_widths_work(mg(1)%cells(1)), z_cell_widths_work(mg(1)%cells(2)))

         ! Allocate arrays for the coordinates of the centers of the finer and coarser cells
         ALLOCATE (y_finer_cell_cent_coords(mg(1)%cells(1)), &
                   z_finer_cell_cent_coords(mg(1)%cells(2)), &
                   y_coarser_cell_cent_coords(mg(1)%cells(1)), &
                   z_coarser_cell_cent_coords(mg(1)%cells(2)))

         ! Compute the widths of the fine-grid cells
         DO i = 1, mg(1)%cells(1)
            y_cell_widths(i) = y_cell_edge_coords(i + 1) - y_cell_edge_coords(i)
            IF (y_cell_widths(i) <= 0.0D0) THEN
               PRINT *, 'Error in gmg_init: Non-positive y_cell_width at index ', i
               STOP
            END IF
         END DO
         DO i = 1, mg(1)%cells(2)
            z_cell_widths(i) = z_cell_edge_coords(i + 1) - z_cell_edge_coords(i)
            IF (z_cell_widths(i) <= 0.0D0) THEN
               PRINT *, 'Error in gmg_init: Non-positive z_cell_width at index ', i
               STOP
            END IF
         END DO

         ! Compute the coordinates of the centers of the fine-grid cells
         DO i = 1, mg(1)%cells(1)
            y_finer_cell_cent_coords(i) = 0.5D0*(y_cell_edge_coords(i) + y_cell_edge_coords(i + 1))
         END DO
         DO i = 1, mg(1)%cells(2)
            z_finer_cell_cent_coords(i) = 0.5D0*(z_cell_edge_coords(i) + z_cell_edge_coords(i + 1))
         END DO

         ! Loop through each grid level to set up the multigrid structure
         DO lvl = 1, lvls
            ! Call the discretization routine to set up the PDE matrix
            CALL discretize(mg(lvl), opts, y_cell_edge_coords, p_cell_edge, y_cell_widths, z_cell_widths, &
                            y_finer_cell_cent_coords, mu_discr, rho_discr, mu_is_const, rad_reduc_factor)
            IF (opts%debug .AND. lvl == 1) THEN
               ! Initialize max and min diagonal values.
               min_diag_val = 1.0D+200
               max_diag_val = -1.0D+200

               ! Find the min and max diagonal values.
               DO i = 1, mg(1)%cells(3)
                  IF (mg(lvl)%A_diag_val_dble(i) < min_diag_val) THEN
                     min_diag_val = mg(lvl)%A_diag_val_dble(i)
                  END IF
                  IF (mg(lvl)%A_diag_val_dble(i) > max_diag_val) THEN
                     max_diag_val = mg(lvl)%A_diag_val_dble(i)
                  END IF
               END DO

               ! Print min and max diagonal values if opts%debugging is TRUE.
               PRINT *, 'Min and max diagonal values on level', lvl, ':', min_diag_val, max_diag_val
            END IF

            ! Convert double precision values to single precision if needed.
            ALLOCATE (mg(lvl)%A_val_sngl(mg(lvl)%A_row_ptr(mg(lvl)%cells(3) + 1) - 1), &
                      mg(lvl)%A_diag_val_sngl(mg(lvl)%cells(3)))
            DO i = 1, mg(lvl)%A_row_ptr(mg(lvl)%cells(3) + 1) - 1
               mg(lvl)%A_val_sngl(i) = REAL(mg(lvl)%A_val_dble(i))
            END DO
            DO i = 1, mg(lvl)%cells(3)
               mg(lvl)%A_diag_val_sngl(i) = REAL(mg(lvl)%A_diag_val_dble(i))
            END DO

            ! Allocate solution and right-hand-side vectors for each grid level
            ALLOCATE (mg(lvl)%sol_dble(mg(lvl)%cells(3)), mg(lvl)%rhs_dble(mg(lvl)%cells(3)))
            ALLOCATE (mg(lvl)%sol_sngl(mg(lvl)%cells(3)), mg(lvl)%rhs_sngl(mg(lvl)%cells(3)))

            ! Allocate buffers for smoother and tridiagonal solver
            ALLOCATE (mg(lvl)%GS_rhs_dble(MAX(mg(lvl)%cells(1), mg(lvl)%cells(2))), &
                      mg(lvl)%GS_up_diag_dble(MAX(mg(lvl)%cells(1), mg(lvl)%cells(2))), &
                      mg(lvl)%res_dble(mg(lvl)%cells(3)))
            ! Allocate buffers for smoother and tridiagonal solver
            ALLOCATE (mg(lvl)%GS_rhs_sngl(MAX(mg(lvl)%cells(1), mg(lvl)%cells(2))), &
                      mg(lvl)%GS_up_diag_sngl(MAX(mg(lvl)%cells(1), mg(lvl)%cells(2))), &
                      mg(lvl)%res_sngl(mg(lvl)%cells(3)))

            ! For every grid except the coarse-grid:
            IF (lvl < lvls) THEN
               ! Call the coarsen routine to allocate necessary buffers and compute values
               CALL coarsen(solver, lvl, y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
                            y_padside, z_padside, y_cell_edge_coords, z_cell_edge_coords, &
                            y_cell_widths, z_cell_widths, p_cell_edge, &
                            y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, &
                            y_cell_widths_work, z_cell_widths_work)

               IF (lvl + 1 < lvls) THEN
                  ! Copy grid center coordinates from coarser grid to finer grid.
                  y_finer_cell_cent_coords(1:mg(lvl)%cells(1)) = y_coarser_cell_cent_coords(1:mg(lvl)%cells(1))
                  z_finer_cell_cent_coords(1:mg(lvl)%cells(2)) = z_coarser_cell_cent_coords(1:mg(lvl)%cells(2))
               END IF
            END IF
         END DO

         ! Allocate additional arrays for residuals and transposes
         ALLOCATE (solver%res_old(mg(1)%cells(3)), solver%rhs_A_trans(mg(1)%cells(3)), &
                   solver%res_A_trans(mg(1)%cells(3)))

         ! Deallocate temporary arrays for cell widths and coordinates
         DEALLOCATE (y_cell_widths, z_cell_widths, &
                     y_cell_widths_work, z_cell_widths_work, &
                     y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
                     y_coarser_cell_cent_coords, z_coarser_cell_cent_coords)
      END ASSOCIATE

   END SUBROUTINE gmg_init

!-----------------------------------------------------------------------
! SUBROUTINE: gmg_solve
! PURPOSE:
!   Computes the multigrid solution for the given parameters.
!
! INPUTS:
!   solver      - TYPE(MultigridSolver), input/output
!   alpha_dble  - DOUBLE PRECISION, reaction coefficient
!   rhs_dble    - DOUBLE PRECISION array, right-hand side vector
!   tol         - DOUBLE PRECISION, tolerance for convergence
!   max_v_cyc   - INTEGER, maximum number of V-cycles
!
! INPUTS/OUTPUTS:
!   sol_dble    - DOUBLE PRECISION array, solution vector
!
! OUTPUTS:
!   res_norm    - DOUBLE PRECISION, final residual norm
!   iter        - INTEGER, number of iterations performed
!   error_code  - INTEGER, status code indicating solver outcome
!
! ERROR CODES:
!   0 - Success: Converged to the desired tolerance `tol` within `max_v_cyc` iterations.
!   1 - Failure: Iterated `max_v_cyc` iterations but did not converge.
!   2 - Failure: Stagnated after consecutive iterations.
!   3 - Failure: The solution appears to be diverging.
!   4 - Error: Invalid input parameters (e.g., non-positive `tol` or `max_v_cyc`).
!-----------------------------------------------------------------------
   SUBROUTINE gmg_solve(alpha_dble, rhs_dble, tol, max_v_cyc, &
                        sol_dble, res_norm, iter, error_code)
      USE gmg_utils, ONLY: trans_mat_vec_mult_crs_dble, resid_dble
      USE gmg_kernel, ONLY: apply_cycle
      IMPLICIT NONE

      ! Declare input variables
      DOUBLE PRECISION, INTENT(IN) :: alpha_dble
      DOUBLE PRECISION, INTENT(IN) :: rhs_dble(solver%multigrid(1)%cells(3))
      DOUBLE PRECISION, INTENT(IN) :: tol
      INTEGER, INTENT(IN) :: max_v_cyc

      ! Declare input/output variables
      DOUBLE PRECISION, INTENT(INOUT) :: sol_dble(solver%multigrid(1)%cells(3))

      ! Declare output variables
      DOUBLE PRECISION, INTENT(OUT) :: res_norm
      INTEGER, INTENT(OUT) :: iter
      INTEGER, INTENT(OUT) :: error_code  ! New output argument for error code

      ! Declare local variables
      INTEGER :: i, lvls, lvl
      REAL :: alpha_sngl
      DOUBLE PRECISION :: rhs_norm, trans_rhs_norm
      DOUBLE PRECISION :: trans_res_norm
      DOUBLE PRECISION :: res_norm_old, trans_res_norm_old
      INTEGER :: clock_rate, start_time, end_time
      DOUBLE PRECISION :: elapsed_time
      DOUBLE PRECISION :: proj

      ! Initialize error_code to Success
      error_code = 0

      ! Error checking
      IF (tol <= 0.0D0 .OR. max_v_cyc <= 0) THEN
         PRINT *, 'Error in gmg_solve: Tolerance and max_v_cyc must be positive.'
         error_code = 4  ! Invalid input parameters
         RETURN
      END IF

      ASSOCIATE (mg => solver%multigrid, opts => solver%options)
         lvls = SIZE(mg)  ! Get the number of levels in the multigrid structure
         alpha_sngl = REAL(alpha_dble)  ! Convert alpha to single precision

         ! Update the discretization matrices on each grid by modifying their diagonal entries
         DO lvl = 1, lvls
            IF (.NOT. opts%mixed_precision .OR. lvl == 1) THEN
               DO i = 1, mg(lvl)%cells(3)
                  mg(lvl)%A_val_dble(mg(lvl)%A_row_ptr(i)) = mg(lvl)%A_diag_val_dble(i) - alpha_dble
               END DO
            END IF
            IF (opts%mixed_precision) THEN
               DO i = 1, mg(lvl)%cells(3)
                  mg(lvl)%A_val_sngl(mg(lvl)%A_row_ptr(i)) = mg(lvl)%A_diag_val_sngl(i) - alpha_sngl
               END DO
            END IF
         END DO

         ! Compute the 2-norm of the right-hand-side
         rhs_norm = SQRT(DOT_PRODUCT(rhs_dble, rhs_dble))
         IF (alpha_dble == 0.0D0) THEN
            CALL trans_mat_vec_mult_crs_dble(mg(1)%cells(3), &
                                             mg(1)%A_row_ptr, &
                                             SIZE(mg(1)%A_col_ind), &
                                             mg(1)%A_col_ind, &
                                             mg(1)%A_val_dble, &
                                             mg(1)%cells(3), &
                                             rhs_dble, solver%rhs_A_trans)
            trans_rhs_norm = SQRT(DOT_PRODUCT(solver%rhs_A_trans, solver%rhs_A_trans))
         ELSE
            trans_rhs_norm = 0.0D0
         END IF

         ! Check if right-hand-side is numerically zero
         IF (rhs_norm == 0.0D0) THEN
            sol_dble = 0.0D0  ! Set solution to zero
            iter = 0
            res_norm = 0.0D0
            RETURN
         END IF

         ! Perform v-cycles until convergence
         iter = 0
         res_norm = 1.0D5  ! Initialize to a large value
         trans_res_norm = 1.0D5
         res_norm_old = res_norm
         trans_res_norm_old = trans_res_norm

         CALL SYSTEM_CLOCK(count_rate=clock_rate)
         CALL SYSTEM_CLOCK(start_time)

         ! Print initial residual information if required
         IF (opts%debug) THEN
            PRINT *, 'alpha = ', alpha_dble, '.'
            PRINT *, 'The residuals in the relative 2-norm were:'
            IF (alpha_dble == 0.0D0) THEN
               PRINT *, '(diff) = the difference (which lies in range(A)) of successive residuals'
               PRINT *, '(A^T) = A^T of the residual (which lies in domain(A))'
            END IF
         END IF

         ! Compute Residual in Double Precision
         CALL resid_dble(mg(1)%cells(3), &
                         mg(1)%A_row_ptr, &
                         SIZE(mg(1)%A_col_ind), &
                         mg(1)%A_col_ind, &
                         mg(1)%A_val_dble, &
                         mg(1)%cells(3), &
                         sol_dble, rhs_dble, mg(1)%res_dble)

         ! Main loop for multigrid iterations
         DO WHILE (.TRUE.)
            ! Iterative refinement potentially with mixed precision
            CALL apply_cycle(solver, sol_dble, alpha_dble, &
                             opts%full_mg .AND. (iter == 0 .OR. alpha_dble == 0.0), &
                             alpha_dble /= 0.0D0 .AND. opts%mixed_precision)

            ! Remove the right nullspace component of the solution if necessary
            IF ((alpha_dble == 0.0D0) .AND. (opts%remov_null == 'cycl')) THEN
               proj = SUM(sol_dble)/mg(1)%cells(3)
               DO i = 1, mg(1)%cells(3)
                  sol_dble(i) = sol_dble(i) - proj
               END DO
            END IF

            iter = iter + 1

            ! Compute the new residual and its relative 2-norm
            CALL resid_dble(mg(1)%cells(3), &
                            mg(1)%A_row_ptr, &
                            SIZE(mg(1)%A_col_ind), &
                            mg(1)%A_col_ind, &
                            mg(1)%A_val_dble, &
                            mg(1)%cells(3), &
                            sol_dble, rhs_dble, mg(1)%res_dble)

            IF ((alpha_dble .NE. 0.0D0) .OR. (iter == 1)) THEN
               res_norm = SQRT(DOT_PRODUCT(mg(1)%res_dble, &
                                           mg(1)%res_dble))/rhs_norm
            ELSE
               DO i = 1, mg(1)%cells(3)
                  solver%res_old(i) = mg(1)%res_dble(i) - solver%res_old(i)
               END DO
               res_norm = SQRT(DOT_PRODUCT(solver%res_old, solver%res_old))/rhs_norm
            END IF

            ! Calculate the transposed residual norm if necessary
            IF (alpha_dble == 0.0D0) THEN
               CALL trans_mat_vec_mult_crs_dble(mg(1)%cells(3), &
                                                mg(1)%A_row_ptr, &
                                                SIZE(mg(1)%A_col_ind), &
                                                mg(1)%A_col_ind, &
                                                mg(1)%A_val_dble, &
                                                mg(1)%cells(3), &
                                                mg(1)%res_dble, &
                                                solver%res_A_trans)
               trans_res_norm = SQRT(DOT_PRODUCT(solver%res_A_trans, solver%res_A_trans))/ &
                                trans_rhs_norm
            END IF

            ! Print the residual after each iteration for debugging
            IF (opts%debug) THEN
               IF (alpha_dble .NE. 0.0D0) THEN
                  PRINT *, 'iter ', iter, ': ', res_norm, '.'
               ELSEIF (iter == 1) THEN
                  PRINT *, 'iter ', iter, ': ', res_norm, '.'
               ELSE
                  PRINT *, 'iter ', iter, ': ', res_norm, ' (diff), ', trans_res_norm, ' (A^T).'
               END IF
            END IF

            ! Check for convergence
            IF (res_norm < tol) THEN
               IF (opts%debug) THEN
                  CALL SYSTEM_CLOCK(end_time)
                  elapsed_time = REAL(end_time - start_time)/REAL(clock_rate)
                  PRINT *, 'GMG reached the desired tolerance after ', elapsed_time, ' seconds.'
               END IF
               IF ((alpha_dble == 0.0D0) .AND. (opts%remov_null == 'mult')) THEN
                  proj = SUM(sol_dble)/mg(1)%cells(3)
                  DO i = 1, mg(1)%cells(3)
                     sol_dble(i) = sol_dble(i) - proj
                  END DO
               END IF
               EXIT
            END IF

            ! Early stopping conditions based on residual growth
            IF ((res_norm > res_norm_old*2.0D0) .AND. ((alpha_dble .NE. 0.0D0) .OR. &
                                                       (trans_res_norm > trans_res_norm_old*2.0D0))) THEN
               IF (opts%debug) THEN
                  CALL SYSTEM_CLOCK(end_time)
                  elapsed_time = REAL(end_time - start_time)/REAL(clock_rate)
                  PRINT *, 'GMG stopped early after ', elapsed_time, ' seconds and appears to be diverging.'
               END IF
               IF ((alpha_dble == 0.0D0) .AND. (opts%remov_null == 'mult')) THEN
                  proj = SUM(sol_dble)/mg(1)%cells(3)
                  DO i = 1, mg(1)%cells(3)
                     sol_dble(i) = sol_dble(i) - proj
                  END DO
               END IF
               error_code = 3  ! Divergence
               EXIT
            ELSE IF (((res_norm/res_norm_old) > 0.999D0) .AND. ((alpha_dble .NE. 0.0D0) .OR. &
                                                                ((trans_res_norm/trans_res_norm_old) > 0.999D0))) THEN
               IF (opts%debug) THEN
                  CALL SYSTEM_CLOCK(end_time)
                  elapsed_time = REAL(end_time - start_time)/REAL(clock_rate)
                  PRINT *, 'GMG stopped early after ', elapsed_time, ' seconds and appears to be stagnating.'
               END IF
               IF ((alpha_dble == 0.0D0) .AND. (opts%remov_null == 'mult')) THEN
                  proj = SUM(sol_dble)/mg(1)%cells(3)
                  DO i = 1, mg(1)%cells(3)
                     sol_dble(i) = sol_dble(i) - proj
                  END DO
               END IF
               error_code = 2  ! Stagnation
               EXIT
            END IF

            res_norm_old = MIN(res_norm, res_norm_old)  ! Update old residual norms
            trans_res_norm_old = MIN(trans_res_norm, trans_res_norm_old)

            ! Store the residual for the next iteration
            IF (alpha_dble == 0.0D0) THEN
               DO i = 1, mg(1)%cells(3)
                  solver%res_old(i) = mg(1)%res_dble(i)
               END DO
            END IF

            ! Check for maximum number of iterations
            IF (iter >= max_v_cyc) THEN
               error_code = 1  ! Did not converge within max_v_cyc
               ! If maximum iterations reached without convergence
               IF (opts%debug) THEN
                  CALL SYSTEM_CLOCK(end_time)
                  elapsed_time = REAL(end_time - start_time)/REAL(clock_rate)
                  PRINT *, 'GMG reached the maximum number of iterations after ', elapsed_time, ' seconds.'
               END IF
               EXIT
            END IF
         END DO

         IF ((alpha_dble == 0.0D0) .AND. (opts%remov_null == 'mult')) THEN
            proj = SUM(sol_dble)/mg(1)%cells(3)
            DO i = 1, mg(1)%cells(3)
               sol_dble(i) = sol_dble(i) - proj
            END DO
         END IF

      END ASSOCIATE

   END SUBROUTINE gmg_solve

   !-----------------------------------------------------------------------
   ! SUBROUTINE: gmg_finalize
   ! PURPOSE:
   !   Frees memory associated with the multigrid data structure.
   !
   ! INPUTS:
   !   solver    - TYPE(MultigridSolver), input/output
   !-----------------------------------------------------------------------
!   SUBROUTINE gmg_finalize(solver)
   SUBROUTINE gmg_finalize()
      IMPLICIT NONE
!      TYPE(MultigridSolver), INTENT(INOUT) :: solver

      ! Declare local variables
      INTEGER :: i

      ! Error checking to ensure grids were allocated
      IF (.NOT. ALLOCATED(solver%multigrid)) THEN
         PRINT *, "Error: Multigrid data structure is not allocated."
         STOP
      END IF

      ASSOCIATE (mg => solver%multigrid, opts => solver%options)
         ! Deallocate pointers in each individual grid variable
         DO i = 1, SIZE(mg)
            IF (ALLOCATED(mg(i)%sol_sngl)) THEN
               DEALLOCATE (mg(i)%sol_sngl, mg(i)%rhs_sngl)
            END IF
            IF (ALLOCATED(mg(i)%sol_dble)) THEN
               DEALLOCATE (mg(i)%sol_dble, mg(i)%rhs_dble)
            END IF
            IF (ALLOCATED(mg(i)%P_row_ptr)) THEN
               DEALLOCATE (mg(i)%P_row_ptr, mg(i)%P_col_ind)
               IF (ALLOCATED(mg(i)%P_val_dble)) THEN
                  DEALLOCATE (mg(i)%P_val_dble)
               END IF
               IF (ALLOCATED(mg(i)%P_val_sngl)) THEN
                  DEALLOCATE (mg(i)%P_val_sngl)
               END IF
            END IF
            IF (ALLOCATED(mg(i)%A_row_ptr)) THEN
               DEALLOCATE (mg(i)%A_row_ptr, mg(i)%A_col_ind)
               IF (ALLOCATED(mg(i)%A_val_dble)) THEN
                  DEALLOCATE (mg(i)%A_val_dble, mg(i)%A_diag_val_dble)
               END IF
               IF (ALLOCATED(mg(i)%A_val_sngl)) THEN
                  DEALLOCATE (mg(i)%A_val_sngl, mg(i)%A_diag_val_sngl)
               END IF
            END IF
            IF (ALLOCATED(mg(i)%GS_rhs_dble)) THEN
               DEALLOCATE (mg(i)%GS_rhs_dble, mg(i)%GS_up_diag_dble)
            END IF
            IF (ALLOCATED(mg(i)%GS_rhs_sngl)) THEN
               DEALLOCATE (mg(i)%GS_rhs_sngl, mg(i)%GS_up_diag_sngl)
            END IF
            IF (ALLOCATED(mg(i)%res_dble)) THEN
               DEALLOCATE (mg(i)%res_dble)
            END IF
            IF (ALLOCATED(mg(i)%res_sngl)) THEN
               DEALLOCATE (mg(i)%res_sngl)
            END IF
         END DO

         ! Deallocate additional arrays for residuals and transposes
         IF (ALLOCATED(solver%res_old)) THEN
            DEALLOCATE (solver%res_old)
         END IF
         IF (ALLOCATED(solver%rhs_A_trans)) THEN
            DEALLOCATE (solver%rhs_A_trans)
         END IF
         IF (ALLOCATED(solver%res_A_trans)) THEN
            DEALLOCATE (solver%res_A_trans)
         END IF
      END ASSOCIATE

      ! Deallocate the grids array
      DEALLOCATE (solver%multigrid)

   END SUBROUTINE gmg_finalize

END MODULE gmg
