#ifdef USE_DOUBLE_PRECISION

#define KIND_generic KIND_dble
#define GS_generic GS_dble
#define resid_generic resid_dble
#define mat_vec_mult_crs_generic mat_vec_mult_crs_dble
#define mat_vec_mult_add_crs_generic mat_vec_mult_add_crs_dble

#define apply_cycle_generic apply_cycle_dble
#define downstroke_generic downstroke_dble
#define upstroke_generic upstroke_dble
#define coarse_solve_generic coarse_solve_dble
#define A_val_generic A_val_dble
#define rhs_generic rhs_dble
#define sol_generic sol_dble
#define GS_rhs_generic GS_rhs_dble
#define GS_up_diag_generic GS_up_diag_dble
#define res_generic res_dble
#define P_val_generic P_val_dble
#define R_val_generic R_val_dble

#elif defined(USE_SINGLE_PRECISION)

#define KIND_generic KIND_sngl
#define GS_generic GS_sngl
#define resid_generic resid_sngl
#define mat_vec_mult_crs_generic mat_vec_mult_crs_sngl
#define mat_vec_mult_add_crs_generic mat_vec_mult_add_crs_sngl

#define apply_cycle_generic apply_cycle_sngl
#define downstroke_generic downstroke_sngl
#define upstroke_generic upstroke_sngl
#define coarse_solve_generic coarse_solve_sngl
#define A_val_generic A_val_sngl
#define rhs_generic rhs_sngl
#define sol_generic sol_sngl
#define GS_rhs_generic GS_rhs_sngl
#define GS_up_diag_generic GS_up_diag_sngl
#define res_generic res_sngl
#define P_val_generic P_val_sngl
#define R_val_generic R_val_sngl
#endif

!-----------------------------------------------------------------------
! SUBROUTINE: apply_cycle_generic
! PURPOSE:
!   Performs a single V-cycle or Full Multigrid (FMG) cycle.
!
! INPUTS:
!   solver     - TYPE(MultigridSolver), input/output
!   alpha      - REAL, reaction coefficient
!   sol_dble   - DOUBLE PRECISION, solution vector in double precision
!   fmg        - LOGICAL. If true, perform FMG cycle.
!-----------------------------------------------------------------------
SUBROUTINE apply_cycle_generic(solver, sol_dble, alpha, fmg)
   USE gmg_data, ONLY: MultigridSolver, Grid
   USE gmg_utils, ONLY: KIND_dble, KIND_sngl
   IMPLICIT NONE

   ! Input/output parameters
   TYPE(MultigridSolver), INTENT(INOUT) :: solver
   REAL(KIND=KIND_generic), INTENT(IN) :: alpha
   REAL(KIND=KIND_dble), INTENT(INOUT) :: sol_dble(:)
   LOGICAL, INTENT(IN) :: fmg

   ! Local variables
   INTEGER :: i

   ASSOCIATE (mg => solver%multigrid)
      ! Solve for correction in chosen precision
      DO i = 1, mg(1)%cells(3)
         ! Convert the residual to the appropriate precision
         mg(1)%rhs_generic(i) = REAL(mg(1)%res_dble(i), KIND=KIND_generic)
         mg(1)%sol_generic(i) = 0.0
      END DO

      ! Down-stroke: Descend the grid hierarchy
      CALL downstroke_generic(solver, alpha)

      ! Coarse-grid solve
      CALL coarse_solve_generic(solver, alpha)

      IF (fmg) THEN
         ! Perform Full Multigrid (FMG) cycle
         DO i = SIZE(mg) - 1, 2, -1
            ! Up-stroke: Ascend the grid hierarchy to level 'i'
            CALL upstroke_generic(solver, alpha, i)
            ! Down-stroke: Descend the grid hierarchy from level 'i'
            CALL downstroke_generic(solver, alpha, i)
            ! Coarse-grid solve
            CALL coarse_solve_generic(solver, alpha)
         END DO
      END IF

      ! Up-stroke: Ascend the grid hierarchy
      CALL upstroke_generic(solver, alpha)

      ! Update double-precision solution
      DO i = 1, mg(1)%cells(3)
         sol_dble(i) = sol_dble(i) + DBLE(mg(1)%sol_generic(i))
      END DO
   END ASSOCIATE
END SUBROUTINE apply_cycle_generic

!-----------------------------------------------------------------------
! SUBROUTINE: downstroke_generic
! PURPOSE:    Performs the down-stroke (restriction phase) of the V-cycle.
!
! INPUTS:
!   solver - MultigridSolver type, input/output
!   alpha      - REAL, reaction term in Helmholtz equation
!   lvl    - INTEGER, starting multigrid level (optional). Default is 1.
!
! OUTPUTS:
!   solver - Updated multigrid data structures after down-stroke
!            (solution vectors, residuals, rhs at coarser grids)
!-----------------------------------------------------------------------
SUBROUTINE downstroke_generic(solver, alpha, lvl_in)
   USE gmg_utils, ONLY: GS_generic, resid_generic, mat_vec_mult_crs_generic, KIND_generic
   IMPLICIT NONE

   ! Input/output parameters
   TYPE(MultigridSolver), INTENT(INOUT) :: solver
   REAL(KIND=KIND_generic), INTENT(IN) :: alpha
   INTEGER, INTENT(IN), OPTIONAL :: lvl_in

   ! Local variables
   INTEGER :: i, lvls, lvl, iter
   CHARACTER(LEN=5) :: GS_type

   IF (.NOT. PRESENT(lvl_in)) THEN
      lvl = 1
   ELSE
      lvl = lvl_in
   END IF

   lvls = SIZE(solver%multigrid)
   ! Down-stroke: Descend the grid hierarchy
   DO i = lvl, lvls - 1
      ASSOCIATE (finer_grid => solver%multigrid(i), coarser_grid => solver%multigrid(i + 1), opts => solver%options)
         IF (alpha == 0.0 .AND. i == 1) THEN
            ! Double the number of iterations for the first level if alpha=0
            iter = opts%pre_smooth_iter*2
         ELSE
            iter = opts%pre_smooth_iter
         END IF

         ! Determine GS_type
         GS_type = determine_GS_type(opts%GS_type, alpha == 0.0, &
                                     solver%multigrid(1)%cells(1), &
                                     solver%multigrid(1)%cells(2), &
                                     finer_grid%cells(1), &
                                     finer_grid%cells(2))

         ! Pre-smoothing on finer grid
         CALL GS_generic(finer_grid%cells(1), finer_grid%cells(2), &
                         finer_grid%A_row_ptr, SIZE(finer_grid%A_col_ind), &
                         finer_grid%A_col_ind, finer_grid%A_val_generic, &
                         finer_grid%rhs_generic, GS_type, iter, &
                         finer_grid%sol_generic, finer_grid%GS_rhs_generic, finer_grid%GS_up_diag_generic)

         ! Compute residual after pre-smoothing
         CALL resid_generic(finer_grid%cells(3), finer_grid%A_row_ptr, &
                            SIZE(finer_grid%A_col_ind), finer_grid%A_col_ind, &
                            finer_grid%A_val_generic, finer_grid%cells(3), &
                            finer_grid%sol_generic, finer_grid%rhs_generic, finer_grid%res_generic)

         ! Restrict residual to coarser grid using the restriction operator
         CALL mat_vec_mult_crs_generic(coarser_grid%cells(3), & ! coarse grid size
                                       finer_grid%R_row_ptr, &   ! Row pointers of R
                                       SIZE(finer_grid%R_col_ind), & ! Number of non-zero entries
                                       finer_grid%R_col_ind, & ! Column indices of R
                                       finer_grid%R_val_generic, & ! Values of R
                                       finer_grid%cells(3), & ! fine grid size
                                       finer_grid%res_generic, & ! residual on fine grid
                                       coarser_grid%rhs_generic) ! rhs on coarse grid

         ! Initialize solution vector on coarser grid to zero
         coarser_grid%sol_generic = 0.0
      END ASSOCIATE
   END DO

END SUBROUTINE downstroke_generic

!-----------------------------------------------------------------------
! SUBROUTINE: coarse_solve_generic
! PURPOSE:    Performs the coarse-grid solution phase of the V-cycle
!
! INPUTS:
!   solver     - MultigridSolver type, input/output
!   alpha      - REAL, reaction term in Helmholtz equation
!
! OUTPUTS:
!   solver     - Updated coarse grid solution
!-----------------------------------------------------------------------
SUBROUTINE coarse_solve_generic(solver, alpha)
   USE gmg_utils, ONLY: GS_generic, KIND_generic
   IMPLICIT NONE

   ! Input/output parameters
   TYPE(MultigridSolver), INTENT(INOUT) :: solver
   REAL(KIND=KIND_generic), INTENT(IN) :: alpha

   ! Local variables
   INTEGER :: iter, lvls
   CHARACTER(LEN=5) :: GS_type

   ! Associate pointers for convenience
   ASSOCIATE (mg => solver%multigrid, opts => solver%options)
      lvls = SIZE(mg)

      IF (alpha == 0.0) THEN
         iter = opts%coarse_solve_iter*2  ! double the number of iterations
      ELSE
         iter = opts%coarse_solve_iter
      END IF
      ! Determine GS_type
      GS_type = determine_GS_type(opts%GS_type, alpha == 0.0, &
                                  mg(1)%cells(1), &
                                  mg(1)%cells(2), &
                                  mg(lvls)%cells(1), &
                                  mg(lvls)%cells(2))

      ! Perform Gauss-Seidel iterations on the coarse grid
      CALL GS_generic(mg(lvls)%cells(1), mg(lvls)%cells(2), &
                      mg(lvls)%A_row_ptr, SIZE(mg(lvls)%A_col_ind), &
                      mg(lvls)%A_col_ind, mg(lvls)%A_val_generic, mg(lvls)%rhs_generic, &
                      GS_type, iter, &
                      mg(lvls)%sol_generic, mg(lvls)%GS_rhs_generic, mg(lvls)%GS_up_diag_generic)

   END ASSOCIATE

END SUBROUTINE coarse_solve_generic

!-----------------------------------------------------------------------
! SUBROUTINE: upstroke_generic
! PURPOSE:    Performs the up-stroke (prolongation phase) of the V-cycle.
!             The solution and rhs vectors are in solver%multigrid(lvl)%sol_<dble|sngl>
!             and solver%multigrid(lvl)%rhs_<dble|sngl>, respectively.
!
! INPUTS:
!   solver - MultigridSolver type, input/output
!   alpha - REAL(KIND=KIND_generic), scaling factor
!   lvl    - INTEGER, starting multigrid level (optional). Default is 1.
!
! OUTPUTS:
!   solver - Updated solution vectors after up-stroke
!-----------------------------------------------------------------------
SUBROUTINE upstroke_generic(solver, alpha, lvl_in)
   USE gmg_utils, ONLY: GS_generic, mat_vec_mult_add_crs_generic, KIND_generic
   IMPLICIT NONE

   ! Input/output parameters
   TYPE(MultigridSolver), INTENT(INOUT) :: solver
   REAL(KIND=KIND_generic), INTENT(IN) :: alpha
   INTEGER, INTENT(IN), OPTIONAL :: lvl_in

   ! Local variables
   INTEGER :: i, lvls, lvl, iter
   REAL(KIND=KIND_generic) :: proj
   CHARACTER(LEN=5) :: GS_type

   IF (.NOT. PRESENT(lvl_in)) THEN
      lvl = 1
   ELSE
      lvl = lvl_in
   END IF

   lvls = SIZE(solver%multigrid)

   DO i = lvls, lvl + 1, -1
      ASSOCIATE (finer_grid => solver%multigrid(i - 1), coarser_grid => solver%multigrid(i), opts => solver%options)
         IF (alpha == 0.0 .AND. i == 2) THEN
            ! Double the number of iterations for the first level if alpha=0
            iter = opts%post_smooth_iter*2
         ELSE
            iter = opts%post_smooth_iter
         END IF

         ! Correct the solution on finer grid after prolongation
         CALL mat_vec_mult_add_crs_generic(finer_grid%cells(3), coarser_grid%P_row_ptr, &
                                           SIZE(coarser_grid%P_col_ind), coarser_grid%P_col_ind, &
                                           coarser_grid%P_val_generic, coarser_grid%cells(3), &
                                           coarser_grid%sol_generic, finer_grid%sol_generic)

         IF ((alpha == 0.0) .AND. (opts%remov_null == 'grid')) THEN
            ! Remove the right nullspace component of the solution
            proj = SUM(finer_grid%sol_generic)/REAL(finer_grid%cells(3), KIND=KIND_generic)
            finer_grid%sol_generic = finer_grid%sol_generic - proj
         END IF

         ! Determine GS_type
         GS_type = determine_GS_type(opts%GS_type, alpha == 0.0, &
                                     solver%multigrid(1)%cells(1), &
                                     solver%multigrid(1)%cells(2), &
                                     finer_grid%cells(1), &
                                     finer_grid%cells(2))

         ! Post-smoothing on finer grid
         CALL GS_generic(finer_grid%cells(1), finer_grid%cells(2), &
                         finer_grid%A_row_ptr, SIZE(finer_grid%A_col_ind), &
                         finer_grid%A_col_ind, finer_grid%A_val_generic, &
                         finer_grid%rhs_generic, GS_type, iter, &
                         finer_grid%sol_generic, finer_grid%GS_rhs_generic, finer_grid%GS_up_diag_generic)
      END ASSOCIATE
   END DO

END SUBROUTINE upstroke_generic

#undef KIND_generic
#undef GS_generic
#undef resid_generic
#undef mat_vec_mult_crs_generic
#undef mat_vec_mult_add_crs_generic

#undef apply_cycle_generic
#undef downstroke_generic
#undef upstroke_generic
#undef coarse_solve_generic
#undef A_val_generic
#undef rhs_generic
#undef sol_generic
#undef GS_rhs_generic
#undef GS_up_diag_generic
#undef res_generic
#undef P_val_generic
#undef R_val_generic