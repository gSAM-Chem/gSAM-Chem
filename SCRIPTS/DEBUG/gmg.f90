! Author: Oliver Yang
! Adapted from code written by Cao Lu.

! Last updated: November 3, 2017

! Please refer to the test scripts included with this file for examples of how to use the library.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Defines a module containing the multigrid data structure to be shared by all subroutines of the Helmoltz solver. All arrays are allocated in pressure_gmg_init and the subroutines called therein. All arrays are deallocated in pressure_gmg_fin.
module pressure_gmg_intern_data
implicit none

! Defines a derived type that contains the data for a snglle grid of multigrid. This includes the number of cells on the grid, the matrix of the rediscretized PDE, the solution vector, the right-hand-side vector, the buffers for the smoother and tridiagonal solver, the residual vector, the restriction operator for the residual on the given grid to the next coarser grid, and the prolongation operator for the solution on the given grid to the next finer grid. The solution vector on the fine-grid will not be allocated, since the array passed into pressure_gmg_solv is used as an initial guess for multigrid and will be overwritten with a new solution iterate in each v-cycle. The right-hand-side on the fine-grid is passed into pressure_gmg_solv as well. The residual on the fine-grid will not be allocated, since it does not need to be computed. The restriction matrix on the coarse-grid and the prolongation matrix on the fine-grid will not be allocated either, since they are not applicable.
type grid
    ! Declares an array that stores the number of cells in the y-direction, the z-direction, and their product (the total number of cells).
    integer, pointer :: cells( :)
    ! Declares the PDE discretization matrix in CRS format.

    ! We store sparse matrices, in which the number of non-zero entries is small compared to the total number of entries in the matrix, in CRS (compressed row storage) format. In CRS, the non-zero values are listed in an array 'A_val' and are given row by row. There is an array 'A_row_ptr' that stores the index in A_val of the first listed entry for each row of the matrix. Our implementation here of CRS creates an additional final entry in A_row_ptr equal to the number of entries in A_val plus 1 (this will easily tell us the number of non-zero entries in the last row of the matrix). Therefore A_row_ptr has length equal to the number of rows plus 1. The entries in A_row_ptr are strictly increasing since the values in A_val are given in row order. There is another array 'A_col_ind' that stores the column indices in the matrix of all the corresponding entries in A_val. Note that the entries in A_col_ind corresponding to a row of the matrix do not have to be increasing (that is, the entries within a row can be listed in any order in CRS format). Both A_val and A_col_ind, of course, have length equal to the number of non-zero entries in the matrix. For more details, see, for example, http://netlib.org/linalg/html_templates/node90.html. Further optimization (not currently implemented) would require the use of block CRS format to be able to take advantage of the architecture-specific block-optimizations in Level 3 BLAS routines for matrix-matrix multiplication (see http://www.netlib.org/blas/).

    ! Recall that the cells on each grid are numbered in the direction of increasing coordinates, in z-major order. That is, consecutive cells run along the z-axis.

    ! The values along the diagonal are copied and stored in another array. The multigrid data structure is set up in pressure_gmg_init for the case when the harmonic alpha = 0. The pressure solver is called in a loop through every non-zero alpha. The diagonal entries of the discretization matrices on each grid are modified with each call. Storing the values for alpha = 0 separately and computing the new values each time directly from the stored values prevents loss of precision.
    integer, pointer :: A_row_ptr( :), A_col_ind( :)
    double precision, pointer :: A_val_dble( :), A_diag_val_dble( :)
    real, pointer :: A_val_sngl( :), A_diag_val_sngl( :)
    ! Declares the solution.
    double precision, pointer :: sol_dble( :)
    real, pointer :: sol_sngl( :)
    ! Declares the right-hand-side.
    double precision, pointer :: rhs_dble( :)
    real, pointer :: rhs_sngl( :)
    ! Declares the buffer for the smoother.
    double precision, pointer :: GS_rhs_dble( :)
    real, pointer :: GS_rhs_sngl( :)
    ! Declares the buffer for the tridiagonal solver.
    double precision, pointer :: GS_up_diag_dble( :)
    real, pointer :: GS_up_diag_sngl( :)
    ! Declares the residual.
    double precision, pointer :: res_dble( :)
    real, pointer :: res_sngl( :)
    ! Declares the restriction matrix.
    integer, pointer :: R_row_ptr( :), R_col_ind( :)
    double precision, pointer :: R_val_dble( :)
    real, pointer :: R_val_sngl( :)
    ! Declares the prolongation matrix.
    integer, pointer :: P_row_ptr( :), P_col_ind( :)
    double precision, pointer :: P_val_dble( :)
    real, pointer :: P_val_sngl( :)
end type grid

! Declares an array of grid variables.
type( grid), allocatable :: multigrid( :)

! We convert the coarse-grid discretization matrix to dense format to avoid using sparse solvers on the coarse-grid. The coarse-grid is small enough that using dense solvers does not have an appreciable effect on efficiency.

! Declares the LU decomposition (with partial pivoting) of the dense coarse-grid discretization matrix, P * A = L * U. L is a lower-triangular matrix with 1s on the diagonal and other entries with magnitude < 1. U is upper-triangular. P is a (row) permutation matrix. We write U and the strictly lower-triangular part of L onto LU and store the permutation as a vector p).
double precision, allocatable :: LU_coarse_dble( :, :)
real, allocatable :: LU_coarse_sngl( :, :)
integer, allocatable :: p_coarse( :)

! Declares the permuted right-hand-side on the coarse-grid.
double precision, allocatable :: p_rhs_coarse_dble( :)
real, allocatable :: p_rhs_coarse_sngl( :)

! The default choice for coarse-grid solver is to use a Krylov subspace method, such as GMRES (generalized minimum residual). When alpha becomes small, however, the discretization matrices become ill-conditioned, and GMRES yields small residuals but large errors. A direct solver such as LU is both more accurate and more stable.

! When the harmonic is 0, the PDEs on each grid become singular Poisson-like and their discretization matrices become singular.

! Declares a copy of the residual on the fine-grid from the previous v_cycle in order to compute the first termination condition.
double precision, allocatable :: res_old( :)
! Declares the transpose of the discretization matrix on the fine-grid applied to the right-hand-side and the residual in order to compute the second termination condition.
double precision, allocatable :: rhs_A_trans( :)
double precision, allocatable :: res_A_trans( :)

! Declares the QR decomposition with column pivoting (QRCP) of the dense coarse-grid discretization matrix.
double precision, allocatable :: Q_coarse_dble( :, :), R_coarse_dble( :, :)
real, allocatable :: Q_coarse_sngl( :, :), R_coarse_sngl( :, :)
integer, allocatable :: e_coarse( :)
double precision, allocatable :: QR_rhs_dble( :)
real, allocatable :: QR_rhs_sngl( :)

! Declares the singular value decomposition (SVD) of the dense coarse-grid discretization matrix, A = U * S * VT. U is the orthogonal matrix of left singular vectors (along the columns). S is a vector with the (non-negative) singular values in decreasing order. VT is the transposed orthogonal matrix of right singular vectors (along the rows).
! double precision, allocatable :: U_coarse( :, :), S_coarse( :), &
!     VT_coarse( :, :)
! real, allocatable :: U_coarse_sngl( :, :), S_coarse_sngl( :), &
!     VT_coarse_sngl( :, :)

end module pressure_gmg_intern_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Defines a module containing the internal subroutines common to both the pressure and stream function solvers.
module gmg_common_intern_proced
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Creates a grid function (e.g. the solution or right-hand-side) for use in testing the multigrid solver and its restriction and prolongation operators. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine create_func( y_points, z_points, y_rel_coords, z_rel_coords, &
    y_range, z_range, func_type, func)

! Declares the input variables.
! Declares the number of points at which the function is to be defined in the y and the z-directions.
integer, intent( in) :: y_points, z_points
! Declares the y and z-coordinates of the points.
double precision, intent( in) :: y_rel_coords( y_points), z_rel_coords( &
    z_points)
! Declares the range in coordinates of the grid in the y and the z-directions.
double precision, intent( in) :: y_range, z_range
! Declares the type of function to create.
character( len = 3), intent( in) :: func_type
! Declares the function. The array will be overwritten.
double precision, intent( out) :: func( y_points * z_points)

! Declares the local variables.
integer :: i, j, ind
double precision, parameter :: pi = 3.141592653589793d+0

select case( func_type)

case( 'cos')
! Computes a product of cosine functions. It satisfies homogeneous Neumann boundary conditions on the grid.
ind = 1
do i = 1, size( y_rel_coords)
    ! The inner loop runs over the z-coordinates.
    do j = 1, size( z_rel_coords)
        func( ind) = cos( pi * y_rel_coords( i) / y_range) * cos( pi * &
            z_rel_coords( j) / z_range)
        ind = ind + 1
    end do
end do

case( 'sin')
! Computes a product of sine functions. It satisfies homogeneous Dirichlet boundary conditions on the grid.
ind = 1
do i = 1, size( y_rel_coords)
    ! The inner loop runs over the z-coordinates.
    do j = 1, size( z_rel_coords)
        func( ind) = sin( pi * y_rel_coords( i) / y_range) * sin( pi * &
            z_rel_coords( j) / z_range)
        ind = ind + 1
    end do
end do

case default
! Computes a bilinear function.
ind = 1
do i = 1, size( y_rel_coords)
    ! The inner loop runs over the z-coordinates.
    do j = 1, size( z_rel_coords)
        func( ind) = (1.0d+0 + y_rel_coords( i) * 2.0d+0 / y_range) * &
            z_rel_coords( j) / z_range
        ind = ind + 1
    end do
end do

end select

end subroutine create_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the diffusion coefficients.
subroutine diff_coeff( length, mu, mu_is_const, rad_reduc_factor)

!use params, only: rad_earth

! Declares the input variables.
! Declares the number of coordinates at which to compute the coefficients.
integer, intent( in) :: length
! Declares the array of coordinates. The routine overwrites it with the coefficients at those coordinates.
double precision, intent( inout) :: mu( length)
! Declares a flag for the coefficient mu. If .false., mu is computed from the standard formula. If .true., all mu are set to 1 in order to compare the GMG solution to the FFT solution.
logical, intent( in) :: mu_is_const
! Declares the reduction factor for the radius of the earth, e.g. if rad_reduc_factor is 2, the radius is halved.
double precision, intent( in) :: rad_reduc_factor

! Declares the local variables.
integer :: i
real(8), parameter :: rad_earth = 6371229.d0
double precision, parameter :: irad = 1.d0/rad_earth
if( mu_is_const) then
    do i = 1, length
        mu( i) = 1.0d+0
    end do
else
    do i = 1, length
        mu( i) = cos( mu( i) * irad * rad_reduc_factor)
    end do
end if


end subroutine diff_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the LU decomposition (with partial pivoting) of a square matrix in-situ.
subroutine LU_part_pivot_dble( A_size, A, p)

! Declares the input variables.
! Declares the number of rows (or columns) of the matrix.
integer, intent( in) :: A_size
! Declares the matrix. It will be overwritten.
double precision, intent( inout) :: A( A_size, A_size)
! Declares the permutation vector. It will be overwritten.
integer, intent( out) :: p( A_size)

! Declares the local variables.
! Declares indices for the matrix.
integer :: mat_row, submat_row, col, pivot_row
! Declares a buffer to perform swaps in the permutation vector.
integer :: p_temp

! Intializes (or resets) the entries of the permutation vector to be equal to the corresponding indices.
do mat_row = 1, A_size
    p( mat_row) = mat_row
end do

! Loops through the rows of the n * n matrix. At the ith iteration, the ith (permuted) row is untouched and the lower right (n - i) * (n - i) block (permuted) is modified. Only n - 1 iterations are necessary, as the last pivot is the last diagonal entry and does not need to be changed.
do mat_row = 1, A_size - 1
    ! Finds the row with the entry of the largest magnitude in the (n - i)th column of the lower right (n - i) * (n - i) block, i.e. the pivot.
    pivot_row = mat_row
    do submat_row = mat_row + 1, A_size
        if( abs( A( p( submat_row), mat_row)) > abs( A( p( pivot_row), &
            mat_row))) then
            pivot_row = submat_row
        end if
    end do
    ! If the row with the pivot is not the next row, swaps those two entries of the permutation vector. Since the rows will henceforth be accessed by the corresponding entry of the permutation vector, we may imagine that the rows of the original matrix have been swapped.
    if( pivot_row > mat_row) then
        p_temp = p( mat_row)
        p( mat_row) = p( pivot_row)
        p( pivot_row) = p_temp
    end if
    ! All the entries in the column below the pivot are divided by the pivot. This gives the multiple of the ith row to be subtracted in Gaussian elimination.
    do submat_row = mat_row + 1, A_size
        A( p( submat_row), mat_row) = A( p( submat_row), mat_row) / A( p( &
            mat_row), mat_row)
    end do
    ! Subtracts multiples of the ith row from all the lower rows.
    do col = mat_row + 1, A_size
        ! Loops through the rows in the inner loop to optimize for cache performance since Fortran stores matrices in column-major order.
        do submat_row = mat_row + 1, A_size
            A( p( submat_row), col) = A( p( submat_row), col) - A( p( &
                submat_row), mat_row) * A( p( mat_row), col)
        end do
    end do
end do

end subroutine LU_part_pivot_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the LU decomposition in single precision.
subroutine LU_part_pivot_sngl( A_size, A, p)

integer, intent( in) :: A_size
real, intent( inout) :: A( A_size, A_size)
integer, intent( out) :: p( A_size)

integer :: mat_row, submat_row, col, pivot_row
integer :: p_temp

do mat_row = 1, A_size
    p( mat_row) = mat_row
end do

do mat_row = 1, A_size - 1
    pivot_row = mat_row
    do submat_row = mat_row + 1, A_size
        if( abs( A( p( submat_row), mat_row)) > abs( A( p( pivot_row), &
            mat_row))) then
            pivot_row = submat_row
        end if
    end do
    if( pivot_row > mat_row) then
        p_temp = p( mat_row)
        p( mat_row) = p( pivot_row)
        p( pivot_row) = p_temp
    end if
    do submat_row = mat_row + 1, A_size
        A( p( submat_row), mat_row) = A( p( submat_row), mat_row) / A( p( &
            mat_row), mat_row)
    end do
    do col = mat_row + 1, A_size
        do submat_row = mat_row + 1, A_size
            A( p( submat_row), col) = A( p( submat_row), col) - A( p( &
                submat_row), mat_row) * A( p( mat_row), col)
        end do
    end do
end do

end subroutine LU_part_pivot_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Applies the Gauss-Seidel iteration on a Helmholtz discretization matrix. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order. It also assumes that the entries are given in the order of the center cell, then the cell below (if applicable), then the cell above (if applicable), then the cell to the left (if applicable), and then the cell to the right (if applicable).
subroutine GS_dble( y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, &
    A_val_dble, sol_dble, rhs_dble, GS_type, GS_rhs_dble, &
    GS_up_diag_dble, iter)

! Declares the input variables.
! Declares the number of grid cells in the y and z-directions.
integer, intent( in) :: y_cells, z_cells
! Declares the row pointer array of the matrix.
integer, intent( in) :: A_row_ptr( y_cells * z_cells + 1)
! Declares the number of non-zero entries of the matrix.
integer, intent( in) :: A_entries
! Declares the entry column index array of the matrix.
integer, intent( in) :: A_col_ind( A_entries)
! Declares the entry value array of the matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the solution. The array will be overwritten.
double precision, intent( inout) :: sol_dble( y_cells * z_cells)
! Declares the right-hand-side.
double precision, intent( in) :: rhs_dble( y_cells * z_cells)
! Declares the type of smoother.
character( len = 5), intent( in) :: GS_type
! Declares the buffer for the smoother. The array will be overwritten.
double precision, intent( out) :: GS_rhs_dble( z_cells)
! Declares the buffer for the tridiagonal solver. The array will be overwritten.
double precision, intent( out) :: GS_up_diag_dble( z_cells - 1)
! Declares the number of smoothing iterations to perform.
integer, intent( in) :: iter

! Declares the local variables.
integer :: i, j, k
! Declares the starting index.
integer :: start_ind
! Declares a flag for whether the number of cells in the y-direction is odd in zebra z-line G-S.
logical :: y_cells_is_odd
! Declares a temporary to store the right-hand-side entry in point-wise G-S.
double precision :: temp

select case( GS_type)

case( 'zline')
! Point-wise G-S is ineffective in this case, since anisotropy in the grid cause errors to be geometrically non-smooth. The theoretical reasons for this are discussed in Trottenberg, 5.1. The solution is to group together strongly connected variables. Since the cell widths are much smaller in the z-direction, the coefficients in the z-direction are much larger. Block G-S uses tridiagonal solves to update these unknowns together. That is why cells are numbered in z-major order.
do i = 1, iter
    ! Smooths the first block.
    ! Computes the new right-hand-side. Multiplies the coefficients for the blocks below with the corresponding solution entries and subtracts them from the right-hand-side.
    do j = 1, z_cells
        GS_rhs_dble( j) = rhs_dble( j) - sol_dble( j + z_cells) * &
            A_val_dble( A_row_ptr( j + 1) - 1)
    end do
    call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_dble, z_cells, 1, sol_dble, GS_rhs_dble, GS_up_diag_dble)

    ! Smooths the intermediate blocks.
    do j = 2, y_cells - 1
        start_ind = (j - 1) * z_cells
        ! Multiplies the coefficients for the blocks above and below with the corresponding solution entries.
        do k = 1, z_cells
            GS_rhs_dble( k) = rhs_dble( start_ind + k) - sol_dble( &
                start_ind + k - z_cells) * A_val_dble( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_dble( start_ind + k + &
                z_cells) * A_val_dble( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_dble, z_cells, j, sol_dble, GS_rhs_dble, GS_up_diag_dble)
    end do

    ! Smooths the last block.
    start_ind = (y_cells - 1) * z_cells
    ! Multiplies the coefficients for the blocks above with the corresponding solution entries.
    do j = 1, z_cells
        GS_rhs_dble( j) = rhs_dble( start_ind + j) - sol_dble( start_ind &
            + j - z_cells) * A_val_dble( A_row_ptr( start_ind + j + 1) - 1)
    end do
    call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_dble, z_cells, y_cells, sol_dble, GS_rhs_dble, &
        GS_up_diag_dble)
end do

case( 'zebra')
! Our tests have shown that sequential sweeps (in the y-direction) of z-line G-S is more effective than zebra-lining.

y_cells_is_odd = (mod( y_cells, 2) > 0)

do i = 1, iter
    ! Smooths the first block.
    ! Computes the new right-hand-side. Multiplies the coefficients for the blocks below with the corresponding solution entries and subtracts them from the right-hand-side.
    do j = 1, z_cells
        GS_rhs_dble( j) = rhs_dble( j) - sol_dble( j + z_cells) * &
            A_val_dble( A_row_ptr( j + 1) - 1)
    end do
    call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_dble, z_cells, 1, sol_dble, GS_rhs_dble, GS_up_diag_dble)

    ! Smooths the intermediate odd-numbered blocks.
    do j = 3, y_cells - 1, 2
        start_ind = (j - 1) * z_cells
        ! Multiplies the coefficients for the blocks above and below with the corresponding solution entries.
        do k = 1, z_cells
            GS_rhs_dble( k) = rhs_dble( start_ind + k) - sol_dble( &
                start_ind + k - z_cells) * A_val_dble( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_dble( start_ind + k + &
                z_cells) * A_val_dble( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_dble, z_cells, j, sol_dble, GS_rhs_dble, GS_up_diag_dble)
    end do

    ! Smooths the last block if it is odd-numbered.
    if( y_cells_is_odd) then
        start_ind = (y_cells - 1) * z_cells
        ! Multiplies the coefficients for the blocks above with the corresponding solution entries.
        do j = 1, z_cells
            GS_rhs_dble( j) = rhs_dble( start_ind + j) - sol_dble( &
                start_ind + j - z_cells) * A_val_dble( A_row_ptr( &
                start_ind + j + 1) - 1)
        end do
        call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_dble, z_cells, y_cells, sol_dble, GS_rhs_dble, &
            GS_up_diag_dble)
    end if

    ! Smooths the intermediate even-numbered blocks.
    do j = 2, y_cells - 1, 2
        start_ind = (j - 1) * z_cells
        ! Multiplies the coefficients for the blocks above and below with the corresponding solution entries.
        do k = 1, z_cells
            GS_rhs_dble( k) = rhs_dble( start_ind + k) - sol_dble( &
                start_ind + k - z_cells) * A_val_dble( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_dble( start_ind + k + &
                z_cells) * A_val_dble( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_dble, z_cells, j, sol_dble, GS_rhs_dble, GS_up_diag_dble)
    end do

    ! Smooths the last block if it is even-numbered.
    if( .not.y_cells_is_odd) then
        start_ind = (y_cells - 1) * z_cells
        ! Multiplies the coefficients for the blocks above with the corresponding solution entries.
        do j = 1, z_cells
            GS_rhs_dble( j) = rhs_dble( start_ind + j) - sol_dble( &
                start_ind + j - z_cells) * A_val_dble( A_row_ptr( &
                start_ind + j + 1) - 1)
        end do
        call tridiag_solv_dble( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_dble, z_cells, y_cells, sol_dble, GS_rhs_dble, &
            GS_up_diag_dble)
    end if
end do

case default
! Performs point-wise G-S.
do i = 1, iter
    do j = 1, y_cells * z_cells
        ! Copies the entry of the right-hand-side for the current cell.
        temp = rhs_dble( j)
        ! Multiplies the coefficients for the other cells in the row with the corresponding solution entries and subtracts them from the right-hand-side entry.
        do k = A_row_ptr( j) + 1, A_row_ptr( j + 1) - 1
            temp = temp - A_val_dble( k) * sol_dble( A_col_ind( k))
        end do
        ! Computes the new solution entry.
        sol_dble( j) = temp / A_val_dble( A_row_ptr( j))
    end do
end do

end select

end subroutine GS_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Applies the Gauss-Seidel iteration in single precision.
subroutine GS_sngl( y_cells, z_cells, A_row_ptr, A_entries, A_col_ind, &
    A_val_sngl, sol_sngl, rhs_sngl, GS_type, GS_rhs_sngl, &
    GS_up_diag_sngl, iter)

integer, intent( in) :: y_cells, z_cells
integer, intent( in) :: A_row_ptr( y_cells * z_cells + 1)
integer, intent( in) :: A_entries
integer, intent( in) :: A_col_ind( A_entries)
real, intent( in) :: A_val_sngl( A_entries)
real, intent( inout) :: sol_sngl( y_cells * z_cells)
real, intent( in) :: rhs_sngl( y_cells * z_cells)
character( len = 5), intent( in) :: GS_type
real, intent( out) :: GS_rhs_sngl( z_cells)
real, intent( out) :: GS_up_diag_sngl( z_cells - 1)
integer, intent( in) :: iter

integer :: i, j, k
integer :: start_ind
logical :: y_cells_is_odd
real :: temp

select case( GS_type)

case( 'zline')
do i = 1, iter
    do j = 1, z_cells
        GS_rhs_sngl( j) = rhs_sngl( j) - sol_sngl( j + z_cells) * &
            A_val_sngl( A_row_ptr( j + 1) - 1)
    end do
    call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_sngl, z_cells, 1, sol_sngl, GS_rhs_sngl, GS_up_diag_sngl)

    do j = 2, y_cells - 1
        start_ind = (j - 1) * z_cells
        do k = 1, z_cells
            GS_rhs_sngl( k) = rhs_sngl( start_ind + k) - sol_sngl( &
                start_ind + k - z_cells) * A_val_sngl( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_sngl( start_ind + k + &
                z_cells) * A_val_sngl( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_sngl, z_cells, j, sol_sngl, GS_rhs_sngl, GS_up_diag_sngl)
    end do

    start_ind = (y_cells - 1) * z_cells
    do j = 1, z_cells
        GS_rhs_sngl( j) = rhs_sngl( start_ind + j) - sol_sngl( start_ind &
            + j - z_cells) * A_val_sngl( A_row_ptr( start_ind + j + 1) - 1)
    end do
    call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_sngl, z_cells, y_cells, sol_sngl, GS_rhs_sngl, &
        GS_up_diag_sngl)
end do

case( 'zebra')

y_cells_is_odd = (mod( y_cells, 2) > 0)

do i = 1, iter
    do j = 1, z_cells
        GS_rhs_sngl( j) = rhs_sngl( j) - sol_sngl( j + z_cells) * &
            A_val_sngl( A_row_ptr( j + 1) - 1)
    end do
    call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
        A_val_sngl, z_cells, 1, sol_sngl, GS_rhs_sngl, GS_up_diag_sngl)

    do j = 3, y_cells - 1, 2
        start_ind = (j - 1) * z_cells
        do k = 1, z_cells
            GS_rhs_sngl( k) = rhs_sngl( start_ind + k) - sol_sngl( &
                start_ind + k - z_cells) * A_val_sngl( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_sngl( start_ind + k + &
                z_cells) * A_val_sngl( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_sngl, z_cells, j, sol_sngl, GS_rhs_sngl, GS_up_diag_sngl)
    end do

    if( y_cells_is_odd) then
        start_ind = (y_cells - 1) * z_cells
        do j = 1, z_cells
            GS_rhs_sngl( j) = rhs_sngl( start_ind + j) - sol_sngl( &
                start_ind + j - z_cells) * A_val_sngl( A_row_ptr( &
                start_ind + j + 1) - 1)
        end do
        call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_sngl, z_cells, y_cells, sol_sngl, GS_rhs_sngl, &
            GS_up_diag_sngl)
    end if

    do j = 2, y_cells - 1, 2
        start_ind = (j - 1) * z_cells
        do k = 1, z_cells
            GS_rhs_sngl( k) = rhs_sngl( start_ind + k) - sol_sngl( &
                start_ind + k - z_cells) * A_val_sngl( A_row_ptr( &
                start_ind + k + 1) - 2) - sol_sngl( start_ind + k + &
                z_cells) * A_val_sngl( A_row_ptr( start_ind + k + 1) - 1)
        end do
        call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_sngl, z_cells, j, sol_sngl, GS_rhs_sngl, GS_up_diag_sngl)
    end do

    if( .not.y_cells_is_odd) then
        start_ind = (y_cells - 1) * z_cells
        do j = 1, z_cells
            GS_rhs_sngl( j) = rhs_sngl( start_ind + j) - sol_sngl( &
                start_ind + j - z_cells) * A_val_sngl( A_row_ptr( &
                start_ind + j + 1) - 1)
        end do
        call tridiag_solv_sngl( y_cells * z_cells, A_row_ptr, A_entries, &
            A_val_sngl, z_cells, y_cells, sol_sngl, GS_rhs_sngl, &
            GS_up_diag_sngl)
    end if
end do

case default
do i = 1, iter
    do j = 1, y_cells * z_cells
        temp = rhs_sngl( j)
        do k = A_row_ptr( j) + 1, A_row_ptr( j + 1) - 1
            temp = temp - A_val_sngl( k) * sol_sngl( A_col_ind( k))
        end do
        sol_sngl( j) = temp / A_val_sngl( A_row_ptr( j))
    end do
end do

end select

end subroutine GS_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a tridiagonal linear subsystem in CRS format usngl the Thomas algorithm.
subroutine tridiag_solv_dble( A_size, A_row_ptr, A_entries, A_val_dble, &
    block_size, bloc, sol_dble, block_rhs, block_up_diag)

! Declares the input variables.
! Declares the number of rows (or columns) of the entire matrix.
integer, intent( in) :: A_size
! Declares the row pointer array of the entire matrix. The matrix entry column index array is not passed in as an input. It is assumed that the entries are given in the order of the diagonal entry, then the lower diagonal (if applicable), and then the upper diagonal entry (if applicable).
integer, intent( in) :: A_row_ptr( A_size + 1)
! Declares the number of non-zero entries of the entire matrix.
integer, intent( in) :: A_entries
! Declares the entry value array of the entire matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the number of rows (or columns) in the block.
integer, intent( in) :: block_size
! Declares the number of the block in the matrix.
integer, intent( in) :: bloc
! Declares the solution of the entire system. Entries for the subsystem will be overwritten.
double precision, intent( inout) :: sol_dble( A_size)
! Declares the right-hand-side for the block.
double precision, intent( in) :: block_rhs( block_size)
! Declares the buffer for computing the modified upper diagonal coefficients. It will be overwritten.
double precision, intent( out) :: block_up_diag( block_size - 1)

! Declares the local variables.
! Declares the current row of the entire matrix.
integer :: row
! Declares the current row of the block.
integer :: block_row
! Declares the common factor for computing the new upper diagonal coefficient and right-hand-side entry in each row.
double precision :: factor

row = (bloc - 1) * block_size + 1

! Starts the forward sweep.
! Modifies the upper diagonal coefficient and the right-hand-side entry in the first row. Writes the right-hand-side entries into the solution vector in preparation for back substitution.
block_up_diag( 1) = A_val_dble( A_row_ptr( row) + 1) / A_val_dble( &
    A_row_ptr( row))
sol_dble( row) = block_rhs( 1) / A_val_dble( A_row_ptr( row))

! Loops through each row until the last one.
row = row + 1
do block_row = 2, block_size - 1
    factor = A_val_dble( A_row_ptr( row)) - A_val_dble( A_row_ptr( row) + &
        1) * block_up_diag( block_row - 1)
    block_up_diag( block_row) = A_val_dble( A_row_ptr( row) + 2) / factor
    sol_dble( row) = (block_rhs( block_row) - A_val_dble( A_row_ptr( row) &
        + 1) * sol_dble( row - 1)) / factor
    row = row + 1
end do

! Modifies the right-hand-side entry in the last row.
factor = A_val_dble( A_row_ptr( row)) - A_val_dble( A_row_ptr( row) + 1) &
    * block_up_diag( block_size - 1)
sol_dble( row) = (block_rhs( block_size) - A_val_dble( A_row_ptr( row) + &
    1) * sol_dble( row - 1)) / factor

! Performs back substitution.
row = row - 1
do block_row = block_size - 1, 1, -1
    sol_dble( row) = sol_dble( row) - block_up_diag( block_row) * &
        sol_dble( row + 1)
    row = row - 1
end do

end subroutine tridiag_solv_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a tridiagonal linear subsystem in snglle precision.
subroutine tridiag_solv_sngl( A_size, A_row_ptr, A_entries, A_val_sngl, &
    block_size, bloc, sol_sngl, block_rhs, block_up_diag)

integer, intent( in) :: A_size
integer, intent( in) :: A_row_ptr( A_size + 1)
integer, intent( in) :: A_entries
real, intent( in) :: A_val_sngl( A_entries)
integer, intent( in) :: block_size
integer, intent( in) :: bloc
real, intent( inout) :: sol_sngl( A_size)
real, intent( in) :: block_rhs( block_size)
real, intent( out) :: block_up_diag( block_size - 1)

integer :: row
integer :: block_row
real :: factor

row = (bloc - 1) * block_size + 1

block_up_diag( 1) = A_val_sngl( A_row_ptr( row) + 1) / A_val_sngl( &
    A_row_ptr( row))
sol_sngl( row) = block_rhs( 1) / A_val_sngl( A_row_ptr( row))

row = row + 1
do block_row = 2, block_size - 1
    factor = A_val_sngl( A_row_ptr( row)) - A_val_sngl( A_row_ptr( row) + &
        1) * block_up_diag( block_row - 1)
    block_up_diag( block_row) = A_val_sngl( A_row_ptr( row) + 2) / factor
    sol_sngl( row) = (block_rhs( block_row) - A_val_sngl( A_row_ptr( row) &
        + 1) * sol_sngl( row - 1)) / factor
    row = row + 1
end do

factor = A_val_sngl( A_row_ptr( row)) - A_val_sngl( A_row_ptr( row) + 1) &
    * block_up_diag( block_size - 1)
sol_sngl( row) = (block_rhs( block_size) - A_val_sngl( A_row_ptr( row) + &
    1) * sol_sngl( row - 1)) / factor

row = row - 1
do block_row = block_size - 1, 1, -1
    sol_sngl( row) = sol_sngl( row) - block_up_diag( block_row) * &
        sol_sngl( row + 1)
    row = row - 1
end do

end subroutine tridiag_solv_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the residual.
subroutine resid_dble( A_rows, A_row_ptr, A_entries, A_col_ind, &
    A_val_dble, A_cols, sol_dble, rhs_dble, res_dble)

! Declares the input variables.
! Declares the number of rows of the matrix.
integer, intent( in) :: A_rows
! Declares the row pointer array of the matrix.
integer, intent( in) :: A_row_ptr( A_rows + 1)
! Declares the number of non-zero entries of the matrix.
integer, intent( in) :: A_entries
! Declares the entry column index array of the matrix.
integer, intent( in) :: A_col_ind( A_entries)
! Declares the entry value array of the matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the number of columns of the matrix.
integer, intent( in) :: A_cols
! Declares the solution.
double precision, intent( in) :: sol_dble( A_cols)
! Declares the right-hand-side.
double precision, intent( in) :: rhs_dble( A_rows)
! Declares the the residual. The array will be overwritten.
double precision, intent( out) :: res_dble( A_rows)

! Declares the local variables.
integer :: i

call mat_vec_mult_crs_dble( A_rows, A_row_ptr, A_entries, A_col_ind, &
    A_val_dble, A_cols, sol_dble, res_dble)
do i = 1, A_rows
    res_dble( i) = rhs_dble( i) - res_dble( i)
end do

end subroutine resid_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the residual in snglle precision.
subroutine resid_sngl( A_rows, A_row_ptr, A_entries, A_col_ind, &
    A_val_sngl, A_cols, sol_sngl, rhs_sngl, res_sngl)

integer, intent( in) :: A_rows
integer, intent( in) :: A_row_ptr( A_rows + 1)
integer, intent( in) :: A_entries
integer, intent( in) :: A_col_ind( A_entries)
real, intent( in) :: A_val_sngl( A_entries)
integer, intent( in) :: A_cols
real, intent( in) :: sol_sngl( A_cols)
real, intent( in) :: rhs_sngl( A_rows)
real, intent( out) :: res_sngl( A_rows)

integer :: i

call mat_vec_mult_crs_sngl( A_rows, A_row_ptr, A_entries, A_col_ind, &
    A_val_sngl, A_cols, sol_sngl, res_sngl)
do i = 1, A_rows
    res_sngl( i) = rhs_sngl( i) - res_sngl( i)
end do

end subroutine resid_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs matrix-vector multiplication in CRS format.
! Adapted from crs_prodMatVec.m in NumGeom/Core/MatrixOps/.
subroutine mat_vec_mult_crs_dble( A_rows, A_row_ptr, A_entries, &
    A_col_ind, A_val_dble, A_cols, x, b)

! Declares the input variables.
! Declares the number of rows of the matrix.
integer, intent( in) :: A_rows
! Declares the row pointer array of the matrix.
integer, intent( in) :: A_row_ptr( A_rows + 1)
! Declares the number of non-zero entries of the matrix.
integer, intent( in) :: A_entries
! Declares the entry column index array of the matrix.
integer, intent( in) :: A_col_ind( A_entries)
! Declares the entry value array of the matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the number of columns of the matrix.
integer, intent( in) :: A_cols
! Declares the vector to be multiplied.
double precision, intent( in) :: x( A_cols)
! Declares the product vector. The array will be overwritten.
double precision, intent( out) :: b( A_rows)

! Declares the local variables.
integer :: i, j
! Declares temporaries to store the row pointers in the register. This allows the routine to read from the row pointer array in main memory only once, instead of twice in each loop.
integer :: ind1
integer :: ind2
! Declares a temporary to store the intermediate sums of the dot products of each row in the register. This avoids writing each intermediate sum to the product vector in main memory.
double precision :: temp

ind1 = A_row_ptr( 1)
do i = 1, A_rows
    temp = 0.0d+0
    ind2 = A_row_ptr( i + 1) - 1
    do j = ind1, ind2
        temp = temp + A_val_dble( j) * x( A_col_ind( j))
    end do
    b( i) = temp
    ind1 = ind2 + 1
end do

end subroutine mat_vec_mult_crs_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs matrix-vector multiplication in single precision.
subroutine mat_vec_mult_crs_sngl( A_rows, A_row_ptr, A_entries, &
    A_col_ind, A_val_sngl, A_cols, x, b)

integer, intent( in) :: A_rows
integer, intent( in) :: A_row_ptr( A_rows + 1)
integer, intent( in) :: A_entries
integer, intent( in) :: A_col_ind( A_entries)
real, intent( in) :: A_val_sngl( A_entries)
integer, intent( in) :: A_cols
real, intent( in) :: x( A_cols)
real, intent( out) :: b( A_rows)

integer :: i, j
integer :: ind1
integer :: ind2
real :: temp

ind1 = A_row_ptr( 1)
do i = 1, A_rows
    temp = 0.0e+0
    ind2 = A_row_ptr( i + 1) - 1
    do j = ind1, ind2
        temp = temp + A_val_sngl( j) * x( A_col_ind( j))
    end do
    b( i) = temp
    ind1 = ind2 + 1
end do

end subroutine mat_vec_mult_crs_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs transposed matrix-vector multiplication in CRS format.
! Adapted from crs_prodMatTransVec.m in NumGeom/Core/MatrixOps/.
subroutine trans_mat_vec_mult_crs( A_rows, A_row_ptr, A_entries, &
    A_col_ind, A_val_dble, A_cols, x, b)

! Declares the input variables.
! Declares the number of rows of the matrix.
integer, intent( in) :: A_rows
! Declares the row pointer array of the matrix.
integer, intent( in) :: A_row_ptr( A_rows + 1)
! Declares the number of non-zero entries of the matrix.
integer, intent( in) :: A_entries
! Declares the entry column index array of the matrix.
integer, intent( in) :: A_col_ind( A_entries)
! Declares the entry value array of the matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the number of columns of the matrix.
integer, intent( in) :: A_cols
! Declares the vector to be multiplied.
double precision, intent( in) :: x( A_cols)
! Declares the product vector. The array will be overwritten.
double precision, intent( out) :: b( A_rows)

! Declares the local variables.
integer :: i, j
! Declares temporaries to store the row pointers in the register. This allows the routine to read from the row pointer array in main memory only once, instead of twice in each loop.
integer :: ind1
integer :: ind2

do i = 1, A_rows
    b( i) = 0.0d+0
end do

ind1 = A_row_ptr( 1)
do i = 1, A_rows - 1
    ind2 = A_row_ptr( i + 1) - 1
    do j = ind1, ind2
        b( A_col_ind( j)) = b( A_col_ind( j)) + A_val_dble( j) * x( i)
    end do
    ind1 = ind2 + 1
end do

end subroutine trans_mat_vec_mult_crs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a dense square linear system by forward and backward substitution given its LU factorization with partial pivoting.
subroutine LU_solv_dble( A_size, A, b)

! Declares the input variables.
! Declares the number of rows (or columns) of the matrix.
integer, intent( in) :: A_size
! Declares the LU factorization of the matrix with partial pivoting, with the strictly lower-triangular part of L and the upper-triangular part of U.
double precision, intent( in) :: A( A_size, A_size)
! Declares the permuted right-hand-side. The array will be overwritten with the solution.
double precision, intent( inout) :: b( A_size)

! Performs forward and back substitution with L and U, respectively.
call forw_solv_dble( A_size, A, b)
call back_solv_dble( A_size, A, b)

end subroutine LU_solv_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a dense square linear system in single precision.
subroutine LU_solv_sngl( A_size, A, b)

integer, intent( in) :: A_size
real, intent( in) :: A( A_size, A_size)
real, intent( inout) :: b( A_size)

call forw_solv_sngl( A_size, A, b)
call back_solv_sngl( A_size, A, b)

end subroutine LU_solv_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a lower-triangular linear system with 1s on the diagonal by forward substitution.
! Adapted from forwardsolve_trans.m in NumGeom/Core/MatrixOps/.
subroutine forw_solv_dble( L_size, L, b)

! Declares the input variables.
! Declares the number of rows (or columns) of the matrix.
integer, intent( in) :: L_size
! Declares the strictly lower-triangular part of L.
double precision, intent( in) :: L( L_size, L_size)
! Declares the right-hand-side. The array will be overwritten with the solution.
double precision, intent( inout) :: b( L_size)

! Declares the local variables.
integer :: i, j

! Loops through the columns.
do i = 1, L_size - 1
    ! Dividing by the diagonal entry to compute the corresponding solution entry is unnecessary since they are all 1s.
    ! Loops through the rows in the inner loop to optimize for cache performance since Fortran stores matrices in column-major order.
    do j = i + 1, L_size
        ! Multiplies the solution entry with the matrix entries in the column and subtracts them from the following solution entries.
        b( j) = b( j) - L( j, i) * b( i)
    end do
end do

end subroutine forw_solv_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a lower-triangular linear system in single precision.
subroutine forw_solv_sngl( L_size, L, b)

integer, intent( in) :: L_size
real, intent( in) :: L( L_size, L_size)
real, intent( inout) :: b( L_size)

integer :: i, j

do i = 1, L_size - 1
    do j = i + 1, L_size
        b( j) = b( j) - L( j, i) * b( i)
    end do
end do

end subroutine forw_solv_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves an upper-triangular linear system by backward substitution.
! Adapted from backsolve.m in NumGeom/Core/MatrixOps/.
subroutine back_solv_dble( U_size, U, b)

! Declares the input variables.
! Declares the number of rows (or columns) of the matrix.
integer, intent( in) :: U_size
! Declares the upper-triangular part of U.
double precision, intent( in) :: U( U_size, U_size)
! Declares the right-hand-side. The array will be overwritten with the solution.
double precision, intent( inout) :: b( U_size)

! Declares the local variables.
integer :: i, j

! Loops through the columns in reverse order.
do i = U_size, 1, -1
    ! Computes the corresponding entry of the solution by dividing by the diagonal entry of the matrix.
    b( i) = b( i) / U( i, i)
    ! Multiplies that solution entry with the matrix entries in the column and subtracts them from the preceding solution entries.
    ! Loops through the rows in the inner loop to optimize for cache performance since Fortran stores matrices in column-major order.
    do j = 1, i - 1
        b( j) = b( j) - U( j, i) * b( i)
    end do
end do

end subroutine back_solv_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves an upper-triangular linear system in single precision.
subroutine back_solv_sngl( U_size, U, b)

integer, intent( in) :: U_size
real, intent( in) :: U( U_size, U_size)
real, intent( inout) :: b( U_size)

integer :: i, j

do i = U_size, 1, -1
    b( i) = b( i) / U( i, i)
    do j = 1, i - 1
        b( j) = b( j) - U( j, i) * b( i)
    end do
end do

end subroutine back_solv_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a dense square linear system using its singular value decomposition.
! Adapted from solve_singular_tsvd.m in NumGeom/Core/MatrixOps/.
subroutine SVD_solv_dble( A_size, U, S, VT, zero_sngl_vals, x, b)

! Declares the input variables.
! Declares the number of rows (or columns) of the matrix.
integer, intent( in) :: A_size
! Declares the SVD of the matrix, with matrix U of left singular vectors, vector S of singular values in decreasing order, and transposed matrix VT of right singular vectors.
double precision, intent( in) :: U( A_size, A_size), S( A_size), VT( &
    A_size, A_size)
! Declares the number of zero singular values, i.e. the number of dimensions to truncate.
integer, intent( in) :: zero_sngl_vals
! Declares the solution. The array will be overwritten.
double precision, intent( out) :: x( A_size)
! Declares the right-hand-side.
double precision, intent( in) :: b( A_size)

! Declares the local variables.
integer :: i, j
! Declares a temporary for applying the SVD.
double precision :: temp

do i = 1, A_size - zero_sngl_vals
    do j = 1, A_size
        temp = U( j, i) * b( j)
    end do
    temp = temp / S( i)
    do j = 1, A_size
        x( j) = x( j) + VT( i, j) * temp
    end do
end do

end subroutine SVD_solv_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solves a dense square linear system using its singular value decomposition in single precision.
subroutine SVD_solv_sngl( A_size, U, S, VT, zero_sngl_vals, x, b)

integer, intent( in) :: A_size
real, intent( in) :: U( A_size, A_size), S( A_size), VT( &
    A_size, A_size)
integer, intent( in) :: zero_sngl_vals
real, intent( out) :: x( A_size)
real, intent( in) :: b( A_size)

integer :: i, j
real :: temp( A_size - zero_sngl_vals)

do i = 1, A_size - zero_sngl_vals
    do j = 1, A_size
        temp( i) = U( j, i) * b( j)
    end do
end do
do i = 1, A_size - zero_sngl_vals
    temp( i) = temp( i) / S( i)
end do
do i = 1, A_size
    do j = 1, A_size - zero_sngl_vals
        x( i) = x( i) + VT( j, i) * temp( j)
    end do
end do

end subroutine SVD_solv_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs matrix-vector multiplication, followed by vector addition, in CRS format.
! Adapted from crs_AXPY.m in NumGeom/Core/MatrixOps/.
subroutine mat_vec_mult_add_crs_dble( A_rows, A_row_ptr, A_entries, &
    A_col_ind, A_val_dble, A_cols, x, y)

! Declares the input variables.
! Declares the number of rows of the matrix.
integer, intent( in) :: A_rows
! Declares the row pointer array of the matrix.
integer, intent( in) :: A_row_ptr( A_rows + 1)
! Declares the number of non-zero entries of the matrix.
integer, intent( in) :: A_entries
! Declares the entry column index array of the matrix.
integer, intent( in) :: A_col_ind( A_entries)
! Declares the entry value array of the matrix.
double precision, intent( in) :: A_val_dble( A_entries)
! Declares the number of columns of the matrix.
integer, intent( in) :: A_cols
! Declares the vector to be multiplied.
double precision, intent( in) :: x( A_cols)
! Declares the vector to be added. The array will be overwritten with the sum.
double precision, intent( inout) :: y( A_rows)

! Declares the local variables.
integer :: i, j
! Declares temporaries to store the row pointers in the register. This allows the routine to read from the row pointer array in main memory only once, instead of twice in each loop.
integer :: ind1
integer :: ind2
! Declares a temporary to store the intermediate sums of the dot products of each row in the register. This avoids writing each intermediate sum to the product vector in main memory.
double precision :: temp

ind1 = A_row_ptr( 1)
do i = 1, A_rows
    temp = y( i)
    ind2 = A_row_ptr( i + 1) - 1
    do j = ind1, ind2
        temp = temp + A_val_dble( j) * x( A_col_ind( j))
    end do
    y( i) = temp
    ind1 = ind2 + 1
end do

end subroutine mat_vec_mult_add_crs_dble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs matrix-vector multiplication, followed by vector addition, in single precision.
subroutine mat_vec_mult_add_crs_sngl( A_rows, A_row_ptr, A_entries, &
    A_col_ind, A_val_sngl, A_cols, x, y)

integer, intent( in) :: A_rows
integer, intent( in) :: A_row_ptr( A_rows + 1)
integer, intent( in) :: A_entries
integer, intent( in) :: A_col_ind( A_entries)
real, intent( in) :: A_val_sngl( A_entries)
integer, intent( in) :: A_cols
real, intent( in) :: x( A_cols)
real, intent( inout) :: y( A_rows)

integer :: i, j
integer :: ind1
integer :: ind2
real :: temp

ind1 = A_row_ptr( 1)
do i = 1, A_rows
    temp = y( i)
    ind2 = A_row_ptr( i + 1) - 1
    do j = ind1, ind2
        temp = temp + A_val_sngl( j) * x( A_col_ind( j))
    end do
    y( i) = temp
    ind1 = ind2 + 1
end do

end subroutine mat_vec_mult_add_crs_sngl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module gmg_common_intern_proced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Defines a module containing the internal subroutines of the pressure solver.
module pressure_gmg_intern_proced
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the number of grids.
function comput_lvls( y_fine_cells, z_fine_cells, coarse_size, &
    semi_coarsen) result( lvls)

! Declares the input variables.
! Declares the number of cells on the fine-grid in the y and the z-directions.
integer, intent( in) :: y_fine_cells, z_fine_cells
! Declares the desired number of cells on one side of the coarse-grid.
integer, intent( in) :: coarse_size
! Declares a flag for whether to semi-coarsen. The routine will semi-coarsen until both sides are between coarse_size and 0.5 * coarse_size + 1. The purpose of semi-coarsening is to reduce the linear system on the coarse-grid to a manageable size. It would otherwise be large, since y_fine_cells >> z_fine_cells.
logical, intent( in) :: semi_coarsen

! Declares the local variables.
! Declares the number of levels required to reduce the number of cells in the y and z-direction to the desired coarse-grid size.
integer :: y_lvls, z_lvls
! Declares a buffer.
integer :: reduc_cells
! Declares the number of grids of multigrid.
integer :: lvls

y_lvls = 1
reduc_cells = y_fine_cells
do while( reduc_cells > coarse_size)
    ! If the number of cells is even, divides by 2.
    if( mod( reduc_cells, 2) < 1) then
        reduc_cells = reduc_cells / 2
    ! Otherwise, divides by 2 (Fortran rounds toward 0) and adds 1.
    else
        reduc_cells = reduc_cells / 2 + 1
    end if
    y_lvls = y_lvls + 1
end do
! Does the same in the z-direction.
z_lvls = 1
reduc_cells = z_fine_cells
do while( reduc_cells > coarse_size)
    if( mod( reduc_cells, 2) < 1) then
        reduc_cells = reduc_cells / 2
    else
        reduc_cells = reduc_cells / 2 + 1
    end if
    z_lvls = z_lvls + 1
end do
! Takes the maximum of the two (uses semi-coarsening).
if( semi_coarsen) then
    lvls = max( y_lvls, z_lvls)
! Takes the minimum of the two (does not use semi-coarsening).
else
    lvls = min( y_lvls, z_lvls)
end if

end function comput_lvls

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the PDE discretization matrix on one grid. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine discr( lvl, y_cell_edge_coords, p_cell_edge, &
    y_cell_widths, z_cell_widths, y_cell_cent_coords, mu_discr, &
    rho_discr, mu_is_const, rad_reduc_factor)
use pressure_gmg_intern_data, only: multigrid
use gmg_common_intern_proced, only: diff_coeff

! Declares the input variables.
! Declares the level of the grid on which to rediscretize.
integer, intent( in) :: lvl
! Declares the y-coordinates of the edges of the cells.
double precision, intent( in) :: y_cell_edge_coords( multigrid( 1) % &
    cells( 1) + 1)
! Declares the background pressures at the edges of the cells in the z-direction.
double precision, intent( in) :: p_cell_edge( multigrid( 1) % cells( 2) + &
    1)
! Declares the widths of the cells in the y and z-directions.
double precision, intent( in) :: y_cell_widths( multigrid( 1) % cells( &
    1)), z_cell_widths( multigrid( 1) % cells( 2))
! Declares the y-coordinates of the centers of the cells.
double precision, intent( in) :: y_cell_cent_coords( multigrid( 1) % &
    cells( 1))
character( len = 1), intent( in) :: mu_discr
character( len = 1), intent( in) :: rho_discr
! Declares a flag for the coefficient mu. If .false., mu is computed from the standard formula. If .true., all mu are set to 1 in order to compare the GMG solution to the FFT solution.
logical, intent( in) :: mu_is_const
! Declares the reduction factor for the radius of the earth, e.g. if rad_reduc_factor is 2, the radius is halved.
double precision, intent( in) :: rad_reduc_factor

! Declares the local variables.
integer :: i, j
! Declares the diffusion coefficients at the centers and edges of the cells in the y-direction.
double precision, allocatable :: mu_cent( :), mu_edge( :)
! Declares the air densities at the centers of the cells in the z-direction.
double precision, allocatable :: rho_cent( :)
! Declares the cell index of the center cell.
integer :: cent_ind
! Declares the matrix entry index.
integer :: ind
! Declares the entry value for the center cell.
double precision :: cent_val
! Declares the entry index for the center cell.
integer :: val_cent_ind

! Allocates an array for the diffusion coefficients at the centers of the cells and computes them.
allocate( mu_cent( multigrid( lvl) % cells( 1)))
do i = 1, multigrid( lvl) % cells( 1)
    mu_cent( i) = y_cell_cent_coords( i)
end do
call diff_coeff( multigrid( lvl) % cells( 1), mu_cent, mu_is_const, &
    rad_reduc_factor)

! Allocates an array for the air densities at the centers of the cells and computes them.
allocate( rho_cent( multigrid( lvl) % cells( 2)))
do i = 1, multigrid( lvl) % cells( 2)
    rho_cent( i) = -( p_cell_edge( i + 1) - p_cell_edge( i)) / &
        z_cell_widths( i)
end do

! Allocates an array for the diffusion coefficients at the edges of the cells and computes them.
allocate( mu_edge( multigrid( lvl) % cells( 1) + 1))
do i = 1, multigrid( lvl) % cells( 1) + 1
    mu_edge( i) = y_cell_edge_coords( i)
end do
call diff_coeff( multigrid( lvl) % cells( 1) + 1, mu_edge, mu_is_const, &
    rad_reduc_factor)

! Allocates the row pointer array. Recall that the cells on each grid are numbered in the direction of increasing coordinates, in z-major order.
allocate( multigrid( lvl) % A_row_ptr( multigrid( lvl) % cells( 3) + 1))
do i = 1, multigrid( lvl) % cells( 3) + 1
    multigrid( lvl) % A_row_ptr( i) = 1
end do
cent_ind = 1
! Loops through the cells.
do i = 1, multigrid( lvl) % cells( 1)
    do j = 1, multigrid( lvl) % cells( 2)
        ! If the cell is not adjacent to the left boundary, increments the row pointer by 1 for the cell to the left.
        if( i > 1) then
            multigrid( lvl) % A_row_ptr( cent_ind + 1) = multigrid( lvl) &
                % A_row_ptr( cent_ind + 1) + 1
        end if
        ! If the cell is not adjacent to the right boundary, increments the row pointer by 1 for the cell to the right.
        if( i < multigrid( lvl) % cells( 1)) then
            multigrid( lvl) % A_row_ptr( cent_ind + 1) = multigrid( lvl) &
                % A_row_ptr( cent_ind + 1) + 1
        end if
        ! If the cell is not adjacent to the lower boundary, increments the row pointer by 1 for the cell below.
        if( j > 1) then
            multigrid( lvl) % A_row_ptr( cent_ind + 1) = multigrid( lvl) &
                % A_row_ptr( cent_ind + 1) + 1
        end if
        ! If the cell is not adjacent to the upper boundary, increments the row pointer by 1 for the cell above.
        if( j < multigrid( lvl) % cells( 2)) then
            multigrid( lvl) % A_row_ptr( cent_ind + 1) = multigrid( lvl) &
                % A_row_ptr( cent_ind + 1) + 1
        end if
        cent_ind = cent_ind + 1
    end do
end do
! Adds each row pointer to the next one, in succession.
do i = 1, multigrid( lvl) % cells( 3)
    multigrid( lvl) % A_row_ptr( i + 1) = multigrid( lvl) % A_row_ptr( i &
        + 1) + multigrid( lvl) % A_row_ptr( i)
end do

! Allocates the arrays for the matrix.
allocate( multigrid( lvl) % A_col_ind( multigrid( lvl) % A_row_ptr( &
    multigrid( lvl) % cells( 3) + 1) - 1))
allocate( multigrid( lvl) % A_val_dble( multigrid( lvl) % A_row_ptr( &
    multigrid( lvl) % cells( 3) + 1) - 1), multigrid( lvl) % A_val_sngl( &
    multigrid( lvl) % A_row_ptr( multigrid( lvl) % cells( 3) + 1) - 1))
allocate( multigrid( lvl) % A_diag_val_dble( multigrid( lvl) % cells( &
    3)), multigrid( lvl) % A_diag_val_sngl( multigrid( lvl) % cells( 3)))

ind = 1
cent_ind = 1

do i = 1, multigrid( lvl) % cells( 1)
    do j = 1, multigrid( lvl) % cells( 2)
        ! Sets the column index for the center cell. Its entry value will be set at the end.
        multigrid( lvl) % A_col_ind( ind) = cent_ind
        cent_val = 0.0d+0
        val_cent_ind = ind
        select case( rho_discr)
        case( 'B')
            ! Sets the column index and entry value for the cell below.
            if( j > 1) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind - 1
                multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
                    mu_cent( i) * (rho_cent( j - 1) * z_cell_widths( j) + rho_cent( j) * z_cell_widths( j - 1)) * 2.0d+0 / &
                    ((z_cell_widths( j - 1) + z_cell_widths( j)) * &
                    z_cell_widths( j) * rho_cent( j) * (z_cell_widths( j &
                    - 1) + z_cell_widths( j)))
                ! Subtracts the entry value for the cell below from the entry value for the center cell.
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
            ! Sets the column index and entry value for the cell above.
            if( j < multigrid( lvl) % cells( 2)) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind + 1
                multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
                    mu_cent( i) * (rho_cent( j) * z_cell_widths( j + 1) + &
                    rho_cent( j + 1) * z_cell_widths( j)) * 2.0d+0 / &
                    ((z_cell_widths( j) + z_cell_widths( j + 1)) * &
                    z_cell_widths( j) * rho_cent( j) * (z_cell_widths( j) &
                    + z_cell_widths( j + 1)))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
        ! case( 'C')
        !     if( j == 1) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind + 1
        !         multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
        !             mu_cent( i) * (rho_cent( j) + rho_cent( j + 1)) * &
        !             2.0d+0 / ((z_cell_widths( j) + z_cell_widths( j + 1)) &
        !             * (z_cell_widths( j) + z_cell_widths( j + 1)) * &
        !             rho_cent( j))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( (j > 1) .and. (j < multigrid( lvl) % cells( 2))) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind - 1
        !         multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
        !             mu_cent( i) * (rho_cent( j - 1) + rho_cent( j)) * &
        !             4.0d+0 / ((z_cell_widths( j - 1) + 2.0d+0 * &
        !             z_cell_widths( j) + z_cell_widths( j + 1)) * &
        !             (z_cell_widths( j - 1) + z_cell_widths( j)) * &
        !             rho_cent( j))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( (j > 1) .and. (j < multigrid( lvl) % cells( 2))) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind + 1
        !         multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
        !             mu_cent( i) * (rho_cent( j) + rho_cent( j + 1)) * &
        !             4.0d+0 / ((z_cell_widths( j - 1) + 2.0d+0 * &
        !             z_cell_widths( j) + z_cell_widths( j + 1)) * &
        !             (z_cell_widths( j) + z_cell_widths( j + 1)) * &
        !             rho_cent( j))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( j == multigrid( lvl) % cells( 2)) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind - 1
        !         multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
        !             mu_cent( i) * (rho_cent( j - 1) + rho_cent( j)) * &
        !             2.0d+0 / ((z_cell_widths( j - 1) + z_cell_widths( j)) &
        !             * (z_cell_widths( j - 1) + z_cell_widths( j)) * &
        !             rho_cent( j))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        case default
            if( j > 1) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind - 1
                multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
                    mu_cent( i) * (rho_cent( j - 1) + rho_cent( j)) / &
                    ((z_cell_widths( j - 1) + z_cell_widths( j)) * &
                    z_cell_widths( j) * rho_cent( j))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
            if( j < multigrid( lvl) % cells( 2)) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind + 1
                multigrid( lvl) % A_val_dble( ind) = mu_cent( i) * &
                    mu_cent( i) * (rho_cent( j) + rho_cent( j + 1)) / &
                    ((z_cell_widths( j) + z_cell_widths( j + 1)) * &
                    z_cell_widths( j) * rho_cent( j))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
        end select
        select case( mu_discr)
        case( 'B')
            ! Sets the column index and entry value for the cell to the left.
            if( i > 1) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind - multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = (mu_cent( i - 1) * &
                    y_cell_widths( i) + mu_cent( i) * y_cell_widths( i - &
                    1)) * mu_cent( i) * 2.0d+0 / ((y_cell_widths( i - 1) &
                    + y_cell_widths( i)) * (y_cell_widths( i - 1) + &
                    y_cell_widths( i)) * y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
            ! Sets the column index and entry value for the cell to the right.
            if( i < multigrid( lvl) % cells( 1)) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind + multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = (mu_cent( i) * &
                    y_cell_widths( i + 1) + mu_cent( i + 1) * &
                    y_cell_widths( i)) * mu_cent( i) * 2.0d+0 / &
                    ((y_cell_widths( i) + y_cell_widths( i + 1)) * &
                    (y_cell_widths( i) + y_cell_widths( i + 1)) * &
                    y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
        ! case( 'C')
        !     if( i == 1) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind + multigrid( &
        !             lvl) % cells( 2)
        !         multigrid( lvl) % A_val_dble( ind) = (mu_cent( i) + &
        !             mu_cent( i + 1)) * mu_cent( i) * 2.0d+0 / &
        !             ((y_cell_widths( i) + y_cell_widths( i + 1)) * &
        !             (y_cell_widths( i) + y_cell_widths( i + 1)))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( (i > 1) .and. (i < multigrid( lvl) % cells( 1))) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind - multigrid( &
        !             lvl) % cells( 2)
        !         multigrid( lvl) % A_val_dble( ind) = (mu_cent( i - 1) + &
        !             mu_cent( i)) * mu_cent( i) * 4.0d+0 / &
        !             ((y_cell_widths( i - 1) + 2.0d+0 * y_cell_widths( i) &
        !             + y_cell_widths( i + 1)) * (y_cell_widths( i - 1) + &
        !             y_cell_widths( i)))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( (i > 1) .and. (i < multigrid( lvl) % cells( 1))) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind + multigrid( &
        !             lvl) % cells( 2)
        !         multigrid( lvl) % A_val_dble( ind) = (mu_cent( i) + &
        !             mu_cent( i + 1)) * mu_cent( i) * 4.0d+0 / &
        !             ((y_cell_widths( i - 1) + 2.0d+0 * y_cell_widths( i) &
        !             + y_cell_widths( i + 1)) * (y_cell_widths( i) + &
        !             y_cell_widths( i + 1)))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        !     if( i == multigrid( lvl) % cells( 1)) then
        !         ind = ind + 1
        !         multigrid( lvl) % A_col_ind( ind) = cent_ind - multigrid( &
        !             lvl) % cells( 2)
        !         multigrid( lvl) % A_val_dble( ind) = (mu_cent( i - 1) + &
        !             mu_cent( i)) * mu_cent( i) * 2.0d+0 / &
        !             ((y_cell_widths( i - 1) + y_cell_widths( i)) * &
        !             (y_cell_widths( i - 1) + y_cell_widths( i)))
        !         cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
        !     end if
        case( 'C')
            if( i > 1) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind - multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = mu_edge( i) * &
                    mu_cent( i) * 2.0d+0 / ((y_cell_widths( i - 1) + &
                    y_cell_widths( i)) * y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
            if( i < multigrid( lvl) % cells( 1)) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind + multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = mu_edge( i + 1) * &
                    mu_cent( i) * 2.0d+0 / ((y_cell_widths( i) + &
                    y_cell_widths( i + 1)) * y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
        case default
            if( i > 1) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind - multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = (mu_cent( i - 1) + &
                    mu_cent( i)) * mu_cent( i) / ((y_cell_widths( i - 1) &
                    + y_cell_widths( i)) * y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
            if( i < multigrid( lvl) % cells( 1)) then
                ind = ind + 1
                multigrid( lvl) % A_col_ind( ind) = cent_ind + multigrid( &
                    lvl) % cells( 2)
                multigrid( lvl) % A_val_dble( ind) = (mu_cent( i) + &
                    mu_cent( i + 1)) * mu_cent( i) / ((y_cell_widths( i) &
                    + y_cell_widths( i + 1)) * y_cell_widths( i))
                cent_val = cent_val - multigrid( lvl) % A_val_dble( ind)
            end if
        end select
        ! Sets the entry value for the center cell.
        multigrid( lvl) % A_val_dble( val_cent_ind) = cent_val
        multigrid( lvl) % A_diag_val_dble( cent_ind) = cent_val
        ind = ind + 1
        cent_ind = cent_ind + 1
    end do
end do

! Typecasts the double arrays into the single ones.
do i = 1, multigrid( lvl) % A_row_ptr( multigrid( lvl) % cells( 3) + 1) - 1
    multigrid( lvl) % A_val_sngl( i) = sngl( multigrid( lvl) % &
        A_val_dble( i))
end do
do i = 1, multigrid( lvl) % cells( 3)
    multigrid( lvl) % A_diag_val_sngl( i) = sngl( multigrid( lvl) % &
        A_diag_val_dble( i))
end do

! Deallocates the arrays for the air densities and diffusion coefficients, since they will not be needed outside the routine.
deallocate( mu_cent, rho_cent, mu_edge)

end subroutine discr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the data associated with one coarser grid.
subroutine coarsen( finer_lvl, y_cell_edge_coords, &
    z_cell_edge_coords, p_cell_edge, y_cell_widths, z_cell_widths, &
    y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
    y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, coarse_size, &
    restr_prec, prol_prec)
use pressure_gmg_intern_data, only: multigrid

! Declares the input variables.
! Declares the level of the old (finer) grid.
integer, intent( in) :: finer_lvl
! Declares the y and z-coordinates of the edges of the old (finer) grid cells. The arrays will be overwritten with the values on the new (coarser) grid.
double precision, intent( inout) :: y_cell_edge_coords( multigrid( 1) % &
    cells( 1) + 1), z_cell_edge_coords( multigrid( 1) % cells( 2) + 1)
! Declares the background pressures at the edges of the finer grid cells in the z-direction. The arrays will be overwritten with the values on the coarser grid.
double precision, intent( inout) :: p_cell_edge( multigrid( 1) % cells( &
    2) + 1)
! Declares the widths of the finer grid cells in the y and z-directions. The arrays will be overwritten with the values on the coarser grid.
double precision, intent( inout) :: y_cell_widths( multigrid( 1) % &
    cells( 1)), z_cell_widths( multigrid( 1) % cells( 2))
! Declares the y and z-coordinates of the centers of the finer and coarser grid cells. The arrays will be overwritten with the values on the new grids.
double precision, intent( inout) :: y_finer_cell_cent_coords( multigrid( &
    1) % cells( 1)), z_finer_cell_cent_coords( multigrid( 1) % cells( &
    2)), y_coarser_cell_cent_coords( multigrid( 1) % cells( 1)), &
    z_coarser_cell_cent_coords( multigrid( 1) % cells( 2))
! Declares the desired number of cells on one side of the coarse-grid. The routine will semi-coarsen until both sides are between coarse_size and 0.5 * coarse_size + 1.
integer, intent( in) :: coarse_size
integer, intent( in) :: restr_prec, prol_prec

! Declares the local variables.
integer :: i
! Declares the coarser level.
integer :: coarser_lvl
! Declares flags for whether the number of finer grid cells in the y and z-directions are odd.
logical :: y_cells_is_odd, z_cells_is_odd
! Declares flags for whether the desired coarse-grid size has already been reached on the finer grid in the y and z-directions.
logical :: y_grid_is_coarse, z_grid_is_coarse

coarser_lvl = finer_lvl + 1

y_cells_is_odd = (mod( multigrid( finer_lvl) % cells( 1), 2) > 0)
z_cells_is_odd = (mod( multigrid( finer_lvl) % cells( 2), 2) > 0)
y_grid_is_coarse = (multigrid( finer_lvl) % cells( 1) < coarse_size + 1)
z_grid_is_coarse = (multigrid( finer_lvl) % cells( 2) < coarse_size + 1)

! Allocates the array for the number of cells on the coarser grid and computes them.
allocate( multigrid( coarser_lvl) % cells( 3))

! If the number of finer grid cells in the y-direction is even and has not reached the desired coarse-grid size, merges the cells in pairs.
if( .not.y_cells_is_odd .and. .not.y_grid_is_coarse) then
    multigrid( coarser_lvl) % cells( 1) = multigrid( finer_lvl) % cells( &
        1) / 2
    do i = 1, multigrid( coarser_lvl) % cells( 1)
        y_cell_widths( i) = y_cell_widths( 2 * i - 1) + y_cell_widths( 2 &
        * i)
    end do
! If the number of cells is odd and has not reached the desired coarse-grid size, merges the cells in pairs and keeps the final cell.
else if( .not.y_grid_is_coarse) then
    ! Divides by 2 (Fortran rounds toward 0) and adds 1.
    multigrid( coarser_lvl) % cells( 1) = multigrid( finer_lvl) % cells( &
        1) / 2 + 1
    do i = 1, multigrid( coarser_lvl) % cells( 1) - 1
        y_cell_widths( i) = y_cell_widths( 2 * i - 1) + y_cell_widths( 2 &
        * i)
    end do
    y_cell_widths( multigrid( coarser_lvl) % cells( 1)) = y_cell_widths( &
        multigrid( finer_lvl) % cells( 1))
! If the desired coarse-grid size has been reached, keeps the finer grid cells.
else
    multigrid( coarser_lvl) % cells( 1) = multigrid( finer_lvl) % cells( 1)
end if

! Does the same in the z-direction.
if( .not.z_cells_is_odd .and. .not.z_grid_is_coarse) then
    multigrid( coarser_lvl) % cells( 2) = multigrid( finer_lvl) % cells( &
        2) / 2
    do i = 1, multigrid( coarser_lvl) % cells( 2)
        z_cell_widths( i) = z_cell_widths( 2 * i - 1) + z_cell_widths( 2 &
            * i)
    end do
else if( .not.z_grid_is_coarse) then
    multigrid( coarser_lvl) % cells( 2) = multigrid( finer_lvl) % cells( &
        2) / 2 + 1
    do i = 1, multigrid( coarser_lvl) % cells( 2) - 1
        z_cell_widths( i) = z_cell_widths( 2 * i - 1) + z_cell_widths( 2 &
            * i)
    end do
    z_cell_widths( multigrid( coarser_lvl) % cells( 2)) = z_cell_widths( &
        multigrid( finer_lvl) % cells( 2))
else
    multigrid( coarser_lvl) % cells( 2) = multigrid( finer_lvl) % cells( 2)
end if

multigrid( coarser_lvl) % cells( 3) = multigrid( coarser_lvl) % cells( 1) &
    * multigrid( coarser_lvl) % cells( 2)

! Computes the coordinates of the edges of the coarser grid cells.
do i = 1, multigrid( coarser_lvl) % cells( 1)
    y_cell_edge_coords( i + 1) = y_cell_edge_coords( i) + y_cell_widths( i)
end do
do i = 1, multigrid( coarser_lvl) % cells( 2)
    z_cell_edge_coords( i + 1) = z_cell_edge_coords( i) + z_cell_widths( i)
end do

! Computes the background pressures at the edges of the coarser grid cells.
if( multigrid( finer_lvl) % cells( 2) > coarse_size) then
    do i = 1, multigrid( coarser_lvl) % cells( 2)
        p_cell_edge( i) = p_cell_edge( 2 * i - 1)
    end do
    p_cell_edge( multigrid( coarser_lvl) % cells( 2) + 1) = p_cell_edge( &
        multigrid( finer_lvl) % cells( 2) + 1)
end if

! Copies the coordinates of the centers of the finer grid cells.
do i = 1, multigrid( finer_lvl) % cells( 1)
    y_finer_cell_cent_coords( i) = y_coarser_cell_cent_coords( i)
end do
do i = 1, multigrid( finer_lvl) % cells( 2)
    z_finer_cell_cent_coords( i) = z_coarser_cell_cent_coords( i)
end do
! Computes the coordinates of the centers of the coarser grid cells.
do i = 1, multigrid( coarser_lvl) % cells( 1)
    y_coarser_cell_cent_coords( i) = 0.5d+0 * (y_cell_edge_coords( i) + &
        y_cell_edge_coords( i + 1))
end do
do i = 1, multigrid( coarser_lvl) % cells( 2)
    z_coarser_cell_cent_coords( i) = 0.5d+0 * (z_cell_edge_coords( i) + &
        z_cell_edge_coords( i + 1))
end do

! Computes the restriction operator.
if( .not.z_grid_is_coarse) then
    call restr( finer_lvl, y_finer_cell_cent_coords, &
        z_finer_cell_cent_coords, y_coarser_cell_cent_coords, &
        z_coarser_cell_cent_coords, restr_prec)
! If the finer grid is already coarse in the z-direction, semi-coarsens in the y-direction.
else
    call semi_restr_y( finer_lvl, y_finer_cell_cent_coords, &
        y_coarser_cell_cent_coords, restr_prec)
end if

! Computes the prolongation operator.
if( .not.z_grid_is_coarse) then
    call prol( coarser_lvl, y_cell_widths, z_cell_widths, &
        y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
        y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, prol_prec)
! If the finer grid is already coarse in the z-direction, semi-prolongates in the y-direction.
else
    call semi_prol_y( coarser_lvl, y_cell_widths, &
        y_finer_cell_cent_coords, y_coarser_cell_cent_coords, prol_prec)
end if

end subroutine coarsen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the restriction operator. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine restr( finer_lvl, y_finer_cell_cent_coords, &
    z_finer_cell_cent_coords, y_coarser_cell_cent_coords, &
    z_coarser_cell_cent_coords, restr_prec)
use pressure_gmg_intern_data, only: multigrid

! Declares the input variables.
! Declares the level of the finer grid, from which to restrict functions.
integer, intent( in) :: finer_lvl
! Declares the y and z-coordinates of the centers of the finer and coarser grid cells.
double precision, intent( in) :: y_finer_cell_cent_coords( multigrid( &
    finer_lvl) % cells( 1)), z_finer_cell_cent_coords( multigrid( &
    finer_lvl) % cells( 2)), y_coarser_cell_cent_coords( multigrid( &
    finer_lvl + 1) % cells( 1)), z_coarser_cell_cent_coords( multigrid( &
    finer_lvl + 1) % cells( 2))
integer, intent( in) :: restr_prec

! Declares the local variables.
integer :: i
! Declares the level of the coarser grid, to which to restrict functions.
integer :: coarser_lvl
logical :: y_finer_cells_is_odd, z_finer_cells_is_odd
integer :: cell_ind
integer :: y_finer_ind, z_finer_ind
integer :: y_coarser_ind, z_coarser_ind
integer :: entry_ind
double precision :: factor

coarser_lvl = finer_lvl + 1

y_finer_cells_is_odd = (mod( multigrid( finer_lvl) % cells( 1), 2) > 0)
z_finer_cells_is_odd = (mod( multigrid( finer_lvl) % cells( 2), 2) > 0)

! Allocates the row pointer array.
allocate( multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) % &
    cells( 3) + 1))
! Sets the first pointer to 1 (the first entry of the matrix).
multigrid( finer_lvl) % R_row_ptr( 1) = 1
! Initializes all the following entries to 4 since all cells on the coarser grid not adjacent to the right or upper boundary have values interpolated from the 4 cells on the finer grid from which they are formed. The exceptions will be addressed below. Note that the pointer for the following row is incremented.
do i = 2, multigrid( coarser_lvl) % cells( 3) + 1
    multigrid( finer_lvl) % R_row_ptr( i) = 4
end do
! If the number of cells in the y-direction on the finer grid is odd, then the cells on the coarser grid adjacent to the right boundary are formed from only 2 cells on the finer grid. The cell at the upper right corner will be addressed below.
cell_ind = (multigrid( coarser_lvl) % cells( 1) - 1) * multigrid( &
    coarser_lvl) % cells( 2) + 1
if( y_finer_cells_is_odd) then
    do i = 1, multigrid( coarser_lvl) % cells( 2) - 1
        multigrid( finer_lvl) % R_row_ptr( cell_ind + i) = 2
    end do
end if
! Similarly, if the number of cells in the z-direction on the finer grid is odd, then the cells on the coarser grid adjacent to the upper boundary are formed from only 2 cells on the finer grid. The cell at the upper right corner will be addressed next.
cell_ind = 1
if( z_finer_cells_is_odd) then
    do i = 1, multigrid( coarser_lvl) % cells( 1) - 1
        cell_ind = cell_ind + multigrid( coarser_lvl) % cells( 2)
        multigrid( finer_lvl) % R_row_ptr( cell_ind) = 2
    end do
end if
! If the number of cells in both the y and z-directions on the finer grid are odd, then the cell on the coarser grid at the upper right corner is formed from only 1 cell on the finer grid.
if( y_finer_cells_is_odd .and. z_finer_cells_is_odd) then
    multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) % cells( &
        3) + 1) = 1
! If the number of cells in only the y or the z-direction on the finer grid is odd, then the cell on the coarser grid at the upper right corner is formed from 2 cells on the finer grid.
elseif( y_finer_cells_is_odd .or. z_finer_cells_is_odd) then
    multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) % cells( &
        3) + 1) = 2
end if
! Adds each row pointer to the next one, in succession.
do i = 1, multigrid( coarser_lvl) % cells( 3)
    multigrid( finer_lvl) % R_row_ptr( i + 1) = multigrid( finer_lvl) % &
        R_row_ptr( i) + multigrid( finer_lvl) % R_row_ptr( i + 1)
end do

! Allocates the entry column index and value arrays.
allocate( multigrid( finer_lvl) % R_col_ind( multigrid( finer_lvl) % &
    R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
allocate( multigrid( finer_lvl) % R_val_dble( multigrid( finer_lvl) % &
    R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
if( restr_prec == 1) then
    allocate( multigrid( finer_lvl) % R_val_sngl( multigrid( finer_lvl) &
        % R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
end if

! Uses bilinear interpolation for the cells on the coarser grid not adjacent to the right or upper boundary.
do y_coarser_ind = 1, multigrid( coarser_lvl) % cells( 1) - 1
    ! Computes the upper index in the y-direction of the four cells on the finer grid from which the current cell on the coarser grid are formed.
    y_finer_ind = 2 * y_coarser_ind
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2) - 1
        ! Computes the entry index for the lower left cell on the finer grid.
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (y_coarser_ind - &
            1) * multigrid( coarser_lvl) % cells( 2) + z_coarser_ind)
        ! Computes the upper index in the z-direction.
        z_finer_ind = 2 * z_coarser_ind
        ! Computes a common factor.
        factor = 1 / ((y_finer_cell_cent_coords( y_finer_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)))
        ! Sets the entry column index and value for the lower left cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 2) &
            * multigrid( finer_lvl) % cells( 2) + z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        ! Sets the entry column index and value for the upper left cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (y_finer_ind &
            - 2) * multigrid( finer_lvl) % cells( 2) + z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        ! Sets the entry column index and value for the lower right cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind + 2) = (y_finer_ind &
            - 1) * multigrid( finer_lvl) % cells( 2) + z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind + 2) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
        ! Sets the entry column index and value for the upper right cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind + 3) = (y_finer_ind &
            - 1) * multigrid( finer_lvl) % cells( 2) + z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 3) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
    end do
end do

! If the number of cells in the y-direction on the finer grid is odd, then uses linear interpolation along the z-direction for the cells on the coarser grid adjacent to the right boundary. The cell at the upper right corner will be addressed below.
if( y_finer_cells_is_odd) then
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2) - 1
        z_finer_ind = 2 * z_coarser_ind
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (multigrid( &
            coarser_lvl) % cells( 1) - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_coarser_ind)
        factor = (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) / &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1))
        ! Sets the entry column index and value for the lower cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (multigrid( &
            finer_lvl) % cells( 1) - 1) * multigrid( finer_lvl) % cells( &
            2) + z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
        ! Sets the entry column index and value for the upper cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (multigrid( &
            finer_lvl) % cells( 1) - 1) * multigrid( finer_lvl) % cells( &
            2) + z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
    end do
! Otherwise, uses bilinear interpolation as before.
else
    y_coarser_ind = multigrid( coarser_lvl) % cells( 1)
    y_finer_ind = 2 * y_coarser_ind
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2) - 1
        z_finer_ind = 2 * z_coarser_ind
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (y_coarser_ind - &
            1) * multigrid( coarser_lvl) % cells( 2) + z_coarser_ind)
        factor = 1.0d+0 / ((z_finer_cell_cent_coords( z_finer_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)))
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 2) &
            * multigrid( finer_lvl) % cells( 2) + z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (y_finer_ind &
            - 2) * multigrid( finer_lvl) % cells( 2) + z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 2) = (y_finer_ind &
            - 1) * multigrid( finer_lvl) % cells( 2) + z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind + 2) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 3) = (y_finer_ind &
            - 1) * multigrid( finer_lvl) % cells( 2) + z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 3) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
    end do
end if

! If the number of cells in the z-direction on the finer grid is odd, then uses linear interpolation along the y-direction for the cells on the coarser grid adjacent to the upper boundary. The cell at the upper right corner will be addressed next.
if( z_finer_cells_is_odd) then
    do y_coarser_ind = 1, multigrid( coarser_lvl) % cells( 1) - 1
        y_finer_ind = 2 * y_coarser_ind
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (y_coarser_ind - &
            1) * multigrid( coarser_lvl) % cells( 2) + multigrid( &
            coarser_lvl) % cells( 2))
        factor = (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) / &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1))
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 1) &
            * multigrid( finer_lvl) % cells( 2)
        multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = y_finer_ind * &
            multigrid( finer_lvl) % cells( 2)
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
    end do
! Otherwise, uses bilinear interpolation as before.
else
    z_coarser_ind = multigrid( coarser_lvl) % cells( 2)
    z_finer_ind = z_coarser_ind * 2
    do y_coarser_ind = 1, multigrid( coarser_lvl) % cells( 1) - 1
        y_finer_ind = y_coarser_ind * 2
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (y_coarser_ind - 1) * multigrid( coarser_lvl) % cells( 2) + z_coarser_ind)
        factor = 1.0d+0 / ((z_finer_cell_cent_coords( z_finer_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)))
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - &
            1) * z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (y_finer_ind &
            - 1) * z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_finer_cell_cent_coords( y_finer_ind) - &
            y_coarser_cell_cent_coords( y_coarser_ind)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 2) = y_finer_ind * &
            z_finer_ind - 1
        multigrid( finer_lvl) % R_val_dble( entry_ind + 2) = &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 3) = y_finer_ind * &
            z_finer_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 3) = &
            (z_coarser_cell_cent_coords( z_coarser_ind) - &
            z_finer_cell_cent_coords( z_finer_ind - 1)) * &
            (y_coarser_cell_cent_coords( y_coarser_ind) - &
            y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
    end do
end if

! If the number of cells in both the y and z-directions on the finer grid are odd, then uses injection for the cells on the coarser grid at the upper right corner.
if( y_finer_cells_is_odd .and. z_finer_cells_is_odd) then
    entry_ind = multigrid( finer_lvl) % R_row_ptr( multigrid( &
        coarser_lvl) % cells( 3))
    multigrid( finer_lvl) % R_col_ind( entry_ind) = multigrid( finer_lvl) &
        % cells( 1) * multigrid( finer_lvl) % cells( 2)
    multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0
! If the number of cells in only the y or z-direction on the finer grid is odd, then uses linear interpolation along the z or y-direction, respectively.
elseif( y_finer_cells_is_odd .and. .not.z_finer_cells_is_odd) then
    z_coarser_ind = multigrid( coarser_lvl) % cells( 2)
    z_finer_ind = 2 * z_coarser_ind
    entry_ind = multigrid( finer_lvl) % R_row_ptr( multigrid( &
        coarser_lvl) % cells( 3))
    factor = (z_coarser_cell_cent_coords( z_coarser_ind) - &
        z_finer_cell_cent_coords( z_finer_ind - 1)) / &
        (z_finer_cell_cent_coords( z_finer_ind) - &
        z_finer_cell_cent_coords( z_finer_ind - 1))
    multigrid( finer_lvl) % R_col_ind( entry_ind) = multigrid( finer_lvl) &
        % cells( 3) - 1
    multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
    multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = multigrid( &
        finer_lvl) % cells( 3)
    multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
elseif( .not.y_finer_cells_is_odd .and. z_finer_cells_is_odd) then
    y_coarser_ind = multigrid( coarser_lvl) % cells( 1)
    y_finer_ind = 2 * y_coarser_ind
    entry_ind = multigrid( finer_lvl) % R_row_ptr( multigrid( &
        coarser_lvl) % cells( 3))
    factor = (y_coarser_cell_cent_coords( y_coarser_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1)) / &
        (y_finer_cell_cent_coords( y_finer_ind) - &
        y_finer_cell_cent_coords( y_finer_ind  - 1))
    multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 1) * &
        multigrid( finer_lvl) % cells( 2)
    multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
    multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = multigrid( &
        finer_lvl) % cells( 3)
    multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
! Otherwise, uses bilinear interpolation.
else
    y_coarser_ind = multigrid( coarser_lvl) % cells( 1)
    z_coarser_ind = multigrid( coarser_lvl) % cells( 2)
    y_finer_ind = multigrid( finer_lvl) % cells( 1)
    z_finer_ind = multigrid( finer_lvl) % cells( 2)
    entry_ind = multigrid( finer_lvl) % R_row_ptr( multigrid( &
        coarser_lvl) % cells( 3))
    factor = 1.0d+0 / ((z_finer_cell_cent_coords( z_finer_ind) - &
        z_finer_cell_cent_coords( z_finer_ind - 1)) * &
        (y_finer_cell_cent_coords( y_finer_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1)))
    multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 1) * &
        multigrid( finer_lvl) % cells( 2) - 1
    multigrid( finer_lvl) % R_val_dble( entry_ind) = &
        (z_finer_cell_cent_coords( z_finer_ind) - &
        z_coarser_cell_cent_coords( z_coarser_ind)) * &
        (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * factor
    multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (y_finer_ind - 1) &
        * multigrid( finer_lvl) % cells( 2)
    multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = &
        (z_coarser_cell_cent_coords( z_coarser_ind) - &
        z_finer_cell_cent_coords( z_finer_ind - 1)) * &
        (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * factor
    multigrid( finer_lvl) % R_col_ind( entry_ind + 2) = y_finer_ind * &
        multigrid( finer_lvl) % cells( 2) - 1
    multigrid( finer_lvl) % R_val_dble( entry_ind + 2) = &
        (z_finer_cell_cent_coords( z_finer_ind) - &
        z_coarser_cell_cent_coords( z_coarser_ind)) * &
        (y_coarser_cell_cent_coords( y_coarser_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
    multigrid( finer_lvl) % R_col_ind( entry_ind + 3) = y_finer_ind * &
        multigrid( finer_lvl) % cells( 2)
    multigrid( finer_lvl) % R_val_dble( entry_ind + 3) = &
        (z_coarser_cell_cent_coords( z_coarser_ind) - &
        z_finer_cell_cent_coords( z_finer_ind - 1)) * &
        (y_coarser_cell_cent_coords( y_coarser_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1)) * factor
end if

if( restr_prec == 1) then
    ! Typecasts the double array into the single one.
    do i = 1, multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) &
        % cells( 3) + 1) - 1
        multigrid( finer_lvl) % R_val_sngl( i) = sngl( multigrid( &
            finer_lvl) % R_val_dble( i))
    end do
    deallocate( multigrid( finer_lvl) % R_val_dble)
end if

end subroutine restr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the restriction operator with semi-coarsening in the y-direction on one grid. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine semi_restr_y( finer_lvl, y_finer_cell_cent_coords, &
    y_coarser_cell_cent_coords, restr_prec)
use pressure_gmg_intern_data, only: multigrid

! Declares the input variables.
! Declares the level of the finer grid, from which to restrict functions.
integer, intent( in) :: finer_lvl
! Declares the y-coordinates of the centers of the finer and coarser grid cells.
double precision, intent( in) :: y_finer_cell_cent_coords( multigrid( &
    finer_lvl) % cells( 1)), y_coarser_cell_cent_coords( multigrid( &
    finer_lvl + 1) % cells( 1))
integer, intent( in) :: restr_prec

! Delcares the local variables.
integer :: i
integer :: coarser_lvl
logical :: y_fine_cells_is_odd
integer :: cell_ind
integer :: y_finer_ind
integer :: y_coarser_ind, z_coarser_ind
integer :: entry_ind
double precision :: factor

coarser_lvl = finer_lvl + 1

y_fine_cells_is_odd = (mod( multigrid( finer_lvl) % cells( 1), 2) > 0)

! Allocates the row pointer array.
allocate( multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) % &
    cells( 3) + 1))
! Sets the first pointer to 1 (the first entry of the matrix).
multigrid( finer_lvl) % R_row_ptr( 1) = 1
! Initializes all the following entries to 2 since all cells on the coarser grid not adjacent to the right boundary have values interpolated from the 2 cells on the finer grid from which they are formed. Note that the pointer for the following row is incremented. The exceptions will be addressed below.
do i = 2, multigrid( coarser_lvl) % cells( 3) + 1
    multigrid( finer_lvl) % R_row_ptr( i) = 2
end do
! If the number of cells in the y-direction on the finer grid is odd, then the cells on the coarser grid adjacent to the right boundary are the same as the ones on the finer grid.
cell_ind = (multigrid( coarser_lvl) % cells( 1) - 1) * multigrid( &
    coarser_lvl) % cells( 2) + 1
if( y_fine_cells_is_odd) then
    do i = 1, multigrid( coarser_lvl) % cells( 2)
        multigrid( finer_lvl) % R_row_ptr( cell_ind + i) = 1
    end do
end if
! Adds each row pointer to the next one, in succession.
do i = 1, multigrid( coarser_lvl) % cells( 3)
    multigrid( finer_lvl) % R_row_ptr( i + 1) = multigrid( finer_lvl) % &
        R_row_ptr( i) + multigrid( finer_lvl) % R_row_ptr( i + 1)
end do

! Allocates the entry column index and value arrays.
allocate( multigrid( finer_lvl) % R_col_ind( multigrid( finer_lvl) % &
    R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
allocate( multigrid( finer_lvl) % R_val_dble( multigrid( finer_lvl) % &
    R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
if( restr_prec == 1) then
    allocate( multigrid( finer_lvl) % R_val_sngl( multigrid( finer_lvl) &
        % R_row_ptr( multigrid( coarser_lvl) % cells( 3) + 1) - 1))
end if

! Uses linear interpolation along the y-direction for the cells on the coarser grid not adjacent to the right boundary.
do y_coarser_ind = 1, multigrid( coarser_lvl) % cells( 1) - 1
    ! Computes the index in the y-direction of the upper of the two cells on the finer grid from which the current cell on the coarser grid are formed.
    y_finer_ind = 2 * y_coarser_ind
    ! Computes a common factor.
    factor = (y_coarser_cell_cent_coords( y_coarser_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1)) / &
        (y_finer_cell_cent_coords( y_finer_ind) - &
        y_finer_cell_cent_coords( y_finer_ind - 1))
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2)
        ! Computes the entry index for the left cell on the finer grid.
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (y_coarser_ind - &
            1) * multigrid( coarser_lvl) % cells( 2) + z_coarser_ind)
        ! Sets the entry column index and value for the left cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (y_finer_ind - 2) &
            * multigrid( finer_lvl) % cells( 2) + z_coarser_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
        ! Sets the entry column index and value for the right cell on the finer grid.
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (y_finer_ind &
            - 1) * multigrid( finer_lvl) % cells( 2) + z_coarser_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
    end do
end do
! If the number of cells in the y-direction on the finer grid is odd, then uses injection for the cells on the coarser grid adjacent to the right boundary.
if( y_fine_cells_is_odd) then
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2)
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (multigrid( &
            coarser_lvl) % cells( 1) - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_coarser_ind)
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (multigrid( &
            finer_lvl) % cells( 1) - 1) * multigrid( finer_lvl) % cells( &
            2) + z_coarser_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0
    end do
! Otherwise, uses linear interpolation as before.
else
    factor = (y_coarser_cell_cent_coords( multigrid( coarser_lvl) % &
        cells( 1)) - y_finer_cell_cent_coords( multigrid( finer_lvl) % &
        cells( 1) - 1)) / (y_finer_cell_cent_coords( multigrid( &
        finer_lvl) % cells( 1)) - y_finer_cell_cent_coords( multigrid( &
        finer_lvl) % cells( 1) - 1))
    do z_coarser_ind = 1, multigrid( coarser_lvl) % cells( 2)
        entry_ind = multigrid( finer_lvl) % R_row_ptr( (multigrid( &
            coarser_lvl) % cells( 1) - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_coarser_ind)
        multigrid( finer_lvl) % R_col_ind( entry_ind) = (multigrid( &
            finer_lvl) % cells( 1) - 2) * multigrid( finer_lvl) % cells( &
            2) + z_coarser_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind) = 1.0d+0 - factor
        multigrid( finer_lvl) % R_col_ind( entry_ind + 1) = (multigrid( &
            finer_lvl) % cells( 1) - 1) * multigrid( finer_lvl) % cells( &
            2) + z_coarser_ind
        multigrid( finer_lvl) % R_val_dble( entry_ind + 1) = factor
    end do
end if

if( restr_prec == 1) then
    ! Typecasts the double array into the single one.
    do i = 1, multigrid( finer_lvl) % R_row_ptr( multigrid( coarser_lvl) &
        % cells( 3) + 1) - 1
        multigrid( finer_lvl) % R_val_sngl( i) = sngl( multigrid( &
            finer_lvl) % R_val_dble( i))
    end do
    deallocate( multigrid( finer_lvl) % R_val_dble)
end if

end subroutine semi_restr_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the prolongation operator on one grid. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine prol( coarser_lvl, y_coarser_cell_widths, &
    z_coarser_cell_widths, y_finer_cell_cent_coords, &
    z_finer_cell_cent_coords, y_coarser_cell_cent_coords, &
    z_coarser_cell_cent_coords, prol_prec)
use pressure_gmg_intern_data, only : multigrid

! Declares the input variables.
! Declares the level of the coarser grid, to which to prolongate functions.
integer, intent( in) :: coarser_lvl
! Declares the widths of the coarser grid cells in the y and z-directions.
double precision, intent( in) :: y_coarser_cell_widths( multigrid( &
    coarser_lvl) % cells( 1)), z_coarser_cell_widths( multigrid( &
    coarser_lvl) % cells( 2))
! Declares the y and z-coordinates of the centers of the finer and coarser grid cells.
double precision, intent( in) :: y_finer_cell_cent_coords( multigrid( &
    coarser_lvl - 1) % cells( 1)), z_finer_cell_cent_coords( multigrid( &
    coarser_lvl - 1) % cells( 2)), y_coarser_cell_cent_coords( multigrid( &
    coarser_lvl) % cells( 1)), z_coarser_cell_cent_coords( multigrid( &
    coarser_lvl) % cells( 2))
integer, intent( in) :: prol_prec

! Declares the local variables.
integer :: i
integer :: finer_lvl
logical :: y_finer_cells_is_even, z_finer_cells_is_even
integer :: cell_ind
integer :: y_finer_ind, z_finer_ind
integer :: y_coarser_ind, z_coarser_ind
integer :: entry_ind
integer :: coarser_ind
integer :: finer_ind
double precision :: factor1, factor2

finer_lvl = coarser_lvl - 1

y_finer_cells_is_even = (mod( multigrid( finer_lvl) % cells( 1), 2) < 1)
z_finer_cells_is_even = (mod( multigrid( finer_lvl) % cells( 2), 2) < 1)

! Allocates the row pointer array.
allocate( multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % &
    cells( 3) + 1))
! Sets the first pointer to 1 (the first entry of the matrix).
multigrid( coarser_lvl) % P_row_ptr( 1) = 1
! Initializes all the following entries to 4 since all cells on the finer grid not adjacent to the boundary have values interpolated from their 4 neighbors on the coarser grid. Note that the pointer for the following row is incremented. The exceptions will be addressed below.
do z_finer_ind = 2, multigrid( finer_lvl) % cells( 3) + 1
    multigrid( coarser_lvl) % P_row_ptr( z_finer_ind) = 4
end do
! The cells on the finer grid adjacent to the corners have only 1 neighbor on the coarser grid.
! Sets the pointer for the lower left corner.
multigrid( coarser_lvl) % P_row_ptr( 2) = 1
! Sets the pointer for the upper left corner.
multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 2) + &
    1) = 1
! Sets the pointer for the lower right corner.
multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 3) - &
    multigrid( finer_lvl) % cells( 2) + 2) = 1
! Sets the pointer for the upper right corner.
multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 3) + &
    1) = 1
! The cells on the finer grid adjacent to the edges but not the corners have 2 neighbors on the coarser grid.
! Sets the pointers for the left and right edges.
cell_ind = multigrid( finer_lvl) % cells( 3) - multigrid( finer_lvl) % &
    cells( 2) + 1
do z_finer_ind = 2, multigrid( finer_lvl) % cells( 2) - 1
    multigrid( coarser_lvl) % P_row_ptr( z_finer_ind + 1) = 2
    multigrid( coarser_lvl) % P_row_ptr( cell_ind + z_finer_ind) = 2
end do
! Sets the pointers for the lower and upper edges.
do z_finer_ind = 2, multigrid( finer_lvl) % cells( 1) - 1
    multigrid( coarser_lvl) % P_row_ptr( (z_finer_ind - 1) * multigrid( &
        finer_lvl) % cells( 2) + 2) = 2
    multigrid( coarser_lvl) % P_row_ptr( z_finer_ind * multigrid( &
        finer_lvl) % cells( 2) + 1) = 2
end do
! Adds each row pointer to the next one, in succession.
do z_finer_ind = 1, multigrid( finer_lvl) % cells( 3)
    multigrid( coarser_lvl) % P_row_ptr( z_finer_ind + 1) = multigrid( &
        coarser_lvl) % P_row_ptr( z_finer_ind) + multigrid( coarser_lvl) &
        % P_row_ptr( z_finer_ind + 1)
end do

! Allocates the entry column index and value arrays.
allocate( multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) - 1))
allocate( multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) - 1))
if( prol_prec == 1) then
    allocate( multigrid( coarser_lvl) % P_val_sngl( multigrid( &
        coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) &
        - 1))
end if

! Uses injection for the cells on the finer grid at the corners.
multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % P_row_ptr( &
    1)) = 1
multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( 1)) = 1.0d+0
multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % P_row_ptr( &
    multigrid( finer_lvl) % cells( 2))) = multigrid( coarser_lvl) % &
    cells( 2)
multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 2))) = 1.0d+0
multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % P_row_ptr( &
    multigrid( finer_lvl) % cells( 3) - multigrid( finer_lvl) % cells( 2) &
    + 1)) = (multigrid( coarser_lvl) % cells( 1) - 1) * multigrid( &
    coarser_lvl) % cells( 2) + 1
multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3) - multigrid( finer_lvl) &
    % cells( 2) + 1)) = 1.0d+0
multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % P_row_ptr( &
    multigrid( finer_lvl) % cells( 3))) = multigrid( coarser_lvl) % &
    cells( 2) * multigrid( coarser_lvl) % cells( 1)
multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3))) = 1.0d+0

! Uses linear interpolation for the cells on the finer grid adjacent to the boundary.
! Loops through the cells in pairs.
do z_finer_ind = 2, multigrid( finer_lvl) % cells( 2) - 3, 2
    ! Computes the index in the z-direction of the lower of the two cells on the coarser grid from which the current pair of cells on the finer grid is formed. (Fortran rounds toward 0.)
    z_coarser_ind = z_finer_ind / 2
    ! Computes a common factor.
    factor1 = (z_finer_cell_cent_coords( z_finer_ind) - &
        z_coarser_cell_cent_coords( z_coarser_ind)) * 2.0d+0 / &
        (z_coarser_cell_widths( z_coarser_ind) + z_coarser_cell_widths( &
        z_coarser_ind + 1))
    ! Sets the entry column indices and values for the lower cell on the finer grid along the left boundary.
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = factor1
    factor2 = (z_coarser_cell_cent_coords( z_coarser_ind + 1) - &
        z_finer_cell_cent_coords( z_finer_ind + 1)) * 2.0d+0 / &
        (z_coarser_cell_widths( z_coarser_ind) + z_coarser_cell_widths( &
        z_coarser_ind + 1))
    ! Sets the entry column indices and values for the upper cell on the finer grid along the left boundary.
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = 1.0d+0 - factor2

    ! Does the same for the cells on the finer grid along the right boundary.
    z_coarser_ind = (multigrid( coarser_lvl) % cells( 1) - 1) * &
        multigrid( coarser_lvl) % cells( 2) + z_coarser_ind
    entry_ind = multigrid( finer_lvl) % cells( 3) - multigrid( finer_lvl) &
        % cells( 2) + z_finer_ind
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind) + 1) = factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind + 1)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind + 1)) = factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind + 1) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( entry_ind + 1) + 1) = 1.0d+0 - factor2
end do

! If the number of cells in the z-direction on the finer grid is even, then does the same for the final pair of cells before the corner.
if( z_finer_cells_is_even) then
    z_finer_ind = multigrid( finer_lvl) % cells( 2) - 2
    z_coarser_ind = z_finer_ind / 2
    factor1 = (z_finer_cell_cent_coords( z_finer_ind) - &
        z_coarser_cell_cent_coords( z_coarser_ind)) * 2.0d+0 / &
        (z_coarser_cell_widths( z_coarser_ind) + z_coarser_cell_widths( &
        z_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = factor1
    factor2 = (z_coarser_cell_cent_coords( z_coarser_ind + 1) - &
        z_finer_cell_cent_coords( z_finer_ind + 1)) * 2.0d+0 / &
        (z_coarser_cell_widths( z_coarser_ind) + z_coarser_cell_widths( &
        z_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = 1.0d+0 - factor2

    z_finer_ind = multigrid( finer_lvl) % cells( 3) - multigrid( &
        finer_lvl) % cells( 2) + z_finer_ind
    z_coarser_ind = (multigrid( coarser_lvl) % cells( 1) - 1) * &
        multigrid( coarser_lvl) % cells( 2) + z_coarser_ind
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1)) = factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind + 1) + 1) = 1.0d+0 - factor2
! Otherwise, uses linear interpolation for the last cell before the corner.
else
    z_finer_ind = multigrid( finer_lvl) % cells( 2) - 1
    z_coarser_ind = z_finer_ind / 2
    factor1 = (z_finer_cell_cent_coords( z_finer_ind) - &
        z_coarser_cell_cent_coords( z_coarser_ind)) * 2.0d+0 / &
        (z_coarser_cell_widths( z_coarser_ind) + z_coarser_cell_widths( &
        z_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = factor1

    z_finer_ind = multigrid( finer_lvl) % cells( 3) - multigrid( &
        finer_lvl) % cells( 2) + z_finer_ind
    z_coarser_ind = (multigrid( coarser_lvl) % cells( 1) - 1) * &
        multigrid( coarser_lvl) % cells( 2) + z_coarser_ind
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = z_coarser_ind + 1
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind) + 1) = factor1
end if

do y_finer_ind = 2, multigrid( finer_lvl) % cells( 1) - 3, 2
    ! Sets the entry column indices and values for the cells on the finer grid along the bottom boundary.
    y_coarser_ind = y_finer_ind / 2
    finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) + 1
    coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % cells( &
        2) + 1
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1
    factor2 = (y_coarser_cell_cent_coords( y_coarser_ind + 1) - &
        y_finer_cell_cent_coords( y_finer_ind + 1)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        coarser_ind + multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        1.0d+0 - factor2

    ! Sets the entry column indices and values for the cells on the finer grid along the upper boundary.
    finer_ind = y_finer_ind * multigrid( finer_lvl) % cells( 2)
    coarser_ind = y_coarser_ind * multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        coarser_ind + multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        1.0d+0 - factor2
end do

! If the number of cells in the y-direction on the finer grid is even, then does the same for the final pair of cells before the corner.
if( y_finer_cells_is_even) then
    y_finer_ind = multigrid( finer_lvl) % cells( 1) - 2
    y_coarser_ind = y_finer_ind / 2
    finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) + 1
    coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % cells( &
        2) + 1
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1
    factor2 =  (y_coarser_cell_cent_coords( y_coarser_ind + 1) - &
        y_finer_cell_cent_coords( y_finer_ind + 1)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        coarser_ind + multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        1.0d+0 - factor2

    finer_ind = y_finer_ind * multigrid( finer_lvl) % cells( 2)
    coarser_ind = y_coarser_ind * multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
        factor2
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        coarser_ind + multigrid( coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + 1) = &
        1.0d+0 - factor2
! Otherwise, uses linear interpolation for the last cell before the corner.
else
    y_finer_ind = multigrid( finer_lvl) % cells( 1) - 1
    y_coarser_ind = y_finer_ind / 2
    coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % cells( &
        2) + 1
    finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) + 1
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1

    coarser_ind = y_coarser_ind * multigrid( coarser_lvl) % cells( 2)
    finer_ind = y_finer_ind * multigrid( finer_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = coarser_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind)) = 1.0d+0 - factor1
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
        coarser_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( finer_ind) + 1) = factor1
end if

! Uses bilinear interpolation for the cells on the coarser grid not adjacent to the boundary.
do y_finer_ind = 2, multigrid( finer_lvl) % cells( 1) - 1
    do z_finer_ind = 2, multigrid( finer_lvl) % cells( 2) - 1
        ! Computes the lower indices in the y and z-directions of the four cells on the coarser grid to which the current cell on the finer grid is adjacent.
        y_coarser_ind = y_finer_ind / 2
        z_coarser_ind = z_finer_ind / 2
        finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) &
            + z_finer_ind
        coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_coarser_ind
        factor1 = 4.0d+0 / ((y_coarser_cell_widths( y_coarser_ind) + &
            y_coarser_cell_widths( y_coarser_ind + 1)) * &
            (z_coarser_cell_widths( z_coarser_ind) + &
            z_coarser_cell_widths( z_coarser_ind + 1)))
        ! Sets the entry column index and value for the lower left cell on the coarser grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = (y_coarser_cell_cent_coords( &
            y_coarser_ind + 1) - y_finer_cell_cent_coords( y_finer_ind)) &
            * (z_coarser_cell_cent_coords( z_coarser_ind + 1) - &
            z_finer_cell_cent_coords( z_finer_ind)) * factor1
        ! Sets the entry column index and value for the upper left cell on the coarser grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = coarser_ind + 1
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = (y_coarser_cell_cent_coords( &
            y_coarser_ind + 1) - y_finer_cell_cent_coords( y_finer_ind)) &
            * (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * factor1
        ! Sets the entry column index and value for the lower right cell on the coarser grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 2) = coarser_ind + multigrid( &
            coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 2) = (y_finer_cell_cent_coords( &
            y_finer_ind) - y_coarser_cell_cent_coords( y_coarser_ind)) * &
            (z_coarser_cell_cent_coords( z_coarser_ind + 1) - &
            z_finer_cell_cent_coords( z_finer_ind)) * factor1
        ! Sets the entry column index and value for the upper right cell on the coarser grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 3) = coarser_ind + multigrid( &
            coarser_lvl) % cells( 2) + 1
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 3) = (y_finer_cell_cent_coords( &
            y_finer_ind) - y_coarser_cell_cent_coords( y_coarser_ind)) * &
            (z_finer_cell_cent_coords( z_finer_ind) - &
            z_coarser_cell_cent_coords( z_coarser_ind)) * factor1
    end do
end do

if( prol_prec == 1) then
    ! Typecasts the double array into the single one.
    do i = 1, multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) &
        % cells( 3) + 1) - 1
        multigrid( coarser_lvl) % P_val_sngl( i) = sngl( multigrid( &
            coarser_lvl) % P_val_dble( i))
    end do
    deallocate( multigrid( coarser_lvl) % P_val_dble)
end if

end subroutine prol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the prolongation operator with semi-prolongation in the y-direction on one grid. The routine assumes that the cells are numbered in the direction of increasing coordinates, in z-major order.
subroutine semi_prol_y( coarser_lvl, y_coarser_cell_widths, &
    y_finer_cell_cent_coords, y_coarser_cell_cent_coords, prol_prec)
use pressure_gmg_intern_data, only: multigrid

! Declares the input variables.
! Declares the level of the coarser grid, to which to prolongate functions.
integer, intent( in) :: coarser_lvl
! Declares the widths of the coarser grid cells in the y-direction.
double precision, intent( in) :: y_coarser_cell_widths( multigrid( &
    coarser_lvl) % cells( 1))
! Declares the y-coordinates of the centers of the finer and coarser grid cells.
double precision, intent( in) :: y_finer_cell_cent_coords( multigrid( &
    coarser_lvl - 1) % cells( 1)), y_coarser_cell_cent_coords( multigrid( &
    coarser_lvl) % cells( 1))
integer, intent( in) :: prol_prec

! Declares the local variables.
integer :: i
integer :: finer_lvl
logical :: y_finer_cells_is_even
integer :: cell_ind
integer :: y_finer_ind, z_finer_ind
integer :: y_coarser_ind
integer :: finer_ind, coarser_ind
double precision :: factor1, factor2

finer_lvl = coarser_lvl - 1

y_finer_cells_is_even = (mod( multigrid( finer_lvl) % cells( 1), 2) < 1)

! Allocates the row pointer array.
allocate( multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1))
! Sets the first pointer to 1 (the first entry of the matrix).
multigrid( coarser_lvl) % P_row_ptr( 1) = 1
! Initializes all the following entries to 2 since all cells on the finer grid not adjacent to the left or right edges have values interpolated from their 2 neighbors on the coarser grid. Note that the pointer for the following row is incremented. The exceptions will be addressed below.
do i = 2, multigrid( finer_lvl) % cells( 3) + 1
    multigrid( coarser_lvl) % P_row_ptr( i) = 2
end do
! The cells on the finer grid adjacent to the left or right edges, including all 4 corners, have only 1 neighbor on the coarser grid.
cell_ind = multigrid( finer_lvl) % cells( 3) - multigrid( finer_lvl) % &
    cells( 2) + 1
do i = 1, multigrid( finer_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_row_ptr( i + 1) = 1
    multigrid( coarser_lvl) % P_row_ptr( cell_ind + i) = 1
end do
! Adds each row pointer to the next one, in succession.
do i = 1, multigrid( finer_lvl) % cells( 3)
    multigrid( coarser_lvl) % P_row_ptr( i + 1) = multigrid( coarser_lvl) &
    % P_row_ptr( i) + multigrid( coarser_lvl) % P_row_ptr( i + 1)
end do

! Allocates the entry column index and value arrays.
allocate( multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) - 1))
allocate( multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
    P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) - 1))
if( prol_prec == 1) then
    allocate( multigrid( coarser_lvl) % P_val_sngl( multigrid( &
        coarser_lvl) % P_row_ptr( multigrid( finer_lvl) % cells( 3) + 1) &
        - 1))
end if

! Uses injection for the cells on the finer grid adjacent to the left and right boundaries.
do z_finer_ind = 1, multigrid( finer_lvl) % cells( 2)
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = z_finer_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( z_finer_ind)) = 1.0d+0
    multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
        P_row_ptr( multigrid( finer_lvl) % cells( 3) - multigrid( &
        finer_lvl) % cells( 2) + z_finer_ind)) = (multigrid( coarser_lvl) &
        % cells( 1) - 1) * multigrid( coarser_lvl) % cells( 2) + &
        z_finer_ind
    multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
        P_row_ptr( multigrid( finer_lvl) % cells( 3) - multigrid( &
        finer_lvl) % cells( 2) + z_finer_ind)) = 1.0d+0
end do

! Uses linear interpolation along the y-direction for the cells on the finer grid not adjacent to the left or right boundaries.
do y_finer_ind = 2, multigrid( finer_lvl) % cells( 1) - 3, 2
    ! Computes the index in the z-direction of the left of the two cells on the coarser grid from which the current pair of cells on the finer grid is formed. Uses MATLAB's bitshift for integer division by 2, which is less expensive. bitshift truncates the fractional part.
    y_coarser_ind = y_finer_ind / 2
    ! Computes common factors.
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    factor2 = (y_coarser_cell_cent_coords( y_coarser_ind + 1) - &
        y_finer_cell_cent_coords( y_finer_ind + 1)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    do z_finer_ind = 1, multigrid( finer_lvl) % cells( 2)
        coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_finer_ind
        finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) &
            + z_finer_ind
        ! Sets the entry column indices and values for the left cell on the finer grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = 1.0d+0 - factor1
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
            coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = factor1
        ! Sets the entry column indices and values for the right cell on the finer grid.
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
            coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
            factor2
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + &
            1) = coarser_ind + multigrid( coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + &
            1) = 1.0d+0 - factor2
    end do
end do

! If the number of cells in the y-direction on the finer grid is even, then does the same for the final pair of cells before the boundary.
if( y_finer_cells_is_even) then
    y_finer_ind = multigrid( finer_lvl) % cells( 1) - 2
    y_coarser_ind = y_finer_ind / 2
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    factor2 = (y_coarser_cell_cent_coords( y_coarser_ind + 1) - &
        y_finer_cell_cent_coords( y_finer_ind + 1)) * 2.0d+0 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    do z_finer_ind = 1, multigrid( finer_lvl) % cells( 2)
        coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_finer_ind
        finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) &
            + z_finer_ind
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = 1.0d+0 - factor1
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
            coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = factor1
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
            coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2))) = &
            factor2
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + &
            1) = coarser_ind + multigrid( coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind + multigrid( finer_lvl) % cells( 2)) + &
            1) = 1.0d+0 - factor2
    end do
! Otherwise, uses linear interpolation for the last cell before the boundary.
else
    y_finer_ind = multigrid( finer_lvl) % cells( 1) - 1
    y_coarser_ind = y_finer_ind / 2
    factor1 = (y_finer_cell_cent_coords( y_finer_ind) - &
        y_coarser_cell_cent_coords( y_coarser_ind)) * 2 / &
        (y_coarser_cell_widths( y_coarser_ind) + y_coarser_cell_widths( &
        y_coarser_ind + 1))
    do z_finer_ind = 1, multigrid( finer_lvl) % cells( 2)
        coarser_ind = (y_coarser_ind - 1) * multigrid( coarser_lvl) % &
            cells( 2) + z_finer_ind
        finer_ind = (y_finer_ind - 1) * multigrid( finer_lvl) % cells( 2) &
            + z_finer_ind
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = coarser_ind
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind)) = 1.0d+0 - factor1
        multigrid( coarser_lvl) % P_col_ind( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = coarser_ind + multigrid( &
            coarser_lvl) % cells( 2)
        multigrid( coarser_lvl) % P_val_dble( multigrid( coarser_lvl) % &
            P_row_ptr( finer_ind) + 1) = factor1
    end do
end if

if( prol_prec == 1) then
    ! Typecasts the double array into the single one.
    do i = 1, multigrid( coarser_lvl) % P_row_ptr( multigrid( finer_lvl) &
        % cells( 3) + 1) - 1
        multigrid( coarser_lvl) % P_val_sngl( i) = sngl( multigrid( &
            coarser_lvl) % P_val_dble( i))
    end do
    deallocate( multigrid( coarser_lvl) % P_val_dble)
end if

end subroutine semi_prol_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Performs a v-cycle.
subroutine v_cyc( alpha_dble, lvl, sol_dble, rhs_dble, GS_type, &
    pre_smooth_iter, post_smooth_iter, remov_null, pre_smooth_prec, &
    resid_prec, restr_prec, coarse_prec, zero_coarse_prec, &
    post_smooth_prec, prol_prec)
use pressure_gmg_intern_data
use gmg_common_intern_proced, only: GS_dble, GS_sngl, resid_dble, &
    resid_sngl, mat_vec_mult_crs_dble, mat_vec_mult_crs_sngl, &
    LU_solv_dble, LU_solv_sngl, back_solv_dble, back_solv_sngl, &
    mat_vec_mult_add_crs_dble, mat_vec_mult_add_crs_sngl

! Declares the input variables.
! Declares the value of the harmonic.
double precision, intent( in) :: alpha_dble
! Declares the grid at which to start the v-cycle. Values not equal to 1 will be used for full multigrid (FMG).
integer, intent( in) :: lvl
! Declares the solution on the fine-grid. The array is passed by reference and will be overwritten with the solution after the v-cycle.
double precision, intent( inout) :: sol_dble( multigrid( lvl) % cells( 3))
! Declares the right-hand-side on the fine-grid.
double precision, intent( in) :: rhs_dble( multigrid( lvl) % cells( 3))
! Declares the type of smoother.
character( len = 5), intent( in) :: GS_type
! Declares the number of smoothing iterations on each grid on the down-stroke.
integer, intent( in) :: pre_smooth_iter
! Declares the number of smoothing iterations on each grid on the up-stroke.
integer, intent( in) :: post_smooth_iter
! Declares the intervals over which to remove the right nullspace.
character( len = 4), intent( in) :: remov_null
! Declares the precisions to be used for pre- and post-smoothing (except on the fine-grid), computation of the residual (except on the fine-grid), transfer operations (restriction and prolongation), and coarse-grid solves in the non-zero and zero alpha cases. If 1, it is single, if 2, it is double.
integer :: pre_smooth_prec, resid_prec, restr_prec, coarse_prec, &
    zero_coarse_prec, post_smooth_prec, prol_prec

! Declares the local variables.
integer :: i, j
! Declares the number of grids of multigrid.
integer :: lvls
! Declares the component of the solution vector in the right nullspace.
double precision :: proj
real :: proj_sngl
! Declares the inputs of the LAPACK functions.
! integer :: INFO

lvls = size( multigrid)

! Descends the grid hierarchy in the down-stroke.
do i = lvl, lvls - 1
    ! The solution and right-hand-side on the fine-grid are passed into  pressure_gmg_solv and are not stored in multigrid.
    if( i == 1) then
        ! Performs pre-smoothing.
        call GS_dble( multigrid( 1) % cells( 1), multigrid( 1) % cells( &
            2), multigrid( 1) % A_row_ptr, size( multigrid( 1) % &
            A_col_ind), multigrid( 1) % A_col_ind, multigrid( 1) % &
            A_val_dble, sol_dble, rhs_dble, GS_type, multigrid( 1) % &
            GS_rhs_dble, multigrid( 1) % GS_up_diag_dble, pre_smooth_iter)
        if( prol_prec == 1) then
            ! Makes a single precision copy of the solution in preparation for the up-stroke.
            do j = 1, multigrid( 1) % cells( 3)
               multigrid( 1) % sol_sngl( j) = sngl( sol_dble( j))
            end do
        end if
        ! Computes the residual of the new solution obtained by smoothing.
        call resid_dble( multigrid( 1) % cells( 3), multigrid( 1) % &
            A_row_ptr, size( multigrid( 1) % A_col_ind), multigrid( 1) % &
            A_col_ind, multigrid( 1) % A_val_dble, multigrid( 1) % cells( &
            3), sol_dble, rhs_dble, multigrid( 1) % res_dble)
        if( ( restr_prec == 1)) then
            ! Makes a single precision copy of the residual in preparation for restriction.
            do j = 1, multigrid( 1) % cells( 3)
                multigrid( 1) % res_sngl( j) = sngl( multigrid( 1) % &
                    res_dble( j))
            end do
        end if
    ! Does the same for all coarser grids except for the coarse-grid.
    else
        if( pre_smooth_prec == 2) then
            call GS_dble( multigrid( i) % cells( 1), multigrid( i) % &
                cells( 2), multigrid( i) % A_row_ptr, size( multigrid( i) &
                % A_col_ind), multigrid( i) % A_col_ind, multigrid( i) % &
                A_val_dble, multigrid( i) % sol_dble, multigrid( i) % &
                rhs_dble, GS_type, multigrid( i) % GS_rhs_dble, &
                multigrid( i) % GS_up_diag_dble, pre_smooth_iter)
            if( (resid_prec == 1) .or. (prol_prec == 1)) then
                do j = 1, multigrid( i) % cells( 3)
                    multigrid( i) % sol_sngl( j) = sngl( multigrid( i) % &
                        sol_dble( j))
                end do
            end if
        else
            call GS_sngl( multigrid( i) % cells( 1), multigrid( i) % &
                cells( 2), multigrid( i) % A_row_ptr, size( multigrid( i) &
                % A_col_ind), multigrid( i) % A_col_ind, multigrid( i) % &
                A_val_sngl, multigrid( i) % sol_sngl, multigrid( i) % &
                rhs_sngl, GS_type, multigrid( i) % GS_rhs_sngl, &
                multigrid( i) % GS_up_diag_sngl, pre_smooth_iter)
            if( (resid_prec == 2) .or. (prol_prec == 2)) then
                do j = 1, multigrid( i) % cells( 3)
                    multigrid( i) % sol_dble( j) = dble( multigrid( i) % &
                        sol_sngl( j))
                end do
            end if
        end if
        if( resid_prec == 2) then
            call resid_dble( multigrid( i) % cells( 3), multigrid( i) % &
                A_row_ptr, size( multigrid( i) % A_col_ind), multigrid( &
                i) % A_col_ind, multigrid( i) % A_val_dble, multigrid( i) &
                % cells( 3), multigrid( i) % sol_dble, multigrid( i) % &
                rhs_dble, multigrid( i) % res_dble)
            if( restr_prec == 1) then
                do j = 1, multigrid( i) % cells( 3)
                    multigrid( i) % res_sngl( j) = sngl( multigrid( i) % &
                        res_dble( j))
                end do
            end if
        else
            call resid_sngl( multigrid( i) % cells( 3), multigrid( i) % &
                A_row_ptr, size( multigrid( i) % A_col_ind), multigrid( &
                i) % A_col_ind, multigrid( i) % A_val_sngl, multigrid( i) &
                % cells( 3), multigrid( i) % sol_sngl, multigrid( i) % &
                rhs_sngl, multigrid( i) % res_sngl)
            if( restr_prec == 2) then
                do j = 1, multigrid( i) % cells( 3)
                    multigrid( i) % res_dble( j) = dble( multigrid( i) % &
                        res_sngl( j))
                end do
            end if
        end if
    end if
    if( restr_prec == 2) then
        ! Restricts the residual to the next coarser grid.
        call mat_vec_mult_crs_dble( multigrid( i + 1) % cells( 3), &
            multigrid( i) % R_row_ptr, size( multigrid( i) % R_col_ind), &
            multigrid( i) % R_col_ind, multigrid( i) % R_val_dble, &
            multigrid( i) % cells( 3), multigrid( i) % res_dble, &
            multigrid( i + 1) % rhs_dble)
        if( (pre_smooth_prec == 1) .or. (post_smooth_prec == 1)) then
            ! Makes a dblele precision copy of the right-hand-side in preparation for smoothing or a coarse-grid solve.
            do j = 1, multigrid( i + 1) % cells( 3)
                multigrid( i + 1) % rhs_sngl( j) = sngl( multigrid( i + &
                    1) % rhs_dble( j))
            end do
        end if
    else
        call mat_vec_mult_crs_sngl( multigrid( i + 1) % cells( 3), &
            multigrid( i) % R_row_ptr, size( multigrid( i) % R_col_ind), &
            multigrid( i) % R_col_ind, multigrid( i) % R_val_sngl, &
            multigrid( i) % cells( 3), multigrid( i) % res_sngl, &
            multigrid( i + 1) % rhs_sngl)
        if( (pre_smooth_prec == 2) .or. (post_smooth_prec == 2)) then
            do j = 1, multigrid( i + 1) % cells( 3)
                multigrid( i + 1) % rhs_dble( j) = dble( multigrid( i + &
                1) % rhs_sngl( j))
            end do
        end if
    end if
    ! Initializes or resets the solution vector on the next coarser grid to 0.
    if( (pre_smooth_prec == 2) .or. (resid_prec == 2) .or. (prol_prec == &
        2) .or. (post_smooth_prec == 2)) then
        do j = 1, multigrid( i + 1) % cells( 3)
            multigrid( i + 1) % sol_dble( j) = 0.0d+0
        end do
    else if( (pre_smooth_prec == 1) .or. (resid_prec == 1) .or. &
        (prol_prec == 1) .or. (post_smooth_prec == 1)) then
        do j = 1, multigrid( i + 1) % cells( 3)
            multigrid( i + 1) % sol_sngl( j) = 0.0e+0
        end do
    end if
end do

if( (restr_prec == 2) .and. ((coarse_prec == 1) .or. (zero_coarse_prec == &
    1))) then
    do i = 1, multigrid( lvls) % cells( 3)
        multigrid( lvls) % rhs_sngl( i) = sngl( multigrid( lvls) % &
            rhs_dble( i))
    end do
else if( (restr_prec == 1) .and. ((coarse_prec == 2) .or. &
    (zero_coarse_prec == 2))) then
    do i = 1, multigrid( lvls) % cells( 3)
        multigrid( lvls) % rhs_dble( i) = dble( multigrid( lvls) % &
            rhs_sngl( i))
    end do
end if
! Computes the coarse-grid solution.
if( alpha_dble /= 0.0d+0) then
    if( coarse_prec == 2) then
        ! Permutes the right-hand-side on the coarse-grid.
        do i = 1, multigrid( lvls) % cells( 3)
            p_rhs_coarse_dble( i) = multigrid( lvls) % rhs_dble( &
                p_coarse( i))
        end do
        ! If alpha is not 0, uses the LU decomposition.
        call LU_solv_dble( multigrid( lvls) % cells( 3), LU_coarse_dble, &
            p_rhs_coarse_dble)
        ! ! Makes a copy of the coarse-grid discretization matrix in dense format.
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     do j = 1, multigrid( lvls) % cells( 3)
        !         LU_coarse_dble( i, j) = 0.0d+0
        !     end do
        ! end do
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
        !         A_row_ptr( i + 1) - 1
        !         LU_coarse_dble( i, multigrid( lvls) % A_col_ind( j)) = &
        !             multigrid( lvls) % A_val_dble( j)
        !     end do
        ! end do
        ! ! Copies the right-hand-side on the coarse-grid.
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     p_rhs_coarse_dble( i) = multigrid( lvls) % rhs_dble( i)
        ! end do
        ! ! Solves the linear system by LU using the LAPACK function dgesv (http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html).
        ! call dgesv( multigrid( lvls) % cells( 3), 1, LU_coarse_dble, &
        !     multigrid( lvls) % cells( 3), p_coarse, p_rhs_coarse_dble, &
        !     multigrid( lvls) % cells( 3), INFO)
        if( prol_prec == 2) then
            do j = 1, multigrid( lvls) % cells( 3)
                multigrid( lvls) % sol_dble( j) = p_rhs_coarse_dble( j)
            end do
        else
            do j = 1, multigrid( lvls) % cells( 3)
                multigrid( lvls) % sol_sngl( j) = sngl( &
                    p_rhs_coarse_dble( j))
            end do
        end if
    else
        do i = 1, multigrid( lvls) % cells( 3)
            p_rhs_coarse_sngl( i) = multigrid( lvls) % rhs_sngl( &
                p_coarse( i))
        end do
        ! If alpha is not 0, uses the LU decomposition.
        call LU_solv_sngl( multigrid( lvls) % cells( 3), LU_coarse_sngl, &
            p_rhs_coarse_sngl)
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     do j = 1, multigrid( lvls) % cells( 3)
        !         LU_coarse_sngl( i, j) = 0.0e+0
        !     end do
        ! end do
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
        !         A_row_ptr( i + 1) - 1
        !         LU_coarse_sngl( i, multigrid( lvls) % A_col_ind( j)) = &
        !             multigrid( lvls) % A_val_sngl( j)
        !     end do
        ! end do
        ! ! Copies the right-hand-side on the coarse-grid.
        ! do i = 1, multigrid( lvls) % cells( 3)
        !     p_rhs_coarse_sngl( i) = multigrid( lvls) % rhs_sngl( i)
        ! end do
        ! ! Solves the linear system by LU using the LAPACK function sgesv (http://www.netlib.org/lapack/explore-html/d0/db8/group__real_g_esolve_ga3b05fb3999b3d7351cb3101a1fd28e78.html).
        ! call sgesv( multigrid( lvls) % cells( 3), 1, LU_coarse_sngl, &
        !     multigrid( lvls) % cells( 3), p_coarse, p_rhs_coarse_sngl, &
        !     multigrid( lvls) % cells( 3), INFO)
        if( prol_prec == 2) then
            do j = 1, multigrid( lvls) % cells( 3)
                multigrid( lvls) % sol_dble( j) = dble( &
                    p_rhs_coarse_sngl( j))
            end do
        else
            do j = 1, multigrid( lvls) % cells( 3)
                multigrid( lvls) % sol_sngl( j) = p_rhs_coarse_sngl( j)
            end do
        end if
    end if
else
    if( zero_coarse_prec == 2) then
        call dgemv( 'T', multigrid( lvls) % cells( 3), multigrid( lvls) % &
            cells( 3) - 1, 1.0d+0, Q_coarse_dble, multigrid( lvls) % &
            cells( 3), multigrid( lvls) % rhs_dble, 1, 0.0d+0, &
            QR_rhs_dble, 1)
        call back_solv_dble( multigrid( lvls) % cells( 3) - 1, &
            R_coarse_dble, QR_rhs_dble)
        do i = 1, multigrid( lvls) % cells( 3) - 1
            multigrid( lvls) % sol_dble( e_coarse( i)) = QR_rhs_dble( i)
        end do
        ! ! If alpha is 0, uses the SVD truncated by 1 dimension, since we know that the right nullspace of our discretization matrices have dimension 1.
        ! call SVD_solv_dble( multigrid( lvls) % cells( 3), U_coarse, &
        !     S_coarse, VT_coarse, 1, multigrid( lvls) % sol_dble, &
        !     multigrid( lvls) % rhs_dble)
        ! if( prol_prec == 1) then
        !     do j = 1, multigrid( lvls) % cells( 3)
        !         multigrid( lvls) % sol_sngl( j) = sngl( multigrid( lvls) &
        !             % sol_dble( j))
        !     end do
        ! end if
    else
        call sgemv( 'T', multigrid( lvls) % cells( 3), multigrid( lvls) % &
            cells( 3) - 1, 1.0e+0, Q_coarse_sngl, multigrid( lvls) % &
            cells( 3), multigrid( lvls) % rhs_sngl, 1, 0.0e+0, &
            QR_rhs_sngl, 1)
        call back_solv_sngl( multigrid( lvls) % cells( 3) - 1, &
            R_coarse_sngl, QR_rhs_sngl)
        do i = 1, multigrid( lvls) % cells( 3) - 1
            multigrid( lvls) % sol_sngl( e_coarse( i)) = QR_rhs_sngl( i)
        end do
        ! call SVD_solv_sngl( multigrid( lvls) % cells( 3), U_coarse_sngl, &
        !     S_coarse_sngl, VT_coarse_sngl, 1, multigrid( lvls) % &
        !     sol_sngl, multigrid( lvls) % rhs_sngl)
        ! if( prol_prec == 2) then
        !     do j = 1, multigrid( lvls) % cells( 3)
        !         multigrid( lvls) % sol_dble( j) = dble( multigrid( lvls) &
        !             % sol_sngl( j))
        !     end do
        ! end if
    end if
end if

! Ascends the grid hierarchy in the up-stroke.
do i = lvls, lvl + 1, -1
    if( i == 2) then
        if( prol_prec == 2) then
            ! Prolongates the error to the next finer grid and corrects the solution.
            call mat_vec_mult_add_crs_dble( multigrid( 1) % cells( 3), &
                multigrid( 2) % P_row_ptr, size( multigrid( 2) % &
                P_col_ind), multigrid( 2) % P_col_ind, multigrid( 2) % &
                P_val_dble, multigrid( 2) % cells( 3), multigrid( 2) % &
                sol_dble, sol_dble)
        else
            call mat_vec_mult_add_crs_sngl( multigrid( 1) % cells( 3), &
                multigrid( 2) % P_row_ptr, size( multigrid( 2) % &
                P_col_ind), multigrid( 2) % P_col_ind, multigrid( 2) % &
                P_val_sngl, multigrid( 2) % cells( 3), multigrid( 2) % sol_sngl, multigrid( 1) % sol_sngl)
            ! Makes a dblele precision copy of the solution in preparation for smoothing.
            do j = 1, multigrid( 1) % cells( 3)
                sol_dble( j) = dble( multigrid( 1) % sol_sngl( j))
            end do
        end if
        if( (alpha_dble == 0.0d+0) .and. (remov_null == 'grid')) then
            ! Removes the right nullspace component of the solution.
            proj = sum( sol_dble) / multigrid( 1) % cells( 3)
            do j = 1, multigrid( 1) % cells( 3)
                sol_dble( j) = sol_dble( j) - proj
            end do
        end if
        ! Performs post-smoothing.
        call GS_dble( multigrid( 1) % cells( 1), multigrid( 1) % cells( &
            2), multigrid( 1) % A_row_ptr, size( multigrid( 1) % &
            A_col_ind), multigrid( 1) % A_col_ind, multigrid( 1) % &
            A_val_dble, sol_dble, rhs_dble, GS_type, multigrid( 1) % &
            GS_rhs_dble, multigrid( 1) % GS_up_diag_dble, post_smooth_iter)
    else
        if( prol_prec == 2) then
            call mat_vec_mult_add_crs_dble( multigrid( i - 1) % cells( &
                3), multigrid( i) % P_row_ptr, size( multigrid( i) % &
                P_col_ind), multigrid( i) % P_col_ind, multigrid( i) % &
                P_val_dble, multigrid( i) % cells( 3), multigrid( i) % &
                sol_dble, multigrid( i - 1) % sol_dble)
            if( post_smooth_prec == 1) then
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_sngl( j) = sngl( multigrid( i &
                        - 1) % sol_dble( j))
                end do
            end if
        else
            call mat_vec_mult_add_crs_sngl( multigrid( i - 1) % cells( &
                3), multigrid( i) % P_row_ptr, size( multigrid( i) % &
                P_col_ind), multigrid( i) % P_col_ind, multigrid( i) % &
                P_val_sngl, multigrid( i) % cells( 3), multigrid( i) % &
                sol_sngl, multigrid( i - 1) % sol_sngl)
            if( post_smooth_prec == 2) then
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_dble( j) = dble( multigrid( i &
                        - 1) % sol_sngl( j))
                end do
            end if
        end if
        if( (alpha_dble == 0.0d+0) .and. (remov_null == 'grid')) then
            if( post_smooth_prec == 2) then
                proj = sum( multigrid( i - 1) % sol_dble) / multigrid( i &
                    - 1) % cells( 3)
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_dble( j) = multigrid( i - 1) &
                        % sol_dble( j) - proj
                end do
            else if( post_smooth_prec == 1) then
                proj_sngl = sum( multigrid( i - 1) % sol_sngl) / &
                    multigrid( i - 1) % cells( 3)
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_sngl( j) = multigrid( i - 1) &
                        % sol_sngl( j) - proj_sngl
                end do
            end if
        end if
        if( post_smooth_prec == 2) then
            call GS_dble( multigrid( i - 1) % cells( 1), multigrid( i - &
                1) % cells( 2), multigrid( i - 1) % A_row_ptr, size( &
                multigrid( i - 1) % A_col_ind), multigrid( i - 1) % &
                A_col_ind, multigrid( i - 1) % A_val_dble, multigrid( i - &
                1) % sol_dble, multigrid( i - 1) % rhs_dble, GS_type, &
                multigrid( i - 1) % GS_rhs_dble, multigrid( i - 1) % &
                GS_up_diag_dble, pre_smooth_iter)
            if( prol_prec == 1) then
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_sngl( j) = sngl( multigrid( i &
                        - 1) % sol_dble( j))
                end do
            end if
        else
            call GS_sngl( multigrid( i - 1) % cells( 1), multigrid( i - &
                1) % cells( 2), multigrid( i - 1) % A_row_ptr, size( &
                multigrid( i - 1) % A_col_ind), multigrid( i - 1) % &
                A_col_ind, multigrid( i - 1) % A_val_sngl, multigrid( i - &
                1) % sol_sngl, multigrid( i - 1) % rhs_sngl, GS_type, &
                multigrid( i - 1) % GS_rhs_sngl, multigrid( i - 1) % &
                GS_up_diag_sngl, pre_smooth_iter)
            if( prol_prec == 2) then
                do j = 1, multigrid( i - 1) % cells( 3)
                    multigrid( i - 1) % sol_dble( j) = dble( multigrid( i &
                        - 1) % sol_sngl( j))
                end do
            end if
        end if
    end if
end do

end subroutine v_cyc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pressure_gmg_intern_proced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Implements the public interface of the pressure solver.
module gmg
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sets up the multigrid data structure for alpha = 0 for the specified grid parameters and equation coefficients.
subroutine pressure_gmg_init( y_fine_cells, z_fine_cells, &
    y_cell_edge_coords, z_cell_edge_coords, p_cell_edge, mu_discr, &
    rho_discr, mu_is_const, rad_reduc_factor)
use pressure_gmg_intern_data
use pressure_gmg_intern_proced, only: comput_lvls, discr, coarsen

! Declares the input variables.
! Declares the number of cells on the fine-grid in the y and the z-directions.
integer, intent( in) :: y_fine_cells, z_fine_cells
! Declares the y and z-coordinates of the edges of the fine-grid cells. The arrays will be overwritten with the values on coarser grids.
double precision, intent( inout) :: y_cell_edge_coords( y_fine_cells + &
    1), z_cell_edge_coords( z_fine_cells + 1)
! Declares the background pressures at the edges of the fine-grid cells in the z-direction. The arrays will be overwritten with the values on coarser grids.
double precision, intent( inout) :: p_cell_edge( z_fine_cells + 1)
! Declares a flag for the coefficient mu. If .false., mu is computed from the standard formula. If .true., all mu are set to 1 in order to compare the GMG solution to the FFT solution.
character( len = 1), intent( in) :: mu_discr
character( len = 1), intent( in) :: rho_discr
logical, intent( in) :: mu_is_const
! Declares the reduction factor for the radius of the earth, e.g. if rad_reduc_factor is 2, the radius is halved.
double precision, intent( in) :: rad_reduc_factor

! Declares the local variables.
integer :: i, j
! Declares the number of grids of multigrid.
integer :: lvls
! Declares the level of the current grid.
integer :: lvl
! Declares the precisions to be used for pre- and post-smoothing (except on the fine-grid), computation of the residual (except on the fine-grid), transfer operations (restriction and prolongation), and coarse-grid solves in the non-zero and zero alpha cases. If 1, it is single, if 2, it is double.
integer :: pre_smooth_prec, resid_prec, restr_prec, coarse_prec, &
    zero_coarse_prec, post_smooth_prec, prol_prec
! Declares the desired number of cells on one side of the coarse-grid.
integer :: coarse_size
! Declares a flag for whether to semi-coarsen. The routine will semi-coarsen until both sides are between coarse_size and 0.5 * coarse_size + 1. The purpose of semi-coarsening is to reduce the linear system on the coarse-grid to a manageable size. It would otherwise be large, since y_fine_cells >> z_fine_cells.
logical :: semi_coarsen
! Declares the widths of the fine-grid cells in the y and z-directions. The arrays will be overwritten with the values on coarser grids.
double precision, allocatable :: y_cell_widths( :), z_cell_widths( :)
! Declares the coordinates of the finer and coarser grid cells in the y and z-directions. The arrays will be overwritten with the values on coarser grids.
double precision, allocatable :: y_finer_cell_cent_coords( :), &
    z_finer_cell_cent_coords( :), y_coarser_cell_cent_coords( :), &
    z_coarser_cell_cent_coords( :)
! Declares the coarse-grid discretization matrix in dense format.
double precision, allocatable :: QR_coarse( :, :), tau_coarse( :), &
    T_coarse( :, :)
real, allocatable :: QR_coarse_sngl( :, :), tau_coarse_sngl( :), &
    T_coarse_sngl( :, :)
! double precision, allocatable :: A_coarse( :, :)
! real, allocatable :: A_coarse_sngl( :, :)
! Declares the inputs of the LAPACK functions.
double precision, allocatable :: DGEQP3_WORK( :), DLARFB_WORK( :, :)
! double precision, allocatable :: DGESVD_WORK( :)
! double precision, allocatable :: DGEJSV_WORK( :)
real, allocatable :: SGEQP3_WORK( :), SLARFB_WORK( :, :)
! real, allocatable :: SGESVD_WORK( :)
! real, allocatable :: SGEJSV_WORK( :)
! integer, allocatable :: GEJSV_WORK( :)
integer :: INFO

pre_smooth_prec = 2
resid_prec = 2
restr_prec = 2
coarse_prec = 2
zero_coarse_prec = 2
post_smooth_prec = 2
prol_prec = 2

coarse_size = 8
semi_coarsen = .true.
! Computes the number of grids.
lvls = comput_lvls( y_fine_cells, z_fine_cells, coarse_size, semi_coarsen)

! Allocates the array of grid variables.
allocate( multigrid( lvls))

! Allocates the array for the number of cells on the fine-grid and computes them.
allocate( multigrid( 1) % cells( 3))
multigrid( 1) % cells( 1) = y_fine_cells
multigrid( 1) % cells( 2) = z_fine_cells
multigrid( 1) % cells( 3) = multigrid( 1) % cells( 1) * multigrid( 1) % &
    cells( 2)

! Allocates arrays for the widths of the cells in the y and z-directions.
allocate( y_cell_widths( multigrid( 1) % cells( 1)), &
    z_cell_widths( multigrid( 1) % cells( 2)))
! Allocates arrays for the coordinates of the centers of the finer and coarser cells in the y and z-directions.
allocate( y_finer_cell_cent_coords( multigrid( 1) % cells( 1)), &
    z_finer_cell_cent_coords( multigrid( 1) % cells( 2)), &
    y_coarser_cell_cent_coords( multigrid( 1) % cells( 1)), &
    z_coarser_cell_cent_coords( multigrid( 1) % cells( 2)))

! Computes the widths of the fine-grid cells.
do i = 1, multigrid( 1) % cells( 1)
    y_cell_widths( i) = y_cell_edge_coords( i + 1) - y_cell_edge_coords( i)
end do
do i = 1, multigrid( 1) % cells( 2)
    z_cell_widths( i) = z_cell_edge_coords( i + 1) - z_cell_edge_coords( i)
end do
! Computes the coordinates of the centers of the fine-grid cells.
do i = 1, multigrid( 1) % cells( 1)
    y_coarser_cell_cent_coords( i) = 0.5d+0 * (y_cell_edge_coords( i) + &
    y_cell_edge_coords( i + 1))
end do
do i = 1, multigrid( 1) % cells( 2)
    z_coarser_cell_cent_coords( i) = 0.5d+0 * (z_cell_edge_coords( i) + &
    z_cell_edge_coords( i + 1))
end do

! Loops through each grid.
do lvl = 1, lvls
    ! Allocates arrays (in the routine) for the PDE discretization matrix for alpha = 0 and computes them.
    call discr( lvl, y_cell_edge_coords, p_cell_edge, &
        y_cell_widths, z_cell_widths, y_coarser_cell_cent_coords, &
        mu_discr, rho_discr, mu_is_const, rad_reduc_factor)
    ! Allocates the solution and right-hand-side vectors, except on the fine-grid.
    if( lvl == 1) then
        allocate( multigrid( lvl) % sol_sngl( multigrid( lvl) % cells( 3)))
    else
        allocate( multigrid( lvl) % sol_dble( multigrid( lvl) % cells( &
            3)), multigrid( lvl) % sol_sngl( multigrid( lvl) % cells( &
            3)), multigrid( lvl) % rhs_dble( multigrid( lvl) % cells( &
            3)), multigrid( lvl) % rhs_sngl( multigrid( lvl) % cells( 3)))
    end if
    ! For every grid except the coarse-grid:
    if( lvl < lvls) then
        ! Allocates the buffers for the smoother and the tridiagonal solver.
        allocate( multigrid( lvl) % GS_rhs_dble( multigrid( lvl) % cells( &
            2)), multigrid( lvl) % GS_rhs_sngl( multigrid( lvl) % cells( &
            2)), multigrid( lvl) % GS_up_diag_dble( multigrid( lvl) % &
            cells( 2) - 1), multigrid( lvl) % GS_up_diag_sngl( multigrid( &
            lvl) % cells( 2) - 1))
        ! Allocates the residual vector.
        allocate( multigrid( lvl) % res_dble( multigrid( lvl) % cells( &
            3)), multigrid( lvl) % res_sngl( multigrid( lvl) % cells( 3)))
        ! Allocates (in the routine) the array for the number of cells on the next coarser grid and computes it.
        ! Computes the widths of the coarser grid cells in the y and z-directions.
        ! Computes the coordinates of the edges of the coarser grid cells in the y and z-directions.
        ! Computes the background pressures at the edges of the coarser grid cells in the z-direction.
        ! Computes the coordinates of the centers of the coarser grid cells in the y and z-directions.
        ! Allocates (in the routine) the restriction operator to the next coarser grid and computes it.
        ! Allocates (in the routine) the prolongation operator from the next coarser grid and computes it.
        call coarsen( lvl, y_cell_edge_coords, z_cell_edge_coords, &
            p_cell_edge, y_cell_widths, z_cell_widths, &
            y_finer_cell_cent_coords, z_finer_cell_cent_coords, &
            y_coarser_cell_cent_coords, z_coarser_cell_cent_coords, &
            coarse_size, restr_prec, prol_prec)
    end if
end do

if( coarse_prec == 2) then
    ! Allocates the LU decomposition (with partial pivoting) of the dense coarse-grid discretization matrix.
    allocate( LU_coarse_dble( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)), p_coarse( multigrid( lvls) % cells( 3)))
    ! Allocates the permuted right-hand-side on the coarse-grid.
    allocate( p_rhs_coarse_dble( multigrid( lvls) % cells( 3)))
else
    allocate( LU_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)), p_coarse( multigrid( lvls) % cells( 3)))
    allocate( p_rhs_coarse_sngl( multigrid( lvls) % cells( 3)))
end if
allocate( res_old( multigrid( 1) % cells( 3)), rhs_A_trans( multigrid( 1) &
    % cells( 3)), res_A_trans( multigrid( 1) % cells( 3)))
if( zero_coarse_prec == 2) then
    allocate( QR_coarse( multigrid( lvls) % cells( 3), multigrid( lvls) % &
        cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            QR_coarse( i, j) = 0.0d+0
        end do
    end do
    do i = 1, multigrid( lvls) % cells( 3)
        do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
            A_row_ptr( i + 1) - 1
            QR_coarse( i, multigrid( lvls) % A_col_ind( j)) = multigrid( &
                lvls) % A_val_dble( j)
        end do
    end do
    allocate( e_coarse( multigrid( lvls) % cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        e_coarse( i) = 0
    end do
    allocate( tau_coarse( multigrid( lvls) % cells( 3)), DGEQP3_WORK( 6 * &
        multigrid( lvls) % cells( 3)))
    call dgeqp3( multigrid( lvls) % cells( 3), multigrid( lvls) % cells( &
        3), QR_coarse, multigrid( lvls) % cells( 3), e_coarse, &
        tau_coarse, DGEQP3_WORK, 6 * multigrid( lvls) % cells( 3), INFO)
    allocate( R_coarse_dble( multigrid( lvls) % cells( 3) - 1, multigrid( &
        lvls) % cells( 3) - 1))
    do i = 1, multigrid( lvls) % cells( 3) - 1
        do j = i, multigrid( lvls) % cells( 3) - 1
            R_coarse_dble( i, j) = QR_coarse( i, j)
        end do
    end do
    allocate( T_coarse( multigrid( lvls) % cells( 3), multigrid( lvls) % &
        cells( 3)))
    call dlarft( 'F', 'C', multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3), QR_coarse, multigrid( lvls) % cells( 3), &
        tau_coarse, T_coarse, multigrid( lvls) % cells( 3))
    allocate( Q_coarse_dble( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            if( i == j) then
                Q_coarse_dble( i, j) = 1.0d+0
            else
                Q_coarse_dble( i, j) = 0.0d+0
            end if
        end do
    end do
    allocate( DLARFB_WORK( multigrid( lvls) % cells( 3), multigrid( lvls) &
        % cells( 3)))
    call dlarfb( 'L', 'N', 'F', 'C', multigrid( lvls) % cells( 3), &
        multigrid( lvls) % cells( 3), multigrid( lvls) % cells( 3), &
        QR_coarse, multigrid( lvls) % cells( 3), T_coarse, multigrid( &
        lvls) % cells( 3), Q_coarse_dble, multigrid( lvls) % cells( 3), &
        DLARFB_WORK, multigrid( lvls) % cells( 3))
    allocate( QR_rhs_dble( multigrid( lvls) % cells( 3) - 1))
    ! ! Allocates the SVD of the dense coarse-grid discretization matrix.
    ! allocate( U_coarse( multigrid( lvls) % cells( 3), multigrid( lvls) % &
    !     cells( 3)), S_coarse( multigrid( lvls) % cells( 3)), VT_coarse( &
    !     multigrid( lvls) % cells( 3), multigrid( lvls) % cells( 3)))
    ! ! Makes a copy of the coarse-grid discretization matrix in dense format.
    ! allocate( A_coarse( multigrid( lvls) % cells( 3), multigrid( lvls) % &
    !     cells( 3)))
    ! do i = 1, multigrid( lvls) % cells( 3)
    !     do j = 1, multigrid( lvls) % cells( 3)
    !         A_coarse( i, j) = 0.0d+0
    !     end do
    ! end do
    ! do i = 1, multigrid( lvls) % cells( 3)
    !     do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
    !         A_row_ptr( i + 1) - 1
    !         A_coarse( i, multigrid( lvls) % A_col_ind( j)) = multigrid( &
    !             lvls) % A_val_dble( j)
    !     end do
    ! end do
    ! ! Allocates the input arrays to d/sgesvd.
    ! allocate( DGESVD_WORK( 10 * multigrid( lvls) % cells( 3)))
    ! ! Computes the SVD using the LAPACK function dgesvd (http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html).
    ! call dgesvd( 'A', 'A', multigrid( lvls) % cells( 3), multigrid( lvls) &
    !     % cells( 3), A_coarse, multigrid( lvls) % cells( 3), S_coarse, &
    !     U_coarse, multigrid( lvls) % cells( 3), VT_coarse, multigrid( &
    !     lvls) % cells( 3), DGESVD_WORK, 10 * multigrid( lvls) % cells( &
    !     3), INFO)
    ! allocate( DGEJSV_WORK( 6 * multigrid( lvls) % cells( 3) * multigrid( &
    !     lvls) % cells( 3)))
    ! allocate( GEJSV_WORK( 4 * multigrid( lvls) % cells( 3)))
    ! call dgejsv( 'C', 'U', 'V', 'N', 'T', 'N', multigrid( lvls) % cells( &
    !     3), multigrid( lvls) % cells( 3), A_coarse, multigrid( lvls) % &
    !     cells( 3), S_coarse, U_coarse, multigrid( lvls) % cells( 3), &
    !     VT_coarse, multigrid( lvls) % cells( 3), DGEJSV_WORK, 6 * &
    !     multigrid( lvls) % cells( 3) * multigrid( lvls) % cells( 3), &
    !     GEJSV_WORK, INFO)
    ! VT_coarse = transpose( VT_coarse)
else
    allocate( QR_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            QR_coarse_sngl( i, j) = 0.0e+0
        end do
    end do
    do i = 1, multigrid( lvls) % cells( 3)
        do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
            A_row_ptr( i + 1) - 1
            QR_coarse_sngl( i, multigrid( lvls) % A_col_ind( j)) = &
                multigrid( lvls) % A_val_sngl( j)
        end do
    end do
    allocate( e_coarse( multigrid( lvls) % cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        e_coarse( i) = 0
    end do
    allocate( tau_coarse_sngl( multigrid( lvls) % cells( 3)), &
        SGEQP3_WORK( 6 * multigrid( lvls) % cells( 3)))
    call sgeqp3( multigrid( lvls) % cells( 3), multigrid( lvls) % cells( &
        3), QR_coarse_sngl, multigrid( lvls) % cells( 3), e_coarse, &
        tau_coarse_sngl, SGEQP3_WORK, 6 * multigrid( lvls) % cells( 3), &
        INFO)
    allocate( R_coarse_sngl( multigrid( lvls) % cells( 3) - 1, multigrid( &
        lvls) % cells( 3) - 1))
    do i = 1, multigrid( lvls) % cells( 3) - 1
        do j = i, multigrid( lvls) % cells( 3) - 1
            R_coarse_sngl( i, j) = QR_coarse_sngl( i, j)
        end do
    end do
    allocate( T_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)))
    call slarft( 'F', 'C', multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3), QR_coarse_sngl, multigrid( lvls) % cells( &
        3), tau_coarse_sngl, T_coarse_sngl, multigrid( lvls) % cells( 3))
    allocate( Q_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
        lvls) % cells( 3)))
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            if( i == j) then
                Q_coarse_sngl( i, j) = 1.0e+0
            else
                Q_coarse_sngl( i, j) = 0.0e+0
            end if
        end do
    end do
    allocate( SLARFB_WORK( multigrid( lvls) % cells( 3), multigrid( lvls) &
        % cells( 3)))
    call slarfb( 'L', 'N', 'F', 'C', multigrid( lvls) % cells( 3), &
        multigrid( lvls) % cells( 3), multigrid( lvls) % cells( 3), &
        QR_coarse_sngl, multigrid( lvls) % cells( 3), T_coarse_sngl, &
        multigrid( lvls) % cells( 3), Q_coarse_sngl, multigrid( lvls) % &
        cells( 3), SLARFB_WORK, multigrid( lvls) % cells( 3))
    allocate( QR_rhs_sngl( multigrid( lvls) % cells( 3) - 1))
    ! allocate( U_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
    !     lvls) % cells( 3)), S_coarse_sngl( multigrid( lvls) % cells( 3)), &
    !     VT_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( lvls) % &
    !     cells( 3)))
    ! allocate( A_coarse_sngl( multigrid( lvls) % cells( 3), multigrid( &
    !     lvls) % cells( 3)))
    ! do i = 1, multigrid( lvls) % cells( 3)
    !     do j = 1, multigrid( lvls) % cells( 3)
    !         A_coarse_sngl( i, j) = 0.0e+0
    !     end do
    ! end do
    ! do i = 1, multigrid( lvls) % cells( 3)
    !     do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
    !         A_row_ptr( i + 1) - 1
    !         A_coarse_sngl( i, multigrid( lvls) % A_col_ind( j)) = sngl( &
    !             multigrid( lvls) % A_val_sngl( j))
    !     end do
    ! end do
    ! allocate( SGESVD_WORK( 10 * multigrid( lvls) % cells( 3)))
    ! ! Computes the SVD using the LAPACK function sgesvd (http://www.netlib.org/lapack/explore-html/d4/dca/group__real_g_esing_gaf03d06284b1bfabd3d6c0f6955960533.html).
    ! call sgesvd( 'A', 'A', multigrid( lvls) % cells( 3), multigrid( lvls) &
    !     % cells( 3), A_coarse_sngl, multigrid( lvls) % cells( 3), &
    !     S_coarse_sngl, U_coarse_sngl, multigrid( lvls) % cells( 3), &
    !     VT_coarse_sngl, multigrid( lvls) % cells( 3), SGESVD_WORK, 10 * multigrid( lvls) % cells( 3), INFO)
    ! allocate( SGEJSV_WORK( 6 * multigrid( lvls) % cells( 3) * multigrid( &
    !     lvls) % cells( 3)))
    ! allocate( GEJSV_WORK( 4 * multigrid( lvls) % cells( 3)))
    ! call sgejsv( 'C', 'U', 'V', 'N', 'T', 'N', multigrid( lvls) % cells( &
    !     3), multigrid( lvls) % cells( 3), A_coarse_sngl, multigrid( lvls) &
    !     % cells( 3), S_coarse_sngl, U_coarse_sngl, multigrid( lvls) % &
    !     cells( 3), VT_coarse_sngl, multigrid( lvls) % cells( 3), &
    !     SGEJSV_WORK, 6 * multigrid( lvls) % cells( 3) * multigrid( lvls) &
    !     % cells( 3), GEJSV_WORK, INFO)
    ! VT_coarse_sngl = transpose( VT_coarse_sngl)
end if

! Deallocates the arrays for the widths of the cells, the coordinates of the centers of the cells, and the input to d/sgesvd, since they will not be needed after the setup stage.
deallocate( y_cell_widths, z_cell_widths, y_finer_cell_cent_coords, &
    z_finer_cell_cent_coords, y_coarser_cell_cent_coords, &
    z_coarser_cell_cent_coords)
if( zero_coarse_prec == 2) then
    deallocate( QR_coarse, tau_coarse, T_coarse, DGEQP3_WORK, DLARFB_WORK)
    ! deallocate( A_coarse, DGESVD_WORK)
    ! deallocate( A_coarse, DGEJSV_WORK, GEJSV_WORK)
else
    deallocate( QR_coarse_sngl, tau_coarse_sngl, T_coarse_sngl, &
        SGEQP3_WORK, SLARFB_WORK)
    ! deallocate( A_coarse_sngl, SGESVD_WORK)
    ! deallocate( A_coarse_sngl, SGEJSV_WORK, GEJSV_WORK)
end if

end subroutine pressure_gmg_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes the multigrid solution for the given parameters.
subroutine pressure_gmg_solv( alpha_dble, rhs_dble, tol, max_v_cyc, &
    sol_dble, res_norm, iter, print_res)
use pressure_gmg_intern_data
use gmg_common_intern_proced, only: trans_mat_vec_mult_crs, resid_dble
use pressure_gmg_intern_proced, only: v_cyc

! Declares the input variables.
! Declares the new value of the harmonic.
double precision, intent( in) :: alpha_dble
! Declares the right-hand-side.
double precision, intent( in) :: rhs_dble( multigrid( 1) % cells( 3))
! Declares the tolerance for the residual in the relative 2-norm.
double precision, intent( in) :: tol
! Declares the maximum number of iterations (v-cycles) for multigrid to use.
integer, intent( in) :: max_v_cyc
! Declares the initial guess for the solution. The array will be overwritten with the multigrid solution.
double precision, intent( inout) :: sol_dble( multigrid( 1) % cells( 3))
! Declares the relative 2-norm of the residual of the multigrid solution. It will be overwritten.
double precision, intent( out) :: res_norm
! Declares the number of iterations (v-cycles) performed by multigrid. It will be overwritten.
integer, intent( out) :: iter
! Declares a flag for whether to print the residual after each iteration (v-cycle).
logical :: print_res

! Declares the local variables.
integer :: i, j
! Declares the number of grids of multigrid.
integer :: lvls
! Declares the level of the current grid.
integer :: lvl
! Declares the new value of the harmonic.
real :: alpha_sngl
! Declares the type of smoother: 'point', 'z-line', 'zebra z-line'.
character( len = 5) :: GS_type
! Declares the number of smoothing iterations on each grid on the down-stroke.
integer :: pre_smooth_iter
! Declares the number of smoothing iterations on each grid on the up-stroke.
integer :: post_smooth_iter
! Declares the intervals over which to remove the right nullspace.
character( len = 4) :: remov_null
! Declares the precisions to be used for pre- and post-smoothing (except on the fine-grid), computation of the residual (except on the fine-grid), transfer operations (restriction and prolongation), and coarse-grid solves in the non-zero and zero alpha cases. If 1, it is single, if 2, it is double.
integer :: pre_smooth_prec, resid_prec, restr_prec, coarse_prec, &
    zero_coarse_prec, post_smooth_prec, prol_prec
! Declares the 2-norm of the right-hand-side.
double precision :: rhs_norm, trans_rhs_norm
double precision :: trans_res_norm
! Declares the minimum value of the relative 2-norm of the residual attained in previous v-cycles.
double precision :: res_norm_old, trans_res_norm_old
! Declares the count rate for the system clock, a starting time, and an ending time.
integer :: clock_rate, start_time, end_time
! Declares the elapsed time.
double precision :: elapsed_time
! Declares the component of the solution vector in the right nullspace.
double precision :: proj
! Declares the inputs of d/sgesv.
integer :: INFO

lvls = size( multigrid)

! GS_type = 'point'
GS_type = 'zline'
! GS_type = 'zebra'

pre_smooth_iter = 1
pre_smooth_iter = 1

! remov_null = 'grid'
remov_null = 'cycl'
! remov_null = 'mult'
! remov_null = 'none'

pre_smooth_prec = 2
resid_prec = 2
restr_prec = 2
coarse_prec = 2
zero_coarse_prec = 2
post_smooth_prec = 2
prol_prec = 2

alpha_sngl = sngl( alpha_dble)
! Updates the discretization matrices on each grid by modifying their diagonal entries.
do lvl = 1, lvls
    do i = 1, multigrid( lvl) % cells( 3)
        multigrid( lvl) % A_val_dble( multigrid( lvl) % A_row_ptr( i)) = &
            multigrid( lvl) % A_diag_val_dble( i) - alpha_dble
    end do
    do i = 1, multigrid( lvl) % cells( 3)
        multigrid( lvl) % A_val_sngl( multigrid( lvl) % A_row_ptr( i)) = &
            multigrid( lvl) % A_diag_val_sngl( i) - alpha_sngl
    end do
end do
if( coarse_prec == 2) then
    ! Makes a copy of the coarse-grid discretization matrix in dense format.
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            LU_coarse_dble( i, j) = 0.0d+0
        end do
    end do
    do i = 1, multigrid( lvls) % cells( 3)
        do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
            A_row_ptr( i + 1) - 1
            LU_coarse_dble( i, multigrid( lvls) % A_col_ind( j)) = &
                multigrid( lvls) % A_val_dble( j)
        end do
    end do
    ! Computes the LU decomposition (with partial pivoting) in-situ.
    ! call LU_part_pivot_dble( multigrid( lvls) % cells( 3), LU_coarse, &
    !     p_coarse)
    call dgetf2( multigrid( lvls) % cells( 3), multigrid( lvls) % cells( &
        3), LU_coarse_dble, multigrid( lvls) % cells( 3), p_coarse, INFO)
else
    do i = 1, multigrid( lvls) % cells( 3)
        do j = 1, multigrid( lvls) % cells( 3)
            LU_coarse_sngl( i, j) = 0.0e+0
        end do
    end do
    do i = 1, multigrid( lvls) % cells( 3)
        do j = multigrid( lvls) % A_row_ptr( i), multigrid( lvls) % &
            A_row_ptr( i + 1) - 1
            LU_coarse_sngl( i, multigrid( lvls) % A_col_ind( j)) = &
                multigrid( lvls) % A_val_sngl( j)
        end do
    end do
    ! call LU_part_pivot_sngl( multigrid( lvls) % cells( 3), &
    !     LU_coarse_sngl, p_coarse)
    call sgetf2( multigrid( lvls) % cells( 3), multigrid( lvls) % cells( &
        3), LU_coarse_sngl, multigrid( lvls) % cells( 3), p_coarse, INFO)
end if

do i = 1, multigrid( 1) % cells( 3)
    multigrid( 1) % sol_sngl( i) = sngl( sol_dble( i))
end do

! Computes the 2-norm of the right-hand-side.
rhs_norm = sqrt( dot_product( rhs_dble, rhs_dble))
if( alpha_dble == 0.0d+0) then
    call trans_mat_vec_mult_crs( multigrid( 1) % cells( 3), &
        multigrid( 1) % A_row_ptr, size( multigrid( 1) % A_col_ind), &
        multigrid( 1) % A_col_ind, multigrid( 1) % A_val_dble, multigrid( &
        1) % cells( 3), rhs_dble, rhs_A_trans)
end if
trans_rhs_norm = sqrt( dot_product( rhs_A_trans, rhs_A_trans))

! If the right-hand-side is numerically 0, then the routine sets the solution as the zero vector and returns.
if( rhs_norm == 0.0d+0) then
    do i = 1, multigrid( 1) % cells( 3)
        sol_dble( i) = 0.0d+0
    end do
    iter = 0
    res_norm = 0.0d+0
    return
end if

! Performs v-cycles.
iter = 0
res_norm = 1.0d+5
trans_res_norm = 1.0d+5
res_norm_old = res_norm
trans_res_norm_old = trans_res_norm

call system_clock( count_rate = clock_rate)
call system_clock( start_time)

if( print_res) then
    print 1, alpha_dble
    1 format( 'alpha = ', es22.15, '.')
    print *, 'The residuals in the relative 2-norm were:'
    if( alpha_dble == 0.0d+0) then
        print *, '(diff) = the difference (which lies in range( A)) &
            &of successive residuals'
        print *, '(A^T) = A^T of the residual (which lies in domain( &
            &A))'
    end if
end if

do while( iter < max_v_cyc)
    call v_cyc( alpha_dble, 1, sol_dble, rhs_dble, GS_type, &
        pre_smooth_iter, post_smooth_iter, remov_null, pre_smooth_prec, &
        resid_prec, restr_prec, coarse_prec, zero_coarse_prec, &
        post_smooth_prec, prol_prec)

    if( (alpha_dble == 0.0d+0) .and. (remov_null == 'cycl')) then
        ! Removes the right nullspace component of the solution.
        proj = sum( sol_dble) / multigrid( 1) % cells( 3)
        do i = 1, multigrid( 1) % cells( 3)
            sol_dble( i) = sol_dble( i) - proj
        end do
    end if

    iter = iter + 1

    ! Computes the new residual and its relative 2-norm.
    call resid_dble( multigrid( 1) % cells( 3), multigrid( 1) % &
        A_row_ptr, size( multigrid( 1) % A_col_ind), multigrid( 1) % &
        A_col_ind, multigrid( 1) % A_val_dble, multigrid( 1) % cells( 3), &
        sol_dble, rhs_dble, multigrid( 1) % res_dble)
    if( (alpha_dble /= 0.0d+0) .or. (iter == 1)) then
        res_norm = sqrt( dot_product( multigrid( 1) % res_dble, &
            multigrid( 1) % res_dble)) / rhs_norm
    else
        do i = 1, multigrid( 1) % cells( 3)
            res_old( i) = multigrid( 1) % res_dble( i) - res_old( i)
        end do
        res_norm = sqrt( dot_product( res_old, res_old)) / rhs_norm
    end if
    if( alpha_dble == 0.0d+0) then
        call trans_mat_vec_mult_crs( multigrid( 1) % cells( 3), &
            multigrid( 1) % A_row_ptr, size( multigrid( 1) % A_col_ind), &
            multigrid( 1) % A_col_ind, multigrid( 1) % A_val_dble, &
            multigrid( 1) % cells( 3), multigrid( 1) % res_dble, &
            res_A_trans)
        trans_res_norm = sqrt( dot_product( res_A_trans, res_A_trans)) / &
            trans_rhs_norm
    end if

    ! Prints the residual after each iteration for testing purposes.
    if( print_res) then
        if( alpha_dble /= 0.0d+0) then
            print 2, iter, res_norm
            2 format( 'iter ', i3, ': ', es22.15, '.')
        ! The transposed residual is not printed for the first iteration in the case of zero alpha, since otherwise optimization by both ifort and gfortran will increase the number of iterations for all alpha. The reason is unknown.
        elseif( iter == 1) then
            print 3, iter, res_norm
            3 format( 'iter ', i3, ': ', es22.15, '.')
        else
            print 4, iter, res_norm, trans_res_norm
            4 format( 'iter ', i3, ': ', es22.15, ' (diff), ', es22.15, &
                ' (A^T).')
        end if
    end if

    if( res_norm < tol) then
        if( print_res) then
            call system_clock( end_time)
            elapsed_time = end_time - start_time
            print 5, elapsed_time / clock_rate
            5 format( 'GMG reached the desired tolerance after ', es9.2, &
                ' seconds.')
        end if
        if( (alpha_dble == 0.0d+0) .and. (remov_null == 'mult')) then
            proj = sum( sol_dble) / multigrid( 1) % cells( 3)
            do i = 1, multigrid( 1) % cells( 3)
                sol_dble( i) = sol_dble( i) - proj
            end do
        end if
        return
    end if

    if( (res_norm > res_norm_old * 2) .and. ((alpha_dble /= 0.0d+0) .or. &
        (trans_res_norm > trans_res_norm_old * 2))) then
        if( print_res) then
            call system_clock( end_time)
            elapsed_time = end_time - start_time
            print 6, elapsed_time / clock_rate
            6 format( 'GMG stopped early after ', es9.2, ' seconds and &
                &appears to be diverging.')
        end if
        if( (alpha_dble == 0.0d+0) .and. (remov_null == 'mult')) then
            proj = sum( sol_dble) / multigrid( 1) % cells( 3)
            do i = 1, multigrid( 1) % cells( 3)
                sol_dble( i) = sol_dble( i) - proj
            end do
        end if
        return
    else if( ((res_norm / res_norm_old) > 0.92) .and. ((alpha_dble /= &
        0.0d+0) .or. ((trans_res_norm / trans_res_norm_old) > 0.92))) then
        if( print_res) then
            call system_clock( end_time)
            elapsed_time = end_time - start_time
            print 7, elapsed_time / clock_rate
            7 format( 'GMG stopped early after ', es9.2, ' seconds and &
                &appears to be stagnating.')
        end if
        if( (alpha_dble == 0.0d+0) .and. (remov_null == 'mult')) then
            proj = sum( sol_dble) / multigrid( 1) % cells( 3)
            do i = 1, multigrid( 1) % cells( 3)
                sol_dble( i) = sol_dble( i) - proj
            end do
        end if
        return
    end if

    res_norm_old = min( res_norm, res_norm_old)
    trans_res_norm_old = min( trans_res_norm, trans_res_norm_old)
    if( alpha_dble == 0.0d+0) then
        do i = 1, multigrid( 1) % cells( 3)
            res_old( i) = multigrid( 1) % res_dble( i)
        end do
    end if
end do

if( print_res) then
    call system_clock( end_time)
    elapsed_time = end_time - start_time
    print 8, elapsed_time / clock_rate
    8 format( 'GMG reached the maximum number of iterations after ', &
        es9.2, ' seconds.')
end if
if( (alpha_dble == 0.0d+0) .and. (remov_null == 'mult')) then
    proj = sum( sol_dble) / multigrid( 1) % cells( 3)
    do i = 1, multigrid( 1) % cells( 3)
        sol_dble( i) = sol_dble( i) - proj
    end do
end if

end subroutine pressure_gmg_solv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Frees memory associated with the multigrid data structure.
subroutine pressure_gmg_fin
use pressure_gmg_intern_data

! Declares the local variables.
integer :: i
! Declares the precisions to be used for pre- and post-smoothing (except on the fine-grid), computation of the residual (except on the fine-grid), transfer operations (restriction and prolongation), and coarse-grid solves in the non-zero and zero alpha cases. If 1, it is single, if 2, it is double.
integer :: pre_smooth_prec, resid_prec, restr_prec, coarse_prec, &
    zero_coarse_prec, post_smooth_prec, prol_prec

pre_smooth_prec = 2
resid_prec = 2
restr_prec = 2
coarse_prec = 2
zero_coarse_prec = 2
post_smooth_prec = 2
prol_prec = 2

! First deallocates pointers in each individual grid variable.
do i = 1, size( multigrid)
    if( i == 1) then
        deallocate( multigrid( i) % sol_sngl)
    end if
    if( i > 1) then
        deallocate( multigrid( i) % sol_dble, multigrid( i) % sol_sngl, &
            multigrid( i) % rhs_dble, multigrid( i) % rhs_sngl, &
            multigrid( i) % P_row_ptr, multigrid( i) % P_col_ind)
        if( prol_prec == 2) then
            deallocate( multigrid( i) % P_val_dble)
        else
            deallocate( multigrid( i) % P_val_sngl)
        end if
    end if
    deallocate( multigrid( i) % cells, multigrid( i) % A_row_ptr, &
        multigrid( i) % A_col_ind, multigrid( i) % A_val_dble, multigrid( &
        i) % A_val_sngl, multigrid( i) % A_diag_val_dble, multigrid( i) % &
        A_diag_val_sngl)
    if( i < size( multigrid)) then
        deallocate( multigrid( i) % GS_rhs_dble, multigrid( i) % &
            GS_rhs_sngl, multigrid( i) % GS_up_diag_dble, multigrid( i) % &
            GS_up_diag_sngl, multigrid( i) % res_dble, multigrid( i) % &
            res_sngl, multigrid( i) % R_row_ptr, multigrid( i) % &
            R_col_ind)
        if( restr_prec == 2) then
            deallocate( multigrid( i) % R_val_dble)
        else
            deallocate( multigrid( i) % R_val_sngl)
        end if
    end if
end do

! Then deallocates the grids array.
deallocate( multigrid)

if( coarse_prec == 2) then
    deallocate( LU_coarse_dble, p_coarse, p_rhs_coarse_dble)
else
    deallocate( LU_coarse_sngl, p_coarse, p_rhs_coarse_sngl)
end if
deallocate( res_old, rhs_A_trans, res_A_trans)
if( zero_coarse_prec == 2) then
    deallocate( Q_coarse_dble, R_coarse_dble, e_coarse, QR_rhs_dble)
    ! deallocate( U_coarse, S_coarse, VT_coarse)
else
    deallocate( Q_coarse_sngl, R_coarse_sngl, e_coarse, QR_rhs_sngl)
    ! deallocate( U_coarse_sngl, S_coarse_sngl, VT_coarse_sngl)
end if

end subroutine pressure_gmg_fin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module gmg
