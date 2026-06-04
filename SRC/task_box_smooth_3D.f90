subroutine box_smooth_3d(A,TERR,nx,ny,nzm,nsx,nsy,M,N,comm,kblock)

!-----------------------------------------------------------------------
!  SUBROUTINE: box_smooth_3d
!
!  PURPOSE
!  -------
!  Perform exact horizontal box averaging of a 3-D field distributed
!  over an MPI-decomposed 2-D horizontal domain, with optional masking
!  (e.g. terrain). The smoothing is applied independently at each
!  vertical level. The algorithm supports arbitrary smoothing window
!  sizes that may exceed individual subdomain dimensions, without using
!  halo exchanges or global gathers.
!
!   works only for MPI
!
!  ALGORITHM
!  ---------
!  1) For each vertical k-block:
!     a) Construct local prefix sums in x of the masked field:
!        X = A * TERR, W = TERR
!     b) Use MPI_Exscan along processor rows to form global x-prefix sums.
!     c) Construct local prefix sums in y to form the local summed-area
!        tables (SAT) SX and SW.
!     d) Use MPI_Exscan along processor columns to complete the global SAT.
!
!  2) Determine the minimal rectangular neighborhood of MPI subdomains
!     that contains all SAT points required to evaluate the box average
!     for every grid point in the local subdomain.
!
!  3) Exchange only the necessary SAT tiles with neighboring ranks
!     (using point-to-point communication), and assemble a local buffer
!     containing the required SAT values.
!
!  4) For each grid point (i,j,k), compute the box-averaged value using
!     four SAT lookups:
!
!        sumX = SX(I2,J2) - SX(I1-1,J2) - SX(I2,J1-1) + SX(I1-1,J1-1)
!        sumW = SW(I2,J2) - SW(I1-1,J2) - SW(I2,J1-1) + SW(I1-1,J1-1)
!
!     where I1,I2,J1,J2 define the horizontal averaging window. The output
!     is B = sumX / sumW for sumW > 0, otherwise zero.
!
!  5) Repeat for all k-blocks until nzm levels are processed.
!
!  This algorithm is exact, independent of subdomain size, and scales
!  efficiently when the smoothing window is small compared to the global
!  domain (M << nxg, N << nyg).
!
!  INTERFACE
!  ---------
!  Input/Output:
!    A(nx,ny,nzm)        : field to be smoothed (will be overwitten!)
!  Input:
!    TERR(nx,ny,nzm)     : Mask field; 1.0 for valid points, 0.0 for
!                          excluded points (e.g. inside terrain).
!    nx, ny              : Local subdomain dimensions.
!    nzm                 : Number of vertical levels.
!    nsx, nsy             : Number of subdomains in x and y directions.
!    M, N                 : Horizontal box size in x and y directions.
!    comm                 : MPI communicator spanning all subdomains.
!    kblock               : Number of vertical levels processed per block
!                          (controls memory usage and message sizes).
!
!  NOTES
!  -----
!  * The global domain size is nxg = nx*nsx, nyg = ny*nsy.
!  * Subdomains are assumed to be mapped line-by-line: rank increases
!    eastward first, then northward.
!  * The routine is independent of the relative sizes of M,N and nx,ny,
!    and does not rely on halo exchanges.
!  * All MPI datatypes must match the real kind used (e.g. MPI_REAL vs
!    MPI_DOUBLE_PRECISION when using compiler options such as -r8).
!
! (C) Marat Khairoutdinov 12/2025
!-----------------------------------------------------------------------

  use iso_fortran_env, only: real32, real64
  use mpi
  implicit none
  integer, intent(in) :: nx,ny,nzm,nsx,nsy,M,N,comm,kblock
  real   , intent(inout)  :: A(nx,ny,nzm)
  real   , intent(in)  :: TERR(nx,ny,nzm)

  integer :: ierr, rank, nprocs
  integer :: rx, ry, nxg, nyg, hx, hy
  integer :: kb, k1, k2, kk, k
  integer :: i,j
  integer :: row_comm, col_comm, color, key

  ! SAT for a k-block
  real, allocatable :: pxX(:,:,:), pxW(:,:,:)
  real, allocatable :: SX(:,:,:),  SW(:,:,:)

  ! Exscan buffers
  real, allocatable :: send_west_X(:,:), recv_west_X(:,:)
  real, allocatable :: send_west_W(:,:), recv_west_W(:,:)
  real, allocatable :: send_south_X(:,:), recv_south_X(:,:)
  real, allocatable :: send_south_W(:,:), recv_south_W(:,:)

  ! Neighborhood tiles of SAT (concatenated into big local arrays)
  integer :: rxL, rxR, ryB, ryT
  integer :: nRx, nRy
  integer :: nxN, nyN
  real, allocatable :: SXN(:,:,:), SWN(:,:,:)   ! (nxN, nyN, kb)

  integer :: Ig0, Jg0, I1, I2, J1, J2
  real :: sumX, sumW

  integer, parameter :: wp = kind(1.0)   ! will be 4 for normal real size , 8 with compiler ptomotion to 8
  integer  mpi_wp

  call t_startf ('task_box_smooth')


  if (wp == real32) then
    mpi_wp = MPI_REAL
  else if (wp == real64) then
    mpi_wp = MPI_REAL8
  else
    ! very rare; if you ever do real128 you need MPI_TYPE_MATCH_SIZE
    call MPI_Abort(comm, 99, ierr)
  end if

  kb = max(1, kblock)

  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  if (nprocs /= nsx*nsy) then
    if (rank==0) write(*,*) 'ERROR: nprocs must equal nsx*nsy'
    call MPI_Abort(comm, 1, ierr)
  end if

  rx = mod(rank, nsx)
  ry = rank / nsx

  nxg = nx*nsx
  nyg = ny*nsy

  hx = M/2
  hy = N/2

  ! local SAT arrays for one k-block
  allocate(pxX(nx,ny,kb), pxW(nx,ny,kb))
  allocate(SX(nx,ny,kb),  SW(nx,ny,kb))

  ! row/col communicators (once)
  color = ry; key = rx
  call MPI_Comm_split(comm, color, key, row_comm, ierr)
  color = rx; key = ry
  call MPI_Comm_split(comm, color, key, col_comm, ierr)

  k1 = 1
  do while (k1 <= nzm)
    k2 = min(nzm, k1+kb-1)

    !--------------------------------------------
    ! 1) local x-prefix for masked X and W
    !--------------------------------------------
    do kk = 1, (k2-k1+1)
      k = k1 + kk - 1
      do j = 1, ny
        if (TERR(1,j,k) >= 0.5) then
          pxX(1,j,kk) = A(1,j,k)
          pxW(1,j,kk) = 1.0
        else
          pxX(1,j,kk) = 0.0
          pxW(1,j,kk) = 0.0
        end if
        do i = 2, nx
          if (TERR(i,j,k) >= 0.5) then
            pxX(i,j,kk) = pxX(i-1,j,kk) + A(i,j,k)
            pxW(i,j,kk) = pxW(i-1,j,kk) + 1.0
          else
            pxX(i,j,kk) = pxX(i-1,j,kk)
            pxW(i,j,kk) = pxW(i-1,j,kk)
          end if
        end do
      end do
    end do

    !--------------------------------------------
    ! 2) stitch x with Exscan in row_comm
    !--------------------------------------------
    allocate(send_west_X(ny,(k2-k1+1)), recv_west_X(ny,(k2-k1+1)))
    allocate(send_west_W(ny,(k2-k1+1)), recv_west_W(ny,(k2-k1+1)))
    do kk = 1, (k2-k1+1)
      do j = 1, ny
        send_west_X(j,kk) = pxX(nx,j,kk)
        send_west_W(j,kk) = pxW(nx,j,kk)
      end do
    end do
    call MPI_Exscan(send_west_X, recv_west_X, ny*(k2-k1+1), mpi_wp, MPI_SUM, row_comm, ierr)
    call MPI_Exscan(send_west_W, recv_west_W, ny*(k2-k1+1), mpi_wp, MPI_SUM, row_comm, ierr)
    if (rx == 0) then
      recv_west_X(:,:) = 0.0
      recv_west_W(:,:) = 0.0
    end if
    do kk = 1, (k2-k1+1)
      do j = 1, ny
        do i = 1, nx
          pxX(i,j,kk) = pxX(i,j,kk) + recv_west_X(j,kk)
          pxW(i,j,kk) = pxW(i,j,kk) + recv_west_W(j,kk)
        end do
      end do
    end do
    deallocate(send_west_X, recv_west_X, send_west_W, recv_west_W)

    !--------------------------------------------
    ! 3) local y-prefix => SX,SW
    !--------------------------------------------
    do kk = 1, (k2-k1+1)
      do i = 1, nx
        SX(i,1,kk) = pxX(i,1,kk)
        SW(i,1,kk) = pxW(i,1,kk)
        do j = 2, ny
          SX(i,j,kk) = SX(i,j-1,kk) + pxX(i,j,kk)
          SW(i,j,kk) = SW(i,j-1,kk) + pxW(i,j,kk)
        end do
      end do
    end do

    !--------------------------------------------
    ! 4) stitch y with Exscan in col_comm
    !--------------------------------------------
    allocate(send_south_X(nx,(k2-k1+1)), recv_south_X(nx,(k2-k1+1)))
    allocate(send_south_W(nx,(k2-k1+1)), recv_south_W(nx,(k2-k1+1)))
    do kk = 1, (k2-k1+1)
      do i = 1, nx
        send_south_X(i,kk) = SX(i,ny,kk)
        send_south_W(i,kk) = SW(i,ny,kk)
      end do
    end do
    call MPI_Exscan(send_south_X, recv_south_X, nx*(k2-k1+1), mpi_wp, MPI_SUM, col_comm, ierr)
    call MPI_Exscan(send_south_W, recv_south_W, nx*(k2-k1+1), mpi_wp, MPI_SUM, col_comm, ierr)
    if (ry == 0) then
      recv_south_X(:,:) = 0.0
      recv_south_W(:,:) = 0.0
    end if
    do kk = 1, (k2-k1+1)
      do i = 1, nx
        do j = 1, ny
          SX(i,j,kk) = SX(i,j,kk) + recv_south_X(i,kk)
          SW(i,j,kk) = SW(i,j,kk) + recv_south_W(i,kk)
        end do
      end do
    end do
    deallocate(send_south_X, recv_south_X, send_south_W, recv_south_W)

    !--------------------------------------------
    ! 5) Fetch just the SAT tiles needed for this tile’s window
    !    General: get rank rectangle that covers [Igmin-hx-1 .. Igmax+hx] etc.
    !--------------------------------------------
    call sat_neighborhood_bounds(rx,ry,nx,ny,nxg,nyg,nsx,nsy,hx,hy, rxL,rxR,ryB,ryT)

    nRx = rxR - rxL + 1
    nRy = ryT - ryB + 1
    nxN = nRx * nx
    nyN = nRy * ny

    allocate(SXN(nxN,nyN,(k2-k1+1)))
    allocate(SWN(nxN,nyN,(k2-k1+1)))
    SXN(:,:,:) = 0.0
    SWN(:,:,:) = 0.0

    call gather_sat_tiles(SX, SXN, nx,ny,(k2-k1+1), rx,ry, rxL,rxR,ryB,ryT, nsx,nsy, comm)
    call gather_sat_tiles(SW, SWN, nx,ny,(k2-k1+1), rx,ry, rxL,rxR,ryB,ryT, nsx,nsy, comm)

    !--------------------------------------------
    ! 6) Compute box mean for this k-block
    !--------------------------------------------
    do kk = 1, (k2-k1+1)
      k = k1 + kk - 1
      do j = 1, ny
        do i = 1, nx
          Ig0 = rx*nx + i
          Jg0 = ry*ny + j

          I1 = max(1, Ig0-hx)
          I2 = min(nxg, Ig0+hx)
          J1 = max(1, Jg0-hy)
          J2 = min(nyg, Jg0+hy)

          sumX = sat_from_neigh(SXN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I2,  J2,  kk) &
               - sat_from_neigh(SXN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I1-1,J2,  kk) &
               - sat_from_neigh(SXN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I2,  J1-1,kk) &
               + sat_from_neigh(SXN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I1-1,J1-1,kk)

          sumW = sat_from_neigh(SWN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I2,  J2,  kk) &
               - sat_from_neigh(SWN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I1-1,J2,  kk) &
               - sat_from_neigh(SWN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I2,  J1-1,kk) &
               + sat_from_neigh(SWN,nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, I1-1,J1-1,kk)

          if (sumW > 0.5) then
            A(i,j,k) = sumX / sumW
          else
            A(i,j,k) = 0.0
          end if
        end do
      end do
    end do

    deallocate(SXN,SWN)

    k1 = k2 + 1
  end do

  call MPI_Comm_free(row_comm, ierr)
  call MPI_Comm_free(col_comm, ierr)

  deallocate(pxX,pxW,SX,SW)

  call t_stopf ('task_box_smooth')


contains

  subroutine sat_neighborhood_bounds(rx,ry,nx,ny,nxg,nyg,nsx,nsy,hx,hy, rxL,rxR,ryB,ryT)
    implicit none
    integer, intent(in) :: rx,ry,nx,ny,nxg,nyg,nsx,nsy,hx,hy
    integer, intent(out):: rxL,rxR,ryB,ryT
    integer :: Igmin, Igmax, Jgmin, Jgmax
    integer :: Ineed_min, Ineed_max, Jneed_min, Jneed_max

    Igmin = rx*nx + 1
    Igmax = rx*nx + nx
    Jgmin = ry*ny + 1
    Jgmax = ry*ny + ny

    ! Need I2 and (I1-1) => cover [Igmin-hx-1 .. Igmax+hx]
    Ineed_min = max(1, Igmin - hx - 1)
    Ineed_max = min(nxg, Igmax + hx)

    Jneed_min = max(1, Jgmin - hy - 1)
    Jneed_max = min(nyg, Jgmax + hy)

    rxL = (Ineed_min-1)/nx
    rxR = (Ineed_max-1)/nx
    ryB = (Jneed_min-1)/ny
    ryT = (Jneed_max-1)/ny

    if (rxL < 0) rxL = 0
    if (ryB < 0) ryB = 0
    if (rxR > nsx-1) rxR = nsx-1
    if (ryT > nsy-1) ryT = nsy-1
  end subroutine sat_neighborhood_bounds

  subroutine gather_sat_tiles(Slocal, Sbig, nx,ny,nk, rx,ry, rxL,rxR,ryB,ryT, nsx,nsy, comm)
    use mpi
    implicit none
    integer, intent(in) :: nx,ny,nk, rx,ry,rxL,rxR,ryB,ryT, nsx,nsy, comm
    real   , intent(in)  :: Slocal(nx,ny,nk)
    real   , intent(inout):: Sbig((rxR-rxL+1)*nx, (ryT-ryB+1)*ny, nk)

    integer :: ierr, rank, p
    integer :: rxi, ryi, rrank
    integer :: nRx, nRy
    integer :: bx, by
    integer :: count
    real, allocatable :: sbuf(:), rbuf(:)

    call MPI_Comm_rank(comm, rank, ierr)

    nRx = rxR - rxL + 1
    nRy = ryT - ryB + 1
    count = nx*ny*nk

    allocate(sbuf(count))
    call pack_tile(Slocal, nx,ny,nk, sbuf)

    ! Put my own tile in place
    bx = (rx - rxL)*nx
    by = (ry - ryB)*ny
    call unpack_tile(Sbig, nRx*nx, nRy*ny, nk, bx,by, nx,ny, sbuf)

    allocate(rbuf(count))

    ! Loop over needed neighbor tiles and fetch them with Sendrecv (simple, bounded)
    do ryi = ryB, ryT
      do rxi = rxL, rxR
        rrank = ryi*nsx + rxi
        if (rrank == rank) cycle

        call MPI_Sendrecv(sbuf, count, mpi_wp, rrank, 777, &
                          rbuf, count, mpi_wp, rrank, 777, comm, MPI_STATUS_IGNORE, ierr)

        bx = (rxi - rxL)*nx
        by = (ryi - ryB)*ny
        call unpack_tile(Sbig, nRx*nx, nRy*ny, nk, bx,by, nx,ny, rbuf)
      end do
    end do

    deallocate(sbuf,rbuf)
  end subroutine gather_sat_tiles

  subroutine pack_tile(S, nx,ny,nk, buf)
    implicit none
    integer, intent(in) :: nx,ny,nk
    real   , intent(in) :: S(nx,ny,nk)
    real   , intent(out):: buf(nx*ny*nk)
    integer :: i,j,k,idx
    idx = 0
    do k = 1, nk
      do j = 1, ny
        do i = 1, nx
          idx = idx + 1
          buf(idx) = S(i,j,k)
        end do
      end do
    end do
  end subroutine pack_tile

  subroutine unpack_tile(Sbig, NXB,NYB,nk, bx,by, nx,ny, buf)
    implicit none
    integer, intent(in) :: NXB,NYB,nk, bx,by, nx,ny
    real   , intent(inout) :: Sbig(NXB,NYB,nk)
    real   , intent(in)    :: buf(nx*ny*nk)
    integer :: i,j,k,idx
    idx = 0
    do k = 1, nk
      do j = 1, ny
        do i = 1, nx
          idx = idx + 1
          Sbig(bx+i, by+j, k) = buf(idx)
        end do
      end do
    end do
  end subroutine unpack_tile

  pure real function sat_from_neigh(Sbig, nx,ny,nRx,nRy,rxL,ryB,nxg,nyg, Ig,Jg, kk)
    implicit none
    integer, intent(in) :: nx,ny,nRx,nRy,rxL,ryB,nxg,nyg,Ig,Jg,kk
    real   , intent(in) :: Sbig(nRx*nx, nRy*ny, *)
    integer :: rxi, ryi, iloc, jloc, bx, by

    if (Ig < 1 .or. Jg < 1) then
      sat_from_neigh = 0.0
      return
    end if
    if (Ig > nxg .or. Jg > nyg) then
      sat_from_neigh = 0.0
      return
    end if

    rxi = (Ig-1)/nx
    ryi = (Jg-1)/ny
    iloc = Ig - rxi*nx
    jloc = Jg - ryi*ny

    bx = (rxi - rxL)*nx
    by = (ryi - ryB)*ny
    sat_from_neigh = Sbig(bx+iloc, by+jloc, kk)
  end function sat_from_neigh

end subroutine box_smooth_3d

