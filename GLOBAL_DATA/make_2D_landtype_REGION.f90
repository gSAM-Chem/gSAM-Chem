implicit none
!======================================================================
! high-resolution global IGBP landtype dataset:
! type 0 - water
! types from 1 to 16 - land
! binary dataset:
character(100):: datafile="./DATA/landtype_16200x8100_global_IGDP_lonflip.bin"
integer nlon_in, nlat_in ! number of longitudes and latitudes
real(8), allocatable :: lat_in(:), lon_in(:) ! lat/lon coordinates
integer, allocatable :: land_in(:,:) ! global data
!=======================================================
! Target dataset variables on regional grid:
integer nlon, nlat ! number of longitudes and latitudes
real, allocatable :: lat(:), lon(:) ! lat/lon coordinates
integer, allocatable :: land(:,:) ! target array
integer, allocatable :: terr(:,:) ! target terrain (given)
real, allocatable :: lakemask(:,:) ! target inland lake mask (given; 0 - land, 1 - lake)

!  local varaibles:
integer :: i, j, ii, jj, k, kk
integer :: ilo, ihi, jlo, jhi
integer :: counts(0:16)
integer :: lt, best_lt, best_cnt
real(8) :: dlon_in, dlat_in, lon0_in, lat0_in
real(8) :: lon_w, lon_e, lat_s, lat_n
logical :: lat_in_increasing
integer :: n_water, n_land


!============================================================================
! Read global data:

open(1,file=datafile,form="unformatted")
read(1) nlon_in
read(1) nlat_in
print*,"original grid size:",nlon_in,nlat_in
allocate(lat_in(nlat_in),lon_in(nlon_in),land_in(nlon_in,nlat_in))
read(1) lon_in(:)
read(1) lat_in(:)
print*,"lon_in:",minval(lon_in),maxval(lon_in)
print*,"lat_in:",minval(lat_in),maxval(lat_in)
read(1) land_in(1:nlon_in,1:nlat_in)
print*,"land_in:",minval(land_in),maxval(land_in)


!===================================================================
! Read terrain and lake mask file that also a source of target grid:
open(4, file="./BIN_D/out_terr.bin",form="unformatted")
read(4) nlon
read(4) nlat
allocate (land(nlon,nlat))
allocate (terr(nlon,nlat))
allocate (lakemask(nlon,nlat))
allocate (lat(nlat), lon(nlon)) 
print*,"target grid size:",nlon,nlat
read(4) lon(1:nlon)
read(4) lat(1:nlat)
print*,"lon:",minval(lon(1:nlon)),maxval(lon(1:nlon))
print*,"lat:",minval(lat),maxval(lat)
read(4) terr(1:nlon,1:nlat)
read(4) 
read(4) lakemask(1:nlon,1:nlat)
print*,'lake mask: ',minval(lakemask),maxval(lakemask)

!=================================================================
!=================================================================
! Majority (most-common) landtype remapping from fine global grid
! to coarser regional grid

! --- infer global grid spacing and orientation
dlon_in = lon_in(2) - lon_in(1)
dlat_in = lat_in(2) - lat_in(1)
lon0_in = lon_in(1)
lat0_in = lat_in(1)
lat_in_increasing = (dlat_in > 0.0_8)   ! not used below, but OK to keep

land(:,:) = 0

do j = 1, nlat
  do i = 1, nlon

    ! --- lakes enforced first
    if (lakemask(i,j) == 1) then
      land(i,j) = 0
      cycle
    end if

    ! --- ocean enforced from terrain (as you stated: terr=0 for ocean)
    if (terr(i,j) <= 0) then
      land(i,j) = 0
      cycle
    end if

    ! --- target cell bounds (lon)
    lon_w = lon(i) - 0.5d0 * (lon(min(i+1,nlon)) - lon(max(i-1,1)))
    lon_e = lon(i) + 0.5d0 * (lon(min(i+1,nlon)) - lon(max(i-1,1)))

    ! --- target cell bounds (lat) (make it explicit and consistent)
    if (j == 1) then
      lat_s = lat(j) - 0.5d0 * (lat(j+1) - lat(j))
      lat_n = 0.5d0 * (lat(j) + lat(j+1))
    else if (j == nlat) then
      lat_s = 0.5d0 * (lat(j) + lat(j-1))
      lat_n = lat(j) + 0.5d0 * (lat(j) - lat(j-1))
    else
      lat_s = 0.5d0 * (lat(j) + lat(j-1))
      lat_n = 0.5d0 * (lat(j) + lat(j+1))
    end if

    ! --- map bounds to global grid indices
    ilo = floor((lon_w - lon0_in) / dlon_in) + 1
    ihi = floor((lon_e - lon0_in) / dlon_in) + 1

    jlo = floor((lat_s - lat0_in) / dlat_in) + 1
    jhi = floor((lat_n - lat0_in) / dlat_in) + 1

    if (jlo > jhi) then
      k = jlo
      jlo = jhi
      jhi = k
    end if

    jlo = max(1, jlo)
    jhi = min(nlat_in, jhi)

    counts(:) = 0

    ! --- count landtypes inside target cell
    do jj = jlo, jhi
      do ii = ilo, ihi
        k = modulo(ii-1, nlon_in) + 1
        lt = land_in(k,jj)
        if (lt >= 0 .and. lt <= 16) counts(lt) = counts(lt) + 1
      end do
    end do

    n_water = counts(0)
    n_land  = sum(counts(1:16))

    ! --------------------------------------------------------------
    ! FIX: allow water to win by majority even if terr>0
    ! --------------------------------------------------------------
    if (n_land == 0) then
      land(i,j) = 0

    else if (n_water > n_land) then
      land(i,j) = 0

    else
      ! --- choose most common landtype among 1..16
      best_lt  = 0
      best_cnt = 0
      do lt = 1, 16
        if (counts(lt) > best_cnt) then
          best_cnt = counts(lt)
          best_lt  = lt
        end if
      end do

      if (best_cnt > 0) then
        land(i,j) = best_lt
      else
        ! fallback: extremely rare case (should not happen if n_land>0)
        land(i,j) = 16  ! bare soil
      end if
    end if

  end do
end do

print*,"land (target):",minval(land),maxval(land)

!===================================================================
! Write output file:
open(3, file="./BIN_D/out_landtype.bin",form="unformatted")
write(3) nlon
write(3) nlat
write(3) lon
write(3) lat
write(3) land(1:nlon,1:nlat)

end
