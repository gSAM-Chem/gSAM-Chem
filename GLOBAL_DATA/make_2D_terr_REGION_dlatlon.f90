implicit none
! interpolate from high-resolution terrain to low resolution terrain
! Marat K, Sep 2020

! The utility produces the files:
! out_terr.bin ! terrain (it cannot be directly used for model input)
! out_latlon.txt ! longitudes and latitudes of the regional domain


!======================================================================
! Edit starting here:

! target grid parameters:

! USA 
real, parameter :: lon0 = -95.+360. ! longitude (degrees east) of domain center
real, parameter :: lat0 = 36. ! latitude of domain's center
integer, parameter :: nlon = 384
integer, parameter :: nlat = 192
real, parameter :: dlon = 0.25 ! resolution in x
real, parameter :: dlat = 0.25 ! resolution in y

integer, parameter :: m_sm = 1  !  number of grid points to smooth terrain

! output directory:
character(100):: outdir="./BIN_D"

! stop editing here:
!======================================================================
! high-resolution terrain dataset:
!character(100):: datafile="./DATA/ETOPO1_latlon_flipped.bin"
character(100):: datafile="./DATA/orog30s_hydroflat_new.bin"
!======================================================================
real, allocatable :: lat_in(:), lon_in(:) ! origonal grid
real, allocatable :: terr_in(:,:) ! original terrain
real terr(-nlon+1:2*nlon,nlat) ! interpolatted terrain
real lat(nlat), lon(nlon) ! interpolate-to grid
real, allocatable :: tmp(:)
real tmpp(1:nlon)
integer nlon_in, nlat_in
integer i,j,jb,jc,m,i1,i2,j1,j2
real lat_min, lat_max, lon_min, lon_max, pi
real summ(-nlon:2*nlon)

! read origonal terrain
open(1,file=datafile,form="unformatted")
read(1) nlon_in
read(1) nlat_in
print*,"original grid size:",nlon_in,nlat_in
allocate(lat_in(nlat_in),lon_in(-nlon_in+1:2*nlon_in),terr_in(-nlon_in+1:2*nlon_in,nlat_in))
allocate(tmp(nlon_in))
read(1) lon_in(1:nlon_in)
read(1) lat_in(:)
read(1) terr_in(1:nlon_in,1:nlat_in)
print*,"lon_in:",minval(lon_in(1:nlon_in)),maxval(lon_in(1:nlon_in))
print*,"lat_in:",minval(lat_in),maxval(lat_in)
print*,"lon(1):",lon_in(1)
print*,"lat(1):",lat_in(1)
print*,"terr_in:",minval(terr_in),maxval(terr_in)


lon_in(1-nlon_in:0) = lon_in(1:nlon_in)-360.
lon_in(nlon_in+1:nlon_in+nlon_in) = lon_in(1:nlon_in)+360.
do j=1,nlat_in
  terr_in(1-nlon_in:0,j) = terr_in(1:nlon_in,j)
  terr_in(nlon_in+1:nlon_in+nlon_in,j) = terr_in(1:nlon_in,j)
end do


!===============================================================
! make the target grid:

pi = acos(-1.)
print*,"dlat= ",dlat
print*,"dlon= ",dlon
do i=1,nlon
 lon(i) = lon0 + (i-nlon/2-0.5)*dlon
end do
do j=1,nlat
 lat(j) = lat0 + (j-nlat/2-0.5)*dlat
end do
print*,"nlon=",nlon
print*,"nlat=",nlat
print*,"lon:",minval(lon),maxval(lon)
print*,"lonW:",minval(lon)-360.,maxval(lon)-360.
print*,"lat:",minval(lat),maxval(lat)


!===============================================================

do j=1,nlat
 jb = max(1,j-1)
 jc = min(nlat,j+1)
 lat_min = 0.5*(lat(j)+lat(jb))
 lat_max = 0.5*(lat(j)+lat(jc))

 j1 = 1
 do m=2,nlat_in 
   if(lat_in(m).ge.lat_min) then 
     j1 = m-1 
     exit
   end if
 end do 
 j2 = nlat_in
 do m=nlat_in-1,1,-1 
   if(lat_in(m).le.lat_max) then 
     j2 = m+1 
     exit
   end if
 end do 

 do i=1,nlon

  lon_min = lon(i)-0.5*dlon
  lon_max = lon(i)+0.5*dlon

  i1=nint(lon_min/(360./nlon_in))+1
  i2=nint(lon_max/(360./nlon_in))+1
  terr(i,j) = sum(terr_in(i1:i2,j1:j2))/real((i2-i1+1)*(j2-j1+1))

 end do

end do


print*,"terr:",minval(terr(1:nlon,1:nlat)),maxval(terr(1:nlon,1:nlat))

open(3, file=trim(outdir)//"/out_terr.bin",form="unformatted")
write(3) nlon
write(3) nlat
write(3) lon
write(3) lat
write(3) terr(1:nlon,1:nlat)

! smooth terrain 

terr(-nlon+1:0,:) = terr(1:nlon,:)
terr(nlon+1:2*nlon,:) = terr(1:nlon,:)

pi = acos(-1.)
summ = 1.
do j=1,nlat
  j1 = max(1,j-m_sm/2)
  j2 = min(nlat,j+m_sm/2)
  do i=1,nlon
    i1 = i-m_sm/2
    i2 = i+m_sm/2
    tmpp(i) = sum(terr(i1:i2,j1:j2))/real((i2-i1+1)*(j2-j1+1)) 
  end do
!  print*,j,m, maxval(terr(1:nlon,j)),maxval(tmpp(1:nlon))
  terr(1:nlon,j) = tmpp(1:nlon)
end do
write(3) terr(1:nlon,1:nlat)

print*,"terr (smoothed):",minval(terr(1:nlon,1:nlat)),maxval(terr(1:nlon,1:nlat))

close(3)
open(3, file=trim(outdir)//"/out_latlon.txt",form="formatted")
write(3,*) lon(1:nlon)
write(3,*) lat(1:nlat)
end


