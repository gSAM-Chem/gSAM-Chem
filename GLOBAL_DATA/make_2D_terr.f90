implicit none
! interpolate from high-resolution terrain to target resolution terrain
! Marat K, Sep 2020

!======================================================================
! Edit starting here:

integer, parameter :: nlon = 12800, nlat = 6400
!integer, parameter :: nlon = 1440, nlat = 720
!integer, parameter :: nlon = 2048, nlat = 512
!integer, parameter :: nlon = 3840, nlat = 1920
!integer, parameter :: nlon = 9216, nlat = 128
!integer, parameter :: nlon = 9216, nlat = 4608
!integer, parameter :: nlon = 1152, nlat = 576

! model grid to interpolate to:
character(100):: gridfile="./../GRIDS/lat_6400.txt"
!character(100):: gridfile="./../GRIDS/lat_1440x720.txt"
!character(100):: gridfile="./../GRIDS/lat_512x2048.txt"
!character(100):: gridfile="./../GRIDS/lat_2000x500_dyvar.txt"
!character(100):: gridfile="./../GRIDS/lat_3840x1920_dyvar.txt"
!character(100):: gridfile="./../GRIDS/lat_128_T85.txt"
!character(100):: gridfile="./../GRIDS/lat_4608.txt"
!character(100):: gridfile="./../GRIDS/lat_576_dyvar.txt"

! output directory:
character(100):: outdir="./BIN_D"

! minimum number of grid points in longitudinal dorection to smooth terrain.

integer:: m_min = max(1,nlon/(2*nlat))

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
print*,"terr_in:",minval(terr_in),maxval(terr_in)


lon_in(1-nlon_in:0) = lon_in(1:nlon_in)-360.
lon_in(nlon_in+1:nlon_in+nlon_in) = lon_in(1:nlon_in)+360.
do j=1,nlat_in
  terr_in(1-nlon_in:0,j) = terr_in(1:nlon_in,j)
  terr_in(nlon_in+1:nlon_in+nlon_in,j) = terr_in(1:nlon_in,j)
end do

open(2,file=gridfile,form="formatted")
read(2,*) lat(:)
do i=1,nlon
 lon(i) = (i-1)*360./real(nlon)
end do
print*,"nlon nlat:",nlon,nlat
print*,"lat:",minval(lat),maxval(lat)
print*,"lon:",minval(lon),maxval(lon)

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

  lon_min = (i-0.5)*(lon(2)-lon(1))
  lon_max = (i+0.5)*(lon(2)-lon(1))

  i1=nint(lon_min/(360./nlon_in))+1
  i2=nint(lon_max/(360./nlon_in))+1
  terr(i,j) = sum(terr_in(i1:i2,j1:j2))/real((i2-i1+1)*(j2-j1+1))

 end do

end do


print*,"terr:",minval(terr),maxval(terr)

open(3, file=trim(outdir)//"/out_terr.bin",form="unformatted")
write(3) nlon
write(3) nlat
write(3) lon
write(3) lat
write(3) terr(1:nlon,1:nlat)

! smooth terrain in longitudal direction

terr(-nlon+1:0,:) = terr(1:nlon,:)
terr(nlon+1:2*nlon,:) = terr(1:nlon,:)
print*,'m_min=',m_min

pi = acos(-1.)
summ = 1.
do j=1,nlat
! if(lat(j).gt.-60) exit
 m = max(m_min,nint((cos(lat(1)*pi/180.)/cos(lat(j)*pi/180.))**2*nlon))
! if(mod(m,2).ne.0) m = m+1
 print*,lat(m),m
 if(m.gt.0) then
  do i=1,nlon
    tmpp(i) = sum(terr(i-m/2:i+m/2,j))/sum(summ(i-m/2:i+m/2)) 
  end do
  print*,j,m, maxval(terr(1:nlon,j)),maxval(tmpp(1:nlon))
  terr(1:nlon,j) = tmpp(1:nlon)
 end if
end do
write(3) terr(1:nlon,1:nlat)

print*,"terr (smoothed):",minval(terr),maxval(terr)

end


