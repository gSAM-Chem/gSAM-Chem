implicit none
! interpolate from high-resolution land IGBP classification to low resolution 
! Marat L, Sep 2020

!!! IMPORTANT !!!!
! this code uses interpolated and smoothed terrain made by make_2D_terr.f90 !
! compile and run that code first to produce 

!!! IMPORTANT !!!!
! if crashing with seg fault, it could be because you need more memory to run
! as the input dataset can be very large.

!======================================================================
! Edit starting here:

character(100):: inoutdir="./BIN_D"  ! directory where terrain file produced by 
                                     ! make_2D_terr.f90 is and where landtype
                                     ! file will be written
! stop editing here:
!======================================================================
! landtype high-resolution dataset:
! high-resolution landtype dataset:
character(100):: datafile="./DATA/landtype_16200x8100_global_IGDP_lonflip.bin"
!======================================================================
real(8), allocatable :: lat_in(:), lon_in(:) ! origonal grid
integer, allocatable :: land_in(:,:) ! original data
integer, allocatable :: land(:,:) ! interpolatted 
integer, allocatable :: land1(:,:) ! interpolatted fixed
integer, allocatable :: terr(:,:) ! interpolatted terrain
real, allocatable :: lat(:), lon(:) ! interpolate-to grid
real, allocatable :: tmp(:)
integer nlon_in, nlat_in
integer nlon, nlat
integer i,j,jb,jc,m,n,k,i1,i2,j1,j2
real(8) lat_min, lat_max, lon_min, lon_max, pi
integer types(17), ii,jj
integer ncorrected

! read terrain file:
open(4, file=trim(inoutdir)//"/out_terr.bin",form="unformatted")
read(4) nlon
read(4) nlat
print*,'interpolate land to to grid (nlon  nlat): ',nlon,nlat
allocate (land(-nlon+1:2*nlon,nlat)) 
allocate (land1(nlon,nlat)) 
allocate (terr(nlon,nlat))
allocate (lat(nlat), lon(nlon))
read(4) lon(1:nlon)
read(4) lat(1:nlat)
read(4) 
read(4) terr(1:nlon,1:nlat)
print*,"lon:",minval(lon),maxval(lon)
print*,"lat:",minval(lat),maxval(lat)

! read original data
open(1,file=datafile,form="unformatted")
read(1) nlon_in
read(1) nlat_in
print*,"original grid size:",nlon_in,nlat_in
allocate(lat_in(nlat_in),lon_in(-nlon_in+1:2*nlon_in),land_in(-nlon_in+1:2*nlon_in,nlat_in))
allocate(tmp(nlon_in))
read(1) lon_in(1:nlon_in)
read(1) lat_in(:)
print*,"lon_in:",minval(lon_in(1:nlon_in)),maxval(lon_in(1:nlon_in))
print*,"lat_in:",minval(lat_in),maxval(lat_in)
read(1) land_in(1:nlon_in,1:nlat_in)
print*,"land_in:",minval(land_in),maxval(land_in)

lon_in(1-nlon_in:0) = lon_in(1:nlon_in)-360.
lon_in(nlon_in+1:nlon_in+nlon_in) = lon_in(1:nlon_in)+360.
do j=1,nlat_in
  land_in(1-nlon_in:0,j) = land_in(1:nlon_in,j)
  land_in(nlon_in+1:nlon_in+nlon_in,j) = land_in(1:nlon_in,j)
end do
print*,"land_in:",minval(land_in),maxval(land_in)

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

  
if(j2.lt.j1) then
 print*,'>>>>>>',j1,j2,lat_min,lat_max,lat_in(j1),lat_in(j2)
end if

!  if(i.eq.nlon/2) then
!   print*,'>>>',lon_min,lon_max,lat_min,lat_max,i1,i2,j1,j2
!  end if

  types = 0
  ! find dominating land type and assign it:
  do jj=j1,j2
   do ii=i1,i2
     types(land_in(ii,jj)+1) =  types(land_in(ii,jj)+1) +1
   end do
  end do
  ii = sum(maxloc(types))-1
  land(i,j) = ii
 end do

end do

print*,"land:",minval(land(1:nlon,1:nlat)),maxval(land(1:nlon,1:nlat))

open(3, file=trim(inoutdir)//"/out_landtype.bin",form="unformatted")
write(3) nlon
write(3) nlat
write(3) lon
write(3) lat
write(3) land(1:nlon,1:nlat)



!pi = cos(01.)

!do j=1,nlat
!  land(1-nlon:0,j) = land(1:nlon,j)
!  land(nlon+1:nlon+nlon,j) = land(1:nlon,j)
!end do
!land1(1:nlon,1:nlat) = land(1:nlon,1:nlat)
!
!ncorrected = 0
!do j=1,nlat
! k = 1./cos(lat(j)*pi/180.)
! do i=1,nlon
!! fix land by making sure that you have some nearest landtype > 0 (not sea) where there is terrain > 0
!  if(terr(i,j).gt.1..and.land(i,j).eq.0) then
!   ! assign nearest land type > 0
!  ! find dominating land type and assign it:
!    do m=1,20
!     n = m*k
!     types = 0
!     do jj=max(1,j-m),min(nlat,j+m)
!      do ii=i-n,i+n
!        if(land(ii,jj).gt.0) types(land(ii,jj)) =  types(land(ii,jj)) +1
!      end do
!     end do
!     if(maxval(types).gt.0) then
!       ii = sum(maxloc(types(1:16)))
!       print*,'corrected:',i,j,ii
!       land1(i,j) = ii
!       ncorrected = ncorrected + 1
!       exit
!     end if
!    end do
!  end if
! end do
!end do
!print*,'corected',ncorrected,'points',real(ncorrected)/(nlat*nlon)*100.,'%'
!
!write(3) land1(1:nlon,1:nlat)


end


