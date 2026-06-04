! This Fortran program creates a variable grid in latitude.
! The grid is isotropic, meaning the resolution is approximately 
! the same in both horizontal directions (in meters),
! from the equator extending polarward up to a specific latitude, 
! which is also computed within the program.
! Beyond this latitude, the resolution in latitude gradually coarsens, 
! and at the poles, it becomes approximately equal to a specified resolution.
! It is implied that the resolution in longitude (in degrees) is uniform.
! Coded by Marat Khairoutdinov, 2018

! Input parameters specified by the user are:
!   nx - number of grid points in longitude.
!   ny - number of grid points in latitude.
!   dy0 - grid spacing in latutude at the pole (in degrees).

! Output:
!   Two files are generated: grid_out and grid_out.txt
!   These files can be renamed by user and serve the following purposes:
!   - grid_out: Used by gSAM to form the lat/lon grid.
!     first column: grid latitudes;
!     second column: grid index;
!     third column: resolution in longitude (in km)
!     fourth column: local ratio of resolution in lat (km) to resolution in lon (km) 
!   - grid_out.txt: Used by the NCL interpolation scripts.

implicit none
!====================================================
!integer, parameter :: ny=128
!integer, parameter :: nx=256
!integer, parameter :: nx=768
!integer, parameter :: ny=384
!integer, parameter :: ny=360
!integer, parameter :: nx=720
!integer, parameter :: ny=4096
!integer, parameter :: nx=8192
!integer, parameter :: ny=1152
!integer, parameter :: nx=4608
!integer, parameter :: ny=4608
!integer, parameter :: ny=2304
!integer, parameter :: nx=9216
!integer, parameter :: ny=512
!integer, parameter :: nx=2000
!integer, parameter :: ny=500
!integer, parameter :: nx=2048
!integer, parameter :: ny=256
integer, parameter :: ny=6400
integer, parameter :: nx=12800
!integer, parameter :: ny=720
!integer, parameter :: nx=1440
!integer, parameter :: ny=256
!integer, parameter :: nx=1024
!integer, parameter :: ny=576
!integer, parameter :: nx=1152
!integer, parameter :: ny=1440
!integer, parameter :: nx=2880
!real(8), parameter :: dy0 = 0.5 
!real(8), parameter :: dy0 = 1. 
real(8), parameter :: dy0 = 2. 
!====================================================
!====================================================
! don't edit below this line.

real(8), parameter :: dx = 360./nx
integer, parameter :: ny0 = ny/2+1

real(8) lat(ny0),lat1(ny0-1),error, r, dlat 
real(8) latout(ny)
real(8):: deg2rad = acos(-1.)/180.
integer niter, n, k

niter = 0.
error = 1.e5
r = 0.5

do while(error.gt.1.e-3)
 niter = niter + 1
1 lat(1) = 89.9
 lat(2) = lat(1) - dy0
 dlat = dy0
 n=3
 do while (.true.)
   dlat = dlat*r
   lat(n) = lat(n-1)-dlat
   if(dlat.lt.dx*cos(deg2rad*0.5*(lat(n)+lat(n-1)))) exit 
!   print*,n,lat(n)
   n = n+1
   if(n.gt.ny0) then
     n= n-1
!     print*,'n exceeded ny0!',r,lat(n),dlat,dx,dx*cos(deg2rad*lat(n))
     if(lat(n).gt.0.01) then
       r = r*1.00001
       goto 1
     else if (lat(n).lt.0.01) then
       r = r*0.99999
       goto 1
     end if
     goto 2
   end if
 end do 
 do while(n.le.ny0)
  lat(n) = lat(n-1)-dx*cos(deg2rad*lat(n-1))
  lat(n) = lat(n-1)-dx*cos(deg2rad*0.5*(lat(n-1)+lat(n)))
  lat(n) = lat(n-1)-dx*cos(deg2rad*0.5*(lat(n-1)+lat(n)))
  lat(n) = lat(n-1)-dx*cos(deg2rad*0.5*(lat(n-1)+lat(n)))
!  print*,n,lat(n)
  n = n+1 
 end do
 n = n-1
 if(mod(niter,1000).eq.0) print*,'niter = ',niter,n, 'r=',r,'lat(n)=',lat(n)
 error = abs(lat(n))
 if(lat(n).lt.0.) then
  r = r*0.99999
 else
  r = r*1.00001
 end if
end do
2 continue
print*,'dx=',dx
k=0
do n=1,ny0-1
 dlat = lat(n+1)-lat(n)
 lat1(n) = 0.5*(lat(n)+lat(n+1))
 k= k+1
 write(*,'(f15.5,i10,2f15.5)') -lat1(n),k,-dlat*20000./180.,-dlat/(dx*cos(deg2rad*lat(n)))
 latout(k) = -lat1(n)
end do
do n=ny0-1,1,-1
 k= k+1
 dlat = lat(n+1)-lat(n)
 write(*,'(f15.5,i10,2f15.5)') lat1(n),k,-dlat*20000./180.,-dlat/(dx*cos(deg2rad*lat(n)))
 latout(k) = lat1(n)
end do

open(1,file='grid_out',form='formatted')
do n=1,ny
 dlat = latout(min(ny,n+1))-latout(n)
 write(1,'(f15.5,i10,2f15.5)') latout(n),n,dlat*20000./180.,40000./nx*cos(deg2rad*latout(n))
end do
open(2,file='grid_out.txt',form='formatted')
write(2,*) latout

end


