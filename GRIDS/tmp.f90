implicit none

! Should be compiled with -r8 flag! !

integer, parameter :: ny0=360
integer, parameter :: nx=720
real, parameter :: dy0 = 5. ! spacing at the poles
real, parameter :: dx = 360./nx
integer, parameter :: ny = ny0/2+1

real lat(ny),lat1(ny-1),error, r, dlat 
real latout(ny0)
real:: deg2rad = acos(-1.)/180.
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
   if(n.gt.ny) then
     n= n-1
!     print*,'n exceeded ny!',r,lat(n),dlat,dx,dx*cos(deg2rad*lat(n))
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
 do while(n.le.ny)
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
do n=1,ny-1
 dlat = lat(n+1)-lat(n)
 lat1(n) = 0.5*(lat(n)+lat(n+1))
 k= k+1
 write(*,'(f15.5,i10,2f15.5)') -lat1(n),k,-dlat*20000./180.,-dlat/(dx*cos(deg2rad*lat(n)))
 latout(k) = -lat1(n)
end do
do n=ny-1,1,-1
 k= k+1
 dlat = lat(n+1)-lat(n)
 write(*,'(f15.5,i10,2f15.5)') -lat1(n),k,-dlat*20000./180.,-dlat/(dx*cos(deg2rad*lat(n)))
 latout(k) = lat1(n)
end do

open(1,file='grid_out',form='formatted')
do n=1,ny0
 dlat = latout(min(ny0,n+1))-latout(n)
 write(1,*) latout(n),n,dlat*20000./180.
end do
open(2,file='grid_out.txt',form='formatted')
write(2,*) latout

end


