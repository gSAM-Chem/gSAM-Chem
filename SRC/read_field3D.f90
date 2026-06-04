#include "fppmacros"

!===========================================================

subroutine read_field3D (ntape,f,nxr,nyr,nzr,lon_gl,lat_gl,z,nxx,nyy,nzz)

! read arrays from a general 3D array
! and interpolate to gSAM grid using bilinear interpolation

use grid, only: rank, masterproc, pres
use params, only: dofliplon, dowallx, dowally, read_meters
implicit none
! Input:
integer ntape
real, intent(out) ::  f(nxx,nyy,nzz)
integer, intent(in) :: nxr,nyr,nzr ! size of array in file
real(8), intent(in) ::  lon_gl(*), lat_gl(*)
real, intent(in) ::  z(*)
integer, intent(in) :: nxx, nyy, nzz
! Local:
real(4),allocatable :: fld(:,:,:), fld1(:,:,:)
real(4),allocatable :: zr(:), pr(:)
real(4),allocatable :: latr(:),lonr(:),lonr1(:)
real x, y, wgt
integer m,n,k1,k2,k3,kk
integer k,i,j,it,jt,i0
call task_rank_to_index(rank,it,jt)
allocate(fld(0:nxr+1,0:nyr+1,2),latr(0:nyr+1),lonr(0:nxr+1),zr(nzr),pr(nzr))
if(dofliplon) allocate(fld1(0:nxr+1,nyr,2),lonr1(1:nxr))
read(ntape) lonr(1:nxr) 
read(ntape) latr(1:nyr)
read(ntape) zr(1:nzr)
read(ntape) pr(1:nzr)
if(dofliplon) then
  i0=-1
  do i=1,nxr
    if(lonr(i).gt.180.) then
      i0=i
      exit
    end if
  end do
  if(i0.lt.0) then
   if(masterproc) print*,'read_field3D_new: could not flip lons. lonr(1)=',lonr(1)
   call task_abort()
  end if
  lonr1(1:nxr-i0+1) = lonr(i0:nxr)-360.
  lonr1(nxr-i0+2:nxr) = lonr(1:i0-1)
  lonr(1:nxr) = lonr1(1:nxr)
  if(masterproc) print*,rank,'lonr(1)=',lonr(1),'lonr(nxr)=',lonr(nxr)
end if
!if(masterproc) print*,rank,minval(zr),maxval(zr),minval(pr),maxval(pr)
lonr(0)=2*lonr(1)-lonr(2)
lonr(nxr+1)=2*lonr(nxr)-lonr(nxr-1)

if(lon_gl(1).lt.lonr(0)) then
   if(masterproc) then
     if(read_meters) then
       print*,'read_field3D_new: cannot interpolate: x_gl(1)=',lon_gl(1),'x_data(0)=',lonr(0)
     else
       print*,'read_field3D_new: cannot interpolate: lon_gl(1)=',lon_gl(1),'lonr(0)=',lonr(0)
     end if
   end if
   call task_abort()
end if

!if(masterproc)print*,'lonr:',lonr(1:nxr)
!if(masterproc)print*,'latr:',latr(1:nyr)
!if(masterproc)print*,'zr:',zr(1:nzr)
!if(masterproc)print*,'pr:',pr(1:nzr)

if(lon_gl(nxx).gt.lonr(nxr+1)) then
  if(masterproc) print*,'read_field3D_new: cannot interpolate: lon_gl(nxx)=',lon_gl(nxx), &
                                                             'lonr(nxr+1)=',lonr(nxr+1) 
  call task_abort()
end if
latr(0)=2*latr(1)-latr(2)
latr(nyr+1)=2*latr(nyr)-latr(nyr-1)
if(lat_gl(1).lt.latr(0)) then
  if(masterproc) print*,'read_field3D_new: cannot interpolate: lat_gl(1)=',lat_gl(1), &
                                                   'is smaller than  latr(0)=',latr(0)
  call task_abort()
end if
if(lat_gl(nyy).gt.latr(nyr+1)) then
  if(masterproc) print*,'read_field3D_new: cannot interpolate: lat_gl(nyy)=',lat_gl(nyy), &
                                                    'is greater than latr(nyr+1)=',latr(nyr+1) 
  call task_abort()
end if
k1=1
k2=2
read(ntape) fld(1:nxr,1:nyr,k1)
read(ntape) fld(1:nxr,1:nyr,k2)
if(dofliplon) then
  fld1(1:nxr-i0+1,1:nyr,:) = fld(i0:nxr,1:nyr,:)
  fld1(nxr-i0+2:nxr,1:nyr,:) = fld(1:i0-1,1:nyr,:)
  fld(1:nxr,1:nyr,:) = fld1(1:nxr,1:nyr,:)
end if
if(dowallx) then
 fld(0,1:nyr,k1) = fld(1,1:nyr,k1)
 fld(nxr+1,1:nyr,k1) = fld(nxr,1:nyr,k1)
 fld(0,1:nyr,k2) = fld(1,1:nyr,k2)
 fld(nxr+1,1:nyr,k2) = fld(nxr,1:nyr,k2)
! fld(0,1:nyr,k1) = 2*fld(1,1:nyr,k1)-fld(2,1:nyr,k1)
! fld(nxr+1,1:nyr,k1) = 2*fld(nxr,1:nyr,k1)-fld(nxr-1,1:nyr,k1)
! fld(0,1:nyr,k2) = 2*fld(1,1:nyr,k2)-fld(2,1:nyr,k2)
! fld(nxr+1,1:nyr,k2) = 2*fld(nxr,1:nyr,k2)-fld(nxr-1,1:nyr,k2)
else
 fld(0,1:nyr,k1) = fld(nxr,1:nyr,k1)
 fld(nxr+1,1:nyr,k1) = fld(1,1:nyr,k1)
 fld(0,1:nyr,k2) = fld(nxr,1:nyr,k2)
 fld(nxr+1,1:nyr,k2) = fld(1,1:nyr,k2)
end if
if(dowally) then
 fld(0:nxr+1,0,k1) = fld(0:nxr+1,1,k1)
 fld(0:nxr+1,nyr+1,k1) = fld(0:nxr+1,nyr,k1)
 fld(0:nxr+1,0,k2) = fld(0:nxr+1,1,k2)
 fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,nyr,k2)
! fld(0:nxr+1,0,k1) = 2*fld(0:nxr+1,1,k1)-fld(0:nxr+1,2,k1)
! fld(0:nxr+1,nyr+1,k1) = 2*fld(0:nxr+1,nyr,k1)-fld(0:nxr+1,nyr-1,k1)
! fld(0:nxr+1,0,k2) = 2*fld(0:nxr+1,1,k2)-fld(0:nxr+1,2,k2)
! fld(0:nxr+1,nyr+1,k2) = 2*fld(0:nxr+1,nyr,k2)-fld(0:nxr+1,nyr-1,k2)
else
 fld(0:nxr+1,0,k1) = fld(0:nxr+1,nyr,k1)
 fld(0:nxr+1,nyr+1,k1) = fld(0:nxr+1,1,k1)
 fld(0:nxr+1,0,k2) = fld(0:nxr+1,nyr,k2)
 fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,1,k2)
end if

kk=1
do k=1,nzz
 if(zr(1).ne.-999.) then
  do while(z(k).gt.zr(kk+1))
   k3=k1
   k1=k2
   k2=k3
   kk=kk+1
   read(ntape) fld(1:nxr,1:nyr,k2)
   if(dofliplon) then
     fld1(1:nxr-i0+1,:,k2) = fld(i0:nxr,:,k2)
     fld1(nxr-i0+2:nxr,:,k2) = fld(1:i0-1,:,k2)
     fld(1:nxr,:,k2) = fld1(1:nxr,:,k2)
   end if
   if(dowallx) then
    fld(0,1:nyr,k2) = fld(1,1:nyr,k2)
    fld(nxr+1,1:nyr,k2) = fld(nxr,1:nyr,k2)
!    fld(0,1:nyr,k2) = 2*fld(1,1:nyr,k2)-fld(2,1:nyr,k2)
!    fld(nxr+1,1:nyr,k2) = 2*fld(nxr,1:nyr,k2)-fld(nxr-1,1:nyr,k2)
   else
    fld(0,1:nyr,k2) = fld(nxr,1:nyr,k2)
    fld(nxr+1,1:nyr,k2) = fld(1,1:nyr,k2)
   end if
   if(dowally) then
    fld(0:nxr+1,0,k2) = fld(0:nxr+1,1,k2)
    fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,nyr,k2)
!    fld(0:nxr+1,0,k2) = 2*fld(0:nxr+1,1,k2)-fld(0:nxr+1,2,k2)
!    fld(0:nxr+1,nyr+1,k2) = 2*fld(0:nxr+1,nyr,k2)-fld(0:nxr+1,nyr-1,k2)
   else
    fld(0:nxr+1,0,k2) = fld(0:nxr+1,nyr,k2)
    fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,1,k2)
   end if
  end do
  wgt = max(0.,min(1.,(z(k)-zr(kk))/(zr(kk+1)-zr(kk)))) ! dont introduce new max/min - Marat(02.2020)
 else
  do while(pres(k).lt.pr(kk+1))
   k3=k1
   k1=k2
   k2=k3
   kk=kk+1
   read(ntape) fld(1:nxr,1:nyr,k2)
   if(dofliplon) then
     fld1(1:nxr-i0+1,:,k2) = fld(i0:nxr,:,k2)
     fld1(nxr-i0+2:nxr,:,k2) = fld(1:i0-1,:,k2)
     fld(1:nxr,:,k2) = fld1(1:nxr,:,k2)
   end if
   if(dowallx) then
    fld(0,1:nyr,k2) = fld(1,1:nyr,k2)
    fld(nxr+1,1:nyr,k2) = fld(nxr,1:nyr,k2)
!    fld(0,1:nyr,k2) = 2*fld(1,1:nyr,k2)-fld(2,1:nyr,k2)
!    fld(nxr+1,1:nyr,k2) = 2*fld(nxr,1:nyr,k2)-fld(nxr-1,1:nyr,k2)
   else
    fld(0,1:nyr,k2) = fld(nxr,1:nyr,k2)
    fld(nxr+1,1:nyr,k2) = fld(1,1:nyr,k2)
   end if
   if(dowally) then
    fld(0:nxr+1,0,k2) = fld(0:nxr+1,1,k2)
    fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,nyr,k2)
!    fld(0:nxr+1,0,k2) = 2*fld(0:nxr+1,1,k2)-fld(0:nxr+1,2,k2)
!    fld(0:nxr+1,nyr+1,k2) = 2*fld(0:nxr+1,nyr,k2)-fld(0:nxr+1,nyr-1,k2)
   else
    fld(0:nxr+1,0,k2) = fld(0:nxr+1,nyr,k2)
    fld(0:nxr+1,nyr+1,k2) = fld(0:nxr+1,1,k2)
   end if
  end do
! wgt = max(0.,min(1.,pr(kk+1)*(pr(kk)-pres(k))/(pres(k)*(pr(kk)-pr(kk+1))))) ! mass weigting -MK 2023
  wgt = (log(pres(k))-log(pr(kk)))/(log(pr(kk+1))-log(pr(kk)))
 end if
 n = 0
 do j=1,nyy
  do while(lat_gl(j+jt).gt.latr(n+1)) 
    n = n+1 
  end do
  y = (lat_gl(j+jt)-latr(n))/(latr(n+1)-latr(n))
  m = 0
  do i=1,nxx
   do while(lon_gl(i+it).gt.lonr(m+1)) 
     m = m+1
   end do
   x = (lon_gl(i+it)-lonr(m))/(lonr(m+1)-lonr(m))
   f(i,j,k) = (1.-wgt)*((1.-x)*(1.-y)*fld(m,n,k1)+x*(1.-y)*fld(m+1,n,k1)+ &
                       (1.-x)*y*fld(m,n+1,k1)+x*y*fld(m+1,n+1,k1)) &
                + wgt*((1.-x)*(1.-y)*fld(m,n,k2)+x*(1.-y)*fld(m+1,n,k2)+ &
                       (1.-x)*y*fld(m,n+1,k2)+x*y*fld(m+1,n+1,k2))  
   end do ! i
  end do ! j
end do ! k
do k=kk+2,nzr
 read(ntape) 
end do
deallocate(fld)
if(dofliplon) deallocate(fld1)
deallocate(lonr)
deallocate(latr)
deallocate(zr)
deallocate(pr)

end

