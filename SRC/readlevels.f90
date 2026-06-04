#include "fppmacros"

subroutine readlevels(filename,field)

! read layered scalar data

use grid
use params, only: read_meters
implicit none
character*(*), intent(in) ::  filename
real, intent(out) ::  field(nx,ny,nzm)   ! output field

integer i,k,nn,nobs,nx1,ny1,nz1
real fld(nx,ny,nzm)
real(4), allocatable :: days(:)
real coef, dayy

if(masterproc) print*,'reading data from ',trim(filename)
open(88,file=trim(filename),status='old',form='unformatted', &
        BUFFEREDYES ACTION='READ')
read(88) nx1,ny1,nz1
read(88) nobs
if(masterproc)print*,'nx ny nz:',nx1,ny1,nz1
allocate (days(nobs))
read(88) days(1:nobs)
if(days(1).gt.365.) days(1:nobs) = days(1:nobs) - (year0-1)*365
nn=1
do i=1,nobs-1
 if(day.gt.days(i)) then
   nn=i
 endif
end do
do i=1,nn-1
 read(88)
 read(88)
 read(88)
 read(88)
 do k=1,nz1
   read(88)
 end do
end do
call read_field3D (88,field,nx1,ny1,nz1,lon_gl,lat_gl,z,nx,ny,nzm)
dayy = day 
if(nobs.gt.1) then
 if(read_meters) then
  call read_field3D (88,fld,nx1,ny1,nz1,x_gl,y_gl,z,nx,ny,nzm)
 else
  call read_field3D (88,fld,nx1,ny1,nz1,lon_gl,lat_gl,z,nx,ny,nzm)
 end if
 if(masterproc) print*,'interpolating: day=',dayy,'day1=',days(nn),'day2=',days(nn+1)
 coef=(dayy-days(nn))/(days(nn+1)-days(nn))
 field(:,:,:) = field(:,:,:)+coef*(fld(:,:,:)-field(:,:,:))
end if
close (88)
call fminmax_print('after interpolating:',field,1,nx,1,ny,nzm)
deallocate(days)
end




