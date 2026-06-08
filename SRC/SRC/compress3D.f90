subroutine compress3D (f,nx,ny,nz,name, long_name, units, &
       	               savebin, dompi, rank, nsubdomains, nfiles, gtype)


! Compress (or not) a given 3D array into the 2-byte integer array
! and writes  into a file, or write directly to netcdf file.

use grid, only: masterproc
use pnetcdf_stuff, only: add_var3D_pnetcdf,dopnetcdf

implicit none

! Input:

real(4), intent(in) :: f(nx,ny,nz)
integer, intent(in) :: nx,ny,nz
character(*), intent(in) :: name,long_name,units
logical, intent(in) :: savebin, dompi
integer, intent(in) :: rank,nsubdomains,nfiles
integer gtype ! type of grid: 1 center of cell, 2 z-inteface, 3 y-inteface, 4 x-interface 

! Local:

integer(2) byte(nx*ny*nz)
real(4) byte4(nx*ny*nz)
integer size,count

integer rrr,ttt,irank
integer integer_max, integer_min
parameter (integer_min=-32000, integer_max=32000)
!	parameter (integer_min=-127, integer_max=127)
real(4) f_max(nz),f_min(nz), f_max1(nz), f_min1(nz), scale
integer i,j,k,req,rankw,nrankw
integer varid, it, jt

! Output is in NetCDF format:

if(dopnetcdf) then
  call add_var3D_pnetcdf(f,name,long_name,units,gtype)
  return
end if

! Output is in internal SAM format:

size=nx*ny*nz
count = 0


if(masterproc) then
             if(gtype.ne.1.and.gtype.ne.2.and.gtype.ne.3.and.gtype.ne.4) then
              print*,'wrong gtype in compress3D: gtype =',gtype
              stop
             end if
           !  print*, trim(name)
             write(46) name,' ',long_name,' ',units
             write(46) savebin, gtype
end if 
if(savebin) then	

  do k=1,nz
   do j=1,ny
    do i=1,nx
     count = count+1
     byte4(count) = f(i,j,k)
    end do
   end do
  end do

  rankw = rank-mod(rank,nsubdomains/nfiles)
  nrankw = nsubdomains/nfiles
  if(rank.eq.rankw) write(46) (byte4(k),k=1,count)
  do irank = rankw+1, rankw+nrankw-1
    call task_barrier()
    if(irank.eq.rank) then
      call task_send_float4(rankw,byte4,count,irank,req)
      call task_wait(req,rrr,ttt)
    end if
    if(rank.eq.rankw) then
      call task_receive_float4(byte4,count,req)
      call task_wait(req,rrr,ttt)
      write(46) (byte4(k),k=1,count)
    end if
  end do


else

!  allocate (byte(size))

  do k=1,nz
   f_max(k)=maxval(f(:,:,k))
   f_min(k)=minval(f(:,:,k))
  end do

  if(dompi) then
    f_max1=f_max
    f_min1=f_min
    call task_max_real4(f_max1,f_max,nz)
    call task_min_real4(f_min1,f_min,nz)
  endif

  do k=1,nz

   if(f_max(k)-f_min(k).lt.1.e-10) then
      scale = 0.
   else
      scale = float(integer_max-integer_min)/(f_max(k)-f_min(k))
   end if

   do j=1,ny
    do i=1,nx
      count=count+1
      byte(count)= integer_min+scale*(f(i,j,k)-f_min(k))   
    end do
   end do

  end do ! k
        
  if(masterproc) then
     write(46) f_max,f_min
  end if

  rankw = rank-mod(rank,nsubdomains/nfiles)
  nrankw = nsubdomains/nfiles
  if(rank.eq.rankw) write(46) (byte(k),k=1,count)
  do irank = rankw+1, rankw+nrankw-1
    call task_barrier()
    if(irank.eq.rank) then
      call task_send_integer2(rankw,byte,count,irank,req)
      call task_wait(req,rrr,ttt)
    end if
    if(rank.eq.rankw) then
      call task_receive_integer2(byte,count,req)
      call task_wait(req,rrr,ttt)
      write(46) (byte(k),k=1,count)
    end if
  end do

end if ! savebin

call task_barrier()

end subroutine compress3D
	
