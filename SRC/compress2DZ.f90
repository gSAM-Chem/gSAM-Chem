subroutine compress2DZ (f,ny,nz,name, long_name, units, &
       	               savebin, dompi, rank, nsubdomains, nsubdomains_x,ztype)


! Compress (or not) a given 2DZ array into the 2-byte integer -array
! and write it to a file.

use grid, only: masterproc
use pnetcdf_stuff, only: add_var2DZ_pnetcdf,dopnetcdf

implicit none

! Input:

real(4), intent(in) :: f(ny,nz)
integer, intent(in) :: ny,nz
character(*), intent(in) :: name,long_name,units
integer, intent(in) :: rank,nsubdomains,nsubdomains_x,ztype
logical, intent(in) :: savebin, dompi

! Local:

integer(2), allocatable :: byte(:)
real(4), allocatable :: byte4(:)
integer size,count,rrr,ttt,irank

real(4) value_min(nz), value_max(nz)
integer integer_max, integer_min
parameter (integer_min=-32000, integer_max=32000)
!	parameter (integer_min=-127, integer_max=127)
real(4) f_max(1),f_min(1), f_max1(1), f_min1(1), scale
integer i,j,k,req,rankw
! Output is in NetCDF format:

if(dopnetcdf) then
  call add_var2DZ_pnetcdf(f,name,long_name,units)
  return
end if

! Output is in internal SAM format:

! Allocate byte array:

size=ny*nz
count = 0

if(masterproc) then
             if(ztype.ne.1.and.ztype.ne.2) then
              print*,'wrong ztype in compress3D: ztype =',ztype
              stop
             end if
           !  print*, trim(name)
             write(46) name,' ',long_name,' ',units
             write(46) savebin,real(ztype,4)
end if 
if(savebin) then	

  allocate (byte4(size))
  call task_barrier()

  do k=1,nz
   do j=1,ny
     count = count+1
     byte4(count) = f(j,k)
   end do
  end do


  rankw = 0
  if(rank.eq.rankw) write(46) (byte4(k),k=1,count)
  do irank = 1, nsubdomains-1
    call task_barrier()
    if(mod(irank,nsubdomains_x).eq.0) then
      if(irank.eq.rank) then
        call task_send_float4(rankw,byte4,count,irank,req)
        call task_wait(req,rrr,ttt)
      end if
      if(rank.eq.rankw) then
        call task_receive_float4(byte4,count,req)
        call task_wait(req,rrr,ttt)
        write(46) (byte4(k),k=1,count)
      end if
    end if
  end do

  call task_barrier()
  deallocate(byte4)


else

  allocate (byte(size))

  do k=1,nz

   f_max=-1.e30
   f_min= 1.e30
   do j=1,ny
     f_max(1) = max(f_max(1),f(j,k)) 	    
     f_min(1) = min(f_min(1),f(j,k)) 	    
   end do
   if(dompi) then
     f_max1=f_max
     f_min1=f_min
     call task_max_real4(f_max1(1),f_max(1),1)
     call task_min_real4(f_min1(1),f_min(1),1)
   endif

   value_max(k) = f_max(1)
   value_min(k) = f_min(1)

   if(f_max(1)-f_min(1).lt.1.e-10) then
      scale = 0.
   else
      scale = float(integer_max-integer_min)/(f_max(1)-f_min(1))
   end if

   do j=1,ny
      count=count+1
      byte(count)= integer_min+scale*(f(j,k)-f_min(1))   
   end do

  end do ! k
        
  if(masterproc) then
     write(46) value_max,value_min
  end if

  rankw = 0
  if(rank.eq.rankw) write(46) (byte(k),k=1,count)
  do irank = 1, nsubdomains-1
    call task_barrier()
    if(mod(irank,nsubdomains_x).eq.0) then
      if(irank.eq.rank) then
        call task_send_integer2(rankw,byte,count,irank,req)
        call task_wait(req,rrr,ttt)
      end if
      if(rank.eq.rankw) then
        call task_receive_integer2(byte,count,req)
        call task_wait(req,rrr,ttt)
        write(46) (byte(k),k=1,count)
      end if
    end if
  end do

  call task_barrier()
  deallocate(byte)

end if ! savebin



end subroutine compress2DZ
	
