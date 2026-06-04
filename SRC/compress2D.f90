subroutine compress2D (f, nx, ny, name, long_name, units, &
       	               savebin, dompi, rank, nsubdomains)


! Compress (or not) a given 2D array into the 2byte integer array
! and write to a file, or write netcdf file instead.

use grid, only: nx_gl,ny_gl
use mpi_stuff, only: masterproc
use pnetcdf_stuff, only: add_var2D_pnetcdf,dopnetcdf

implicit none
! Input:

real(4), intent(in) :: f(nx,ny)  
integer, intent(in) ::  nx,ny
character(*), intent(in) :: name,long_name,units
integer, intent(in) :: rank,nsubdomains
logical, intent(in) :: savebin, dompi

! Local:

integer(2), allocatable :: byte(:)
real(4), allocatable :: byte4(:)
integer(2), allocatable :: byte_2D(:)
real(4), allocatable :: byte4_2D(:)

character(7) form
integer integer_max, integer_min
parameter (integer_min=-32000, integer_max=32000)
!	parameter (integer_min=-127, integer_max=127)
real(4) f_max(1),f_min(1), f_max1(1), f_min1(1), scale
integer i,j,m,it,jt,rrr,ttt,irank
integer n_in, count, reqs_in(max(1,nsubdomains-1)), tag
logical flag(nsubdomains-1)

! Output is in NetCDF format:

if(dopnetcdf) then
  call add_var2D_pnetcdf(f,name,long_name,units)
  return
end if

! Output is in internal SAM format:

if(masterproc) then
    write(46) name,' ',long_name,' ',units
    write(46) savebin
end if 

if(savebin) then	

  allocate (byte4(nx*ny))
  allocate(byte4_2D(nx_gl*ny_gl))
  call task_barrier()

   count = 0
   do j=1,ny
    do i=1,nx
     count = count+1
     byte4(count) = f(i,j)
    end do
   end do


  if(.not.dompi) then
    write(46) byte4(:)
  else
    call task_gatherv_float4(byte4,byte4_2D,nx*ny)
    if(masterproc) write(46) byte4_2D(:)
  end if

  call task_barrier()
  deallocate(byte4)
  deallocate(byte4_2D)


else

   allocate (byte(nx*ny))
   allocate(byte_2D(nx_gl*ny_gl))
   call task_barrier()

   f_max(1) = maxval(f)
   f_min(1) = minval(f)
   if(dompi) then
     f_max1=f_max
     f_min1=f_min
     call task_max_real4(f_max1,f_max,1)
     call task_min_real4(f_min1,f_min,1)
   endif

   if(abs(f_max(1)).lt.10..and.abs(f_min(1)).lt.10.) then
          form='(f10.7)'
   else if(abs(f_max(1)).lt.100..and.abs(f_min(1)).lt.100.) then
          form='(f10.6)'
   else if(abs(f_max(1)).lt.1000..and.abs(f_min(1)).lt.1000.) then
          form='(f10.5)'
   else if(abs(f_max(1)).lt.10000..and.abs(f_min(1)).lt.10000.) then
          form='(f10.4)'
   else if(abs(f_max(1)).lt.100000..and.abs(f_min(1)).lt.100000.) then
          form='(f10.3)'
   else if(abs(f_max(1)).lt.1000000..and.abs(f_min(1)).lt.1000000.) then
          form='(f10.2)'
   else if(abs(f_max(1)).lt.10000000..and.abs(f_min(1)).lt.10000000.) then
          form='(f10.1)'
   else if(abs(f_max(1)).lt.100000000..and.abs(f_min(1)).lt.100000000.) then
          form='(f10.0)'
   else
          form='(f10.0)'
!          f_min=-999.
!          f_max= 999.
   end if

   if(f_max(1)-f_min(1).lt.1.e-10) then
      scale = 0.
   else
      scale = float(integer_max-integer_min)/(f_max(1)-f_min(1))
   end if

   count = 0
   do j=1,ny
    do i=1,nx
      count = count+1
      byte(count)= integer_min+scale*(f(i,j)-f_min(1))   
    end do
   end do

  if(masterproc) then
     write(46) real(f_max(1),4),real(f_min(1),4)
    ! print*,'writing ',trim(long_name)
  end if

  if(.not.dompi) then
    write(46) byte(:)
  else
    call task_gatherv_integer2(byte,byte_2D,nx*ny)
    if(masterproc) write(46) byte_2D(:)
  end if

  call task_barrier()
  deallocate(byte)
  deallocate(byte_2D)

end if ! savebin

end subroutine compress2D
	
