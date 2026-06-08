
!=========================================================================
! This is a collection of some utilities.

integer function lenstr (string)
      
! returns string's length ignoring the rightmost blank and null characters

implicit none
character(*) string
integer k
lenstr = 0
do k = 1,len(string)
 if (string(k:k).ne.' '.and.string(k:k).ne.char(0)) then
     lenstr = lenstr+1
 end if
end do
111  return
end

!===========================================================

subroutine averageXY(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)

! average scalar field
use grid, only: nx,ny,nzm,wgt
use terrain, only: terra
implicit none
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm) 
real(8) ff,factor
integer i,j,k	
factor = 1./dble(nx*ny)
do k =1,nzm
 ff = 0.
 do j =1,ny
  do i =1,nx
    ff = ff + f(i,j,k)*wgt(j,k)*terra(i,j,k)
  end do
 end do
 ff = ff*factor
 fm(k) = real(ff)
end do 
end

!===========================================================

subroutine averageXY_SFC(f,dimx1,dimx2,dimy1,dimy2,fm)

! average scalar field
use grid, only: dompi,nsubdomains,nx,ny,wgty
use terrain, only: terra
implicit none
integer dimx1, dimx2, dimy1, dimy2
real f(dimx1:dimx2, dimy1:dimy2),fm
real(8) fm1(1),fm2(1),factor
integer i,j
factor = 1./dble(nx*ny)
fm1(1) = 0.
do j =1,ny
 do i =1,nx
   fm1(1) = fm1(1) + f(i,j)*wgty(j)
 end do
end do
fm1(1) = fm1(1) * factor
if(dompi) then
 fm2(1) = fm1(1)
 call task_sum_real8(fm2,fm1,1)
 fm=real(fm1(1)/dble(nsubdomains))
else
 fm=real(fm1(1))
endif
end


!===========================================================
	
subroutine averageXY_MPI(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)

! average scalar field
use grid, only: dompi,nsubdomains,nx,ny,nzm,wgt
use terrain, only: terra
implicit none
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm)
real(8) fm1(nzm),fm2(nzm),factor
integer i,j,k
factor = 1./dble(nx*ny)
do k =1,nzm
 fm1(k) = 0.
 do j =1,ny
  do i =1,nx
    fm1(k) = fm1(k) + f(i,j,k)*wgt(j,k)*terra(i,j,k)
  end do
 end do
 fm1(k) = fm1(k) * factor
end do
if(dompi) then
 do k =1,nzm
   fm2(k) = fm1(k)
 end do
 call task_sum_real8(fm2,fm1,nzm)
 do k=1,nzm
  fm(k)=real(fm1(k)/dble(nsubdomains))
 end do
else
 do k=1,nzm
  fm(k)=real(fm1(k))
 end do
endif
end


		
!===========================================================
	
	
subroutine fminmax_print(name,f,dimx1,dimx2,dimy1,dimy2,dimz)

use grid
implicit none
integer, intent(in) ::  dimx1, dimx2, dimy1, dimy2, dimz
real, intent(in) :: f(dimx1:dimx2, dimy1:dimy2, 1:dimz)
character(*) name
real fmin(1),fmax(1),fff(1),fsum(1)
integer i,j,k
	
fmin(1) = minval(f(dimx1:dimx2, dimy1:dimy2, 1:dimz))
fmax(1) = maxval(f(dimx1:dimx2, dimy1:dimy2, 1:dimz))
	
if(dompi) then
  fff(1)=fmax(1)
  call task_max_real(fff(1),fmax(1),1)
  fff(1)=fmin(1)
  call task_min_real(fff(1),fmin(1),1)
end if
if(masterproc) print *,name,fmin,fmax 
end

!===========================================================
	
subroutine fminmax_printm(name,f,dimx1,dimx2,dimy1,dimy2,dimz,mask1)

use grid
implicit none
integer, intent(in) ::  dimx1, dimx2, dimy1, dimy2, dimz
real, intent(in) :: f(dimx1:dimx2, dimy1:dimy2, 1:dimz)
logical, intent(in) :: mask1(dimx1:dimx2, dimy1:dimy2, 1:dimz)
character(*) name
real fmin(1),fmax(1),fff(1),fsum(1)
integer i,j,k

fmin(1) = minval(f(dimx1:dimx2, dimy1:dimy2, 1:dimz), mask = mask1)
fmax(1) = maxval(f(dimx1:dimx2, dimy1:dimy2, 1:dimz), mask = mask1)

if(dompi) then
  fff(1)=fmax(1)
  call task_max_real(fff(1),fmax(1),1)
  fff(1)=fmin(1)
  call task_min_real(fff(1),fmin(1),1)
end if
if(masterproc) print *,name,fmin,fmax
end


!===========================================================

! re-implementation of sum() intrinsic function
! which often is not accurate for summing large numbers

real(8) function total_sum(f,dimx1,dimx2,dimy1,dimy2,dimz)
implicit none
integer, intent(in) ::  dimx1, dimx2, dimy1, dimy2, dimz
real, intent(in) :: f(dimx1:dimx2, dimy1:dimy2, 1:dimz)
integer i,j,k
total_sum = 0._8
do k=1,dimz
 do j=dimy1,dimy2
  do i=dimx1,dimx2
   total_sum = total_sum + real(f(i,j,k),8)
  end do
 end do
end do

end

!===========================================================
	
subroutine setvalue(f,n,f0)
implicit none
integer n
real f(n), f0
integer k
do k=1,n
 f(k)=f0
end do
end

!===========================================================

! determine number of byte in a record in direct access files (can be anything, from 1 to 8):
! can't assume 1 as it is compiler and computer dependent
integer function bytes_in_rec()
implicit none
character(8) str
integer n, err
open(1,status ='scratch',access ='direct',recl=1)
do n = 1,8
 write(1,rec=1,iostat=err) str(1:n)
 if (err.ne.0) exit
 bytes_in_rec = n
enddo
close(1,status='delete')
end

!===========================================================
!----------------------------------------
! average a global 3D array alone x direction
! with the result known to each subdomain

subroutine mean_x_3D(in, out)

use grid
implicit none
real, intent(in) ::  in(nx,ny,nzm)
real, intent(out) :: out(ny,nzm)
real, allocatable ::  buf(:,:,:)
integer i,j,k,m
integer n_in, count, reqs_in(max(1,nsubdomains_x-1)), rnk, tag
integer n_out, reqs_out(nsubdomains_x-1), tags(nsubdomains_x-1), ranks(nsubdomains_x-1)
logical flag(nsubdomains_x)

! average locally first:
do k=1,nzm
 do j=1,ny
  out(j,k) = sum(in(1:nx,j,k))
 end do
end do
out(:,:)=out(:,:)/dble(nx)

if(dompi) then


 if(mod(rank,nsubdomains_x).ne.0) then
   call task_bsend_float(rank-mod(rank,nsubdomains_x),out(1:ny,1:nzm),ny*nzm, 31)
   call task_breceive_float(out(1:ny,1:nzm),ny*nzm,rnk,tag)
 else
   allocate(buf(ny,nzm,max(1,nsubdomains_x-1)))
   n_in = 0
   do i=1,nsubdomains_x-1
     n_in = n_in + 1
     call task_receive_float(buf(1:ny,1:nzm,n_in),ny*nzm,reqs_in(n_in))
     flag(n_in) = .false.
   end do
   count = n_in
   do while (count .gt. 0)
    do m = 1,n_in
     if(.not.flag(m)) then
         call task_test(reqs_in(m), flag(m), rnk, tag)
         if(flag(m)) then
           count=count-1
           out(1:ny,1:nzm) = out(1:ny,1:nzm) + buf(1:ny,1:nzm,m)
         endif
     endif
    end do
   end do
   out = out/dble(nsubdomains_x)
   n_out = 0.
   do i=1,nsubdomains_x-1
     n_out=n_out+1
     call task_send_float(rank+i,out(1:ny,1:nzm),ny*nzm,31,reqs_out(n_out))
   end do
   call task_waitall(n_out,reqs_out,tags,ranks)
   deallocate(buf)
 end if

end if

call task_barrier()

end


!===========================================================

!----------------------------------------
! average a global 2D array (x-y) alone x direction
! with the result known to each subdomain

subroutine mean_x_2D(in, out)

use grid
implicit none
real, intent(in) ::  in(nx,ny)
real, intent(out) :: out(ny)
real buf(ny,max(1,nsubdomains_x-1))
integer i,j,k,m
integer n_in, count, reqs_in(max(1,nsubdomains_x-1)), rnk, tag
integer n_out, reqs_out(nsubdomains_x-1), tags(nsubdomains_x-1), ranks(nsubdomains_x-1)
logical flag(nsubdomains_x)

! average locally first:
do j=1,ny
 out(j) = sum(in(1:nx,j))/dble(nx)
end do

if(dompi) then

 if(mod(rank,nsubdomains_x).ne.0) then
   call task_bsend_float(rank-mod(rank,nsubdomains_x),out(1:ny),ny, 31)
   call task_breceive_float(out(1:ny),ny,rnk,tag)
 else
   n_in = 0
   do i=1,nsubdomains_x-1
     n_in = n_in + 1
     call task_receive_float(buf(1:ny,n_in),ny,reqs_in(n_in))
     flag(n_in) = .false.
   end do
   count = n_in
   do while (count .gt. 0)
    do m = 1,n_in
     if(.not.flag(m)) then
         call task_test(reqs_in(m), flag(m), rnk, tag)
         if(flag(m)) then
           count=count-1
           out(1:ny) = out(1:ny) + buf(1:ny,m)
         endif
     endif
    end do
   end do
   out = out/dble(nsubdomains_x)
   n_out=0.
   do i=1,nsubdomains_x-1
     n_out = n_out+1
     call task_send_float(rank+i,out(1:ny),ny,31,reqs_out(n_out))
   end do
   call task_waitall(n_out,reqs_out,tags,ranks)
 end if

end if

call task_barrier()
end

!===========================================================

subroutine systemf(command)
character(*) command
call system(command)
end

!===========================================================

subroutine write_ready_to_restart_file(string)
  implicit none
  character(LEN=*), intent(in) :: string
  logical :: log_Exists

  !-------------------------------------------------------------------
!bloss: Make sure that the run will not restart if it does not stop cleanly
  inquire(file='./RESTART/ReadyForRestart',exist=log_Exists)
  if(.NOT.log_Exists) then
    open(unit=47,file='./RESTART/ReadyForRestart',form='FORMATTED',status='NEW')
    close(47)
  end if

  open(unit=47,file='./RESTART/ReadyForRestart',form='FORMATTED',status='REPLACE')
  rewind(47)
  write(47,*) TRIM(string)
  close(47)

end subroutine write_ready_to_restart_file

!===========================================================

subroutine get_system_time_real8(time_in_seconds)
  implicit none
  real(8), intent(out) :: time_in_seconds
  integer(8) :: count, count_rate
  call system_clock(count, count_rate)
  time_in_seconds = real(count,KIND=8)/real(count_rate,KIND=8)
end subroutine get_system_time_real8
