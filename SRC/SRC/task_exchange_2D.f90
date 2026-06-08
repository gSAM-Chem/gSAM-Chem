! Non-blocking receives before blocking sends for efficiency on IBM SP.



subroutine task_exchange_2D(f, dimx1, dimx2, dimz, i_1, i_2, idir)
	
! sends and receives the boundary messages for 2D ayyats (x-z and y-z)	
	
use grid
use mpi_stuff
use advection, only: NADVS
implicit none
	
integer dimx1, dimx2, dimz
integer i_1, i_2
integer idir ! =1 x-z array, =2 y-z array
real f(dimx1:dimx2,dimz)
	
integer bufflen
parameter (bufflen = max(nx,ny)*(3+NADVS)*nz)
real buff_send(bufflen)	! buff to send data
real buff_recv(bufflen,2)	! buff to receive data
integer reqs_in(2), ids(2), count  
integer i_start(2), i_end(2), ranks(2), nn
logical flag(2)

integer i, k, n, rf, tag, m
integer i1, i2, mask
	
        if(idir.ne.1.and.idir.ne.2) then
         print*,'task_exchange_2D: idir =',idir,' Exitting...'
         stop
        end if
	i1 = i_1 - 1
	i2 = i_2 - 1

        if(idir.eq.1) then
	  ranks(1) = rankee
	  ranks(2) = rankww
	  i_start(1) = nx-i1
	  i_end  (1) = nx
          nn = nx
        else
	  ranks(1) = ranknn
	  ranks(2) = rankss
	  i_start(1) = ny-i1
	  i_end  (1) = ny
          nn = ny
        end if
	i_start(2) = 1
	i_end  (2) = 1+i2

	mask=i1*10000+i2*1000+1*100+1*10
	ids(1) = mask+2
	ids(2) = mask+1
!----------------------------------------------------------------------
!  Send/receive buffs to/from neighbors 
!----------------------------------------------------------------------


! Non-blocking receives first:

	do m = 1,2

	 if(rank.ne.ranks(m)) then 

          call task_receive_float(buff_recv(1,m),bufflen,reqs_in(m))
	  flag(m) = .false.

	 else

          flag(m) = .true.

         end if

	end do ! m

! Blocking sends:

	do m = 1,2

	   n=0
	   do k=1,dimz
	       do i=i_start(m),i_end(m)
	         n = n+1
	         buff_send(n) = f(i,k)
	       end do
           end do

	if(rank.ne.ranks(m)) then 

	  call task_bsend_float(ranks(m),buff_send,n,ids(m))

	 else

          do i=1,n
            buff_recv(i,m) =  buff_send(i)
          end do

	 end if

	end do ! m


! Fill the data from the buffers that have the same rank (were not sent):

	count = 0
	do m = 1,2
	 if(flag(m)) then
           count = count+1
           tag = ids(m)
           call task_assign_bnd_2D(f,dimx1,dimx2,dimz,nn,buff_recv(1,m),tag)
	 end if 
	end do


! Monitor the progress of receiving; fill the data:


        do while (count .lt. 2)
	  do m = 1,2
	   if(.not.flag(m)) then
	    call task_test(reqs_in(m),flag(m),rf,tag)
	    if(flag(m)) then 
	      count=count+1
	      call task_assign_bnd_2D(f,dimx1,dimx2,dimz,nn,buff_recv(1,m),tag)
	    endif   
	   endif
	  end do
	end do

 call task_barrier()

end




subroutine task_assign_bnd_2D(f,dimx1,dimx2,dimz,nn,buff,tag)

! this routine assignes the boundary info after MPI exchange

use grid
implicit none

integer dimx1, dimx2, dimz
real f(dimx1:dimx2, dimz)

real buff(*)    ! buff for sending data
integer nn,tag

integer i, k, n, proc
integer i1, i2

!       The dimensions of the fields in common com3d. Needed by MPI

!       Decode the tag:

i1 = tag/10000
i2 = (tag-i1*10000)/1000
proc = tag-i1*10000-i2*1000-1*100-1*10

          if(proc.eq.1) then

             n=0
             do k=1,dimz
                 do i=nn+1,nn+1+i2
                   n = n+1
                   f(i,k) = buff(n)
                 end do
             end do

          elseif(proc.eq.2) then

             n=0
             do k=1,dimz
                 do i=-i1,0
                   n = n+1
                   f(i,k) = buff(n)
                 end do
             end do


          endif


end

