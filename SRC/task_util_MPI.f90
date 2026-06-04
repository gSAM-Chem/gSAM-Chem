
	subroutine task_start(rank,numtasks)
        use mpi
        use mpi_stuff, only: comm, docommsplit
        use grid, only : nsubdomains
        implicit none
	integer rank,numtasks,rc,ierr
        integer color
	call MPI_INIT(ierr)
        if(ierr .ne. 0) then
        	print *,'Error starting MPI program. Terminating.'
        	call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        	call MPI_FINALIZE(ierr)
        	stop
     	end if		
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
        if(rank.eq.0) print*,'total number of requested cores:',numtasks
        comm = MPI_COMM_WORLD
        if(numtasks.lt.nsubdomains) then
           if(rank.eq.0) print *,'number of processors is smaller than nsubdomains!',&
             '  numtasks=',numtasks,'   nsubdomains=',nsubdomains
           call task_abort()
        endif
        ! split communicator to allow running on smaller number of cores than allocated
        ! useful when number of cores that can be allocated is bigger then you need
        ! MK, May 2020
        if(.not.docommsplit) then
         if(rank.eq.0) print*,'number of processors is larger than nsubdomains!',&
            'set docommsplit to .true. in SRC/mpi_stuff if want that behavior.'
           call task_stop
        end if 
        if(numtasks.gt.nsubdomains) then
          if(rank.lt.nsubdomains) then
            color = 1
          else
            color = 0
          end if
          call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, rank, comm, ierr)
          numtasks = nsubdomains
          if(rank.ge.nsubdomains) call task_stop()
        end if
	return
	end
	
!----------------------------------------------------------------------
	
	subroutine task_abort()
        use mpi
        use mpi_stuff, only: dompi, comm	
        implicit none
	integer ierr, rc
        if(dompi) then
         call MPI_ABORT(comm,ierr,rc) 
        end if
	call task_stop()

	end
!----------------------------------------------------------------------
	
        subroutine task_abort_msg(msg)
        use mpi
        use mpi_stuff, only : rank, dompi, comm        
        implicit none
        character*(*) msg
        integer ierr, rc
        print*,rank,' exitted with msg:',msg
        if(dompi) then
          call MPI_ABORT(comm,ierr,rc)
        end if
        call task_stop()

        end

!----------------------------------------------------------------------
	subroutine task_stop()
        use mpi
        use grid, only: dompi,nstep,nstop,nelapse
        implicit none	
	integer ierr

	if(dompi) then
          call MPI_FINALIZE(ierr)	  
	endif
        stop

!	if(nstep.ge.nstop) then
!          call exit(9) ! avoid resubmission when finished
!        elseif(nelapse.eq.0) then
!           call exit(0) !bloss: clean exit condition for restart
!        else
!           call exit(1) !bloss: avoid resubmission if ending in error
!        end if

	end
!----------------------------------------------------------------------

        subroutine task_barrier()
        use mpi
        use grid, only: dompi
        use mpi_stuff, only : comm
        implicit none
        integer ierr

        if(dompi) then
          call MPI_BARRIER(comm,ierr)
        end if

        return
        end

!----------------------------------------------------------------------

        subroutine task_barrier_all()
        use mpi
        use grid, only: dompi
        use mpi_stuff, only : comm
        implicit none
	integer ierr
        
	if(dompi) then
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if  

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_real(rank_from,buffer,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_from       ! broadcasting task's rank
        real    buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer(1)).eq.sizeof(tmp)) then
         real_size = MPI_REAL
        else
         real_size = MPI_REAL8
        end if

        call MPI_BCAST(buffer,length,real_size,rank_from,comm,ierr)
        call MPI_BARRIER(comm,ierr)  ! sometime collective BCAST can hang up, call to avoid... -MK

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_float4(rank_from,buffer,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_from       ! broadcasting task's rank
        real(4) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_BCAST(buffer,length,MPI_REAL,rank_from,comm,ierr)
        call MPI_BARRIER(comm,ierr)  ! sometime collective BCAST can hang up, call to avoid... -MK

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_real8(rank_from,buffer,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_from       ! broadcasting task's rank
        real(8) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_BCAST(buffer,length,MPI_REAL8,rank_from,comm,ierr)
        call MPI_BARRIER(comm,ierr)  ! sometime collective BCAST can hang up, call to avoid... -MK
        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_integer(rank_from,buffer,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_from       ! broadcasting task's rank
        integer buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_BCAST(buffer,length,MPI_INTEGER,rank_from,comm,ierr)
        call MPI_BARRIER(comm,ierr)  ! sometime collective BCAST can hang up, call to avoid... -MK

        return
        end

!----------------------------------------------------------------------

	subroutine task_bsend_float(rank_to,buffer,length,tag)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer rank_to		! receiving task's rank
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if
	 call MPI_SEND(buffer,length,real_size,rank_to,tag,comm,ierr)
	
	return
	end

!----------------------------------------------------------------------

	subroutine task_bsend_float4(rank_to,buffer,length,tag)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer rank_to		! receiving task's rank
	real(4) buffer(*)		! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer ierr

	call MPI_SEND(buffer,length,MPI_REAL,rank_to,tag,comm,ierr)
	
	return
	end

!----------------------------------------------------------------------

        subroutine task_bsend_float8(rank_to,buffer,length,tag)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_to         ! receiving task's rank
        real(8) buffer(*)               ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer ierr

        call MPI_SEND(buffer,length,MPI_REAL8,rank_to,tag,comm,ierr)

        return
        end


!----------------------------------------------------------------------

        subroutine task_bsend_character(rank_to,buffer,length,tag)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_to         ! receiving task's rank
        character*1 buffer(*)   ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer ierr

        call MPI_SEND(buffer,length,MPI_CHARACTER,rank_to,tag,comm,ierr)

        return
        end
!----------------------------------------------------------------------

        subroutine task_gather_character(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        character*1 buffer_out(*)   ! buffer of data
        character*1 buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr

        call MPI_GATHER(buffer_out,length,MPI_CHARACTER,buffer_in,length,MPI_CHARACTER, &
                        0,comm,ierr)

        return
        end
!----------------------------------------------------------------------

        subroutine task_gatherv_character(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none

        character*1 buffer_out(*)   ! buffer of data
        character*1 buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr
        integer i
        integer :: rcounts(nsubdomains), displs(nsubdomains)
        do i=1,nsubdomains
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_GATHERV(buffer_out,length,MPI_CHARACTER,buffer_in,rcounts,displs, &
                        MPI_CHARACTER,0,comm,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_gather_float4(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(4) buffer_out(*)   ! buffer of data
        real(4) buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr

        call MPI_GATHER(buffer_out,length,MPI_REAL,buffer_in,length,MPI_REAL, &
                        0,comm,ierr)

        return
        end
!----------------------------------------------------------------------

        subroutine task_gatherv_float(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none

        real buffer_out(*)   ! buffer of data
        real buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr
        integer :: rcounts(nsubdomains), displs(nsubdomains)
        integer i
        integer :: real_size, real_size1
        real(4) tmp
        if(sizeof(buffer_out(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if
        if(sizeof(buffer_in(1)).eq.sizeof(tmp)) then
         real_size1=MPI_REAL
        else
         real_size1=MPI_REAL8
        end if
        do i=1,nsubdomains
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_GATHERV(buffer_out,length,real_size,buffer_in,rcounts,displs, &
                        real_size1,0,comm,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_gatherv_float4(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none

        real(4) buffer_out(*)   ! buffer of data
        real(4) buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer i,ierr
        integer :: rcounts(nsubdomains), displs(nsubdomains)
        do i=1,nsubdomains
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_GATHERV(buffer_out,length,MPI_REAL,buffer_in,rcounts,displs, &
                        MPI_REAL,0,comm,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_allgatherv_float(buffer_in,buffer_out,nx,ny,nsubsx,nsubsy)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real,intent(out):: buffer_out(nx*nsubsx,ny*nsubsy)  
        real,intent(in):: buffer_in(*)   ! buffer of data 
        integer,intent(in)::  nx, ny, nsubsx, nsubsy
        integer length          ! buffers' length
        integer ierr
        integer :: real_size, real_size1
        integer i,j,m,it,jt
        integer :: rcounts(nsubsx*nsubsy), displs(nsubsx*nsubsy)
        real(4), allocatable :: tmp(:)
        real(4) ttt
        length = nx*ny
        if(sizeof(buffer_out(1,1)).eq.sizeof(ttt)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if
        if(sizeof(buffer_in(1)).eq.sizeof(ttt)) then
         real_size1=MPI_REAL
        else
         real_size1=MPI_REAL8
        end if
        do i=1,nsubsx*nsubsy
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do
        allocate(tmp(nx*nsubsx*ny*nsubsy),stat=i); if(i.gt.0) call task_abort_msg("task_allgatherv_float4: alloc tmp failed!")

        call MPI_ALLGATHERV(buffer_in,length,real_size1,tmp,rcounts,displs, &
                        real_size,comm,ierr)
        m=0
        do jt=0,(nsubsy-1)*ny,ny
         do it=0,(nsubsx-1)*nx,nx
          do j=1,ny
           do i=1,nx
            m=m+1
            buffer_out(it+i,jt+j)= tmp(m)
           end do
          end do
         end do
        end do
        deallocate(tmp)

        return
        end

!----------------------------------------------------------------------

        subroutine task_scatterv_float(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none
        real buffer_out(*)   ! buffer of data (sender only, that is rank=0)
        real buffer_in(*)   ! buffer of data (local data)
        integer length          ! buffers' length
        integer ierr
        integer i
        integer :: rcounts, displs(nsubdomains), sendcounts(nsubdomains)
        integer :: real_size, real_size1
        real(4) tmp
        if(sizeof(buffer_out(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if
        if(sizeof(buffer_in(1)).eq.4) then
         real_size1=MPI_REAL
        else
         real_size1=MPI_REAL8
        end if
        rcounts = length
        do i=1,nsubdomains
          sendcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_SCATTERV(buffer_out,sendcounts,displs,real_size,buffer_in,rcounts, &
                        real_size1,0,comm,ierr)

        return
        end


!----------------------------------------------------------------------

        subroutine task_gather_integer(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer buffer_out(*)   ! buffer of data
        integer buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr

        call MPI_GATHER(buffer_out,length,MPI_INTEGER,buffer_in,length,MPI_INTEGER, &
                        0,comm,ierr)

        return
        end
!----------------------------------------------------------------------

        subroutine task_gatherv_integer(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none

        integer buffer_out(*)   ! buffer of data
        integer buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr
        integer i
        integer :: rcounts(nsubdomains), displs(nsubdomains)
        do i=1,nsubdomains
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_GATHERV(buffer_out,length,MPI_INTEGER,buffer_in,rcounts,displs, &
                        MPI_INTEGER,0,comm,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_gatherv_integer2(buffer_out,buffer_in,length)
        use mpi
        use mpi_stuff, only : comm
        use grid, only: nsubdomains
        implicit none

        integer(2) buffer_out(*)   ! buffer of data
        integer(2) buffer_in(*)   ! buffer of data (receiver only, that is rank=0)
        integer length          ! buffers' length
        integer ierr
        integer i
        integer :: rcounts(nsubdomains), displs(nsubdomains)
        do i=1,nsubdomains
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do

        call MPI_GATHERV(buffer_out,length,MPI_INTEGER2,buffer_in,rcounts,displs, &
                        MPI_INTEGER2,0,comm,ierr)

        return
        end

!----------------------------------------------------------------------
        subroutine task_allgatherv_integer(buffer_in,buffer_out,nx,ny,nsubsx,nsubsy)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer,intent(out):: buffer_out(nx*nsubsx,ny*nsubsy)
        integer,intent(in):: buffer_in(*)   ! buffer of data 
        integer,intent(in)::  nx, ny, nsubsx, nsubsy
        integer length          ! buffers' length
        integer ierr
        integer i,j,m,it,jt
        integer :: rcounts(nsubsx*nsubsy), displs(nsubsx*nsubsy)
        integer,allocatable ::  tmp(:)
        length = nx*ny
        do i=1,nsubsx*nsubsy
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do
        allocate(tmp(nx*nsubsx*ny*nsubsy),stat=i); if(i.gt.0) call task_abort_msg("task_allgatherv_integer: alloc tmp failed!")

        call MPI_ALLGATHERV(buffer_in,length,MPI_INTEGER,tmp,rcounts,displs, &
                        MPI_INTEGER,comm,ierr)
        m=0
        do jt=0,(nsubsy-1)*ny,ny
         do it=0,(nsubsx-1)*nx,nx
          do j=1,ny
           do i=1,nx
            m=m+1
            buffer_out(it+i,jt+j)= tmp(m)
           end do
          end do
         end do
        end do
        deallocate(tmp)

        return
        end
!----------------------------------------------------------------------
        subroutine task_allgatherv_integer2(buffer_in,buffer_out,nx,ny,nsubsx,nsubsy)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer(2),intent(out):: buffer_out(nx*nsubsx,ny*nsubsy)
        integer(2),intent(in):: buffer_in(*)   ! buffer of data
        integer,intent(in)::  nx, ny, nsubsx, nsubsy
        integer length          ! buffers' length
        integer ierr
        integer i,j,m,it,jt
        integer :: rcounts(nsubsx*nsubsy), displs(nsubsx*nsubsy)
        integer(2),allocatable ::  tmp(:)
        length = nx*ny
        do i=1,nsubsx*nsubsy
          rcounts(i) = length
          displs(i) = length*(i-1)
        end do
        allocate(tmp(nx*nsubsx*ny*nsubsy),stat=i); if(i.gt.0) call task_abort_msg("task_allgatherv_integer2: alloc tmp failed!")

        call MPI_ALLGATHERV(buffer_in,length,MPI_INTEGER2,tmp,rcounts,displs, &
                        MPI_INTEGER2,comm,ierr)
        m=0
        do jt=0,(nsubsy-1)*ny,ny
         do it=0,(nsubsx-1)*nx,nx
          do j=1,ny
           do i=1,nx
            m=m+1
            buffer_out(it+i,jt+j)= tmp(m)
           end do
          end do
         end do
        end do
        deallocate(tmp)

        return
        end

!----------------------------------------------------------------------

	subroutine task_send_float(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer rank_to		! receiving task's rank
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request		! request id
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_ISEND(buffer,length,real_size,rank_to,tag,comm,request,ierr)

	
	return
	end
!----------------------------------------------------------------------

        subroutine task_send_float4(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_to         ! receiving task's rank
        real(4) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request         ! request id
        integer ierr
        integer :: real_size

        call MPI_ISEND(buffer,length,MPI_REAL,rank_to,tag,comm,request,ierr)


        return
        end

!----------------------------------------------------------------------

        subroutine task_send_float8(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_to         ! receiving task's rank
        real(8) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request         ! request id
        integer ierr
        integer :: real_size

        call MPI_ISEND(buffer,length,MPI_REAL8,rank_to,tag,comm,request,ierr)


        return
        end

!----------------------------------------------------------------------

	subroutine task_send_integer(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer rank_to		! receiving task's rank
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	integer ierr

	call MPI_ISEND(buffer,length,MPI_INTEGER,rank_to,tag, &
					comm,request,ierr)

	return
	end
	
!----------------------------------------------------------------------

        subroutine task_send_integer2(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer rank_to         ! receiving task's rank
        integer(2) buffer(*)       ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request
        integer ierr

        call MPI_ISEND(buffer,length,MPI_INTEGER2,rank_to,tag, &
                                        comm,request,ierr)

        return
        end

!----------------------------------------------------------------------

	subroutine task_send_character(rank_to,buffer,length,tag,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer rank_to		! receiving task's rank
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	integer ierr

	call MPI_ISEND(buffer,length,MPI_CHARACTER,rank_to,tag, &
					comm,request,ierr)

	return
	end
	
!----------------------------------------------------------------------

        subroutine task_breceive_float(buffer,length,rank,tag)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer status(MPI_STATUS_SIZE)
	integer rank, tag
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

	call MPI_RECV(buffer,length,real_size,MPI_ANY_SOURCE, &
		MPI_ANY_TAG,comm,status,ierr)
	rank = status(MPI_SOURCE)
	tag = status(MPI_TAG)
	return
	end
!----------------------------------------------------------------------

        subroutine task_breceive_float4(buffer,length,rank,tag)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(4) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer status(MPI_STATUS_SIZE)
        integer rank, tag
        integer ierr

        call MPI_RECV(buffer,length,MPI_REAL,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,comm,status,ierr)
        rank = status(MPI_SOURCE)
        tag = status(MPI_TAG)
        return
        end

!----------------------------------------------------------------------

        subroutine task_breceive_float8(buffer,length,rank,tag)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(8) buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer status(MPI_STATUS_SIZE)
        integer rank, tag
        integer ierr

        call MPI_RECV(buffer,length,MPI_REAL8,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,comm,status,ierr)
        rank = status(MPI_SOURCE)
        tag = status(MPI_TAG)
        return
        end

!----------------------------------------------------------------------

        subroutine task_breceive_character(buffer,length,rank,tag)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        character*1 buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer status(MPI_STATUS_SIZE)
        integer rank, tag       
        integer ierr
        
        call MPI_RECV(buffer,length,MPI_CHARACTER,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,comm,status,ierr)
        rank = status(MPI_SOURCE)
        tag = status(MPI_TAG)
        return
        end

!----------------------------------------------------------------------

        subroutine task_receive_float(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer request
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

	call MPI_IRECV(buffer,length,real_size,MPI_ANY_SOURCE, &
		MPI_ANY_TAG,comm,request,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_receive_float4(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real(4) buffer(*)		! buffer of data
	integer length		! buffers' length
	integer request
	integer ierr

	call MPI_IRECV(buffer,length,MPI_REAL,MPI_ANY_SOURCE, &
		MPI_ANY_TAG,comm,request,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_receive_float8(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(8) buffer(*)               ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr

        call MPI_IRECV(buffer,length,MPI_REAL8,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,comm,request,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_receive_integer(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	integer ierr

	call MPI_IRECV(buffer,length,MPI_INTEGER,MPI_ANY_SOURCE, &
		MPI_ANY_TAG,comm,request,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_receive_integer2(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer(2) buffer(*)       ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr

        call MPI_IRECV(buffer,length,MPI_INTEGER2,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,comm,request,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_receive_character(buffer,length,request)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	integer ierr

	call MPI_IRECV(buffer,length,MPI_CHARACTER,MPI_ANY_SOURCE, &
		MPI_ANY_TAG,comm,request,ierr)

	return
	end

!----------------------------------------------------------------------
        subroutine task_wait(request,rank,tag)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	integer status(MPI_STATUS_SIZE),request
	integer rank, tag
	integer ierr
	call MPI_WAIT(request,status,ierr) 
	rank = status(MPI_SOURCE)
	tag = status(MPI_TAG)

	return
	end

!----------------------------------------------------------------------
        
        subroutine task_waitall(count,reqs,ranks,tags)
        use mpi
	use grid, only: dompi
	implicit none
 	integer count,reqs(count)
	integer stats(MPI_STATUS_SIZE,1000),ranks(count),tags(count)
	integer ierr, i
	if(dompi) then
	call MPI_WAITALL(count,reqs,stats,ierr)
        if(count.gt.1000) then
            print*,'task_waitall: count > 1000 !'
	    call task_abort()
	end if
	do i = 1,count
	  ranks(i) = stats(MPI_SOURCE,i)
	  tags(i) = stats(MPI_TAG,i)
	end do
	end if

	return
	end

!----------------------------------------------------------------------
        subroutine task_test(request,flag,rank,tag)
        use mpi
	implicit none
	integer status(MPI_STATUS_SIZE),request
	integer rank, tag
	logical flag
	integer ierr
	call MPI_TEST(request,flag,status,ierr)
	if(flag) then 
	  rank = status(MPI_SOURCE)
	  tag = status(MPI_TAG)
	endif

	return
	end

!----------------------------------------------------------------------

        subroutine task_sum_real(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real buffer_in(*)	! buffer of data
	real buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer_in(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

	call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                           real_size,MPI_SUM,comm,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_sum_real8(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real(8) buffer_in(*)	! buffer of data
	real(8) buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                         MPI_REAL8,MPI_SUM,comm,ierr)

	return
	end
!----------------------------------------------------------------------

        subroutine task_sum_integer(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer buffer_in(*)	! buffer of data
	integer buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                        MPI_INTEGER,MPI_SUM,comm,ierr)

	return
	end

!----------------------------------------------------------------------
        subroutine task_sum_integer1(int_in,int_out)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer int_in 
        integer int_out
        integer buffer_in(1)    ! buffer of data
        integer buffer_out(1)   ! buffer of data
        integer length          ! buffers' length
        integer ierr
        buffer_in(1) = int_in
        length = 1

        call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                        MPI_INTEGER,MPI_SUM,comm,ierr)
        int_out = buffer_out(1)

        return
        end

!----------------------------------------------------------------------

        subroutine task_max_real(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real buffer_in(*)	! buffer of data
	real buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer_in(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                          length,real_size,MPI_MAX,comm,ierr)

	return
        end
!----------------------------------------------------------------------

        subroutine task_max_real4(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real(4) buffer_in(*)	! buffer of data
	real(4) buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                          length,MPI_REAL,MPI_MAX,comm,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_max_real8(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(8) buffer_in(*)    ! buffer of data
        real(8) buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_ALLREDUCE(buffer_in,buffer_out, &
                          length,MPI_REAL,MPI_MAX,comm,ierr)

        return
        end

!----------------------------------------------------------------------

        subroutine task_max_integer(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer buffer_in(*)	! buffer of data
	integer buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                        length,MPI_INTEGER,MPI_MAX,comm,ierr)

	return
	end

!----------------------------------------------------------------------

        subroutine task_min_real(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real buffer_in(*)	! buffer of data
	real buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr
        integer :: real_size
        real(4) tmp

        if(sizeof(buffer_in(1)).eq.sizeof(tmp)) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                            length,real_size,MPI_MIN,comm,ierr)
	return
	end
!----------------------------------------------------------------------

        subroutine task_min_real8(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        real(8) buffer_in(*)    ! buffer of data
        real(8) buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_ALLREDUCE(buffer_in,buffer_out, &
                            length,MPI_REAL8,MPI_MIN,comm,ierr)
        return
        end

!----------------------------------------------------------------------

        subroutine task_min_real4(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	real(4) buffer_in(*)	! buffer of data
	real(4) buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                            length,MPI_REAL,MPI_MIN,comm,ierr)
	return
	end
!----------------------------------------------------------------------

        subroutine task_min_integer(buffer_in,buffer_out,length)
        use mpi
        use mpi_stuff, only : comm
	implicit none
	
	integer buffer_in(*)	! buffer of data
	integer buffer_out(*)	! buffer of data
	integer length		! buffers' length
	integer ierr

	call MPI_ALLREDUCE(buffer_in,buffer_out, &
                  length,MPI_INTEGER,MPI_MIN,comm,ierr)

	return
	end
!----------------------------------------------------------------------

!----------------------------------------------------------------------

        subroutine task_bcast_fourdim_array_real(rank_from,array,n1,n2,n3,n4)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real, intent(inout) :: array(n1,n2,n3,n4) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4        ! dimension lengths
        integer ierr

        real, allocatable :: rtmp(:)
        integer :: count, nn, myrank, i,j,k,m
        real :: tmp1, tmp2

        nn = n1*n2*n3*n4
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_fourdim_array_real: alloc rtmp failed!")
        if(ierr.ne.0) then
          write(*,*) 'Error in allocating array in task_bcast_fourdim_array'
          call task_abort()
        end if

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do m = 1,n4
            do k = 1,n3
              do j = 1,n2
                do i = 1,n1
                  rtmp(count) = array(i,j,k,m)
                  count = count + 1
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do m = 1,n4
            do k = 1,n3
              do j = 1,n2
                do i = 1,n1
                  array(i,j,k,m) = rtmp(count)
                  count = count + 1
                end do
              end do
            end do
          end do
        end if

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,991) myrank, tmp1, tmp2
!!$        991 format('Consistency check in task_bcast_fourdim: rank/sum/sum2 = ',I4,2E16.8)


        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_fourdim_array_real8(rank_from,array,n1,n2,n3,n4)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i,j,k,m
        real(8) :: tmp1, tmp2

        nn = n1*n2*n3*n4
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_fourdim_array_real8: alloc rtmp failed!")
        if(ierr.ne.0) then
          write(*,*) 'Error in allocating array in task_bcast_fourdim_array'
          call task_abort()
        end if

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do m = 1,n4
            do k = 1,n3
              do j = 1,n2
                do i = 1,n1
                  rtmp(count) = array(i,j,k,m)
                  count = count + 1
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do m = 1,n4
            do k = 1,n3
              do j = 1,n2
                do i = 1,n1
                  array(i,j,k,m) = rtmp(count)
                  count = count + 1
                end do
              end do
            end do
          end do
        end if

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,992) myrank, tmp1, tmp2
!!$        992 format('Consistency check in task_bcast_fourdim: rank/sum/sum2 = ',I4,2E16.8)


        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_fivedim_array_real8(rank_from,array,n1,n2,n3,n4,n5)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4,n5) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4,n5        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i1,i2,i3,i4,i5
        real :: tmp1, tmp2

        nn = n1*n2*n3*n4*n5
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_fivedim_array_real8: alloc rtmp failed!")
        if(ierr.ne.0) then
          write(*,*) 'Error in allocating array in task_bcast_fivedim_array_real8'
          call task_abort()
        end if

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do i5 = 1,n5
            do i4 = 1,n4
              do i3 = 1,n3
                do i2 = 1,n2
                  do i1 = 1,n1
                    rtmp(count) = array(i1,i2,i3,i4,i5)
                    count = count + 1
                  end do
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do i5 = 1,n5
            do i4 = 1,n4
              do i3 = 1,n3
                do i2 = 1,n2
                  do i1 = 1,n1
                    array(i1,i2,i3,i4,i5) = rtmp(count)
                    count = count + 1
                  end do
                end do
              end do
            end do
          end do
        end if

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,993) myrank, tmp1, tmp2
!!$        993 format('Consistency check in task_bcast_fivedim: rank/sum/sum2 = ',I4,2E16.8)

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_sixdim_array_real8(rank_from,array,n1,n2,n3,n4,n5,n6)
        use mpi
        use mpi_stuff, only : comm
        implicit none

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4,n5,n6) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4,n5,n6        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i1,i2,i3,i4,i5,i6
        real(8) :: tmp1, tmp2

        nn = n1*n2*n3*n4*n5*n6
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_sixdim_array_real8: alloc rtmp failed!")
        if(ierr.ne.0) then
          write(*,*) 'Error in allocating array in task_bcast_sixdim_array_real8'
          call task_abort()
        end if

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do i6 = 1,n6
            do i5 = 1,n5
              do i4 = 1,n4
                do i3 = 1,n3
                  do i2 = 1,n2
                    do i1 = 1,n1
                      rtmp(count) = array(i1,i2,i3,i4,i5,i6)
                      count = count + 1
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do i6 = 1,n6
            do i5 = 1,n5
              do i4 = 1,n4
                do i3 = 1,n3
                  do i2 = 1,n2
                    do i1 = 1,n1
                      array(i1,i2,i3,i4,i5,i6) = rtmp(count)
                      count = count + 1
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if


!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,994) myrank, tmp1, tmp2
!!$        994 format('Consistency check in task_bcast_sixdim: rank/sum/sum2 = ',I4,2E16.8)

        return
        end


