	
	
subroutine bound_uv

! Periodic boundary exchange 

use vars
implicit none
         
integer i,j,k

	  do k=1,nzm
	   do j=1,ny
	     u(nxp1,j,k) = u(1,j,k)
	   end do
	  end do

	  if(RUN3D) then

	    do k=1,nzm
	     do i=1,nx
	      v(i,nyp1,k) = v(i,1,k)
	     end do
	    end do

	  endif
	
end subroutine bound_uv
