! Simple randaom number generator in the range [0,1]
! ranset_(iseed) initializes with iseed
! ranf_() returns next random numer




      real(4) function ranf_()
       implicit none
        call random_number(ranf_)
       return
      end


      subroutine ranset_(iseed)
      implicit none
      real(4) ranf_
      integer iseed, i, m, nsteps
      nsteps = iseed*1000
      do i = 1,nsteps
	m = ranf_()
      end do	
      return
      end




