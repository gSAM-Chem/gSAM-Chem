subroutine setsurface()

use vars
use params
use slm_vars, only: slm_init

implicit none
integer i,j,nobs
real field(nx,ny,2)
real days(2)
real(8) day1
real coef

if(readsst) then
    if(dofcast) then
      day1 = day
      day = day0
    end if
    call readsurface(sstfile,sstobs,daysstobs)
    if(nrestart.eq.0) then
      if(daysstobs(2).eq.0) then
        sstxy0(:,:) = sstobs(:,:,1)
      else
        coef=(day-daysstobs(1))/(daysstobs(2)-daysstobs(1))
        sstxy0(:,:) = sstobs(:,:,1)+coef*(sstobs(:,:,2)-sstobs(:,:,1))
      end if
      sstxy(:,:) = sstxy0(:,:) - t00
    end if
    if(dofcast) day = day1
endif

! intialize SLM

call slm_init()

end
