
subroutine radiation()

!	Radiation interface

use grid
use params, only: dosmoke, doradsimple
implicit none

if(doradsimple) then

!  A simple predefined radiation (longwave only)

    if(dosmoke) then
       call rad_simple_smoke()
    else
       call rad_simple()
    end if
	 
else


! Call full radiation package:
 
call t_startf ('radiation')

    call rad_full()	
 
call t_stopf ('radiation')

endif


end


