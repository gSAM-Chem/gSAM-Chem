! Initialize MISC output

subroutine stat_Minit()

use vars
implicit none

     if(mod(nstep,nsave2D).eq.0) then
      cflz_max_xy(:,:) = 0.
      cflh_max_xy(:,:) = 0. 
      cfl_max_xy(:,:) = 0. 
      zcflz_max_xy(:,:) = 0.
      zcflh_max_xy(:,:) = 0.
      zcfl_max_xy(:,:) = 0.
      sst_min_xy(:,:) = 1000.
      sst_max_xy(:,:) = 0.
     end if
     albvis_xy(:,:) = 0.
     albnir_xy(:,:) = 0.

end

