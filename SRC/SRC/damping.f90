
subroutine damping()

!  "Spange"-layer damping at the domain top region and around poles
!   Use implicit scheme.

use vars
use microphysics, only: micro_field, index_water_vapor
use params, only: dodamping, dolatlon, dodamping_poles, pi, dodamping_w, &
                  damping_u_cu, damping_w_cu, ggr, nub, dohs94, dodamping_u, doterrain, doflat
use terrain, only: terrau,terrav,terraw
use check_energy
implicit none

real nu, nustar, tau_max
real tauy(ny),umax(ny),wmax,tau,cu,wmax1,zzz
integer i, j, k, n_damp(nx,ny)


call t_startf ('damping')

call sumenergy(1,'damping')

tau_max = dtn/dt 

!--------------------------------------------------------------------------------
! based on MetOffice Unified Model 6.0 and on Wood et al doi:10.1002/qj.2235
! damp vertical velocity at the top of domain 
taudamp = 0.

if(dodamping) then

 do k=1,nzm
   nu = (zi(k)-zi(1))/(zi(nzm)-zi(1))
   if(nu.gt.nub) then
 !     taudamp(k) = tau_max*sin(0.5*pi*(nu-nub)/(1.-nub))**2
      zzz = 100.*((nu-nub)/(1.-nub))**2
      taudamp(k) = 0.333*tau_max*zzz/(1.+zzz)
   end if
 end do

  do k = 1,nzm
   if(taudamp(k).gt.0.) then
     do j = 1,ny
      do i=1,nx
        w(i,j,k)= w(i,j,k)/(1.+taudamp(k))*terraw(i,j,k)
      end do! i 
     end do! j 
   end if
  end do ! k

end if ! dodamping

if(dodamping_w) then

  do k = 1,nzm
   if(taudamp(k).eq.0.) then
     wmax1 = damping_w_cu*dz*adzw(k)/dtn
     do j=1,ny
       wmax = wmax1 !*dx*mu(j)/(dy*ady(j))
       do i=1,nx
          w(i,j,k)= (w(i,j,k) + min(wmax,max(-wmax,w(i,j,k)))*tau_max) &
                    /(1.+tau_max) * terraw(i,j,k)
       end do! i 
     end do! j
   end if
  end do ! k

end if

!--------------------------------------------------------------------------------
if(dolatlon.and..not.doflat.and.dodamping_poles) then

  do j=1,ny
    tauy(j) = tau_max*(1.-mu(j)**2)**200
    umax(j) = damping_u_cu*dx*mu(j)/dtn
  end do
  do k = 1,nzm
   if(dodamping_u.and.pres(k).lt.70.) then
    tau = tau_max
    do j=1,ny
      do i=1,nx
           u(i,j,k)= (u(i,j,k) + min(umax(j),max(-umax(j),u(i,j,k)))*tau) &
               /(1. + tau) * terrau(i,j,k)
           v(i,j,k)= (v(i,j,k) + min(umax(j),max(-umax(j),v(i,j,k)))*tau) &
               /(1. + tau) * terrav(i,j,k)
      end do
    end do! j
   else
    do j=1,ny
      do i=1,nx
           u(i,j,k)= (u(i,j,k) + min(umax(j),max(-umax(j),u(i,j,k)))*tauy(j)) &
                     /(1. + tauy(j)) * terrau(i,j,k)
           v(i,j,k)= (v(i,j,k) + min(umax(j),max(-umax(j),v(i,j,k)))*tauy(j)) &
                     /(1. + tauy(j)) * terrav(i,j,k)
      end do
    end do
   end if
  end do ! k

end if 

call sumenergy(-1,'damping')

call t_stopf('damping')




end subroutine damping
