
subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor
use sgs, only: setperturb_sgs
use terrain, only: k_terra, terra, terrau, terrav, terraw
use rad, only: qrad

implicit none

integer i,j,k,ptype,it,jt,n
real(4) rrr,ranf_
real xxx,yyy,zzz
real nu, nuv, phi, lamda, phi0, lamda0, phi1, lamda1
real, allocatable :: psi(:,:)

call ranset_(100+3*(rank+nsubdomains*nens))

ptype = perturb_type

call setperturb_sgs(ptype)  ! set sgs fields

call task_rank_to_index(rank,it,jt)

select case (ptype)

  case(-1)

! do nothing

  case(-11)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.le.10) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

 case(-10)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
           u(i,j,k)=u(i,j,k)+1*rrr*terrau(i,j,k)
       end do
      end do
     end do

 case(-2)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.lt.10) then
           u(i,j,k)=u(i,j,k)+0.1*rrr*terrau(i,j,k)
         end if
       end do
      end do
     end do

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.ge.k_terra(i,j).and.k.le.k_terra(i,j)+5) then
            t(i,j,k)=t(i,j,k)+0.4*rrr*terra(i,j,k)*dx/4000.
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.ge.k_terra(i,j).and.k.le.k_terra(i,j)+5) then
            u(i,j,k)=u(i,j,k)+0.1*rrr*terrau(i,j,k)
         endif
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt-ny_gl/2.-0.5)
         do i=1,nx
          xxx = dx*(i+it-nx_gl/2.)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr*terra(i,j,k)
              if(mod(rank,2).eq.0.and.i.eq.(nx/2).and.k.eq.5.and.(j.eq.1.or.j.eq.ny)) print*,rank,j,t(i,j,k)
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr*terra(i,j,k)
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*terra(i,j,k)
         endif
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)*terra(i,j,k)
         endif
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*terra(i,j,k)
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr*terra(i,j,k)
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*terra(i,j,k)
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr*terra(i,j,k)
         endif
       end do
      end do
     end do

  case(7)  ! solid-body rotation on sphere

     allocate(psi(nx+1,ny+1))
     ! new North Pole coorditates:
     phi0 =   0. ! pi/2.
     lamda0 = 0.
     do j=1,ny+1
       phi1 = latv_gl(j+jt)*pi/180.
       do i=1,nx+1
         if(i+it.ne.nx_gl+1) then
           lamda1 = (lon_gl(i+it)-0.5*(lon_gl(1+COL1)-lon_gl(1)))*pi/180.
         else
           lamda1 = (lon_gl(1)-0.5*(lon_gl(1+COL1)-lon_gl(1)))*pi/180.
         end if
         phi = asin(sin(phi1)*sin(phi0)-cos(phi1)*cos(phi0)*cos(lamda1))
         psi(i,j) = sin(phi1)*sin(phi0)-cos(phi1)*cos(phi0)*cos(lamda1)
     !    psi(i,j) = 0.5*cos(2.*phi)
       end do
     end do
     psi = psi * 184350381.222549876563222 / earth_factor
     do j=1,ny
       do i=1,nx
         u(i,j,:) = (psi(i,j+1)-psi(i,j))/(ady(j)*dy)
         qrad(i,j,:) = psi(i,j)
       end do
     end do
     do j=1,ny+1
       if(j+jt.ne.1.and.j+jt.ne.ny_gl+1) then
        do i=1,nx
          v(i,j,:) = -(psi(i+1,j)-psi(i,j))/muv(j)/dx
        end do
       else
          v(:,j,:) = 0.
       end if
     end do
     deallocate(psi)
     do j=1,ny
       do i=1,nx
         phi1 = lat_gl(j+jt)*pi/180.
         lamda1 = lon_gl(i+it)*pi/180.
         phi = asin(sin(phi1)*sin(phi0)-cos(phi1)*cos(phi0)*cos(lamda1))
         lamda = atan2(cos(phi1)*sin(lamda1),sin(phi)*cos(phi0)+cos(phi1)*cos(lamda1)*sin(phi0))
       !  lamda = asin(cos(phi1)*sin(lamda1)/cos(phi))
    !     micro_field(i,j,:,1) = 1.e-3*(0.5+0.5*cos(2.*phi))
    !     micro_field(i,j,:,1) = sin(phi1)**30
         if(abs(lat_gl(j+jt)).lt.10.) then
             micro_field(i,j,:,1) = 1.
         else
             micro_field(i,j,:,1) = 0.
         end if
       end do
     end do


  case(75)  ! Rossby-Haurwitz wave on shere

     n = 4 ! number of waves
     allocate(psi(nx+1,ny+1))
     do j=1,ny+1
       phi = latv_gl(j+jt)*pi/180.
       do i=1,nx+1
         if(i+it.ne.nx_gl+1) then
           lamda = (lon_gl(i+it)-0.5*(lon_gl(1+COL1)-lon_gl(1)))*pi/180.
         else
           lamda = (lon_gl(1)-0.5*(lon_gl(1+COL1)-lon_gl(1)))*pi/180.
         end if
         psi(i,j) = -sin(phi)+cos(phi)**n*sin(phi)*cos(n*lamda)
       end do
     end do
     psi = -psi * (50./n) *rad_earth/earth_factor
     do j=1,ny
       do i=1,nx
         u(i,j,:) = (psi(i,j+1)-psi(i,j))/(ady(j)*dy)
       end do
     end do
     do j=1,ny+1
       if(j+jt.ne.1.and.j+jt.ne.ny_gl+1) then
        do i=1,nx
          v(i,j,:) = -(psi(i+1,j)-psi(i,j))/muv(j)/dx
        end do
       else
          v(:,j,:) = 0.
       end if
     end do
     deallocate(psi)

  case(76)  ! Non-divergent horizontal flow in a box

     allocate(psi(nx+1,ny+1))
     do j=1,ny+1
       do i=1,nx+1
         psi(i,j) = sin(pi*(i+it-1.)/nx_gl)*sin(pi*(j+jt-1.)/ny_gl)
       end do
     end do
     do j=1,ny
       do i=1,nx
         u(i,j,:) = (psi(i,j+1)-psi(i,j))/(ady(j)*dy)
       end do
     end do
     do j=1,ny+1
       if(j+jt.ne.1.and.j+jt.ne.ny_gl+1) then
        do i=1,nx
          v(i,j,:) = -(psi(i+1,j)-psi(i,j))/muv(j)/dx
        end do
       else
          v(:,j,:) = 0.
       end if
     end do
     deallocate(psi)

  case(77)  ! Non-divergent vertical x-z  flow in a box

     allocate(psi(nx+1,nzm+1))
     do k=1,nzm+1
       do i=1,nx+1
         psi(i,k) = 0.2*sin(pi*(i+it-1.)/nx_gl)*sin(pi*(k-1.)/nzm)
       end do
     end do
     do k=1,nzm
       do i=1,nx
         u(i,:,k) = (psi(i,k+1)-psi(i,k))/(adz(k)*dz)
       end do
     end do
     do k=1,nzm+1
        do i=1,nx
          w(i,:,k) = -(psi(i+1,k)-psi(i,k))/dx
        end do
     end do
     deallocate(psi)



  case(8)  ! Jablonowski & Williamson (2006) hydrostatic test

     do k=1,nzm
      nu = pres(k)/1000.
      nuv = pi/2.*(nu-0.252)
      if(nu.ge.0.2) then
         t(:,:,k) = 288.*nu**(287.*0.005/9.80616)
      else
         t(:,:,k) = 288.*nu**(287.*0.005/9.80616)+4.8e5*(0.2-nu)**5
      end if 
      t(1:nx,1:ny,k) = t(1:nx,1:ny,k)+gamaz(k)
      do j=1,ny
       do i=1,nx
           phi = latitude(i,j)*pi/180.
           t(i,j,k) = t(i,j,k) +0.75*nu*pi*35./287.*sin(nuv)*cos(nuv)**0.5 &
             *((-2.*sin(phi)**6*(cos(phi)**2+1./3.)+10./63.)*2.*35.*cos(nuv)**1.5 &
               +(8./5.*cos(phi)**3*(sin(phi)**2+2./3.)-pi/4.)*rad_earth*7.29212e-5)
           u(i,j,k) = 35.*cos(nuv)**1.5*sin(2.*phi)**2
           lamda = longitude(i,j)*pi/180.
!           u(i,j,k) = u(i,j,k) + &
!          exp(-100.*acos(sin(2.*pi/9.)*sin(phi)+cos(2.*pi/9.)*cos(phi)*cos(lamda-pi/9.))**2)
       end do
      end do
     end do
     micro_field(:,:,:,:) = 0.
     call diagnose()
     bet(:) = ggr/tabs0(:)

     case(81) 
! baroclinic test perturbation:
     do k=1,nzm
      do j=1,ny
       do i=1,nx
         lamda = longitude(i,j)*pi/180.
         phi = latitude(i,j)*pi/180.
         u(i,j,k) = u(i,j,k) + &
          exp(-100.*acos(sin(2.*pi/9.)*sin(phi)+cos(2.*pi/9.)*cos(phi)*cos(lamda-pi/9.))**2)
       end do
      end do
     end do

  case(9)  ! point source

     do j=1,ny
      do i=1,nx
        if(i+it.ge.nx_gl/10-2.and.i+it.le.nx_gl/10+2.and. &
           j+jt.ge.ny_gl/2-ny_gl/4.and.j+jt.le.ny_gl/2+ny_gl/4) &
          micro_field(i,j,1:3,index_water_vapor)=micro_field(i,j,1:3,index_water_vapor)+15.e-3
      end do
     end do
     print*,'>>',rank,minval(micro_field),maxval(micro_field)

  case(10)  ! Held and Suarez (1994) test

     do k=1,nzm
      do j=1,ny
       do i=1,nx
           phi = latitude(i,j)*pi/180.
           t(i,j,k) = max(200.,(315.-60.*sin(phi)**2-10.*log(pres(k)/1000.)*cos(phi)**2) &
                              /prespot(k)) + gamaz(k)
           rrr=1.-2.*ranf_()
           if(k.ge.k_terra(i,j).and.k.le.k_terra(i,j)+10) then
              t(i,j,k)=t(i,j,k)+1.*rrr*terra(i,j,k)
           endif
       end do
      end do
     end do

  case(11)  ! DCMIP 2.1 test

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         phi = latitude(i,j)*pi/180.
         u(i,j,k) = u(i,j,k)*cos(phi)
         t(i,j,k) = 300.+gamaz(k)
       end do
      end do
     end do
     call diagnose()
     bet(:) = ggr/tabs0(:)

  case(12)  ! Sheared flow

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.lt.36) then
           u(i,j,k) = u(i,j,k)+0.1*rrr
         end if
       end do
      end do
     end do

  case(13)  ! Rossby Mountain test by Jablonowski

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         phi = latitude(i,j)*pi/180.
         u(i,j,k) = 20.*cos(phi)
       end do
      end do
     end do
     call diagnose()


  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select


end

