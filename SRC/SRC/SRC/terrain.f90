#include "fppmacros"

module terrain

!incorporation of terrain into SAM
!Author: Marat Khairoutdinov, Nov. 2011
!--------------------------------------------------------------------

use grid
use vars, only: elevation, elevationg, landmask, longitude, latitude
use params, only: terrainfile, terrain_type, doterrain, pi,dowallx,dowally, &
                  rad_earth, earth_factor, readterr, alpha_hybrid, doterrepair, docheck

implicit none

! terrain masks: inside terrain? (1 no, 0 yes, in between - kind of)
real:: terrau   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) = 1.
real:: terrav   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) = 1.
real:: terraw   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) = 1.
real:: terra    (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) = 1.

! weighting coefficient for velocity in hybrid advection scheme
real:: alphah   (0:nx+1, 0:ny+1, nzm) = 1.

integer:: k_terra(dimx1_s:dimx2_s, dimy1_s:dimy2_s) = 1 ! vertical grid index of terrain top
integer:: k_terrau(dimx1_u:dimx2_u, dimy1_u:dimy2_u) = 1
integer:: k_terrav(dimx1_v:dimx2_v, dimy1_v:dimy2_v) = 1

integer:: kmax = 1  ! maximum vertical index of terrain in the current subdomain

real, parameter :: alpha_min = 0. ! minimum alpha around terrain 
                             !(overwritten by alpha_hybrid parameter)

contains

subroutine setterrain

!  Random noise

implicit none

integer i,j,k,m,it,jt,nx1,ny1,kb,kc,iter,nrep,nrep1,ii,jj
real xxx,yyy,zzz,rad, aaa, bbb, ccc, xxx0, yyy0, dxxx, dyyy, bbb1,ccc1
real tmp(nx,ny,2)
real days(2)
real(4), external :: ranf_ 

if(masterproc) then
        print*,'doterrain=',doterrain
        if(doterrain) print*,'terrain_type=',terrain_type
        print*,'doterrepair:',doterrepair
end if

call task_rank_to_index(rank,it,jt)

select case (terrain_type)

case(0)

case(1) ! bell-shaped mountain

  ! Bouldar case:
       aaa = 10000. ! half-width
       bbb = 2000. ! height
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           elevation(i,j) = bbb*aaa**2/(aaa**2+(xxx-0.5*nx_gl*dx)**2)
      !     elevation(i,j) = bbb*aaa**2/(aaa**2+(xxx-0.5*nx_gl*dx)**2+&
      !                YES3D*(yyy**2))
           end do
         end do
       
case(11) ! bell-shaped mountain

  ! Bouldar case:
       aaa = 10000. ! half-width
       bbb = 100. ! height
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           elevation(i,j) = bbb*aaa**2/(aaa**2+(xxx-0.5*nx_gl*dx)**2)
      !     elevation(i,j) = bbb*aaa**2/(aaa**2+(xxx-0.5*nx_gl*dx)**2+&
      !                YES3D*(yyy**2))
           end do
         end do

  
case(2) ! Cube

       aaa=16*dx ! width
       bbb=16*dz ! height
       elevation = 0.
       do j=1,ny
        yyy = dy*(j+jt-1)
         do i=1,nx
          xxx = dx*(i+it-1)
           if(abs(xxx-0.5*(nx_gl+0.5)*dx).le.0.5*aaa.and. &
              abs(yyy-0.5*(ny_gl+0.5)*dx).le.0.5*aaa ) then
                elevation(i,j) = bbb
           end if
         end do
       end do

case(21) ! Building from CEVAL A1-1 and A1-5

       aaa= 20. ! length in x
       bbb= 30. ! lenth in y
       ccc= 25. ! height
       elevation = 0.
       xxx0 = (0.5*nx_gl)*dx
       yyy0 = 0.
       do j=1,ny
         do i=1,nx
          xxx = dx*(i+it)
           if(abs(xxx-xxx0).lt.aaa/2..and.abs(y_gl(j+jt)-yyy0).lt.bbb/2.) then
                elevation(i,j) = ccc
           end if
         end do
       end do

case(210) ! rotate Building from CEVAL A1-1 and A1-5 by 45 degrees counter clockwise

       aaa= 20. ! length in x
       bbb= 30. ! lenth in y   ! reduce size a bit because of sperious  size increase due to rotation
       ccc= 25. ! height
       elevation = 0.
       xxx0 = 0.5*nx_gl*dx
       yyy0 = 0.
       bbb1 = (bbb-dx)/sqrt(2.)
       ccc1 = (aaa-dx)/sqrt(2.)
       do j=1,ny
         do i=1,nx
          xxx = dx*(i+it)-xxx0
          yyy = y_gl(j+jt)-yyy0
           if(yyy.lt.xxx+bbb1.and.yyy.gt.xxx-bbb1 &
              .and.yyy.lt.-xxx+ccc1.and.yyy.gt.-xxx-ccc1) elevation(i,j) = ccc
         end do
       end do

case(22) ! Array of Buildings from CEDVAL B1-1

       aaa=20*dz ! length in x
       bbb=30*dz ! lenth in y
       ccc=25*dz ! height
       xxx0 = 0.5*nx_gl*dx
       yyy0 = 0.
       dxxx = 40*dz
       dyyy = 50*dz
       elevation = 0.
       do jj=-1,1
       do ii=-4,2
       do j=1,ny
        yyy = dy*(j+jt-1)
         do i=1,nx
          xxx = dx*(i+it+0.5)
           if(abs(xxx-(xxx0+ii*dxxx)).lt.aaa/2..and.abs(y_gl(j+jt)-(yyy0+jj*dyyy)).lt.bbb/2.) then
                elevation(i,j) = ccc
           end if
         end do
       end do
       end do
       end do

case (26) ! cube rotated by 45 degrees relative to the flow CEDVAL case A1-6

       aaa= 25. ! length in x
       bbb= 25. ! lenth in y   ! reduce size a bit because of sperious  size increase due to rotation
       ccc= 25. ! height
       elevation = 0.
       xxx0 = 0.5*nx_gl
       yyy0 = 0.
       bbb1 = (bbb-dx)/sqrt(2.)
       ccc1 = (aaa-dx)/sqrt(2.)
       do j=1,ny
         do i=1,nx
          xxx = dx*(i+it)-xxx0
          yyy = y_gl(j+jt)-yyy0
           if(yyy.lt.xxx+bbb1.and.yyy.gt.xxx-bbb1 &
              .and.yyy.lt.-xxx+ccc1.and.yyy.gt.-xxx-ccc1) elevation(i,j) = ccc
         end do
       end do

!       aaa=25*dx ! size
!       xxx0 = 0.5*nx_gl*dx
!       yyy0 = 0.
!       elevation = 0.
!       do j=1,ny
!         do i=1,nx
!          xxx = dx*(i+it-1)
!           if(abs(xxx-xxx0).lt.aaa/sqrt(2.).and.abs(y_gl(j+jt)).lt.aaa/sqrt(2.)-abs(xxx-xxx0)) &
!                elevation(i,j) = aaa
!         end do
!       end do
!

case (27) ! vertical cyclinder

       aaa=20.*dx ! radius
       bbb = 40.*dx ! height
       xxx0 = 0.5*nx_gl*dx ! position of center
       yyy0 = 0.
       elevation = 0.
       do j=1,ny
         do i=1,nx
          xxx = dx*(i+it)
           if(sqrt((xxx-xxx0)**2+y_gl(j+jt)**2).lt.aaa) &
                elevation(i,j) = bbb
         end do
       end do


case(28) ! Array of 8 bars in pi-Chamber

       aaa=40./2000.*nx_gl*dx ! width of bar
       ccc=12./1000.*zi(nz) ! height
       dxxx =  nx_gl*dx/(8+1) ! distance between bars
       elevation = 0.
       do ii=1,8
         do i=1,nx
          xxx = dx*(i+it)
           if(abs(xxx-(ii*dxxx)).lt.aaa/2.) then
                elevation(i,:) = ccc
           end if
         end do
       end do


case(3) ! kink
  
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           elevation(i,j)=max(0.,0.95*zi(2)*(1.-abs(xxx-0.5*nx_gl*dx)/(30.*dx)))
         end do
       end do
 

case(4) ! squared-cosine hill
 
       aaa = nx_gl/4*dx ! width
       bbb = 3. ! height
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
          zzz = min(0.25,((xxx-(0.5*nx_gl+0.5)*dx)/aaa)**2 &
                    +((yyy-(0.5*ny_gl+0.5)*dy)/aaa)**2 &
                    )
          elevation(i,j) = bbb*cos(pi*sqrt(zzz))**2
         end do
       end do
  
case(5) ! squared-cosine valley

       aaa = nx_gl/8.*dx ! half-width
       bbb = 10. ! depth
       do j=1,ny
         yyy = y_gl(ny_gl)
         do i=1,nx
          xxx = dx*(i+it-1)
          zzz = ((xxx-0.5*nx_gl*dx)/aaa)**2+(YES3D*(yyy-0.5*yv_gl(ny_gl+1))/aaa)**2
          elevation(i,j) = max(0.,bbb-bbb*cos(real(pi)*sqrt(zzz))**2)
        end do
       end do

     case(6) ! round island
       
       elevation = 0.
       where(landmask.eq.1) elevation = z(3)

case(7) ! random terrain

       elevation = 0.
       call ranset_(5*rank)
       do j=1,ny
        do i=1,nx
         k=ranf_()*4
         elevation(i,j) = zi(k)
        end do
       end do

case(8) ! Schar test (Cartesian)

       aaa = 5000. ! half-width
       bbb = 250. ! height
       ccc = 4000. ! wavelength
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           rad = sqrt((xxx-0.5*nx_gl*dx)**2+(yyy-0.5*yv_gl(ny_gl+1))**2)
           elevation(i,j) = bbb*exp(-rad**2/aaa**2)*cos(pi*rad/ccc)**2
         end do
       end do

case(81) ! Schar ridge (2D; Cartesian)

       aaa = 5000. ! half-width
       bbb = 250. ! height
       ccc = 4000. ! wavelength
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           rad = xxx-0.25*nx_gl*dx
           elevation(i,j) = bbb*exp(-rad**2/aaa**2)*cos(pi*rad/ccc)**2
         end do
       end do

case(82) ! Schar circular ridge (2D; Cartesian)

       aaa = 5000. ! half-width
       bbb = 250. ! height
       ccc = 4000. ! wavelength
       do j=1,ny
         yyy = y_gl(j+jt)
         do i=1,nx
          xxx = dx*(i+it-1)
           rad = sqrt((xxx-0.25*nx_gl*dx)**2+yyy**2)
           elevation(i,j) = bbb*exp(-rad**2/aaa**2)*cos(pi*rad/ccc)**2
         end do
       end do

case(9) ! DCMIP test 2.1

       aaa = 5000. ! half-width
       bbb = 250. ! height
       ccc = 4000. ! wavelength
       do j=1,ny
         do i=1,nx
          rad = rad_earth/earth_factor*acos(cos(pi/180.*latitude(i,j)) &
                *cos(pi/180.*longitude(i,j)-0.25*pi))
          elevation(i,j) = bbb*exp(-rad**2/aaa**2)*cos(pi*rad/ccc)**2
         end do
       end do

case(10) ! Rossby Mountain test by Jablonowski

       aaa = 1500000. ! half-width
       bbb = 2000. ! height
       do j=1,ny
         do i=1,nx
          rad = rad_earth/earth_factor*acos(sin(pi/180.*latitude(i,j))*sin(pi/6.)+ &
                cos(pi/180.*latitude(i,j))*cos(pi/6.)*cos(pi/180.*longitude(i,j)-0.5*pi))
          elevation(i,j) = bbb*exp(-rad**2/aaa**2)
         end do
       end do

case default

       if(masterproc) print*,'terrain_type ',terrain_type, &
                     ' is not defined in setperturb(). Exitting...'
       call task_abort()

end select

if(readterr) then
    call readsurface(terrainfile,tmp(:,:,1:2),days)
    elevation(:,:) = tmp(:,:,1)
    where(elevation.lt.0.) elevation=0.
end if

if(docheck) return

do k=1,nzm
  do j=1,ny
    do i=1,nx
        if(elevation(i,j).ge.z(k)) then 
         terra(i,j,k) = 0.
         kmax = max(kmax,k) 
        end if
    end do
  end do
end do

if(dompi) then
 call task_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
else
 call bound_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
endif


if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,0
         terra(i,j,k) = terra(1,j,k)
       end do
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=nx+1,dimx2_s
         terra(i,j,k) = terra(nx,j,k)
       end do
     end do
    end do
  end if
end if

if(dowally) then
  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
     do j=dimy1_s,0
       do i=dimx1_s,dimx2_s
         terra(i,j,k) = terra(i,1,k)
       end do
     end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
     do j=ny+1,dimy2_s
       do i=dimx1_s,dimx2_s
         terra(i,j,k) = terra(i,ny,k)
       end do
     end do
    end do
  end if
end if

do k=1,nzm
  do j=dimy1_s,dimy2_s
    do i=dimx1_s,dimx2_s
       if(terra(i,j,k).lt.1.) k_terra(i,j) = k+1
    end do
  end do
end do

if(doterrepair) then
! repair terrain by removing pits/holes with 1 grid-cell width
! MK March 2020
! remove pits
iter = 0
nrep1 = 0
do while(iter.lt.20)
  iter = iter+1
  nrep = 0
  ! first repair in x direction:
  do j=1,ny
      do i=1,nx
        if(k_terra(i,j).lt.k_terra(i-1,j).and.k_terra(i,j).lt.k_terra(i+1,j)) then 
         nrep = nrep + 1
         if(k_terra(i-1,j).lt.k_terra(i+1,j)) then
            terra(i,j,:) = terra(i-1,j,:)  
            k_terra(i,j) = k_terra(i-1,j)  
         else
            terra(i,j,:) = terra(i+1,j,:)  
            k_terra(i,j) = k_terra(i+1,j)  
         end if
        end if
      end do
    end do
    do j=1,ny
      do i=1,nx
        if(k_terra(i,j).lt.k_terra(i,j-1).and.k_terra(i,j).lt.k_terra(i,j+1)) then
         nrep = nrep + 1
         if(k_terra(i,j-1).lt.k_terra(i,j+1)) then
            terra(i,j,:) = terra(i,j-1,:)
            k_terra(i,j) = k_terra(i,j-1)
         else
            terra(i,j,:) = terra(i,j+1,:)
            k_terra(i,j) = k_terra(i,j+1)
         end if
        end if
      end do
    end do
  call task_barrier()
  if(dompi) then
   call task_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  else
   call bound_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  endif
  if(dowallx) then
   if(mod(rank,nsubdomains_x).eq.0) then
     terra(dimx1_s:0,:,:) = terra(-dimx1_s+1:1:-1,:,:)
   end if
   if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
     terra(nx+1:dimx2_s,:,:) =  terra(nx:2*nx-dimx2_s+1:-1,:,:)
   end if
  end if
  if(RUN3D.and.dowally) then
   if(rank.lt.nsubdomains_x) then
     terra(:,dimy1_s:0,:) = terra(:,-dimy1_s+1:1:-1,:)
   end if
   if(rank.gt.nsubdomains-nsubdomains_x-1) then
     terra(:,ny+1:dimy2_s,:) = terra(:,ny:2*ny-dimy2_s+1:-1,:)
   end if
  end if
  if(dowallx) then
    if(mod(rank,nsubdomains_x).eq.0) then
      do k=1,nzm
       do j=dimy1_s,dimy2_s
         do i=dimx1_s,0
           terra(i,j,k) = terra(1,j,k)
         end do
       end do
      end do
    end if
    if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
      do k=1,nzm
       do j=dimy1_s,dimy2_s
         do i=nx+1,dimx2_s
           terra(i,j,k) = terra(nx,j,k)
         end do
       end do
      end do
    end if
  end if
  if(dowally) then
    if(rank.lt.nsubdomains_x) then
      do k=1,nzm
       do j=dimy1_s,0
         do i=dimx1_s,dimx2_s
           terra(i,j,k) = terra(i,1,k)
         end do
       end do
      end do
    end if
    if(rank.gt.nsubdomains-nsubdomains_x-1) then
      do k=1,nzm
       do j=ny+1,dimy2_s
         do i=dimx1_s,dimx2_s
           terra(i,j,k) = terra(i,ny,k)
         end do
       end do
      end do
    end if
  end if
  do k=1,nzm
    do j=dimy1_s,dimy2_s
      do i=dimx1_s,dimx2_s
         if(terra(i,j,k).lt.1.) k_terra(i,j) = k+1
      end do
    end do
  end do
  nrep1 = nrep1+nrep
end do
if(nrep.ne.0) then
  print*,'terrain pits could not be repaired. rank=',rank,' nrep =',nrep
  stop
end if
if(dompi) then
  call task_sum_integer1(nrep1,nrep)
else
  nrep = nrep1
end if
if(masterproc) print*,'terrain pits repaired. Number of repaires:',nrep
end if

if(doterrepair) then
! repair terrain by removing peaks with 1 grid-cell width
iter = 0
nrep1 = 0
do while(iter.lt.20)
  iter = iter+1
  nrep = 0
  ! first repair in x direction:
  do j=1,ny
      do i=1,nx
        if(k_terra(i,j).gt.k_terra(i-1,j).and.k_terra(i,j).gt.k_terra(i+1,j)) then 
         nrep = nrep + 1
         if(k_terra(i-1,j).gt.k_terra(i+1,j)) then
            terra(i,j,:) = terra(i-1,j,:)  
            k_terra(i,j) = k_terra(i-1,j)  
         else
            terra(i,j,:) = terra(i+1,j,:)  
            k_terra(i,j) = k_terra(i+1,j)  
         end if
        end if
      end do
    end do
    do j=1,ny
      do i=1,nx
        if(k_terra(i,j).gt.k_terra(i,j-1).and.k_terra(i,j).gt.k_terra(i,j+1)) then
         nrep = nrep + 1
         if(k_terra(i,j-1).gt.k_terra(i,j+1)) then
            terra(i,j,:) = terra(i,j-1,:)
            k_terra(i,j) = k_terra(i,j-1)
         else
            terra(i,j,:) = terra(i,j+1,:)
            k_terra(i,j) = k_terra(i,j+1)
         end if
        end if
      end do
    end do
  call task_barrier()
  if(dompi) then
   call task_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  else
   call bound_exchange(terra,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  endif
  if(dowallx) then
   if(mod(rank,nsubdomains_x).eq.0) then
     terra(dimx1_s:0,:,:) = terra(-dimx1_s+1:1:-1,:,:)
   end if
   if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
     terra(nx+1:dimx2_s,:,:) =  terra(nx:2*nx-dimx2_s+1:-1,:,:)
   end if
  end if
  if(RUN3D.and.dowally) then
   if(rank.lt.nsubdomains_x) then
     terra(:,dimy1_s:0,:) = terra(:,-dimy1_s+1:1:-1,:)
   end if
   if(rank.gt.nsubdomains-nsubdomains_x-1) then
     terra(:,ny+1:dimy2_s,:) = terra(:,ny:2*ny-dimy2_s+1:-1,:)
   end if
  end if
  if(dowallx) then
    if(mod(rank,nsubdomains_x).eq.0) then
      do k=1,nzm
       do j=dimy1_s,dimy2_s
         do i=dimx1_s,0
           terra(i,j,k) = terra(1,j,k)
         end do
       end do
      end do
    end if
    if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
      do k=1,nzm
       do j=dimy1_s,dimy2_s
         do i=nx+1,dimx2_s
           terra(i,j,k) = terra(nx,j,k)
         end do
       end do
      end do
    end if
  end if
  if(dowally) then
    if(rank.lt.nsubdomains_x) then
      do k=1,nzm
       do j=dimy1_s,0
         do i=dimx1_s,dimx2_s
           terra(i,j,k) = terra(i,1,k)
         end do
       end do
      end do
    end if
    if(rank.gt.nsubdomains-nsubdomains_x-1) then
      do k=1,nzm
       do j=ny+1,dimy2_s
         do i=dimx1_s,dimx2_s
           terra(i,j,k) = terra(i,ny,k)
         end do
       end do
      end do
    end if
  end if
  do k=1,nzm
    do j=dimy1_s,dimy2_s
      do i=dimx1_s,dimx2_s
         if(terra(i,j,k).lt.1.) k_terra(i,j) = k+1
      end do
    end do
  end do
  nrep1 = nrep1+nrep
end do
if(nrep.ne.0) then
  print*,'terrain pits could not be repaired. rank=',rank,' nrep =',nrep
  stop
end if
if(dompi) then
  call task_sum_integer1(nrep1,nrep)
else
  nrep = nrep1
end if
if(masterproc) print*,'terrain peaks repaired. Number of repaires:',nrep
end if


do k=1,nzm
  do j=1,ny
    do i=1,nx
        if(terra(i,j,k).lt.1.) then 
         elevationg(i,j) = zi(k+1)
        end if
    end do
  end do
end do

do k=1,nzm
  do j=1,ny
    do i=1,nx
        if(terra(i,j,k).lt.1) then
         terraw(i,j,k) = 0
         terraw(i,j,k+1) = 0
        end if
    end do
  end do
end do


do k=1,nzm
  do j=1,ny
    do i=1,nx
       terrau(i,j,k) = min(terra(i-1,j,k),terra(i,j,k))
       terrav(i,j,k) = min(terra(i,j-YES3D,k),terra(i,j,k))
    end do
  end do
end do

call task_barrier()

if(dompi) then
 call task_exchange(terrau,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV)
 call task_exchange(terrav,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3)
 call task_exchange(terraw,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV)
else
 call bound_exchange(terrau,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV)
 call bound_exchange(terrav,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3)
 call bound_exchange(terraw,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV)
endif

call task_barrier()

call boundaries(1)
call boundaries(2)

do k=1,nzm,1
  do j=dimy1_u,dimy2_u
    do i=dimx1_u,dimx2_u
       if(terrau(i,j,k).lt.1.) k_terrau(i,j) = k+1
    end do
  end do
end do
do k=1,nzm,1
  do j=dimy1_v,dimy2_v
    do i=dimx1_v,dimx2_v
       if(terrav(i,j,k).lt.1.) k_terrav(i,j) = k+1
    end do
  end do
end do

call fminmax_print('elevation (grid):',elevationg(:,:),1,nx,1,ny,1)

call make_alpha_hybrid()

end subroutine setterrain


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine make_alpha_hybrid()

! calculate the hybrid 2nd/3rd order advection scheme's alpha coefficeints
! around terrain to reduce noise
use vars, only: t
use params, only: doterrhyb

real(8) a,b,c,d,e
integer i,j,k,kb,kc,iter
integer, parameter :: niter = 200 ! determines extend of penetration outside the terrain
real alpha(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real tmp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

if(masterproc) then
 print*,'doterrhyb:',doterrhyb
 print*,'alpha_hybrid = ', alpha_hybrid
end if

if(doterrhyb) then
 where(terra.eq.0.) 
  tmp = min(alpha_hybrid,alpha_min)
 elsewhere
  tmp = alpha_hybrid
 end where
! also add some damping at the loweat model levels:
tmp(:,:,1) = min(alpha_hybrid,alpha_min) 


 do iter=1,niter
  do k=2,nzm-1
   kc = min(k+1,nzm)
   kb = max(1,k-1)
   do j=1,ny
    a = mu(j)*muv(j+1)*(dx/dy)**2/(ady(j)*adyv(j+1))
    b = mu(j)*muv(j)*(dx/dy)**2/(ady(j)*adyv(j))
    c = mu(j)*mu(j)/(adz(k)*adzw(k+1))
    d = mu(j)*mu(j)/(adz(k)*adzw(k))
 !   a = mu(j)*muv(j+1)*(dx/dy)**2/(ady(j)*adyv(j+1))
 !   b = mu(j)*muv(j)*(dx/dy)**2/(ady(j)*adyv(j))
 !   c = mu(j)*mu(j)*(dx/dz)**2/(adz(k)*adzw(k+1))
 !   d = mu(j)*mu(j)*(dx/dz)**2/(adz(k)*adzw(k))
    e = 1._8/(2._8+a+b+c+d)
    do i=1,nx
      if(terra(i,j,k).eq.1.) tmp(i,j,k) = &
                 (tmp(i+1,j,k)+tmp(i-1,j,k)+a*tmp(i,j+1,k)+ &
                 b*tmp(i,j-1,k)+c*tmp(i,j,kc)+d*tmp(i,j,kb))*e
    end do
   end do
  end do
  if(dompi) then
    call task_exchange(tmp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  else
    call bound_exchange(tmp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS)
  endif
 end do
 
 alphah(0:nx+1,0:ny+1,1:nzm) = tmp(0:nx+1,0:ny+1,1:nzm)
end if

call fminmax_print('alphah:',alphah(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)

if(dowally) then
  if(rank.lt.nsubdomains_x) then
    alphah(:,0:3,:) = 1.
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    alphah(:,ny-2:ny+1,:) = 1.
  end if
end if


end subroutine make_alpha_hybrid

end module terrain

