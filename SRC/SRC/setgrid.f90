subroutine setgrid

! Initialize vertical grid

use vars	
use params
use terrain, only: setterrain, terra, terrau, terrav, terraw

implicit none
	
real(8) latit, long
integer i,j,it,jt,k, kmax
real(8) sums(nzm),sumu(nzm),sumv(nzm),sumw(nzm),sumss(nzm)
real(8) coef, buffer(nzm,5), buffer1(nzm,5)
real(8) lat(dimy1_s:ny_gl+dimy2_s-ny) 
real(8) adyy(dimy1_s:ny_gl+dimy2_s-ny) 
real field(nx,ny,2)
real days(2)


if(dooceanonly) goto 555

 if(doconstdz) then

    zi(1) = 0.
    do k=2,nz
     zi(k)=dz*(k-1)
    end do
    do k=1,nzm
      adz(k) = (zi(k+1)-zi(k))/dz
      z(k) = 0.5*(zi(k+1)+zi(k))
    end do
    do k=2,nzm
      adzw(k) = (z(k)-z(k-1))/dz
    end do
    adzw(1) = 1.
    adzw(nz) = adzw(nzm)

 else

    open(8,file='./CASES/'//trim(case)//'/grd',status='old',form='formatted') 
    read(8,*) z(1)

    if(z(1).ne.0.) then 

! old (original) style of specifying grid levels in grd as height of scalar levels:

      do k=2,nzm     
        read(8,fmt=*,end=111) z(k)
        kmax=k
      end do
      goto 222
111   do k=kmax+1,nzm
       z(k)=z(k-1)+(z(k-1)-z(k-2))
      end do
222   continue
      close (8)
      dz = 0.5*(z(1)+z(2))
      do k=2,nzm
        adzw(k) = (z(k)-z(k-1))/dz
      end do
      adzw(1) = 1.
      adzw(nz) = adzw(nzm)
      adz(1) = 1.
      do k=2,nzm-1
        adz(k) = 0.5*(z(k+1)-z(k-1))/dz
      end do
      adz(nzm) = adzw(nzm)
      zi(1) = 0.
      do k=2,nz
        zi(k) = zi(k-1) + adz(k-1)*dz
      end do
      dz = zi(2)-zi(1)
      do k=1,nzm
        adz(k) = (zi(k+1)-zi(k))/dz
        z(k) = 0.5*(zi(k+1)+zi(k))
      end do
      do k=2,nzm
        adzw(k) = (z(k)-z(k-1))/dz
      end do
      adzw(1) = 1.
      adzw(nz) = adzw(nzm)

    else 

! new style with grd specifying height of grid interfaces:

      zi(1) = z(1)
      do k=2,nz
        read(8,fmt=*,end=333) zi(k)
        kmax=k
      end do
      goto 444
333   do k=kmax+1,nz
       zi(k)=zi(k-1)+(zi(k-1)-zi(k-2))
      end do
444   continue
      close (8)
      dz = zi(2)-zi(1)
      do k=1,nzm
        adz(k) = (zi(k+1)-zi(k))/dz
        z(k) = 0.5*(zi(k+1)+zi(k))
      end do
      do k=2,nzm
        adzw(k) = (z(k)-z(k-1))/dz
      end do
      adzw(1) = 1.
      adzw(nz) = adzw(nzm)

    end if
 end if

if(masterproc) then
   print*,'nzm=',nzm
   print*,'z=',z
   print*,'zi=',zi
   print*,'nx_gl=',nx_gl
   print*,'ny_gl=',ny_gl
   print*,'nx=',nx
   print*,'ny=',ny
end if

do k=1,nzm
  gamaz(k)=ggr/cp*z(k)
end do


if(z(nzm).gt.40000.) then
  if(masterproc) print*,'the domain top is too high: z(nzm) = ',z(nzm), '   exit...'
  call task_abort()
end if

555 continue

call task_rank_to_index(rank,it,jt)

if(dolatlon) then
  if(readlatlon.or.readlat) then
    if(doflat) then
      if(masterproc) print*,'doflat cannot be used with readlat or readlatlon set to True'
      call task_abort()
    end if
    if(masterproc) print*,'reading lat/lon coords from file',trim(latlonfile)
    if(readlatlon) then
       if(latlonfilebin) then
         open(1, file=trim(latlonfile), form='unformatted', status='old')
         read(1) lon_gl(1:nx_gl)
         read(1) lat(1:ny_gl)
         close(1)
       else
         open(1, file=trim(latlonfile), form='formatted', status='old')
         read(1,*) lon_gl(1:nx_gl)
         read(1,*) lat(1:ny_gl)
         close(1)
       end if
       if(.not.doglobal) then
         latitude0=0.5*(lat(1)+lat(ny_gl))
         longitude0=lon_gl(1)
       ! note that for dolatlon=.true., dx is always horizontal spacing at the equator, 
       ! even when domain is away from the equator.
         dx = (lon_gl(1+COL1)-lon_gl(1))*deg2rad*rad_earth/earth_factor
       end if
    else
       if(latlonfilebin) then
         open(1, file=trim(latlonfile), form='unformatted', status='old')
         read(1) lat(1:ny_gl)
         close(1)
       else
         open(1, file=trim(latlonfile), form='formatted', status='old')
         do j=1,ny_gl
           read(1,*) lat(j)
         end do
         close(1)
       end if
       if(.not.doglobal) latitude0=0.5*(lat(1)+lat(ny_gl))
    end if
    doyvar = .true.
    if(dy.eq.0.) dy = minval(lat(2:ny_gl)-lat(1:ny_gl-1))*deg2rad*rad_earth/earth_factor
  else 
     do j=1,ny_gl+YES3D
      yv_gl(j) = (j-ny_gl/(1.+YES3D)-1.)*dy  
     end do
     y_gl(1:ny_gl) = 0.5*(yv_gl(1+YES3D:ny_gl+YES3D)+yv_gl(1:ny_gl))
     do j=1,ny_gl
       lat(j)=latitude0+y_gl(j)/rad_earth*earth_factor*rad2deg      
       if(lat(j).gt.90..or.lat(j).lt.-90.) then
        if(masterproc) print*,'setgrid: dy on lat/lon grid is too large. exit...', &
            ' j=',j,' lat=',lat(j)
        call task_abort()
       end if
     end do
  end if
  do j=dimy1_s,0
    lat(j) = lat(1)
  end do
  do j=ny_gl+1,ny_gl+dimy2_s-ny
    lat(j) = lat(ny_gl)
  end do
  lat_gl(1:ny_gl) = lat(1:ny_gl)
  do j=1,ny_gl+1
     if(j.eq.1) then
      if(doglobal) then
       latv_gl(j) = -90.
      else
       latv_gl(j) = lat(1)-0.5*(lat(1+YES3D)-lat(1))
      end if
     else if(j.eq.ny_gl+1) then
      if(doglobal) then
       latv_gl(j) = 90.
      else
       latv_gl(ny_gl+1) = lat(ny_gl)+0.5*(lat(ny_gl)-lat(ny_gl-YES3D))
      end if
     else
       latv_gl(j) = 0.5*(lat(j)+lat(j-1))
     end if
     if(latv_gl(j).gt.90..or.latv_gl(j).lt.-90.) then
       if(masterproc) print*,'setgrid: latv_gl on lat/lon grid is > 90. exit...', &
           ' j=',j,' latv_gl=',latv_gl(j)
       call task_abort()
     end if
  end do
  if(readlatlon.or.readlat) yv_gl(1:ny_gl+1) = latv_gl(1:ny_gl+1)*deg2rad*rad_earth/earth_factor
  y_gl(1:ny_gl) = 0.5*(yv_gl(2:ny_gl+1)+yv_gl(1:ny_gl))
  mu_gl(1:ny_gl) = cos(y_gl(1:ny_gl)*(1.d0/dble(rad_earth))*dble(earth_factor))
  do j=dimy1_s,dimy2_s
    if(j+jt.ge.1.and.j+jt.le.ny_gl) then
      mu(j) = mu_gl(j+jt)
    else
      mu(j) = cos(deg2rad*lat(j+jt))
    end if
    if(lat(j+jt).gt.0.) then
      tanr(j) = sqrt(1.-mu(j)**2)/mu(j)/(rad_earth/earth_factor)
    else
      tanr(j) = -sqrt(1.-mu(j)**2)/mu(j)/(rad_earth/earth_factor)
    end if
  !  tanr(j) = sin(deg2rad*lat(j+jt))/cos(deg2rad*lat(j+jt))/(rad_earth/earth_factor)
  end do
  dy = y_gl(ny_gl/(1+YES3D)+YES3D)-y_gl(ny_gl/(1+YES3D))
  do j=1,ny_gl
   adyy(j) = (yv_gl(j+1)-yv_gl(j))/dy
  end do
  do j=dimy1_s,0
    adyy(j) = adyy(1)
  end do
  do j=ny_gl+1,ny_gl+dimy2_s-ny
    adyy(j) = adyy(ny_gl)
  end do
  do j=dimy1_s,dimy2_s
    ady(j) = adyy(j+jt)
  end do
  do j=dimy1_v,dimy2_v
    adyv(j) = 0.5*(ady(j)+ady(j-1))
  end do
  do j=dimy1_v,dimy2_v
     muv(j) = (ady(j-1)*mu(j)+ady(j)*mu(j-1))/(ady(j-1)+ady(j))
  end do
  imu = 1./mu
  imuv = 1./muv
  do j=2,ny_gl
     muv_gl(j) = (adyy(j-1)*mu_gl(j)+adyy(j)*mu_gl(j-1))/(adyy(j-1)+adyy(j))
  end do
  muv_gl(1) = cos(deg2rad*latv_gl(1))
  muv_gl(ny_gl+1) = cos(deg2rad*latv_gl(ny_gl+1))
  if(doflat) then
   mu = 1.
   muv = 1.
   imu = 1.
   imuv = 1.
   mu_gl = 1.
   muv_gl = 1.
   tanr = 0.
  end if
  if(masterproc) then
     print*,'readlatlon=',readlatlon
!     print*,'lat(dimy1_s:ny_gl+dimy2_s-ny):',lat
!     print*,'lat_gl(1:ny_gl):',lat_gl
!     print*,'latv_gl(1:ny_gl+1):',latv_gl
  end if

else  ! not dolatlon

  do j=1,ny_gl
   y_gl(j) = (j-1)*dy
  end do
  yv_gl(1:ny_gl) = y_gl(1:ny_gl) - 0.5*dy
  yv_gl(ny_gl+1) = yv_gl(ny_gl)  + dy
  do j=1,ny_gl
    lat(j)=latitude0+(y_gl(j)-yv_gl(ny_gl/2+1))/rad_earth*earth_factor*rad2deg
    if(lat(j).gt.90..or.lat(j).lt.-90.) then
     if(masterproc) print*,'setgrid: dy on lat/lon grid is too large. exit...', &
         ' j=',j,' lat=',lat(j)
     call task_abort()
    end if
  end do
  if(readlatlon) then
   if(latlonfilebin) then
    open(1, file=trim(latlonfile), form='unformatted', status='old')
    read(1) lon_gl(1:nx_gl)
    read(1) lat(1:ny_gl)
    close(1)
   else
    open(1, file=trim(latlonfile), form='formatted', status='old')
    read(1,*) lon_gl(1:nx_gl)
    read(1,*) lat(1:ny_gl)
    close(1)
   end if
   longitude0 = lon_gl(nx_gl/(1+COL1))
   latitude0 = lat(ny_gl/(1+YES3D))
  end if
  lat_gl(1:ny_gl) = lat(1:ny_gl)
  do j=1,ny_gl+1
     if(j.eq.1) then
       latv_gl(j) = lat(1)-0.5*(lat_gl(1+YES3D)-lat_gl(1))
     else if(j.eq.ny_gl+1) then
       latv_gl(ny_gl+1) = lat_gl(ny_gl)+0.5*(lat_gl(ny_gl)-lat_gl(ny_gl-YES3D))
     else
       latv_gl(j) = 0.5*(lat_gl(j)+lat_gl(j-1))
     end if
     if(latv_gl(j).gt.90..or.latv_gl(j).lt.-90.) then
       if(masterproc) print*,'setgrid: latv_gl on Cartesian grid is > 90. exit...', &
           ' j=',j,' latv_gl=',latv_gl(j)
       call task_abort()
     end if
  end do

end if ! dolatlon


if(.not.readlatlon) then
  do i=1,nx_gl
    lon_gl(i) = longitude0+dx/cos(latitude0*deg2rad)*(i-1)/rad_earth*earth_factor*rad2deg
  end do
end if
! longitudes for zonal velocity 
lonu_gl(1:nx_gl) = lon_gl(1:nx_gl)-0.5*(lon_gl(1+COL1)-lon_gl(1))
lonu_gl(nx_gl+1) = lon_gl(nx_gl)+0.5*(lon_gl(1+COL1)-lon_gl(1))
! 0,0 coordinate in Cartesian grid starts at the u-wind component position, not gridcell center
do i=1,nx_gl
 xu_gl(i) = (i-1)*dx
end do
xu_gl(nx_gl+1) = xu_gl(nx_gl)+dx
x_gl(1:nx_gl) = xu_gl(1:nx_gl)+0.5*dx


if (doradlat) then
  do j=0,ny+1
   latitude(:,j) = lat(j+jt)
  end do
else
  latitude(:,:) = latitude0
end if

if (doradlon) then
  do i=1,nx
    longitude(i,:) = lon_gl(i+it)
  end do
else
  longitude(:,:) = longitude0
end if
call fminmax_print('latitude:',real(latitude(1:nx,1:ny)),1,nx,1,ny,1)
call fminmax_print('longitude:',real(longitude(1:nx,1:ny)),1,nx,1,ny,1)

if(masterproc) then
  print*,'dy=',dy
  print*,'dx=',dx
  print*,'lat_gl (min/max):',minval(lat_gl),maxval(lat_gl)
  print*,'latv_gl( min/max)::',minval(latv_gl),maxval(latv_gl)
  print*,'lon_gl (min/max):',minval(lon_gl),maxval(lon_gl)
  print*,'lonu_gl (min/max):',minval(lonu_gl),maxval(lonu_gl)
  print*,'mu_gl (min/max):',minval(mu_gl(1:ny_gl)),maxval(mu_gl(1:ny_gl))
!  print*,'mu_gl:',mu_gl(1:ny_gl)
  print*,'muv_gl (min/max):',minval(muv_gl(1:ny_gl+1)),maxval(muv_gl(1:ny_gl+1))
!  print*,'muv_gl:',muv_gl(1:ny_gl+1)
  print*,'x_gl (min/max):',minval(x_gl(1:nx_gl)),maxval(x_gl(1:nx_gl))
  print*,'xu_gl (min/max):',minval(xu_gl(1:nx_gl+1)),maxval(xu_gl(1:nx_gl+1))
  print*,'y_gl (min/max):',minval(y_gl(1:ny_gl)),maxval(y_gl(1:ny_gl))
  print*,'yv_gl (min/max):',minval(yv_gl(1:ny_gl+1)),maxval(yv_gl(1:ny_gl+1))
  print*,'ady:',minval((yv_gl(2:ny_gl+1)-yv_gl(1:ny_gl))/dy), &
                maxval((yv_gl(2:ny_gl+1)-yv_gl(1:ny_gl))/dy)
end if
!====================================================
! Initialize landmask:

if(OCEAN) then
  landmask(1:nx,1:ny)=0.
else if(LAND) then
  landmask(1:nx,1:ny)=1.
else if(ISLAND) then
  if(.not.readlandmask) then
      if(doroundisland) then
        do j=1,ny
          do i=1,nx
            if((dx*(i+it)-dx*nx_gl/2)**2+YES3D*(dy*(j+jt)-dy*ny_gl/2)**2.lt.island_radius**2) then
               landmask(i,j) = 1.
            else
               landmask(i,j) = 0.
            end if
          end do
        end do
      else
        do j=1,ny
          do i=1,nx
             if((RUN2D.or.j+jt.ge.island_y1*nx_gl.and.j+jt.le.island_y2*nx_gl).and.     &
                i+it.ge.island_x1*nx_gl.and.i+it.le.island_x2*nx_gl) then
                landmask(i,j) = 1.
             else
                landmask(i,j) = 0.
             end if
         end do
        end do
       end if
  else
    call readsurface(landmaskfile,field,days)
    landmask(:,:) = nint(field(:,:,1))
  end if
end if

!=========================================================================
! Initialize terrain

if(doterrain) call setterrain()

!=========================================================================
! compute weights for averaging:
sums = 0.
sumu = 0.
sumv = 0.
sumw = 0.
sumss = 0.
do k=1,nzm
 do j=1,ny
  do i=1,nx
    sums(k) = sums(k) + mu(j)*ady(j)*terra(i,j,k)
    sumu(k) = sumu(k) + mu(j)*ady(j)*terrau(i,j,k)
    sumv(k) = sumv(k) + muv(j)*adyv(j)*terrav(i,j,k)
    sumw(k) = sumw(k) + mu(j)*ady(j)*terraw(i,j,k)
    sumss(k) = sumss(k) + mu(j)*ady(j)
  end do
 end do
end do
if(dompi) then
  buffer(:,1) = sums(:)
  buffer(:,2) = sumu(:)
  buffer(:,3) = sumv(:)
  buffer(:,4) = sumw(:)
  buffer(:,5) = sumss(:)
  call task_sum_real8(buffer,buffer1,nzm*5)
  sums(:) = buffer1(:,1)
  sumu(:) = buffer1(:,2)
  sumv(:) = buffer1(:,3)
  sumw(:) = buffer1(:,4)
  sumss(:) = buffer1(:,5)
end if
if(masterproc.and.(minval(sums).eq.0..or.minval(sumu).eq.0. &
                .or.minval(sumv).eq.0..or.minval(sumw).eq.0.)) then
    print*,'error in computing wgts in setgrid: terrain cannot be everywhere!'
    call task_abort()
end if
if(docheck) then
  coef = dble(nx*ny)
else
  coef = dble(nx_gl*ny_gl)
end if
do k=1,nzm
 do j=1,ny
   wgt(j,k) = mu(j)*ady(j)*coef/sums(k)
   wgts(j,k) = mu(j)*ady(j)*coef/sumss(k)
   wgtu(j,k) = mu(j)*ady(j)*coef/sumu(k)
   wgtv(j,k) = muv(j)*adyv(j)*coef/sumv(k)
   wgtw(j,k) = mu(j)*ady(j)*coef/sumw(k)
 end do
end do
wgtw(:,nz) = wgtw(:,nzm)
do j=1,ny
  wgtxys(:,j) = wgt(j,1)*terra(1:nx,j,1)
  wgtxyt(:,j) = wgt(j,nzm)*terra(1:nx,j,nzm)
  wgty(j) = mu(j)*ady(j)*nx_gl*ny_gl/sumss(1)
  wgtxys(:,j) = wgty(j)
  wgtxyt(:,j) = wgty(j)
end do

  
if(docoriolis) then
 if(dofplane) then
   if(fcor.eq.-999.) then
     fcor = 4*pi/86400.*sin(latitude0*deg2rad)
     fcorz = 0.
     if(docoriolisz) fcorz =  4*pi/86400.*cos(latitude0*deg2rad) 
   else
     fcorz = sqrt((4*pi/86400.)**2-fcor**2)
   endif
   fcory(:) = fcor
   fcorzy(:) = fcorz
   longitude(:,:) = longitude0
 else
  do j=0,ny
     fcory(j)= 4.*pi/86400.*sin(lat(j+jt)*deg2rad)*earth_factor
  end do
  if(docoriolisz) then
    do j=1,ny
      fcorzy(j) = 4.*pi/86400.*cos(lat(j+jt)*deg2rad)*earth_factor
    end do
  end if
 end if ! dofplane
else
  fcor = 0.
  fcorz = 0.
  fcorzy = 0.
  fcory = 0.
end if

if(masterproc) then
 print*,'lat_gl(min/max):',minval(lat_gl),maxval(lat_gl)
 print*,'latv_gl(min/max):',minval(latv_gl),maxval(latv_gl)
 print*,'lon_gl(min/max):',minval(lon_gl),maxval(lon_gl)
end if

if(maxval(abs(latitude)).ge.90.) then
  print*,'rank=',rank,': error in setting latitude: max(abs(lat))=',maxval(abs(latitude))
  stop
end if
if(maxval(abs(longitude)).gt.360.) then
  print*,'rank=',rank,': error in setting longitude: max(lon)=',maxval(abs(longitude))
  stop
end if

end
