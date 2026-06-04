#include "fppmacros"
     
subroutine write_fields3D
	
use vars
use rad, only: qrad
use params
use sgs, only: tke, dodns, dosmagor,  tkh
use microphysics, only: flag_micro3Dout, index_water_vapor, &
                        total_micro_3Dout,micro_write_fields3D, GET_reffc, Get_reffi
use tracers
use chemistry, only: chem_write_fields3D, nchem_fields_3Dsave
use buildings
use pnetcdf_stuff, only: open_file3D_pnetcdf, close_file_pnetcdf, dopnetcdf
use terrain

implicit none

character(200)filename
character(80) long_name
character(8)  name
character(10) timechar
character(5)  sepchar
character(3)  filetype
character(4)  rankchar
character(10) units
integer i,j,k,n,nfields,nfields1
real(4) tmp(nx,ny,nzm)
real rdx, rdy, rdz
integer, external :: lenstr
real coef
character(9) fileaction
logical flag
integer nsteplast

if(docheck) return

dopnetcdf = save3Dnetcdf

nfields=6 ! number of 3D fields to save
if(dosgs.and..not.dodns) then
  nfields=nfields+1
  if(.not.dosmagor) nfields=nfields+1
end if
if((dolongwave.or.doshortwave).and..not.doradhomo) nfields=nfields+1
!bloss: add 3D outputs for microphysical fields specified by flag_micro3Dout
!       except for water vapor (already output as a SAM default).
if(docloud) nfields=nfields+total_micro_3Dout-flag_micro3Dout(index_water_vapor)
if(dochem) nfields = nfields + nchem_fields_3Dsave
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(dotracers) nfields=nfields+SUM(flag_tracer3Dout)
if(doregion.and.donudge3D) nfields = nfields+5
if(dobuildings) nfields = nfields+6
if(dobuildings.and.doshadows) nfields = nfields+2
if(dobuildings.and.(doirgroundfromwalls &
  .or.doirwallsfromground.or.doirwallsfromwalls)) nfields = nfields+3

nfields1=0

if(.not.dopnetcdf) then

 if(rank.eq.rank-mod(rank,nsubdomains/nfiles3D)) then

  if(nfiles3D.eq.1) then
     sepchar=""
  else
     write(rankchar,'(i4)') rank/(nsubdomains/nfiles3D)
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  end if
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  filetype = '.3D'

  if(save3Dsep) then
     fileaction = 'WRITE'
     filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        trim(date_pr)//'_'//timechar(1:10)//filetype//sepchar
  else
     fileaction = 'READWRITE'
     filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//filetype//sepchar
  end if

  open(46,file=filename,status='unknown', &
          form='unformatted',BUFFEREDYES ACTION=trim(fileaction))
  if(masterproc) print*, 'Writting to file: '//trim(filename)

  nsteplast = 0
  if(.not.save3Dsep) then
       if(nfiles3D.gt.1) then
         if(masterproc) print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         if(masterproc) print*,'save3Dsep requires nfiles3D=1. exit...'
         if(masterproc) print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         call task_abort()
       end if
       do while(.true.)
        if(masterproc) then
         read(46,end=222,err=222)  nsteplast
         if(nsteplast.ge.nstep) then
           backspace(46)
           goto 222   ! yeh, I know, it's bad ....
         end if
         read(46)
         read(46)
         read(46)
         read(46)
         read(46)
         read(46)
        end if
        do i=1,nfields
           read(46)
           read(46) flag
           if(masterproc.and..not.flag) read(46)
           do n=1,nsubdomains 
            read(46)
           end do
        end do
       end do
 222   continue
       if(masterproc) print*,'nsteplast=',nsteplast
  end if

  if(masterproc) then

      write(46) nstep, time
      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields,nfiles3D
      write(46) real(dx,4),real(dy,4)
      write(46) real(dble(nstep)*dt/(3600._8*24._8)+day0,8),datechar,timeUTsec
      write(46) dolatlon
      write(46) 'atm'
      write(46) real(z(1:nzm),4),real(pres(1:nzm),4), &
                lat_gl,lon_gl,latv_gl(1:ny_gl), &
                real(dz*adz(1:nzm),4),real(rho(1:nzm),4),real(rhow(1:nzm),4), &
                real(presi(1:nzm),4),wgty

  end if ! masterproc
 
 end if ! rank.... 

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
       trim(date_pr)//'_'//timechar(1:10)//'.3D_atm.nc'

 call open_file3D_pnetcdf(filename)

 if(masterproc) print*, 'Writting to file: '//trim(filename)
end if ! .not.dopnetcdf



  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=u(i,j,k) + ug
    end do
   end do
  end do
  name='U'
  long_name='X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,4)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=v(i,j,k) + vg
    end do
   end do
  end do
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,3)

!if(doregion.and.donudge3D) then
!
!  nfields1=nfields1+1
!  call box_smooth_3d(u(1:nx,1:ny,1:nzm),terrau(1:nx,1:ny,1:nzm),tmp,nx,ny,nzm, &
!                     nsubdomains_x,nsubdomains_y,10,10,comm,8)
!  name='U_SM'
!  long_name='X Wind Component (SMOOTHED)'   
!  units='m/s'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,4)
!
!  nfields1=nfields1+1
!  call box_smooth_3d(v(1:nx,1:ny,1:nzm),terrav(1:nx,1:ny,1:nzm),tmp,nx,ny,nzm, &
!                     nsubdomains_x,nsubdomains_y,10,10,comm,8)
!  name='V_SM'
!  long_name='Y Wind Component (SMOOTHED)'
!  units='m/s'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,3)
!
!  nfields1=nfields1+1
!  call box_smooth_3d_real8(tabs(1:nx,1:ny,1:nzm),terra(1:nx,1:ny,1:nzm),tmp,nx,ny,nzm, &
!                     nsubdomains_x,nsubdomains_y,10,10,comm,8)
!  name='TABS_SM'
!  long_name='Temperature (SMOOTHED)'
!  units='K'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
!end if


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=w(i,j,k)
    end do
   end do
  end do
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,2)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=p(i,j,k,nb)*rho(k)
    end do
   end do
  end do
  name='P'
  long_name='Pressure Perturbation'
  units='Pa'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tabs(i,j,k)
!      tmp(i,j,k)=t(i,j,k)-t0(k)
    end do
   end do
  end do
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

!  nfields1=nfields1+1
!  tmp = 0.
!  do k=1,nzm
!   coef = rho(k)*adz(k)*dz
!   rdz = 1./coef
!   do j=1,ny
!     rdx = imu(j)/dx
!     rdy = imu(j)/(dy*ady(j))
!     do i=1,nx
!      tmp(i,j,k) = (u(i+1,j,k)-u(i,j,k))*rdx + (muv(j+YES3D)*v(i,j+YES3D,k)-muv(j)*v(i,j,k))*rdy + &
!                   (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
!     end do
!   end do
!  end do
!  name='DIV'
!  long_name='Divergence'
!  units='1/s'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
!
!  nfields1=nfields1+1
!  tmp = 0.
!  do k=1,nzm
!   do j=1,ny
!    do i=1,nx
!      tmp(i,j,k)=imuv(j)*((v(i,j,k)-v(i-1,j,k))/dx-(mu(j)*u(i,j,k)-mu(j-YES3D)*u(i,j-YES3D,k))/(dy*adyv(j)))
!    end do
!   end do
!  end do
!  name='VOR'
!  long_name='Vorticity'
!  units='1/s'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

if(doregion.and.donudge3D) then

  if(dayfld3Dobs(1).gt.0.) then
    coef=(day-dayfld3Dobs(1))/(dayfld3Dobs(2)-dayfld3Dobs(1))
  else
    coef = 0.
  end if

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=uobs(i,j,k,1)+(uobs(i,j,k,2)-uobs(i,j,k,1))*coef
    end do
   end do
  end do
  name='UOBS'
  long_name='Observed X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)= vobs(i,j,k,1)+(vobs(i,j,k,2)-vobs(i,j,k,1))*coef
    end do
   end do
  end do
  name='VOBS'
  long_name='Observed Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)= wobs(i,j,k,1)+(wobs(i,j,k,2)-wobs(i,j,k,1))*coef
    end do
   end do
  end do
  name='WOBS'
  long_name='Observed Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,2)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tobs(i,j,k,1)+(tobs(i,j,k,2)-tobs(i,j,k,1))*coef
    end do
   end do
  end do
  name='TOBS'
  long_name='Observed Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=1.e3*(qobs(i,j,k,1)+(qobs(i,j,k,2)-qobs(i,j,k,1))*coef)
    end do
   end do
  end do
  name='QOBS'
  long_name='Observed Total Water'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

end if ! doregion...

if(dosgs.and..not.dodns) then
   nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=tkh(i,j,k)
     end do
    end do
   end do
   name='TKH'
   long_name='SGS Eddy Diffusivity'
   units='m2/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
   if(.not.dosmagor) then
    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=tke(i,j,k)
      end do
     end do
    end do
    name='STKE'
    long_name='SGS TKE'
    units='m2/s2'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
   end if
end if


if((dolongwave.or.doshortwave).and..not.doradhomo) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qrad(i,j,k)*86400.
    end do
   end do
  end do
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/day'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
end if

if(dobuildings) then

  tmp = 0.

  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=cell_wall(i,j,k)
    end do
   end do
  end do
  name='CELLWALL'
  long_name='Wall Cell indicator'
  units=' '
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=t_face(i,j,k)
    end do
   end do
  end do
  name='T_FACE'
  long_name='Wall skin temperature'
  units='C'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=t_wall(i,j,k)
    end do
   end do
  end do
  name='T_CORE'
  long_name='Wall core temperature'
  units='C'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_wall(i,j,k)
    end do
   end do
  end do
  name='FLXWALL'
  long_name='flux from outside to wall'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_face(i,j,k)
    end do
   end do
  end do
  name='FLXAIR'
  long_name='flux from walls to air'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_bld(i,j,k)
    end do
   end do
  end do
  name='FLXBLD'
  long_name='flux from walls to building'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  if(doshadows) then

    nfields1=nfields1+1
    do k=1,k_face_max
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=cell_sunny(i,j,k)
      end do
     end do
    end do
    name='CELLSUN'
    long_name='Sunny cells'
    units=' '
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                   save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
    nfields1=nfields1+1
    do k=1,k_face_max
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=flx_sw(i,j,k)
      end do
     end do
    end do
    name='FLXSW'
    long_name='Building Solar Flux'
    units='W/m2'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                   save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  end if ! doshadows
  
  if(doirgroundfromwalls.or.doirwallsfromground.or.doirwallsfromwalls) then

  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_sfc_wall(i,j,k)
    end do
   end do
  end do
  name='FLXWS'
  long_name='flux to wall from surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_wall_wall(i,j,k)
    end do
   end do
  end do
  name='FLXWW'
  long_name='flux walls from walls'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,k_face_max
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=flx_lw(i,j,k)
    end do
   end do
  end do
  name='FLXLW'
  long_name='Building IR flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  end if
end if

if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffc()
  name='REL'
  long_name='Effective Radius for Cloud Liquid Water'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
end if
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffi()
  name='REI'
  long_name='Effective Radius for Cloud Ice'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
end if


!  nfields1=nfields1+1
!  do k=1,nzm
!   do j=1,ny
!    do i=1,nx
!      tmp(i,j,k)=t(i,j,k)
!    end do
!   end do
!  end do
!  name='TT'
!  long_name='Liquid-Ice water static energy'
!  units='K'
!  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qv(i,j,k)*1.e3
    end do
   end do
  end do
  name='QV'
  long_name='Water Vapor'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)


!-------------------------------------------------------------------
! output microphysics fields:

 call micro_write_fields3D(nfields1)

!-----------------------------------------------------------------
! output chemical fields:

if(dochem) call chem_write_fields3D(nfields1)

!-----------------------------------------------------------------


 do n = 1,ntracers
   if(dotracers.AND.flag_tracer3Dout(n).gt.0) then
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,ny
            do i=1,nx
               tmp(i,j,k)=tracer(i,j,k,n)
            end do
         end do
         ! remove factor of rho from number, if this field is a number concentration
      end do
      name=TRIM(tracername(n))
      long_name=TRIM(tracername(n))
      units=TRIM(tracerunits(n))
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           .true.,dompi,rank,nsubdomains,nfiles3D,1)
   end if
 end do

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields=',nfields,'  nfields1=',nfields1
    call task_abort()
  end if

  if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
  else
   if(masterproc) then
         close(46)
         if(save2Dsep.and.dogzip3D) call systemf('gzip -f '//filename)
   end if
  end if
  if(masterproc) print*, 'Done. ', nfields1,'fields'
 
end
