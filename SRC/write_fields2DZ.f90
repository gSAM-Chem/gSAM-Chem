#include "fppmacros"
     
subroutine write_fields2DZ
	
use vars
use rad, only: qrad
use params
use terrain, only: terrau,terrav,terraw,terra
use sgs, only: tkh
use microphysics, only: nmicro_fields, micro_field, flag_number, &
     flag_micro2DZout, mkname, mklongname, mkunits, mkoutputscale, &
     index_water_vapor, GET_reffc, Get_reffi
use pnetcdf_stuff, only: open_file2DZ_pnetcdf, close_file_pnetcdf, dopnetcdf

implicit none

character *200 filename
character *80 long_name
character *8 name, name1
character *10 timechar
character *4 filetype
character *10 units
integer i,j,k,n,nfields,nfields1,nsteplast
real tmp(nx,ny,nzm),tmpz(ny,nzm),wgtx(ny,nzm)
real(4) tmpzz(ny,nzm)
character*7 filestatus
integer, external :: lenstr
logical:: save2DZbin =.true.
logical flag

if(docheck) return

if(RUN2D) return

dopnetcdf = save2Dnetcdf

nfields=10 ! number of 2DZ fields to save
if(.not.docloud) nfields=nfields-1
if(.not.doprecip) nfields=nfields-1
if(docloud) nfields=nfields+SUM(flag_micro2DZout)-flag_micro2DZout(index_water_vapor)
if((dolongwave.or.doshortwave).and..not.doradhomo) nfields=nfields+1
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(doterrain) nfields=nfields+1

nfields1=0

if(.not.dopnetcdf) then

 if(masterproc) then
 
  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do

! Make sure that the new run doesn't overwrite the file from the old run

    filetype = '.2DZ'
   
    if(save2DZsep) then
      filename='./OUT_2DZ/'//trim(case)//'_'//trim(caseid)//'_'// &
             trim(date_pr)//'_'//timechar(1:10)//filetype
      open(46,file=filename,form='unformatted',BUFFEREDYES ACTION='WRITE')
      if(masterproc) print*, 'Writting to file: '//trim(filename)
    else
      filename='./OUT_2DZ/'//trim(case)//'_'//trim(caseid)//filetype
      open(46,file=filename,form='unformatted',BUFFEREDYES ACTION='READWRITE')
      if(masterproc) print*, 'Writting to file: '//trim(filename)
      do while(.true.)
         read(46,end=222)  nsteplast
         if(nsteplast.ge.nstep) then
           backspace(46)
           goto 222   ! yeh, I know, it's bad ....
         end if
         read(46)
         read(46)
         read(46)
         do i=1,nfields
           read(46) name1
           read(46) flag
           if(.not.flag) read(46)
           do j=0,nsubdomains-1
            if(mod(j,nsubdomains_x).eq.0) then
               read(46)
            end if
           end do
         end do
      end do
 222  continue
      print*,'nsteplast=',nsteplast
    end if

    write(46) nstep, time,datechar,timeUTsec
    write(46) 'atm'
    write(46) nx,ny,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dy,4),lat_gl,dble(nstep)*dt/(3600._8*24._8)+day0,&
              (real(z(k),4),k=1,nzm),(real(pres(k),4),k=1,nzm), &
               real(dz*adz(1:nzm),4),real(rho(1:nzm),4),real(rhow(1:nzm),4), &
               latv_gl,real(y_gl,4),real(yv_gl,4),wgty

 end if ! masterproc

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_2DZ/'//trim(case)//'_'//trim(caseid)//'_'// &
       trim(date_pr)//'_'//timechar(1:10)//'.2DZ_atm.nc'

 call open_file2DZ_pnetcdf(filename)

 if(masterproc) print*, 'Writting to file: '//trim(filename)

end if ! .not.dopnetcdf


!------------------------------
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terrau(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,wgtx)
  do k=1,nzm
   do j=1,ny
     wgtx(j,k)=1./(wgtx(j,k)+1.e-10)
   end do
  end do

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(u(i,j,k)+ug)*wgtx(j,k)*terrau(i,j,k)
    end do
   end do
  end do
  name='U'
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  long_name='X Wind Component'
  units='m/s'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terrav(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,wgtx)
  do k=1,nzm
   do j=1,ny
     wgtx(j,k)=1./(wgtx(j,k)+1.e-10)
   end do
  end do

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(v(i,j,k)+vg)*wgtx(j,k)*terrav(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terraw(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,wgtx)
  do k=1,nzm
   do j=1,ny
     wgtx(j,k)=1./(wgtx(j,k)+1.e-10)
   end do
  end do

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=w(i,j,k)*wgtx(j,k)*terraw(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,2)

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,wgtx)
  do k=1,nzm
   do j=1,ny
     wgtx(j,k)=1./(wgtx(j,k)+1.e-10)
   end do
  end do

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tkh(i,j,k)*wgtx(j,k)*terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='KZH'
  long_name='eddy diffusivity'
  units='m2/s'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=pp(i,j,k)*wgtx(j,k)*terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='PP'
  long_name='Pressure'
  units='hPa'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)


if((dolongwave.or.doshortwave).and..not.doradhomo) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qrad(i,j,k)*wgtx(j,k)*terra(i,j,k)*86400.
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/day'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffc()
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tmp(i,j,k)*wgtx(j,k)*terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='REL'
  long_name='Effective Radius for Cloud Liquid Water'
  units='mkm'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffi()
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tmp(i,j,k)*wgtx(j,k)*terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='REI'
  long_name='Effective Radius for Cloud Ice'
  units='mkm'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if

if(doterrain) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='TERRA'
  long_name='Terrain Mask'
  units=' '
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tabs(i,j,k)*terra(i,j,k)*wgtx(j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=t(i,j,k)*wgtx(j,k)*terra(i,j,k)
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='TT'
  long_name='Liquid-Ice water static energy'
  units='K'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qv(i,j,k)*wgtx(j,k)*terra(i,j,k)*1.e3
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='QV'
  long_name='Water Vapor'
  units='g/kg'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)

if(docloud) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qcl(i,j,k)+qci(i,j,k))*wgtx(j,k)*terra(i,j,k)*1.e3
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='QN'
  long_name='Cloud Water + cloud Ice)'
  units='g/kg'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if


if(doprecip) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qpl(i,j,k)+qpi(i,j,k))*wgtx(j,k)*terra(i,j,k)*1.e3
    end do
   end do
  end do
  call mean_x_3D(tmp,tmpz)
  tmpzz = tmpz
  name='QP'
  long_name='Precipitating Water/Ice'
  units='g/kg'
  call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
                                 save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
end if


do n = 1,nmicro_fields
   if(docloud.AND.flag_micro2DZout(n).gt.0.AND.n.ne.index_water_vapor) then
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,ny
            do i=1,nx
               tmp(i,j,k)=micro_field(i,j,k,n)*wgtx(j,k)*terra(i,j,k)*mkoutputscale(n)
            end do
         end do
         ! remove factor of rho from number, if this field is a number concentration
         if(flag_number(n).gt.0) tmp(:,:,k) = tmp(:,:,k)*rho(k)
      end do
      call mean_x_3D(tmp,tmpz)
      tmpzz = tmpz
      name=TRIM(mkname(n))
      long_name=TRIM(mklongname(n))
      units=TRIM(mkunits(n))
      call compress2DZ(tmpzz,ny,nzm,name,long_name,units, &
           save2DZbin,dompi,rank,nsubdomains,nsubdomains_x,1)
   end if
end do

call task_barrier()

if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields2DZ error: nfields=',nfields,'  nfields1=',nfields1
    call task_abort()
end if

if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
else
  if(masterproc) then
     close(46)
     if(save2DZsep.and.dogzip2DZ) call systemf('gzip -f '//filename)
  endif
end if
if(masterproc) print*, 'Done.', nfields1,'fields'

 
end
