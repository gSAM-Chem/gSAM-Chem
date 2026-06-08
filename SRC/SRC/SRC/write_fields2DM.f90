#include "fppmacros"
     
subroutine write_fields2DM
	
use vars
use params
use pnetcdf_stuff, only: open_file2D_pnetcdf, close_file_pnetcdf, dopnetcdf

implicit none

character(200)filename
character(80) long_name
character(8)  name, name1
character(10) timechar
character(3)  filetype
character(10) units
integer i,j,k,nfields,nfields1,nsteplast
real(4) tmp(nx,ny)
real coef, coefa
integer, external :: lenstr
logical flag

nfields= 12
if(.not.(ISLAND.or.LAND)) nfields=nfields-2
!===================================================================

dopnetcdf = saveMnetcdf

nfields1=0


if(.not.dopnetcdf) then

 if(masterproc) then

  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do

! Make sure that the new run doesn't overwrite the file from the old run 

    filetype = '.2D'
    
    if(saveMsep) then
       filename='./OUT_MISC/'//trim(case)//'_'//trim(caseid)//'_MISC_'// &
          trim(date_pr)//'_'//timechar(1:10)//filetype 
          open(46,file=filename,form='unformatted',BUFFEREDYES ACTION='WRITE')
       print*, 'Writting to file: '//trim(filename)
    else
       filename='./OUT_MISC/'//trim(case)//'_'//trim(caseid)//'_MISC'//filetype
       open(46,file=filename,form='unformatted',BUFFEREDYES ACTION='READWRITE')	
       print*, 'Writting to file: '//trim(filename)
       do while(.true.)
         read(46,end=222)  nsteplast
         if(nsteplast.ge.nstep) then
           backspace(46)
           goto 222   ! yeh, I know, it's bad ....
         end if
         read(46)
         read(46)
         read(46)
         read(46)
         read(46)
         do i=1,nfields
           read(46)
           read(46) flag
           if(.not.flag) read(46)
           read(46)
         end do
       end do
 222    continue
       print*,'nsteplast=',nsteplast
       notopenedM=.false. 
    end if

    write(46) nstep, dolatlon, time, datechar, timeUTsec
    write(46) 'atm'
    write(46) nx,ny,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dx,4),lat_gl,latv_gl,real(y_gl,4),real(yv_gl,4),wgty
    write(46) real(dy,4),lon_gl,lonu_gl
    write(46) dble(nstep)*dt/(3600._8*24._8)+day0

 end if! masterproc

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_MISC/'//trim(case)//'_'//trim(caseid)//'_MISC_'// &
       trim(date_pr)//'_'//timechar(1:10)//'.2D_atm.nc'

 call open_file2D_pnetcdf(filename)

 if(masterproc) print*, 'Writting to file: '//trim(filename)

end if ! .not.dopnetcdf



if(.not.nstep.eq.1.and.saveMavg) then
   coef = 1./float(nsaveM)
else
   coef = 1.
end if
coefa = 1./float(nsaveM)


! MISC fields:


if(ISLAND.or.LAND) then

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=albvis_xy(i,j)*coef
    end do
   end do
  name='ALBVIS'
  long_name='Visible Albedo (diffuse)'
  units=' '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=albnir_xy(i,j)*coef
    end do
   end do
  name='ALBNIR'
  long_name='Near-IR Albedo (diffuse)'
  units=' '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

end if

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= alb_xy(i,j)
    end do
   end do
  name='ALBZ'
  long_name='visible albedo for cos(z)=1'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= albc_xy(i,j)
    end do
   end do
  name='ALBCZ'
  long_name='clear-sky visible albedo for cos(z)=1'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= sst_min_xy(i,j)-273.16
    end do
   end do
  name='TMIN'
  long_name='Minumum Surface Temperature'
  units='C'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= sst_max_xy(i,j)-273.16
    end do
   end do
  name='TMAX'
  long_name='Maximum Surface Temperature'
  units='C'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= cflz_max_xy(i,j)
    end do
   end do
  name='CFLZ'
  long_name='Maximum CFL in Z'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= cflh_max_xy(i,j)
    end do
   end do
  name='CFLH'
  long_name='Maximum CFL in xy'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= cfl_max_xy(i,j)
    end do
   end do
  name='CFL'
  long_name='Maximum CFL'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= zcflz_max_xy(i,j)
    end do
   end do
  name='ZCFLZ'
  long_name='Height of Maximum CFL in Z'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= zcflh_max_xy(i,j)
    end do
   end do
  name='ZCFLH'
  long_name='Height of Maximum CFL in xy'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)= zcfl_max_xy(i,j)
    end do
   end do
  name='ZCFL'
  long_name='Height of Maximum CFL'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               saveMbin,dompi,rank,nsubdomains)

!=====================================================

call task_barrier()

!===================================================================


if(nfields.ne.nfields1) then
  if(masterproc) print*,'write_fields2DM: error in nfields!!',' nfields=',nfields,'nfields1=',nfields1
  call task_abort()
end if

if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
else
  if(masterproc) then
     close(46)
     if(saveMsep.and.dogzipM) call systemf('gzip -f '//filename)
  endif
end if
if(masterproc) print*, 'Done.', nfields1,'fields'


end
