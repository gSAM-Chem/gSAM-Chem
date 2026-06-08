
!================================================================
! write terrain masks to invariant directory:
! Feb 2021, MK

subroutine write_terrain()

use terrain
use vars, only: rho, rhow
use grid, only: nfiles3D, datechar, timeUTsec
use params, only: dolatlon, dobuildings, docheck
use buildings, only: normal_face1, normal_face2, k_face_max
use pnetcdf_stuff, only: open_file3D_pnetcdf, close_file_pnetcdf, dopnetcdf

implicit none

character *200 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *3 filetype
character *10 units
integer i,j,k,n,nfields,nf
real(4) tmp(nx,ny,nzm)
integer, external :: lenstr
logical savebin

if(docheck) return

dopnetcdf = save3Dnetcdf

nfields = 5  !number of 3D fields to save
if(dobuildings) nfields = nfields+2
savebin = .true.

if(.not.dopnetcdf) then

 nf = 1 ! nfiles3D
 filetype = '.3D'
 filename='./OUT_INV/'//trim(case)//'_'//trim(caseid)//'_TERR_MASKS'//filetype
 open(46,file=filename,status='unknown',form='unformatted',ACTION='WRITE')
 if(masterproc) print*,'saving terrain masks to ',trim(filename)
 
 if(masterproc) then

      write(46) nstep, time, datechar, timeUTsec
      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields,nfiles3D
      write(46) real(dx,4),real(dy,4)
      write(46) real(dble(nstep)*dt/(3600._8*24._8)+day0,8),datechar,0._8
      write(46) dolatlon
      write(46) 'atm'
      write(46) real(z(1:nzm),4),real(pres(1:nzm),4), &
                lat_gl,lon_gl,latv_gl(1:ny_gl), &
                real(dz*adz(1:nzm),4),real(rho(1:nzm),4),real(rhow(1:nzm),4), &
                real(presi(1:nzm),4),wgty

 end if ! masterproc


else

! open netcdf file:

 write(rankchar,'(i4)') nsubdomains
 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_INV/'//trim(case)//'_'//trim(caseid)//'_'// &
       rankchar(5-lenstr(rankchar):4)//'_TERR_MASKS.3D_atm.nc'

 if(masterproc) print*,'saving terrain masks to ',trim(filename)
 call open_file3D_pnetcdf(filename)

end if ! .not.dopnetcdf


  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terra(i,j,k)
    end do
   end do
  end do
  name='TERRA'
  long_name='TERR MASK SCALARS'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 savebin,dompi,rank,nsubdomains,nf,1)

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terrau(i,j,k)
    end do
   end do
  end do
  name='TERRAU'
  long_name='TERR MASK U-WIND'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 savebin,dompi,rank,nsubdomains,nf,4)

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terrav(i,j,k)
    end do
   end do
  end do
  name='TERRAV'
  long_name='TERR MASK V-WIND'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 savebin,dompi,rank,nsubdomains,nf,3)

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=terraw(i,j,k)
    end do
   end do
  end do
  name='TERRAW'
  long_name='TERR MASK W-WIND'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 savebin,dompi,rank,nsubdomains,nf,2)

  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=alphah(i,j,k)
    end do
   end do
  end do
  name='ALPHA'
  long_name='ALPHA HYBRID'
  units=''
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 savebin,dompi,rank,nsubdomains,nf,1)

  if(dobuildings) then

   tmp=0.

   do k=1,k_face_max
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=normal_face1(i,j,k)
     end do
    end do
   end do
   name='NORMALX'
   long_name='x component of normal'
   units=''
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  savebin,dompi,rank,nsubdomains,nf,1)

   do k=1,k_face_max
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=normal_face2(i,j,k)
     end do
    end do
   end do
   name='NORMALY'
   long_name='y component of normal'
   units=''
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  savebin,dompi,rank,nsubdomains,nf,1)
  end if

if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
else
  close (46)
end if


end subroutine write_terrain

