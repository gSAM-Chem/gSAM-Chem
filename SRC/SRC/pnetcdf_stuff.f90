#include "fppmacros"

#ifdef pnetcdf_out

! module to handle parallel NetCDF output using pNetCDF library.
! Marat Khairoutdinov 2023

module pnetcdf_stuff

use mpi
use pnetcdf
use grid

implicit none

integer ncid      ! NetCDF file ID 
integer dimids2D(3) ! Dimension IDs 
integer dimids3DS(4) ! Dimension IDs 
integer dimids3DU ! Dimension IDs 
integer dimids3DV ! Dimension IDs 
integer dimids3DW ! Dimension IDs 
integer dimids2DZ(3)! Dimension IDs 
integer(kind=MPI_OFFSET_KIND) :: dim_sizes(4) ! Sizes of the dimensions
logical:: dopnetcdf = .false.

CONTAINS

!================================================================================================

SUBROUTINE open_file3D_pnetcdf(filename)

  use vars, only: rho, pres, rhow,presi
  use params, only: dolatlon
  character(*), intent(in) :: filename ! Name of the NetCDF file
  real x(nx_gl), y(ny_gl)
  integer ierr, i, j
  integer(kind=MPI_OFFSET_KIND) start(1), count(1)
  integer vdimids(4) ! Dimension variable IDs for scalar (cell center)
  integer vdimidsU(4) ! Dimension variable IDs for u
  integer vdimidsV(4) ! Dimension variable IDs for v
  integer vdimidsW(4) ! Dimension variable IDs for w
  integer dayid,rhoid,pid,dzid,rhowid,pwid
  real(8) tmp(1:nx_gl)
  real dzz(nzm)

  ! Open the NetCDF file in parallel mode
  ierr =  nf90mpi_create(comm, filename, IOR(NF90_CLOBBER, NF90_64BIT_DATA), MPI_INFO_NULL, ncid)
  IF (ierr .ne. NF90_NOERR) THEN
    print *, 'Error opening NetCDF file', trim(filename),ierr,NF90_NOERR,nf90mpi_strerror(ierr)
    call task_abort()
  END IF

  ! Define the dimensions
  
  if(dolatlon) then
   ierr =  nf90mpi_def_dim(ncid, 'lon', int(nx_gl,MPI_OFFSET_KIND), dimids3DS(1))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining lon dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'lat', int(ny_gl,MPI_OFFSET_KIND), dimids3DS(2))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining lat dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'lonu', int(nx_gl,MPI_OFFSET_KIND), dimids3DU)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining lonu dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'latv', int(ny_gl,MPI_OFFSET_KIND), dimids3DV)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining latv dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
  else
   ierr =  nf90mpi_def_dim(ncid, 'x', int(nx_gl,MPI_OFFSET_KIND), dimids3DS(1))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining x dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'y', int(ny_gl,MPI_OFFSET_KIND), dimids3DS(2))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining y dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'xu', int(nx_gl,MPI_OFFSET_KIND), dimids3DU)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining lonu dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'yv', int(ny_gl,MPI_OFFSET_KIND), dimids3DV)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining yv dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
  end if
  if(dolatlon) then
   ierr =  nf90mpi_def_dim(ncid, 'lev', int(nzm,MPI_OFFSET_KIND), dimids3DS(3))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining z dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'levi', int(nzm,MPI_OFFSET_KIND), dimids3DW)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining zi dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
  else
   ierr =  nf90mpi_def_dim(ncid, 'z', int(nzm,MPI_OFFSET_KIND), dimids3DS(3))
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining z dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
   ierr =  nf90mpi_def_dim(ncid, 'zi', int(nzm,MPI_OFFSET_KIND), dimids3DW)
   IF (ierr .ne. NF90_NOERR) THEN
       print *, 'Error defining zi dimension in NetCDF file:',nf90mpi_strerror(ierr)
       call task_abort()
   END IF
  end if
  ierr =  nf90mpi_def_dim(ncid, 'time', int(NF90_UNLIMITED,MPI_OFFSET_KIND), dimids3DS(4))
  IF (ierr .ne. NF90_NOERR) THEN
      print *, 'Error defining time dimension in NetCDF file:',nf90mpi_strerror(ierr)
      call task_abort()
  END IF

  if(dolatlon) then
    ierr = nf90mpi_def_var(ncid, 'lon',NF90_DOUBLE,dimids3DS(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','longitude')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','degrees_east')
    ierr = nf90mpi_def_var(ncid, 'lat',NF90_DOUBLE,dimids3DS(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','latitude')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','degrees_north')
    ierr = nf90mpi_def_var(ncid, 'lonu',NF90_DOUBLE,dimids3DU,vdimidsU(1))
    ierr = nf90mpi_put_att(ncid, vdimidsU(1),'long_name','u-longitude')
    ierr = nf90mpi_put_att(ncid, vdimidsU(1),'units','degrees_east')
    ierr = nf90mpi_def_var(ncid, 'latv',NF90_DOUBLE,dimids3DV,vdimidsV(2))
    ierr = nf90mpi_put_att(ncid, vdimidsV(2),'long_name','v-latitude')
    ierr = nf90mpi_put_att(ncid, vdimidsV(2),'units','degrees_north')
    vdimidsU(2) = vdimids(2)
    vdimidsV(1) = vdimids(1)
    vdimidsW(1:2) = vdimids(1:2)
  else
    ierr = nf90mpi_def_var(ncid, 'x',NF90_FLOAT,dimids3DS(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','distance in y')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','m')
    ierr = nf90mpi_def_var(ncid, 'y',NF90_FLOAT,dimids3DS(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','distance in x')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','m')
    ierr = nf90mpi_def_var(ncid, 'xu',NF90_FLOAT,dimids3DU,vdimidsU(1))
    ierr = nf90mpi_put_att(ncid, vdimidsU(1),'long_name','distance in y')
    ierr = nf90mpi_put_att(ncid, vdimidsU(1),'units','m')
    ierr = nf90mpi_def_var(ncid, 'yv',NF90_FLOAT,dimids3DV,vdimidsV(2))
    ierr = nf90mpi_put_att(ncid, vdimidsV(2),'long_name','distance in x')
    ierr = nf90mpi_put_att(ncid, vdimidsV(2),'units','m')
    vdimidsU(2) = vdimids(2)
    vdimidsV(1) = vdimids(1)
    vdimidsW(1:2) = vdimids(1:2)
  end if
  if(dolatlon) then
    ierr = nf90mpi_def_var(ncid, 'lev',NF90_FLOAT,dimids3DS(3),vdimids(3))
    ierr = nf90mpi_put_att(ncid, vdimids(3),'long_name','level pressure')
    ierr = nf90mpi_put_att(ncid, vdimids(3),'units','m')
    ierr = nf90mpi_def_var(ncid, 'levi',NF90_FLOAT,dimids3DW,vdimidsW(3))
    ierr = nf90mpi_put_att(ncid, vdimidsW(3),'long_name','interface height')
    ierr = nf90mpi_put_att(ncid, vdimidsW(3),'units','m')
  else
    ierr = nf90mpi_def_var(ncid, 'z',NF90_FLOAT,dimids3DS(3),vdimids(3))
    ierr = nf90mpi_put_att(ncid, vdimids(3),'long_name','level height')
    ierr = nf90mpi_put_att(ncid, vdimids(3),'units','m')
    ierr = nf90mpi_def_var(ncid, 'zi',NF90_FLOAT,dimids3DW,vdimidsW(3))
    ierr = nf90mpi_put_att(ncid, vdimidsW(3),'long_name','interface height')
    ierr = nf90mpi_put_att(ncid, vdimidsW(3),'units','m')
  end if
  ierr = nf90mpi_def_var(ncid, 'time',NF90_DOUBLE,dimids3DS(4),vdimids(4))
  ierr = nf90mpi_put_att(ncid, vdimids(4),'long_name','time')
  if(timeUTsec.gt.0._8) then
   ierr = nf90mpi_put_att(ncid, vdimids(4),'units', &
                        'seconds since 1900-01-01 00:00:00.0')
   ierr = nf90mpi_put_att(ncid, vdimids(4),'calendar','gregorian')
   ierr = nf90mpi_def_var(ncid, 'day',NF90_DOUBLE,dimids3DS(4),dayid)
   ierr = nf90mpi_put_att(ncid, dayid,'units','day')
   ierr = nf90mpi_put_att(ncid, dayid,'long_name','day')
  else
   ierr = nf90mpi_put_att(ncid, vdimids(4),'units','day')
  end if
  vdimidsU(3:4) = vdimids(3:4)
  vdimidsV(3:4) = vdimids(3:4)
  vdimidsW(4) = vdimids(4)
  ierr = nf90mpi_def_var(ncid, 'dz',NF90_FLOAT,dimids3DS(3),dzid)
  ierr = nf90mpi_put_att(ncid, dzid,'long_name','layer thickness')
  ierr = nf90mpi_put_att(ncid, dzid,'units','m')
  ierr = nf90mpi_def_var(ncid, 'rho',NF90_FLOAT,dimids3DS(3),rhoid)
  ierr = nf90mpi_put_att(ncid, rhoid,'long_name','air density')
  ierr = nf90mpi_put_att(ncid, rhoid,'units','kg/m3')
  ierr = nf90mpi_def_var(ncid, 'rhoi',NF90_FLOAT,dimids3DW,rhowid)
  ierr = nf90mpi_put_att(ncid, rhowid,'long_name','air density at interface')
  ierr = nf90mpi_put_att(ncid, rhowid,'units','kg/m3')
  if(dolatlon) then
   ierr = nf90mpi_def_var(ncid, 'z',NF90_FLOAT,dimids3DS(3),pid)
   ierr = nf90mpi_put_att(ncid, pid,'long_name','level height')
   ierr = nf90mpi_put_att(ncid, pid,'units','m')
   ierr = nf90mpi_def_var(ncid, 'zi',NF90_FLOAT,dimids3DW,pwid)
   ierr = nf90mpi_put_att(ncid, pwid,'long_name','interface height')
   ierr = nf90mpi_put_att(ncid, pwid,'units','m')
  else
   ierr = nf90mpi_def_var(ncid, 'p',NF90_FLOAT,dimids3DS(3),pid)
   ierr = nf90mpi_put_att(ncid, pid,'long_name','level pressure')
   ierr = nf90mpi_put_att(ncid, pid,'units','mb')
   ierr = nf90mpi_def_var(ncid, 'pi',NF90_FLOAT,dimids3DW,pwid)
   ierr = nf90mpi_put_att(ncid, pwid,'long_name','pressure at interface')
   ierr = nf90mpi_put_att(ncid, pwid,'units','mb')
  end if

  ierr = nf90mpi_enddef(ncid)

  ierr = nf90mpi_begin_indep_data(ncid)
  if(dolatlon) then
      ierr = nf90mpi_put_var(ncid,vdimids(1),lon_gl(1:nx_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),lat_gl(1:ny_gl))
      tmp(1:nx_gl) = lon_gl(1:nx_gl)-0.5*(lon_gl(2)-lon_gl(1))
      ierr = nf90mpi_put_var(ncid,vdimidsU(1),tmp)
      ierr = nf90mpi_put_var(ncid,vdimidsV(2),latv_gl(1:ny_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(3),pres(1:nzm))
      ierr = nf90mpi_put_var(ncid,vdimidsW(3),presi(1:nzm))
      ierr = nf90mpi_put_var(ncid,pid,z(1:nzm))
      ierr = nf90mpi_put_var(ncid,pwid,zi(1:nzm))
  else
      do i=1,nx_gl
        x(i) = (i-1)*dx
      end do
      do j=1,ny_gl
        y(j) = (j-1)*dy
      end do
      ierr = nf90mpi_put_var(ncid,vdimids(1),x(1:nx_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),y(1:ny_gl))
      do i=1,nx_gl
        x(i) = (i-1.5)*dx
      end do
      do j=1,ny_gl
        y(j) = (j-1.5)*dy
      end do
      ierr = nf90mpi_put_var(ncid,vdimidsU(1),x(1:nx_gl))
      ierr = nf90mpi_put_var(ncid,vdimidsV(2),y(1:ny_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(3),z(1:nzm))
      ierr = nf90mpi_put_var(ncid,vdimidsW(3),zi(1:nzm))
      ierr = nf90mpi_put_var(ncid,pid,pres(1:nzm))
      ierr = nf90mpi_put_var(ncid,pwid,presi(1:nzm))
  end if
  if(timeUTsec.gt.0._8) then
    ierr = nf90mpi_put_var(ncid,vdimids(4),timeUTsec)
    ierr = nf90mpi_put_var(ncid,dayid,day)
  else
    ierr = nf90mpi_put_var(ncid,vdimids(4),day)
  end if
  dzz(:) = dz*adz(:)
  ierr = nf90mpi_put_var(ncid,dzid,dzz(1:nzm))
  ierr = nf90mpi_put_var(ncid,rhoid,rho(1:nzm))
  ierr = nf90mpi_put_var(ncid,rhowid,rhow(1:nzm))

  ierr = nf90mpi_end_indep_data(ncid)

  call task_barrier()

END SUBROUTINE open_file3D_pnetcdf

!================================================================================================

SUBROUTINE open_file2D_pnetcdf(filename)

  use params, only: dolatlon
  character(*), intent(in) :: filename ! Name of the NetCDF file

  character(10) dimname(3)
  real x(nx_gl), y(ny_gl)
  integer ierr, i, j
  integer(kind=MPI_OFFSET_KIND) start(1), count(1)
  integer dayid, vdimids(3) ! Dimension variable IDs 

  ! Open the NetCDF file in parallel mode
  ierr =  nf90mpi_create(comm, filename, IOR(NF90_CLOBBER, NF90_64BIT_DATA), MPI_INFO_NULL, ncid)
  IF (ierr .ne. NF90_NOERR) THEN
    print *, 'Error opening NetCDF file', trim(filename), ierr,NF90_NOERR,nf90mpi_strerror(ierr)
    call task_abort()
  END IF

  ! Define the dimensions
  
  if(dolatlon) then
    dimname(1) = 'lon'
    dimname(2) = 'lat'
  else
    dimname(1) = 'x'
    dimname(2) = 'y'
  end if
  dimname(3) = 'time'
  dim_sizes(1) = nx_gl
  dim_sizes(2) = ny_gl
  dim_sizes(3) = NF90_UNLIMITED

  DO i = 1, 3
    ierr =  nf90mpi_def_dim(ncid, trim(dimname(i)), dim_sizes(i), dimids2D(i))
    IF (ierr .ne. NF90_NOERR) THEN
      print *, 'Error defining dimension in NetCDF file:', nf90mpi_strerror(ierr)
      call MPI_Abort(comm, ierr, ierr)
      RETURN
    END IF
  ENDDO

  if(dolatlon) then
    ierr = nf90mpi_def_var(ncid, 'lon',NF90_DOUBLE,dimids2D(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','longitude')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','degrees_east')
    ierr = nf90mpi_def_var(ncid, 'lat',NF90_DOUBLE,dimids2D(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','latitude')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','degrees_north')
  else
    ierr = nf90mpi_def_var(ncid, 'x',NF90_FLOAT,dimids2D(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','distance in y')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','m')
    ierr = nf90mpi_def_var(ncid, 'y',NF90_FLOAT,dimids2D(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','distance in x')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','m')
  end if
  ierr = nf90mpi_def_var(ncid, 'time',NF90_DOUBLE,dimids2D(3),vdimids(3))
  ierr = nf90mpi_put_att(ncid, vdimids(3),'long_name','time')
  if(timeUTsec.gt.0._8) then
   ierr = nf90mpi_put_att(ncid, vdimids(3),'units', &
                        'seconds since 1900-01-01 00:00:00.0')
   ierr = nf90mpi_put_att(ncid, vdimids(3),'calendar','gregorian')
   ierr = nf90mpi_def_var(ncid, 'day',NF90_DOUBLE,dimids2D(3),dayid)
   ierr = nf90mpi_put_att(ncid, dayid,'units','day')
   ierr = nf90mpi_put_att(ncid, dayid,'long_name','day')
  else
   ierr = nf90mpi_put_att(ncid, vdimids(3),'units','day')
  end if

  ierr = nf90mpi_enddef(ncid)

  ierr = nf90mpi_begin_indep_data(ncid)
  if(dolatlon) then
      ierr = nf90mpi_put_var(ncid,vdimids(1),lon_gl(1:nx_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),lat_gl(1:ny_gl))
  else
      do i=1,nx_gl
        x(i) = (i-1)*dx
      end do
      do j=1,ny_gl
        y(j) = (j-1)*dy
      end do
      ierr = nf90mpi_put_var(ncid,vdimids(1),x(1:nx_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),y(1:ny_gl))
  end if
  if(timeUTsec.gt.0._8) then
    ierr = nf90mpi_put_var(ncid,vdimids(3),timeUTsec)
    ierr = nf90mpi_put_var(ncid,dayid,day)
  else
    ierr = nf90mpi_put_var(ncid,vdimids(3),day)
  end if
  ierr = nf90mpi_end_indep_data(ncid)

  call task_barrier()

END SUBROUTINE open_file2D_pnetcdf

!======================================================================================================

SUBROUTINE open_file2DZ_pnetcdf(filename)

  use vars, only: rho
  use params, only: dolatlon
  character(*), intent(in) :: filename ! Name of the NetCDF file

  character(10) dimname(3)
  real y(ny_gl)
  integer ierr, i, j
  integer(kind=MPI_OFFSET_KIND) start(1), count(1)
  integer dayid, vdimids(3) ! Dimension variable IDs 
  integer dzid, pid, rhoid
  real dzz(nzm)

  ! Open the NetCDF file in parallel mode
  ierr =  nf90mpi_create(comm, filename, IOR(NF90_CLOBBER, NF90_64BIT_DATA), MPI_INFO_NULL, ncid)
  IF (ierr .ne. NF90_NOERR) THEN
    print *, 'Error opening NetCDF file', trim(filename),ierr,NF90_NOERR, nf90mpi_strerror(ierr)
    call task_abort()
  END IF

  ! Define the dimensions
  
  if(dolatlon) then
    dimname(1) = 'lat'
    dimname(2) = 'lev'
  else
    dimname(1) = 'y'
    dimname(2) = 'z'
  end if
  dimname(3) = 'time'
  dim_sizes(1) = ny_gl
  dim_sizes(2) = nzm
  dim_sizes(3) = NF90_UNLIMITED

  DO i = 1, 3
    ierr =  nf90mpi_def_dim(ncid, trim(dimname(i)), dim_sizes(i), dimids2DZ(i))
    IF (ierr .ne. NF90_NOERR) THEN
      print *, 'Error defining dimension in NetCDF file:', nf90mpi_strerror(ierr)
      call MPI_Abort(comm, ierr, ierr)
      RETURN
    END IF
  ENDDO

  if(dolatlon) then
    ierr = nf90mpi_def_var(ncid, 'lat',NF90_DOUBLE,dimids2DZ(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','latitude')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','degrees_north')
    ierr = nf90mpi_def_var(ncid, 'lev',NF90_FLOAT,dimids2DZ(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','pressure')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','mb')
  else
    ierr = nf90mpi_def_var(ncid, 'y',NF90_FLOAT,dimids2DZ(1),vdimids(1))
    ierr = nf90mpi_put_att(ncid, vdimids(1),'long_name','distance in x')
    ierr = nf90mpi_put_att(ncid, vdimids(1),'units','m')
    ierr = nf90mpi_def_var(ncid, 'z',NF90_FLOAT,dimids2DZ(2),vdimids(2))
    ierr = nf90mpi_put_att(ncid, vdimids(2),'long_name','height')
    ierr = nf90mpi_put_att(ncid, vdimids(2),'units','m')
  end if
  ierr = nf90mpi_def_var(ncid, 'time',NF90_DOUBLE,dimids2DZ(3),vdimids(3))
  ierr = nf90mpi_put_att(ncid, vdimids(3),'long_name','time')
  if(timeUTsec.gt.0._8) then
   ierr = nf90mpi_put_att(ncid, vdimids(3),'units', &
                        'seconds since 1900-01-01 00:00:00.0')
   ierr = nf90mpi_put_att(ncid, vdimids(3),'calendar','gregorian')
   ierr = nf90mpi_def_var(ncid, 'day',NF90_DOUBLE,dimids2DZ(3),dayid)
   ierr = nf90mpi_put_att(ncid, dayid,'units','day')
   ierr = nf90mpi_put_att(ncid, dayid,'long_name','day')
  else
   ierr = nf90mpi_put_att(ncid, vdimids(3),'units','day')
  end if
  ierr = nf90mpi_def_var(ncid, 'dz',NF90_FLOAT,dimids2DZ(2),dzid)
  ierr = nf90mpi_put_att(ncid, dzid,'long_name','layer thickness')
  ierr = nf90mpi_put_att(ncid, dzid,'units','m')
  ierr = nf90mpi_def_var(ncid, 'rho',NF90_FLOAT,dimids2DZ(2),rhoid)
  ierr = nf90mpi_put_att(ncid, rhoid,'long_name','air density')
  ierr = nf90mpi_put_att(ncid, rhoid,'units','kg/m3')
  if(dolatlon) then
   ierr = nf90mpi_def_var(ncid, 'z',NF90_FLOAT,dimids2DZ(2),pid)
   ierr = nf90mpi_put_att(ncid, pid,'long_name','level height')
   ierr = nf90mpi_put_att(ncid, pid,'units','m')
  else
   ierr = nf90mpi_def_var(ncid, 'p',NF90_FLOAT,dimids2DZ(2),pid)
   ierr = nf90mpi_put_att(ncid, pid,'long_name','level pressure')
   ierr = nf90mpi_put_att(ncid, pid,'units','mb')
  end if

  ierr = nf90mpi_enddef(ncid)

  ierr = nf90mpi_begin_indep_data(ncid)
  if(dolatlon) then
      ierr = nf90mpi_put_var(ncid,vdimids(1),lat_gl(1:ny_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),pres(1:nzm))
      ierr = nf90mpi_put_var(ncid,pid,z(1:nzm))
  else
      do j=1,ny_gl
        y(j) = (j-1)*dy
      end do
      ierr = nf90mpi_put_var(ncid,vdimids(1),y(1:ny_gl))
      ierr = nf90mpi_put_var(ncid,vdimids(2),z(1:nzm))
      ierr = nf90mpi_put_var(ncid,pid,pres(1:nzm))
  end if
  if(timeUTsec.gt.0._8) then
    ierr = nf90mpi_put_var(ncid,vdimids(3),timeUTsec)
    ierr = nf90mpi_put_var(ncid,dayid,day)
  else
    ierr = nf90mpi_put_var(ncid,vdimids(3),day)
  end if
  dzz(:) = dz*adz(:)
  ierr = nf90mpi_put_var(ncid,dzid,dzz(1:nzm))
  ierr = nf90mpi_put_var(ncid,rhoid,rho(1:nzm))

  ierr = nf90mpi_end_indep_data(ncid)

  call task_barrier()

END SUBROUTINE open_file2DZ_pnetcdf

!======================================================================================================

SUBROUTINE close_file_pnetcdf()
 integer ierr
 call task_barrier()
 ierr = nf90mpi_close(ncid)
 if(masterproc) print*,'closed netcdf file'
END SUBROUTINE close_file_pnetcdf

!======================================================================================================

SUBROUTINE add_var3D_pnetcdf(f,name,long_name,units,gtype)
 real(4) f(:,:,:)
 character(*) name,long_name,units
 integer gtype
 integer dimids3D(4) ! Dimension IDs 
 integer ierr, varid, it, jt
 integer(kind=MPI_OFFSET_KIND) start(4), count(4)

 ierr = nf90mpi_redef(ncid)
 select case (gtype)
  case(1)
    dimids3D = dimids3DS
  case(2)  
    dimids3D = dimids3DS
    dimids3D(3) = dimids3DW
  case(3) 
    dimids3D = dimids3DS
    dimids3D(2) = dimids3DV
  case(4) 
    dimids3D = dimids3DS
    dimids3D(1) = dimids3DU
  case default
    if(masterproc) print*,'wrong gtype for netcdf varibale ',trim(name)
 end select
 ierr = nf90mpi_def_var(ncid,trim(name),NF90_FLOAT,dimids3D,varid)
 ierr = nf90mpi_put_att(ncid, varid,'long_name',trim(long_name)) 
 ierr = nf90mpi_put_att(ncid, varid,'units',trim(units)) 
 ierr = nf90mpi_enddef(ncid)
 call task_rank_to_index(rank,it,jt)
 start(1) = it+1 
 start(2) = jt+1 
 start(3) = 1 
 start(4) = 1 
 count(1) = nx
 count(2) = ny
 count(3) = nzm
 count(4) = 1
 ierr = nf90mpi_put_var_all(ncid,varid,values=f,start=start,count=count)
 call task_barrier()
 
END SUBROUTINE add_var3D_pnetcdf

!======================================================================================================

SUBROUTINE add_var2D_pnetcdf(f,name,long_name,units)
 real(4) f(:,:)
 character(*) name,long_name,units
 integer ierr, varid, it, jt
 integer(kind=MPI_OFFSET_KIND) start(3), count(3)

 ierr = nf90mpi_redef(ncid)
 ierr = nf90mpi_def_var(ncid,trim(name),NF90_FLOAT,dimids2D,varid)
 ierr = nf90mpi_put_att(ncid, varid,'long_name',trim(long_name))
 ierr = nf90mpi_put_att(ncid, varid,'units',trim(units))
 ierr = nf90mpi_enddef(ncid)
 call task_rank_to_index(rank,it,jt)
 start(1) = it+1
 start(2) = jt+1
 start(3) = 1
 count(1) = nx
 count(2) = ny
 count(3) = 1
 ierr = nf90mpi_put_var_all(ncid,varid,values=f,start=start,count=count)
 call task_barrier()
 if(masterproc) print*,'wrote field: ',trim(name)

END SUBROUTINE add_var2D_pnetcdf

!======================================================================================================

SUBROUTINE add_var2DZ_pnetcdf(f,name,long_name,units)
 real(4) f(:,:)
 character(*) name,long_name,units
 integer ierr, varid, it, jt
 integer(kind=MPI_OFFSET_KIND) start(3), count(3)

 ierr = nf90mpi_redef(ncid)
 ierr = nf90mpi_def_var(ncid,trim(name),NF90_FLOAT,dimids2DZ,varid)
 ierr = nf90mpi_put_att(ncid, varid,'long_name',trim(long_name))
 ierr = nf90mpi_put_att(ncid, varid,'units',trim(units))
 ierr = nf90mpi_enddef(ncid)
 call task_rank_to_index(rank,it,jt)
 start(1) = jt+1 
 start(2) = 1 
 start(3) = 1  
 count(1) = ny
 count(2) = nzm
 count(3) = 1
 ierr = nf90mpi_put_var_all(ncid,varid,values=f,start=start,count=count)
 call task_barrier()
 
END SUBROUTINE add_var2DZ_pnetcdf


end module pnetcdf_stuff

!-----------------------------------------------------------

#else

module pnetcdf_stuff

logical:: dopnetcdf = .false.

contains

! dummy section for the case when pnetcdf ;ibrary is not used

SUBROUTINE open_file3D_pnetcdf(filename)
character(*) filename
end subroutine open_file3D_pnetcdf

SUBROUTINE open_file2D_pnetcdf(filename)
character(*) filename
end subroutine open_file2D_pnetcdf

SUBROUTINE open_file2DZ_pnetcdf(filename)
character(*) filename
end subroutine open_file2DZ_pnetcdf

SUBROUTINE close_file_pnetcdf()
end subroutine close_file_pnetcdf

SUBROUTINE add_var3D_pnetcdf(f,name,long_name,units,gtype)
real(4) f(:,:,:)
character(*) name,long_name,units
integer gtype
end subroutine add_var3D_pnetcdf

SUBROUTINE add_var2D_pnetcdf(f,name,long_name,units)
real(4) f(:,:)
character(*) name,long_name,units
end subroutine add_var2D_pnetcdf

SUBROUTINE add_var2DZ_pnetcdf(f,name,long_name,units)
real(4) f(:,:)
character(*) name,long_name,units
end subroutine add_var2DZ_pnetcdf

end module pnetcdf_stuff

#endif
