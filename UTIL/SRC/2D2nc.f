c
c (C) 2000 Marat Khairoutdinov
c

	implicit none
	include 'netcdf.inc'

c---------------------------------------------------------------
c variables:

	character(200) filename,filename1
	character(80) long_name
	character(10) units
	character(8)  name
	character(1)  blank
	character(3) nm
	character(3) comp
        character(6) gridtype
        character(14) datechar

	integer(2), allocatable :: byte(:)
	real(4), allocatable :: byte4(:)
	real(4), allocatable :: fld(:,:)
        real(8), allocatable :: lat(:), lon(:), latv(:), lonu(:), wgt(:)
	real(4) fmax,fmin,day0
	real(4) dx,dy
	real(8) time_old, pi, time, mean,timesec, timeUTsec
	real(4), allocatable ::  x(:),y(:),xv(:),yv(:)
	logical condition
	integer(4) nsubs,nsubsx,nsubsy,nx,ny,nz,nfields,nstep
	integer(4) i,j,k,k1,k2,n,m,i0,j0,nx_gl,ny_gl,ifields,it,jt

	integer(4) vdimids(3),start(3),count(3),ndimids
	integer(4) ncid,err,yid,xid,dayid,timeid,ntime,varrr,oldmode
        logical(4) isbin, dolatlon
c External functions:

	integer(4) iargc,strlen1
	external iargc,strlen1

c---------------------------------------------------------------
c---------------------------------------------------------------
c
c Read the file-name from the comman line:
c
	i=COMMAND_ARGUMENT_COUNT()
	if(i.eq.0) then
	  print*,'no input-file name is specified.'
	  print*,'Format: 2D2nc input.2D (optional)latlon'
	  stop
	end if
	call getarg(1,filename)
        if(i.eq.2) then
          call getarg(2,gridtype)
          if(gridtype.ne.'latlon') then
            print*,'Optional argument not latlon.'
	    print*,'Format: 2D2nc input.2D latlon'
            stop
          end if
        end if
c---------------------------------------------------------------
c Read files; merge data from different subdomains;
c save as a netcdf file.
c
	open(1,file=filename,status='old',form='unformatted',action='read')


	ntime=1
        time_old = -1.

      do while(.true.) ! infinite loop 

c
c Check input file's extension:

        condition = ntime.eq.1
        if(condition) then
 	  filename1=filename
          do i=1,190
            if(filename1(i:i+2).eq.'.2D') then
	      EXIT
	    else if(i.eq.190) then
	      print*,'wrong file name extension!'
	      stop
            endif
            k = i
          end do
        end if

	read(1,end=3333,err=3333) nstep, dolatlon, timesec, datechar, timeUTsec
        print*,'NSTEP=',nstep
        print*,'DATE=',trim(datechar)
        print*,'timeUTsec=',timeUTsec
        print*,'timesec=',timesec
        if(gridtype.eq.'latlon') dolatlon = .true.
        read(1) comp
	read(1) nx,ny,nz,nsubs,nsubsx,nsubsy,nfields
	print*,'nx,ny,nz,nsubs,nsubsx,nsubsy,nfields:'
	print*,nx,ny,nz,nsubs,nsubsx,nsubsy,nfields
        if(.not.allocated(lat)) then
         allocate(lat(ny*nsubsy),lon(nx*nsubsx),latv(ny*nsubsy+1),lonu(nx*nsubsx),wgt(ny*nsubsy))
         allocate(x(nx*nsubsx),xv(nx*nsubsx),y(ny*nsubsy),yv(ny*nsubsy+1))
        end if
        read(1) dx,lat,latv,y,yv
        read(1) dy,lon,lonu
        print*,'lat:',minval(lat),maxval(lat)
        print*,'lon:',minval(lon),maxval(lon)
        print*,'latv:',minval(latv),maxval(latv)
        print*,'lonu:',minval(lonu),maxval(lonu)
        print*,'y:',minval(y),maxval(y)
        print*,'vv:',minval(yv),maxval(yv)
	read(1) time
	print*,nstep,time,dx,dy

        if(condition) then
          filename1=filename
          do i=1,190
            if(filename1(i:i+2).eq.'.2D') then
              write(filename1(i+3:i+9),'(a7)') '_'//comp//'.nc'
              print*,filename1
              EXIT
            endif
            k = i
          end do
          if(trim(filename).eq.trim(filename1)) then
               print*,'attempt to ovewrite binary file!'
               print*,'filename=',trim(filename)
               print*,'filename1=',trim(filename1)
               stop
          end if
        end if

	
	nx_gl=nx*nsubsx
	ny_gl=ny*nsubsy
	if(condition) then
	 print*,'nx_gl=',nx_gl, '    dx=',dx
	 print*,'ny_gl=',ny_gl, '    dy=',dy
	end if

        pi = acos(-1.d0)
        if(dolatlon) then
          mean = 0.
	  do j=1,ny_gl
           wgt(j) = (latv(j+1)-latv(j))*cos(pi/180.d0*lat(j))
           mean = mean + wgt(j)
	  end do
          wgt = wgt/mean*ny_gl
        else
	  do i=1,nx_gl
           x(i) = dx*(i-1)
	  end do
          xv(:) = x(:)-0.5*(x(2)-x(1))
          wgt = 1.d0
        end if

	if(ntime.eq.1) then
	  day0=time
          allocate(byte(nx_gl*ny_gl))
          allocate(byte4(nx_gl*ny_gl))
	  allocate(fld(nx_gl,ny_gl))
        end if
c
c Initialize netcdf stuff, define variables,etc.
c
       if(condition) then



        err = NF_CREATE(filename1,
c     &                  IOR(NF_CLOBBER, NF_NETCDF4), ncid)
     &          IOR(NF_CLOBBER, NF_64BIT_DATA), ncid)
c     &                  ior(nf_clobber,nf_64bit_offset), ncid)
c	err = NF_CREATE(filename1, NF_CLOBBER, ncid)
c     Disable fill (faster writes for huge vars)
       err = NF_SET_FILL(ncid, NF_NOFILL, oldmode)
c	err = NF_REDEF(ncid)

        if(dolatlon) then
	  err = NF_DEF_DIM(ncid, 'lon', nx_gl, xid)
	  if(ny_gl.ne.1)err = NF_DEF_DIM(ncid, 'lat', ny_gl, yid)
	  err = NF_DEF_DIM(ncid, 'time', NF_UNLIMITED, timeid)
          err = NF_DEF_VAR(ncid, 'lon', NF_DOUBLE, 1, xid, varrr)
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'units',12,'degrees_east')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',9,'longitude')
	  if(ny_gl.ne.1) then
           err = NF_DEF_VAR(ncid, 'lat', NF_DOUBLE, 1, yid, varrr)
	   err = NF_PUT_ATT_TEXT(ncid,varrr,'units',13,'degrees_north')
           err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',8,'latitude')
           err = NF_DEF_VAR(ncid, 'wgt', NF_DOUBLE, 1, yid, varrr)
           err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',17,'averaging wieghts')
	  endif
         else
          err = NF_DEF_DIM(ncid, 'x', nx_gl, xid)
          if(ny_gl.ne.1)err = NF_DEF_DIM(ncid, 'y', ny_gl, yid)
          err = NF_DEF_DIM(ncid, 'time', NF_UNLIMITED, timeid)
          err = NF_DEF_VAR(ncid, 'x', NF_FLOAT, 1, xid, varrr)
          err = NF_PUT_ATT_TEXT(ncid,varrr,'units',1,'m')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',1,'x')
          if(ny_gl.ne.1) then
           err = NF_DEF_VAR(ncid, 'y', NF_FLOAT, 1, yid, varrr)
           err = NF_PUT_ATT_TEXT(ncid,varrr,'units',1,'m')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',1,'m')
          endif
         end if
        err = NF_DEF_VAR(ncid, 'time', NF_DOUBLE, 1, timeid, varrr)
        if(timeUTsec.gt.0._8) then
          err = NF_PUT_ATT_TEXT(ncid,varrr,'units',35,'seconds since 1900-01-01 00:00:00.0')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',4,'time')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'calendar',9,'gregorian')
          err = NF_DEF_VAR(ncid,'day',NF_DOUBLE, 1, timeid, varrr)
          err = NF_PUT_ATT_TEXT(ncid,varrr,'units',3,'day')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',3,'day')
        else
          err = NF_PUT_ATT_TEXT(ncid,varrr,'units',3,'day')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',4,'time')
        end if
        err = NF_DEF_VAR(ncid,'timesec',NF_DOUBLE, 1, timeid, varrr)
        err = NF_PUT_ATT_TEXT(ncid,varrr,'units',1,'s')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',15,'run time in sec')


	err = NF_ENDDEF(ncid)

        if(dolatlon) then
	  err = NF_INQ_VARID(ncid,'lon',varrr)
	  err = NF_PUT_VAR_DOUBLE(ncid, varrr, lon)
	  if(ny_gl.ne.1) then
	  err = NF_INQ_VARID(ncid,'lat',varrr)
	  err = NF_PUT_VAR_DOUBLE(ncid, varrr, lat)
	  err = NF_INQ_VARID(ncid,'wgt',varrr)
	  err = NF_PUT_VAR_DOUBLE(ncid, varrr, wgt)
	  endif
        else
          err = NF_INQ_VARID(ncid,'x',varrr)
          err = NF_PUT_VAR_REAL(ncid, varrr, x)
          if(ny_gl.ne.1) then
          err = NF_INQ_VARID(ncid,'y',varrr)
          err = NF_PUT_VAR_REAL(ncid, varrr, y)
          endif
        end if


	end if ! condition


	if(ny_gl.ne.1) then
	 ndimids=3
	 vdimids(1) = xid
	 vdimids(2) = yid
	 vdimids(3) = timeid
	 start(1) = 1
         start(2) = 1
         start(3) = ntime 
	 count(1) = nx_gl
         count(2) = ny_gl
         count(3) = 1 
	else
	 ndimids=2
	 vdimids(1) = xid
	 vdimids(2) = timeid
	 start(1) = 1
         start(2) = ntime 
	 count(1) = nx_gl
         count(2) = 1 
	endif

!       it is faster with large netcdf file to first define all variables
!       before starting putting their values into the files

        if(condition) then

        do ifields=1,nfields

         read(1) name,blank,long_name,blank,units

         err = NF_REDEF(ncid)
         err = NF_DEF_VAR(ncid,name,NF_FLOAT,
     &                           ndimids,vdimids,varrr)
         err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',
     &          strlen1(trim(long_name)),trim(long_name))
         err = NF_PUT_ATT_TEXT(ncid,varrr,'units',
     &          strlen1(trim(units)),trim(units))
         err = NF_ENDDEF(ncid)
!         print*,'Done defining name ',name
         read(1) isbin
         if(isbin) then
           read(1) 
         else
           read(1) 
           read(1) 
         end if
        end do ! ifields

!       rewind file:
        rewind(1)
        read(1)
        read(1)
        read(1)
        read(1)
        read(1)
        read(1)

        end if ! condition
        print*,'                Max             Min           Mean '

	do ifields=1,nfields
	
         read(1) name,blank,long_name,blank,units
         mean = 0.
         read(1) isbin
         if(isbin) then
     	   read(1) byte4(1:nx_gl*ny_gl)
           m=0
           do jt=0,(nsubsy-1)*ny,ny
            do it=0,(nsubsx-1)*nx,nx
             do j=1,ny
              do i=1,nx
               m=m+1
	       fld(it+i,jt+j)= byte4(m)
              end do
             end do
            end do
           end do
         else
           read(1) fmax,fmin
           read(1) byte(1:nx_gl*ny_gl)
           m=0
           do jt=0,(nsubsy-1)*ny,ny
            do it=0,(nsubsx-1)*nx,nx
             do j=1,ny
              do i=1,nx
               m=m+1
               fld(it+i,jt+j)= fmin+real(byte(m)+32000)*(fmax-fmin)/64000.
              end do
             end do
            end do
           end do
         end if
         do j=1,ny_gl
          do i=1,nx_gl
           mean=mean+fld(i,j)*wgt(j)
          end do
         end do
         mean = mean/float(nx_gl*ny_gl)
         print*,name,'   ',minval(fld),maxval(fld),real(mean,4)
	
         if(time.gt.time_old) then
	   err = NF_INQ_VARID(ncid,name,varrr)
           err = NF_PUT_VARA_REAL(ncid,varrr,start,count,fld)
         end if

	end do ! ifields

        if(time.gt.time_old) then
           print*,'time=',time,' time_old=',time_old
	   err = NF_INQ_VARID(ncid,'time',varrr)
           if(timeUTsec.gt.0._8) then
            err = NF_PUT_VAR1_DOUBLE(ncid,varrr,ntime,timeUTsec)
	    err = NF_INQ_VARID(ncid,'day',varrr)
            err = NF_PUT_VAR1_DOUBLE(ncid,varrr,ntime,time)
           else
            err = NF_PUT_VAR1_DOUBLE(ncid,varrr,ntime,time)
           end if
           err = NF_INQ_VARID(ncid,'timesec',varrr)
           err = NF_PUT_VAR1_DOUBLE(ncid,varrr,ntime,timesec)
           time_old = time
           ntime=ntime+1
        end if    

      end do ! while


 3333	continue

	err = NF_CLOSE(ncid)

	end
	




	integer function strlen1(str)
	character*(*) str
	strlen1=len(str)
	do i=len(str),1,-1
	  if(str(i:i).ne.' ') then
	    strlen1=i
	    return
	  endif 
	end do
        return
	end
