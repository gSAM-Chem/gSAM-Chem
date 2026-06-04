c
c (C) 2000 Marat Khairoutdinov
c

	implicit none
	include 'netcdf.inc'

c---------------------------------------------------------------
c variables:

	character(200) filename
	character(80) long_name
	character(10) units
	character(8)  name
	character(1)  blank
	character(3) nm
	character(3) comp

	real(4), allocatable :: byte4(:)
	real(4), allocatable :: fld(:,:)
	real(4) fmax,fmin,day0
	real(4) dx,dy,z(500),p(500),x(100000)
	real(8) y(100000),time, time_old
	logical condition
	integer(4) nsubs,nsubsx,nsubsy,nx,ny,nz,nfields,nstep
	integer(4) i,j,k,k1,k2,n,i0,j0,nx_gl,ny_gl,length,ifields

	integer(4) vdimids(3),start(3),count(3),ndimids
	integer(4) ncid,err,yid,zid,timeid,ntime,varrr
c External functions:

	integer(4) iargc,strlen1
	external iargc,strlen1

	real fldmin, fldmax

c---------------------------------------------------------------
c---------------------------------------------------------------
c
c Read the file-name from the comman line:
c
	i=COMMAND_ARGUMENT_COUNT()
	if(i.eq.0) then
	  print*,'no input-file name is specified.'
	  print*,'Format: 2Dcom2nc input.com2DZ'
	  stop
	end if
	call getarg(1,filename)

c---------------------------------------------------------------
c Read files; merge data from different subdomains;
c save as a netcdf file.
c
	open(1,file=filename,status='old',form='unformatted',action='read')


	ntime=1
        time_old = -1.

      do while(.true.) ! infinit loop 

c
c The input file's extension:

        condition=ntime.eq.1

        if(condition) then
          do i=1,190
            if(filename(i:i+3).eq.'.2DZ') then
	      EXIT
	    else if(i.eq.190) then
	      print*,'wrong file name extension!'
	      stop
            endif
          end do
        end if

	read(1,end=3333,err=3333) nstep
        read(1) comp
        print*,'component:',comp
	read(1) nx,ny,nz,nsubs,nsubsx,nsubsy,nfields
	ny_gl=ny*nsubsy
	read(1) dy,y(1:ny_gl),time,z(1:nz),p(1:nz)
        print*,'time=',time
        print*,'ny_gl=',y(ny_gl)
        print*,'lat:',y(1:ny_gl)
        print*,'z:',z(1:nz)
        print*,'p:',p(1:nz)

c
c The output filename:

        if(condition) then
          do i=1,190
            if(filename(i:i+3).eq.'.2DZ') then
              write(filename(i+4:i+10),'(a7)') '_'//comp//'.nc'
              print*,filename
              EXIT
            endif
          end do
        end if


	print*,'nx,ny,nz,nsubs,nsubsx,nsubsy,nfields:'
	print*,nx,ny,nz,nsubs,nsubsx,nsubsy,nfields
	
	if(condition) then
	  day0=time
          allocate(byte4(nz*ny))
	  allocate(fld(ny_gl,nz))
        end if
c
c Initialize netcdf stuff, define variables,etc.
c
       if(condition) then

	err = NF_CREATE(filename, NF_CLOBBER, ncid)
	err = NF_REDEF(ncid)

	err = NF_DEF_DIM(ncid, 'lev', nz, zid)
	err = NF_DEF_DIM(ncid, 'lat', ny_gl, yid)
	err = NF_DEF_DIM(ncid, 'time', NF_UNLIMITED, timeid)

        err = NF_DEF_VAR(ncid, 'lev', NF_FLOAT, 1, zid, varrr)
	err = NF_PUT_ATT_TEXT(ncid,varrr,'units',1,'mb')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',8,'pressure')
         err = NF_DEF_VAR(ncid, 'lat', NF_DOUBLE, 1, yid, varrr)
	 err = NF_PUT_ATT_TEXT(ncid,varrr,'units',13,'degrees_north')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',8,'latitude')
        err = NF_DEF_VAR(ncid, 'time', NF_DOUBLE, 1, timeid, varrr)
	err = NF_PUT_ATT_TEXT(ncid,varrr,'units',3,'day')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',4,'time')
        err = NF_DEF_VAR(ncid, 'z', NF_FLOAT, 1, zid,varrr)
        err = NF_PUT_ATT_TEXT(ncid,varrr,'units',1,'m')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',6,'height')

	err = NF_ENDDEF(ncid)

	err = NF_INQ_VARID(ncid,'lev',varrr)
	err = NF_PUT_VAR_REAL(ncid, varrr, p(nz:1:-1))
	err = NF_INQ_VARID(ncid,'lat',varrr)
	err = NF_PUT_VAR_DOUBLE(ncid, varrr, y)
        err = NF_INQ_VARID(ncid,'z',varrr)
        err = NF_PUT_VAR_REAL(ncid, varrr, z(nz:1:-1))

	end if ! condition

	 ndimids=3
	 vdimids(1) = yid
	 vdimids(2) = zid
	 vdimids(3) = timeid
	 start(1) = 1
         start(2) = 1
         start(3) = ntime 
         count(1) = ny_gl
         count(2) = nz
         count(3) = 1 

	

	do ifields=1,nfields
	
         read(1) name,blank,long_name,blank,units
         print*,trim(name),trim(long_name),trim(units)
         read(1)

	  do n=0,nsubs-1,nsubsx
	  
            read(1) (byte4(k),k=1,ny*nz)
            j0 = n/nsubsx
            j0 = j0 * (ny_gl/nsubsy)
            length=0
            do k=1,nz
             do j=1+j0,ny+j0
                length=length+1
                fld(j,k)=byte4(length)
             end do
            end do

	  end do ! n

	 if(condition) then
	  err = NF_REDEF(ncid)
          err = NF_DEF_VAR(ncid,name,NF_FLOAT,
     &                           ndimids,vdimids,varrr)
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',
     &		strlen1(long_name),long_name(1:strlen1(long_name)))
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'units',
     &		strlen1(units),units(1:strlen1(units)))
	  err = NF_ENDDEF(ncid)
	 end if
	  
	
         if(time.gt.time_old) then
	   err = NF_INQ_VARID(ncid,name,varrr)
           err = NF_PUT_VARA_REAL(ncid,varrr,start,count,fld(:,nz:1:-1))
         end if

	end do

        if(time.gt.time_old) then
	   err = NF_INQ_VARID(ncid,'time',varrr)
           err = NF_PUT_VAR1_DOUBLE(ncid,varrr,ntime,time)
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
