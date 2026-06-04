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

	integer(2), allocatable :: byte_large(:,:,:),byte(:,:)
	real(4), allocatable :: byte4_large(:,:,:),byte4(:,:)
	real(4), allocatable :: fld(:,:)
        real(4), allocatable ::  lats(:,:),lons(:,:)
	real(4) fmax,fmin,day0,mean
	real(4) dx,dy,time_old, pi,time
	logical condition
	integer(4) nsubs,nsubsx,nsubsy,nx,ny,nz,nfields,nstep,ny_max
	integer(4) i,j,k,k1,k2,n,m,mn,i0,j0,nx_gl,ny_gl,ifields
        integer(4) nys(10000)
        real(4) x1(100000),y1(100000),x(100000),y(100000)

	integer(4) vdimids(3),start(3),count(3),ndimids
	integer(4) ncid,err,yid,xid,timeid,ntime,varrr
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
	  print*,'Format: 2D2nc input.2D'
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

	read(1,end=3333,err=3333) nstep, dolatlon
        read(1) comp
	read(1) nx,ny_max,nz,nsubs,nsubsx,nsubsy,nfields
	print*,'nx,ny_max,nz,nsubs,nsubsx,nsubsy,nfields:'
	print*,nx,ny_max,nz,nsubs,nsubsx,nsubsy,nfields
        read(1) dx
        read(1) dy
	read(1) time
	print*,nstep,time,dx,dy

        j=0
        i=0
        read(1) nys(1:nsubs)
        allocate(lats(ny_max,nsubs),lons(nx,nsubs))
        read(1) lats,lons
        print*,minval(lats),maxval(lats)
        do n=1,nsubs
          ny = nys(n)
          y1(1:ny) = lats(1:ny,n)
          if(mod(n-1,nsubsx).eq.0) then
            y(j+1:j+ny) = y1(1:ny)
            j = j+ny
          end if
        !  if(n.le.nsubsx-1) then
        !    x(i+1:1+nx) = x1(1:nx)
        !    i=i+nx
        !  end if
        end do
        nx_gl = nx*nsubsx
        do i=1,nsubsx*nx
         x(i) = lons(1,1)+(i-1)*(lons(2,1)-lons(1,1))
        end do
        print*,'nx_gl=',nx_gl
        ny_gl = j
        print*,'ny_gl=',ny_gl
        print*,'lat:',y(1:ny_gl)
        print*,'lon:',x(1:nx_gl)
        print*,'lat:',minval(y(1:ny_gl)),maxval(y(1:ny_gl))
        print*,'lon:',minval(x(1:nx_gl)),maxval(x(1:nx_gl))

        if(condition) then
          filename1=filename
          do i=1,190
            if(filename1(i:i+3).eq.'.2DC') then
              write(filename1(i+4:i+10),'(a7)') '_'//comp//'.nc'
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

	if(ntime.eq.1) then
	  day0=time
          allocate(byte_large(nx,ny_max,nsubs))
          allocate(byte4_large(nx,ny_max,nsubs))
          allocate(byte(nx,ny_max))
          allocate(byte4(nx,ny_max))
	  allocate(fld(nx_gl,ny_gl))
        end if
c
c Initialize netcdf stuff, define variables,etc.
c
       if(condition) then



        err = NF_CREATE(filename1,
     &                  ior(nf_clobber,nf_64bit_offset), ncid)
!	err = NF_CREATE(filename1, NF_CLOBBER, ncid)
	err = NF_REDEF(ncid)

        if(dolatlon) then
	  err = NF_DEF_DIM(ncid, 'lon', nx_gl, xid)
	  if(ny_gl.ne.1)err = NF_DEF_DIM(ncid, 'lat', ny_gl, yid)
	  err = NF_DEF_DIM(ncid, 'time', NF_UNLIMITED, timeid)
          err = NF_DEF_VAR(ncid, 'lon', NF_FLOAT, 1, xid, varrr)
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'units',12,'degrees_east')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',9,'longitude')
	  if(ny_gl.ne.1) then
           err = NF_DEF_VAR(ncid, 'lat', NF_FLOAT, 1, yid, varrr)
	   err = NF_PUT_ATT_TEXT(ncid,varrr,'units',13,'degrees_north')
          err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',8,'latitude')
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
        err = NF_DEF_VAR(ncid, 'time', NF_FLOAT, 1, timeid, varrr)
	err = NF_PUT_ATT_TEXT(ncid,varrr,'units',3,'day')
        err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',4,'time')

	err = NF_ENDDEF(ncid)

        if(dolatlon) then
	  err = NF_INQ_VARID(ncid,'lon',varrr)
	  err = NF_PUT_VAR_REAL(ncid, varrr, x)
	  if(ny_gl.ne.1) then
	  err = NF_INQ_VARID(ncid,'lat',varrr)
	  err = NF_PUT_VAR_REAL(ncid, varrr, y)
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

	

	err = NF_REDEF(ncid)
	do ifields=1,nfields
	
         read(1) name,blank,long_name,blank,units
         read(1) isbin
         if(.not.isbin) read(1) 
     	 read(1) 

	 if(condition) then
          err = NF_DEF_VAR(ncid,name,NF_FLOAT,
     &                           ndimids,vdimids,varrr)
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'long_name',
     &		strlen1(long_name),long_name(1:strlen1(long_name)))
	  err = NF_PUT_ATT_TEXT(ncid,varrr,'units',
     &		strlen1(units),units(1:strlen1(units)))
	  end if
         end do 
         print*,'Done defining netcdf'
	 err = NF_ENDDEF(ncid)
         rewind(1)
         read(1)
         read(1)
         read(1)
         read(1)
         read(1)
         read(1)
         read(1)
         read(1)

	do ifields=1,nfields
         mean = 0.
         read(1) name,blank,long_name,blank,units
         read(1) isbin
         print*,isbin
         if(isbin) then
          read(1) byte4_large(:,:,:)
         else
          read(1) fmax,fmin
          read(1) byte_large(:,:,:)
         end if
         mn=0
         j0=0
         do n=0,nsubsy-1
          do m=0,nsubsx-1
             mn = mn+1
             ny = nys(mn)
             i0 = m*(nx_gl/nsubsx)
             if(isbin) then
               byte4(:,:) = byte4_large(:,:,mn)
               do j=1,ny
                 do i=1,nx
                   fld(i0+i,j0+j)=byte4(i,j)
                 end do
               end do
             else
                byte(:,:) = byte_large(:,:,mn)
                 do j=1,ny
                  do i=1,nx
                    fld(i0+i,j0+j)=fmin+(byte(i,j)+32000.)*(fmax-fmin)/64000.
                  end do
                 end do
             end if
           end do
           j0=j0+ny
         end do ! n


         print*,name,strlen1(long_name),'   ',minval(fld),maxval(fld)
	
         if(time.gt.time_old) then
	   err = NF_INQ_VARID(ncid,name,varrr)
           err = NF_PUT_VARA_REAL(ncid,varrr,start,count,fld)
         end if

	end do

        if(time.gt.time_old) then
	   err = NF_INQ_VARID(ncid,'time',varrr)
           err = NF_PUT_VAR1_REAL(ncid,varrr,ntime,time)
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
