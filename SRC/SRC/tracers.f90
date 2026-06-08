module tracers


! This module serves as a template for adding tracer transport in the model. The tracers can be 
! chemical tracers, or bin microphysics drop/ice categories, etc. 
! The number of tracers is set by the parameter ntracers which is set in domain.f90.
! Also, the logical flag dotracers should be set to .true. in namelist (default is .false.).
! The model will transport the tracers around automatically (advection and SGS diffusion).
! The user must supply the initialization in the subroutine tracers_init() in this module.
! By default, the surface flux of all tracers is zero. Nonzero values can be set in tracers_flux().
! The local sinks/sources of tracers should be supplied in tracers_physics().



 use grid
 implicit none

 real tracer  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers) 
 real fluxbtr (nx, ny, 0:ntracers) ! surface flux of tracers
 real fluxttr (nx, ny, 0:ntracers) ! top boundary flux of tracers
 real trwle(nz,0:ntracers)  ! resolved vertical flux 
 real trwsb(nz,0:ntracers)  ! SGS vertical flux 
 real tradv(nz,0:ntracers)  ! tendency due to vertical advection 
 real trdiff(nz,0:ntracers)  ! tendency due to vertical diffusion 
 real trphys(nz,0:ntracers)  ! tendency due to physics 
 integer flag_tracer3Dout(0:ntracers) 
 integer flag_tracer2DZout(0:ntracers)
 real tr_xy(nx,ny,0:ntracers) ! tracer mass path
 real tr_ac_xy(nx,ny,0:ntracers) ! accumulated tracers on the surface
 real tr_acs_xy(nx,ny,0:ntracers) ! integrated tracers' concentr.  right near the surface
 
 logical:: initialized_list = .false. ! flag to read source file
 integer, parameter :: nlist_max = 100 ! maximum allocated list size
 integer nlist           ! actual number of lines in the list
 real lons(nlist_max), lats(nlist_max) ! lon and lat of sources
 real power(nlist_max) ! power of the sources
 integer i_list(nlist_max), j_list(nlist_max) ! local i,j indexes
 real p_list(nlist_max)  ! emission rate (arbitrary units)

 character *4 tracername(0:ntracers)
 character *10 tracerunits(0:ntracers)

CONTAINS

 subroutine tracers_init()

  integer k,ntr
  character *2 ntrchar
  integer, external :: lenstr

 fluxbtr = 0.
 fluxttr = 0.

! Add your initialization code here. Default is to set to 0 in setdata.f90.

 if(nrestart.eq.0) then

  tracer = 0.

 end if

! Specify te tracers' default names:

   ! Default names are TRACER01, TRACER02, etc:

   flag_tracer3Dout(0) = 0
   flag_tracer2DZout(0) = 0
   do ntr = 1,ntracers
     write(ntrchar,'(i2)') ntr
     do k=1,3-lenstr(ntrchar)-1
        ntrchar(k:k)='0'
     end do
     tracername(ntr) = 'TR'//ntrchar(1:2)
     tracerunits(ntr) = '[TR]'
     flag_tracer3Dout(ntr) = 1
     flag_tracer2DZout(ntr) = 1
  end do

 end subroutine tracers_init

!------------------------------------------------------------------------------------

 subroutine tracers_flux()

 use terrain, only: k_terra
 use vars, only: raf
 use params, only: dotrsfcflux
 integer i,j
 real ramin

! Set surface and top fluxes of tracers. Default is 0 set in setdata.f90

 if(dotrsfcflux) then
  do j=1,ny
   do i=1,nx
     ramin = 10.*dtn/(z(k_terra(i,j))-zi(k_terra(i,j)))
     fluxbtr(i,j,1:ntracers) = -tracer(i,j,k_terra(i,j),1:ntracers)/max(ramin,raf(i,j))
   end do
  end do
 end if

 end subroutine tracers_flux


!------------------------------------------------------------------------------------
 
 subroutine tracers_physics()
 
 integer itr

 ! add here a call to a subroutine that does something to tracers besides advection and diffusion.
 ! The transport is done automatically. 

  trphys = 0. ! Default tendency due to physics. You code should compute this to output statistics.

 end subroutine tracers_physics



!------------------------------------------------------------------------------------


 subroutine tracers_source()

 use params, only: tracer_source_type
 use terrain, only: k_terra, terra
 use vars, only: rho
 integer i,j,it,jt,k,kkk,m
 character(40) ctmp

 call task_rank_to_index(rank,it,jt)
 select case (tracer_source_type)

  case(100)

! read sources from the file
  if(.not.initialized_list) then
    open(55,file='./CASES/'//trim(case)//'/nyc_power_plants_100mi.dat', status='old',form='formatted')
    read(55,*) 
    m=0
    do i=1,nlist_max
     read(55,*,end=333) ctmp,lats(i),lons(i),power(i)
     if(lons(i).lt.0.) lons(i) = lons(i)+360.
     m=i
    end do
    333 continue
    nlist = 0
    do k = 1,m
      do j=1,ny
        do i=1,nx
          if(lats(k).ge.latv_gl(j+jt).and.lats(k).le.latv_gl(j+jt+1).and. &
             lons(k).ge.lonu_gl(i+it).and.lons(k).le.lonu_gl(i+it+1)) then
               nlist = nlist+1
               i_list(nlist) = i 
               j_list(nlist) = j 
               p_list(nlist) = power(k) 
          end if
        end do
      end do 
    end do
    if(nlist.gt.0) then
     print*,'rank=',rank,'  nlist =',nlist
    end if
    initialized_list = .true.
  end if
  m = 1
  do k=1,nlist
    kkk = k_terra(i_list(k),j_list(k))+1
    tracer(i_list(k),j_list(k),kkk,m) = max(0.,tracer(i_list(k),j_list(k),kkk,m) +  p_list(k)*dtn)
  end do

  case(1)
 ! LI Northprt power plant for 486x486 600m grid
   if(nx_gl.ne.486.and.ny_gl.ne.486) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   kkk=1
   do k=1,nzm
    if(zi(k).gt.180.) then
      kkk = k
      exit
    end if
   end do
   m = 1
   do j=1,ny
    do i=1,nx
      if(i+it.eq.238.and.j+jt.eq.250.) tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
    end do
   end do

  case(2)
   ! LI Northprt power plant for 2916x2916 100m grid
   if(nx_gl.ne.2916.and.ny_gl.ne.2916) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   kkk=1
   do k=1,nzm
    if(zi(k).gt.180.) then
      kkk = k
      exit
    end if
   end do
   m = 1
   do j=1,ny
    do i=1,nx
      if(i+it.eq.1462.and.j+jt.eq.1497.) tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      if(i+it.eq.1462.and.j+jt.eq.1496.) tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      if(i+it.eq.1462.and.j+jt.eq.1495.) tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
    end do
   end do

  case(3)

   if(nx_gl.ne.2916.and.ny_gl.ne.2916) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   m = 1
   do j=1,ny
    do i=1,nx
   ! Indian Point  power plant for 2916x2916 100m grid
      if(i+it.eq.946.and.j+jt.eq.1865) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.758.and.j+jt.eq.313) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.2443.and.j+jt.eq.1945) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.1459.and.j+jt.eq.1492) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.841.and.j+jt.eq.1306) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if

    end do
   end do

  case(35)

   if(ntracers.ne.4.and.nx_gl.ne.2916.and.ny_gl.ne.2916) then
     if(masterproc) print*,'wrong point source for given grid.'
     if(masterproc) print*,'ntracers should be 4. Actual value is ',ntracers
     call task_abort()
   end if
   m=1
   do j=1,ny
    do i=1,nx
      if(i+it.eq.1143.and.j+jt.eq.1200) then  ! Rovenswood
          kkk = k_terra(i,j)+4
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  10000.*dtn)
      end if
      if(i+it.eq.1174.and.j+jt.eq.1226) then ! Astoria, Queens, stack heigh (assumed) 80m
          kkk = k_terra(i,j)+4
          tracer(i,j,kkk,m+1) = max(0.,tracer(i,j,kkk,m+1) +  10000.*dtn)
      end if
      if(i+it.eq.877.and.j+jt.eq.1060) then ! Linden, New Jersey, stack heigh (assumed) 80m
          kkk = k_terra(i,j)+4
          tracer(i,j,kkk,m+2) = max(0.,tracer(i,j,kkk,m+2) +  10000.*dtn)
      end if
      if(i+it.eq.1749.and.j+jt.eq.1364) then ! Northport, stack height 180m
          kkk = k_terra(i,j)+7
          tracer(i,j,kkk,m+3) = max(0.,tracer(i,j,kkk,m+3) +  10000.*dtn)
      end if

    end do
   end do

  case(36)

  ! arbitrary tracer sources:
   if(nx_gl.ne.6000.and.ny_gl.ne.6000) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   m = 1
   do j=1,ny
    do i=1,nx
   ! Indian Point  power plant for 2916x2916 100m grid
      if(i+it.eq.1798.and.j+jt.eq.1838) then  
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if
      if(i+it.eq.4933.and.j+jt.eq.2568) then 
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if
      if(i+it.eq.2038.and.j+jt.eq.3623) then 
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if
      if(i+it.eq.638.and.j+jt.eq.1973) then 
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if
      if(i+it.eq.2383.and.j+jt.eq.4858) then 
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if
      if(i+it.eq.589.and.j+jt.eq.3880) then 
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1000.*dtn)
      end if

    end do
   end do

  case(4)

   ! arbitrary tracers sources for 4608x4608 5m grid
   if(nx_gl.ne.4608.and.ny_gl.ne.4608) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   m = 1
   do j=1,ny
    do i=1,nx
      if(i+it.eq.1300.and.j+jt.eq.3700) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.2150.and.j+jt.eq.3225) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.1615.and.j+jt.eq.2260) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.1306.and.j+jt.eq.1390) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.3560.and.j+jt.eq.2863) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.2926.and.j+jt.eq.1048) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
    end do
   end do

  case(5)

   if(nx_gl.ne.2880.and.ny_gl.ne.2880) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   m = 1
   do j=1,ny
    do i=1,nx
   ! Indian Point  
      if(i+it.eq.928.and.j+jt.eq.1847) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.740.and.j+jt.eq.285) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.2425.and.j+jt.eq.1927) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.1441.and.j+jt.eq.1474) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if
      if(i+it.eq.823.and.j+jt.eq.1288) then
          kkk = k_terra(i,j)+1
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  1.e3*dtn)
      end if

    end do
   end do

  case(6)

   ! arbitrary tracers sources for 1296x1296 3m grid 
   if(nx_gl.ne.1296.and.ny_gl.ne.1296) then
     if(masterproc) print*,'wrong point source for given grid.'
     call task_abort()
   end if
   m = 1
   do j=1,ny
    do i=1,nx
!    sources for SC2:
!      if(i+it.eq.484.and.j+jt.eq.577) then  ! SW
!          kkk = k_terra(i,j)+1
!          print*,'release at:',lon_gl(i+it),lat_gl(j+jt),z(kkk)
!          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  4900./(dx*dy*(zi(kkk+1)-zi(kkk)))*dtn)  
!      end if
!      if(i+it.eq.649.and.j+jt.eq.649) then  ! CR
!          kkk = k_terra(i,j)+1
!          print*,'release at:',lon_gl(i+it),lat_gl(j+jt),z(kkk)
!          tracer(i,j,kkk,m+1) = max(0.,tracer(i,j,kkk,m+1) +  4020000./3600./(dx*dy*(zi(kkk+1)-zi(kkk)))*dtn) 
!      end if
   ! sources for BC2:
!      if(i+it.eq.685.and.j+jt.eq.557) then  ! SE
!          kkk = k_terra(i,j)+1
!          print*,'release at:',lon_gl(i+it),lat_gl(j+jt),z(kkk)
!          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) +  4750./1800./(dx*dy*(zi(kkk+1)-zi(kkk)))*dtn) 
!      end if
      if(i+it.eq.732.and.j+jt.eq.708) then  ! ENE
          kkk = k_terra(i,j)
          print*,'release at:',lon_gl(i+it),lat_gl(j+jt),z(kkk)
          tracer(i,j,kkk,m) = max(0.,tracer(i,j,kkk,m) + &
                 4.75/3600./(rho(kkk)*dx*dy*(zi(kkk+1)-zi(kkk)))*dtn*terra(i,j,kkk))  
      end if
      if(i+it.eq.649.and.j+jt.eq.649) then  ! CR
          kkk = k_terra(i,j)
          print*,'release at:',lon_gl(i+it),lat_gl(j+jt),z(kkk)
          tracer(i,j,kkk,m+1) = max(0.,tracer(i,j,kkk,m+1) +  &
                 3.97/3600./(rho(kkk)*dx*dy*(zi(kkk+1)-zi(kkk)))*dtn*terra(i,j,kkk)) 
      end if

    end do
   end do

 end select

 end subroutine tracers_source



!------------------------------------------------------------------------------------


 subroutine tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

! Initialize the list of tracers statistics variables written in statistics.f90

   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr


   do ntr=1,ntracers

     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))
     deflist(count) = trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))
     status(count) = 1
     average_type(count) = 0		
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLX'
     deflist(count) = 'Total flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLXS'
     deflist(count) = 'SGS flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'ADV'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical advection')
     unitlist(count) = trim(tracerunits(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'DIFF'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical SGS transport')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'PHYS'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to physics')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
   end do

 end subroutine tracers_hbuf_init

end module tracers
