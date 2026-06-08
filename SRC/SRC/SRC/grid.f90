! grid.f90
! Coded by Marat Khairoutdinov
!-----------------------------------------------------------------------------------------------
! contains many model parameters and grid-related arrays such as horizontally invariant profiles
! and horizontal grid arrays.
! some parameters are set by the model itself, others by the PARAMETER namelist
!
!-----------------------------------------------------------------------------------------------

module grid

use domain
use advection, only: NADV, NADVS
use mpi_stuff, only: comm, dompi, rank, masterproc

implicit none

character(6), parameter :: version = '1.8.7'
character(8), parameter :: version_date = 'Feb 2026'
        
!----------------------------------------------------------------------------------------
integer, parameter :: nx = nx_gl/nsubdomains_x ! subdomain's number of grid points in x 
integer, parameter :: ny = ny_gl/nsubdomains_y ! subdomain's number of grid points in y
integer, parameter :: nz = nz_gl+1  ! number of intefaces (where w is defined)
integer, parameter :: nzm = nz-1    ! number of scalar levels or mid-levels

!----------------------------------------------------------------------------------------
integer, parameter :: nsubdomains = nsubdomains_x * nsubdomains_y ! number of subdomains

logical, parameter :: RUN3D = ny_gl.gt.1  ! flag for 3D run (.true.)
logical, parameter :: RUN2D = .not.RUN3D  ! flag for 2D run (.true.)

!----------------------------------------------------------------------------------------
! some axiliary indexes:

integer, parameter :: nxp1 = nx + 1
integer, parameter :: nyp1 = ny + 1 * YES3D
integer, parameter :: nxp2 = nx + 2
integer, parameter :: nyp2 = ny + 2 * YES3D
integer, parameter :: nxp3 = nx + 3
integer, parameter :: nyp3 = ny + 3 * YES3D
integer, parameter :: nxp4 = nx + 4
integer, parameter :: nyp4 = ny + 4 * YES3D

!----------------------------------------------------------------------------------------
! array start and end indexes in horizontal (u,v,w - wind, s - scalar)
integer, parameter :: dimx1_u = -1                        !!-1        -1        -1        -1
integer, parameter :: dimx2_u = nxp3                      !!nxp3      nxp3      nxp3      nxp3
integer, parameter :: dimy1_u = 1-(2+NADV)*YES3D          !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
integer, parameter :: dimy2_u = nyp2+NADV*YES3D           !!nyp5      nyp4      nyp3      nyp2
integer, parameter :: dimx1_v = -1-NADV                   !!-4        -3        -2        -1
integer, parameter :: dimx2_v = nxp2+NADV                 !!nxp5      nxp4      nxp3      nxp2
integer, parameter :: dimy1_v = 1-2*YES3D                 !!1-2*YES3D 1-2*YES3D 1-2*YES3D 1-2*YES3D
integer, parameter :: dimy2_v = nyp3                      !!nyp3       nyp3      nyp3      nyp3
integer, parameter :: dimx1_w = -1-NADV                   !!-4        -3        -2        -1
integer, parameter :: dimx2_w = nxp2+NADV                 !!nxp5      nxp4      nxp3      nxp2
integer, parameter :: dimy1_w = 1-(2+NADV)*YES3D          !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
integer, parameter :: dimy2_w = nyp2+NADV*YES3D           !!nyp5      nyp4      nyp3      nyp2
integer, parameter :: dimx1_s = -2-NADVS                  !!-4        -3        -2        -2
integer, parameter :: dimx2_s = nxp3+NADVS                !!nxp5      nxp4      nxp3      nxp3
integer, parameter :: dimy1_s = 1-(3+NADVS)*YES3D         !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-3*YES3D
integer, parameter :: dimy2_s = nyp3+NADVS*YES3D          !!nyp5      nyp4      nyp3      nyp3

!----------------------------------------------------------------------------------------
! Vertical grid parameters:
real z(nzm)      ! height of the pressure levels above surface,m
real pres(nzm)  ! pressure,mb at scalar levels
real zi(nz)     ! height of the interface levels
real presi(nz)  ! pressure,mb at interface levels
real adz(nzm)   ! ratio of the thickness of scalar levels to dz 
real adzw(nz)	! ratio of the thinckness of w levels to dz
real pres0      ! Reference surface pressure, Pa

!----------------------------------------------------------------------------------------
! Horizontal grid
real(8):: ady(dimy1_s:dimy2_s) = 1.  ! ratio of gridstep in y to dy
real(8):: adyv(dimy1_v:dimy2_v) = 1. ! ratio of gridstep in y to dy (for v)
real(8):: y_gl(1:ny_gl) = 0.         ! global grid's y .
real(8):: yv_gl(1:ny_gl+1) = 0.      ! global grid's y at v-points 
real(8):: x_gl(1:nx_gl) = 0.         ! global grid's x .
real(8):: xu_gl(1:nx_gl+1) = 0.      ! global grid's x at v-points 
real(8):: mu(dimy1_s:dimy2_s) = 1.   ! cos(lat) if dolatlon = .true.
real(8):: muv(dimy1_v:dimy2_v) = 1.  ! cos(lat) at v-points if dolatlon = .true.
real(8):: tanr(dimy1_s:dimy2_s) = 0. ! tan(lat)/r if dolatlon = .true.
real(8):: imu(dimy1_s:dimy2_s) = 1.  ! 1./mu
real(8):: imuv(dimy1_v:dimy2_v)= 1.  ! 1/muv
real(8):: lon_gl(1:nx_gl) = 0.       ! global grid's lon if dolatlon = .true.
real(8):: lonu_gl(1:nx_gl+1) = 0.    ! global grid's lon at u-points if dolatlon = .true.
real(8):: lat_gl(1:ny_gl) = 0.       ! global grid's lat if dolatlon = .true.
real(8):: latv_gl(1:ny_gl+1) = 0.    ! global grid's lat at v-points if dolatlon = .true.
real(8):: mu_gl(1:ny_gl) = 1.        ! global grid's cos(lat) if dolatlon = .true.
real(8):: muv_gl(1:ny_gl+1) = 1.     ! global grid's cos(lat) at v-points if dolatlon = .true.
real(8):: wgty(1:ny_gl) = 1.         ! weights for computing mean column integrals
real(8):: wgt(ny,nzm) = 1.           ! averaging weights (should be used together with terra() arrays)
real(8):: wgts(ny,nzm) = 1.          ! averaging weights for simple global mean (no regard for terrain)
real(8):: wgtxyt(nx,ny) = 1.         ! averaging weights at top
real(8):: wgtxys(nx,ny) = 1.         ! averaging weights at surface
real(8):: wgtu(ny,nzm) = 1.          ! averaging weights for u points
real(8):: wgtv(ny,nzm) = 1.          ! averaging weights for v points
real(8):: wgtw(ny,nz) = 1.           ! averaging weights for w points

real(8):: gu(ny,nzm)                 ! Jacobian for u
real(8):: gv(ny,nzm)                 ! Jacobian for v
real(8):: gw(ny,nzm)                 ! Jacobian for w
real(8):: igu(ny,nzm),igv(ny,nzm),igw(ny,nzm) ! inverse Jacobians

!----------------------------------------------------------------------------------------
! internally set parameters

integer:: nadams = 3             ! order of AB scheme, set depending on nadv_mom

character(80) case               ! id-string to identify a case-name(set in CaseName file)
integer:: nstep =0               ! current number of performed time steps 
integer:: nstep_run =0           ! current number of performed time steps of current run
integer  ncycle                  ! number of subcycles over the dynamical timestep
integer icycle                   ! current subcycle 
real cfl_adv                     ! CFL due to advection
real cfl_advh                    ! CFL due to advection (horizontal)
real cfl_advz                    ! CFL due to advection (vertical)
real cfl_sgs                     ! CFL due to SGS
integer:: na=1, nb=2, nc=3       ! indeces for swapping the rhs arrays for AB scheme
real:: at = 0., bt = 0., ct = 0. ! coefficients for the Adams-Bashforth scheme 
real dtn                         ! current dynamical timestep (can be smaller or larger than dt)
real dtp                         ! actual time interval between changes in nstep
real dts                         ! actual time interval between call to collect statistics
real dtst                        ! accumulated time when sampling to collect statistics
real dtpa                        ! actual time interval between changes in nstep for preca_xy
real dtstat                      ! accumulated time between statistics writes
real(8):: dt3(3) = 0. 	         ! dynamical timesteps for three most recent time steps
integer:: ncallocean = 1         ! number of atmospheric steps to call ocean model
real(8):: time=0.	         ! current time in sec. (from beginning of the run)
real(8):: timeUTsec=0.	         ! current UT time in seconds since Jan. 1, 1900 00:00:00
real(8) day	                 ! current day (including fraction)
real(8) dtfactor                 ! dtn/dt
real:: error_max = 0             ! maximum relative error in GMG pressure solver
integer:: niter_max = 0.         ! maximum number of v-vyvles in GMG pressure solver 
logical :: gmg_initialized = .false. ! GMG initialization flag
logical:: doyvar = .false.       ! if true, than variable resolution in y
logical :: collect_coars =.false.! flag to collect coarse-grid statistics at a given time step

logical dostatis                 ! flag to permit the gathering of statistics
logical dostatisrad              ! flag to permit the gathering of radiation statistics
integer nstatis	                 ! the interval between substeps to compute statistics

logical :: compute_reffc = .false. ! flag indicating that effective radius (liquid) is computed
logical :: compute_reffi = .false. ! flag indicating that effective radius (ice) is computed

logical notopenedM               ! flag to see if the MISC output datafile is opened	
logical notopened2D              ! flag to see if the 2D output datafile is opened	
logical notopened2DL             ! flag to see if the 2D output datafile is opened	
logical notopened2DZ             ! flag to see if the 2D zonal-mean output datafile is opened	
logical notopened3D              ! flag to see if the 3D output datafile is opened	
logical notopenedmom             ! flag to see if the statistical moment file is opened

! parameters for safe restart (dosaferestart=.true.)
integer:: restart_number_w = 1 ! number for identifying current RESTART* directory to write restart 
integer:: restart_number_r = 1 ! number for identifying current RESTART* directory to read restart 

!-----------------------------------------
! parameters set by PARAMETER namelist (the rest can be found in params.f90)

real:: dx =0. 	! grid spacing in x direction (at equator)
real:: dy =0.	! grid spacing in y direction
real:: dlon =0. ! grid spacing in degrees in x direction (at equator)
real:: dlat =0.	! grid spacing in degrees in y direction
real:: dz =0.	! constant grid spacing in z direction (when doconstdz =.true.)
logical:: doconstdz = .false.  ! do constant vertical grid spacing set by dz

integer:: nstop =0           ! time step number to stop the integration
integer:: nelapse =999999999 ! time step number to elapse before stoping

real:: dt=0.	        ! timestep (s) interval between steps. The actual timestep can fluctuate.
real(8):: day0=0.	! starting day (including fraction)
integer:: year0=1	! base year 
integer:: year=0	! current year 
integer(8):: date0=0	! starting date in YYYMMDDHHMMSS format (integer)

integer:: nprint =1	! frequency of printing a listing (steps)
integer:: nrestart =0 ! switch to control starting/restarting of the model
integer:: nstat =1	! the interval in time steps to compute statistics
integer:: nstatfrq =1 ! frequency of computing statistics 
logical:: restart_sep =.false.  ! write separate restart files for sub-domains
integer:: nrestart_skip =0 ! number of skips of writing restart (default 0)
integer:: nrestart_steps = 1 ! frequency of writing restart  
logical:: dosaferestart = .false. ! if .true., write restart in alternating RESTART or RESTART1 directories for given run

character(80):: caseid =''! id-string to identify a run	
character(80):: caseid_restart =''! id-string for branch restart file 
character(80):: case_restart =''! id-string for branch restart file 

logical:: doisccp = .false.
logical:: domodis = .false.
logical:: domisr = .false.
logical:: dosimfilesout = .false.

logical:: doSAMconditionals = .false. !core updraft,downdraft conditional statistics
logical:: dosatupdnconditionals = .false.!cloudy updrafts,downdrafts and cloud-free

logical:: doscamiopdata = .false.! initialize the case from a SCAM IOP netcdf input file
logical:: dozero_out_day0 = .false.
character(len=120):: iopfile=''
logical :: isInitialized_scamiopdata = .false.
logical :: wgls_holds_omega = .false.
character(256):: rundatadir ='./RUNDATA' ! path to data directory

integer:: nsave3D =1     ! frequency of writting 3D fields (steps)
integer:: nsave3Dstart =99999999! timestep to start writting 3D fields
integer:: nsave3Dend  =99999999 ! timestep to end writting 3D fields
logical:: save3Dbin =.false.   ! save 3D data in binary format(no 2-byte compression)
logical:: save3Dsep =.true.   ! use separate file for each time sample
integer:: nfiles3D = 1  ! number of files per time sample
real   :: qnsave3D =0.    !threshold manimum cloud water(kg/kg) to save 3D fields
logical:: dogzip3D =.false.    ! gzip compress a 3D output file   
logical:: rad3Dout = .false. ! output additional 3D radiation foelds (like reff)
logical:: save3Dnetcdf = .false. ! write 3D output in NetCDF format (requires pNetCDF library)

integer:: nsave2D =1     ! frequency of writting 2D fields (steps)
integer:: nsave2Dstart =99999999! timestep to start writting 2D fields
integer:: nsave2Dend =99999999  ! timestep to end writting 2D fields
logical:: save2Dbin =.false.   ! save 2D data in binary format, rather than compressed
logical:: save2Dsep =.true.   ! write separate file for each time point for 2D output
logical:: save2Davg =.false.   ! flag to time-average 2D output fields 
logical:: save2Drada =.false.   ! flag to save radiation averaged radiation fields 
logical:: save2Dradac =.false.   ! flag to save radiation accumulated fields 
logical:: dogzip2D =.false.    ! gzip compress a 2D output file if save2Dsep=.true.   
logical:: save2Dnetcdf = .false. ! write 2D output in NetCDF format (requires pNetCDF library)
logical:: snow2Dout =.false.   ! Output snow 2D fields 

integer:: nsaveM =1     ! frequency of writting MISC fields (steps)
integer:: nsaveMstart =99999999! timestep to start writting MISC fields
integer:: nsaveMend =99999999  ! timestep to end writting MISC fields
logical:: saveMbin =.false.   ! save MISC data in binary format, rather than compressed
logical:: saveMsep =.true.   ! write separate file for each time point for MISC output
logical:: saveMavg =.false.   ! flag to time-average MISC output fields 
logical:: dogzipM =.false.    ! gzip compress a MISC output file if save2Dsep=.true.
logical:: saveMnetcdf = .false. ! write MISC output in NetCDF format (requires pNetCDF library)

integer:: nsave2DL =1     ! frequency of writting 2D fields from LSM (steps)
integer:: nsave2DLstart =99999999! timestep to start writting 2D fields
integer:: nsave2DLend =99999999  ! timestep to end writting 2D fields
logical:: save2DLbin =.false.   ! save 2D data in binary format, rather than compressed
logical:: save2DLsep =.true.   ! write separate file for each time point for 2D output
logical:: save2DLavg =.false.   ! flag to time-average 2D output fields (default .false.)
logical:: dogzip2DL =.false.    ! gzip compress a 2D output file if save2Dsep=.true.   
logical:: save2DLnetcdf = .false. ! write 2D output in NetCDF format (requires pNetCDF library)

integer:: nsave2DZ =1     ! frequency of writting 2D fields (steps)
integer:: nsave2DZstart =99999999! timestep to start writting 2D fields
integer:: nsave2DZend =99999999  ! timestep to end writting 2D fields
logical:: save2DZsep =.true.   ! write separate file for each time point for 2D output
logical:: dogzip2DZ =.false.    ! gzip compress a 2D output file if save2Dsep=.true.
logical:: save2DZnetcdf = .false. ! write 2D output in NetCDF format (requires pNetCDF library)

integer:: nstatmom =1! frequency of writting statistical moment fields (steps)
integer:: nstatmomstart =99999999! timestep to start writting statistical moment fields
integer:: nstatmomend =99999999  ! timestep to end writting statistical moment fields
logical:: savemomsep =.true.! use one file with stat moments for each time point
logical:: savemombin =.false.! save statistical moment data in binary format

integer:: nstatcoars =1! frequency of writting coarse fields (steps)
integer:: nstatcoarsstart =99999999! timestep to start writting coarse fields
integer:: nstatcoarsend =99999999  ! timestep to end writting coarse fields
logical:: savecoarsbin =.false.! save ciarse fields in binary format
integer:: nfilescoars = 1 ! number of files in one time sample

integer:: nmovie =1! frequency of writting movie fields (steps)
integer:: nmoviestart =99999999! timestep to start writting statistical moment fields
integer:: nmovieend =99999999  ! timestep to end writting statistical moment fields

!-----------------------------------------
character(14) datechar  ! current date in YYYMMDDHHMMSS format (character)
character(19) date_pr   ! date in YYYY-MM-DD-HH:MM:SS format

end module grid
