! Set the domain dimensionality, size and number of subdomains.

module domain

       integer, parameter :: COL1 = 1   ! = 0 for single column run, and 1 otherwise
       integer, parameter :: YES3D = 0  ! Domain dimensionality: 1 - 3D, 0 - 2D
!       integer, parameter :: nx_gl = 96 ! Number of grid points in X
!       integer, parameter :: ny_gl = 96 ! Number of grid points in X
!       integer, parameter :: nz_gl = 80 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 8 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 12 ! No of subdomains in y
       integer, parameter :: nx_gl = 64 ! Number of grid points in X
       integer, parameter :: ny_gl = 1 ! Number of grid points in X
       integer, parameter :: nz_gl = 256 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 64 ! No of subdomains in x
       integer, parameter :: nsubdomains_y  = 1 ! No of subdomains in y
!       integer, parameter :: nx_gl = 512 ! Number of grid points in X
!       integer, parameter :: ny_gl = 256 ! Number of grid points in X
!       integer, parameter :: nz_gl = 83 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 16 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 8 ! No of subdomains in y
!       integer, parameter :: nx_gl = 768 ! Number of grid points in X
!       integer, parameter :: ny_gl = 768 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 24 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 32 ! No of subdomains in y
!       integer, parameter :: nx_gl = 2048 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1024 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 64 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 32 ! No of subdomains in y
!       integer, parameter :: nx_gl = 3456 ! Number of grid points in X
!       integer, parameter :: ny_gl = 3456 ! Number of grid points in X
!       integer, parameter :: nz_gl = 83 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 54 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 64 ! No of subdomains in y
!       integer, parameter :: nx_gl = 1536 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1536 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 32 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 48 ! No of subdomains in y
!       integer, parameter :: nx_gl = 3072 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1536 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 64 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 48 ! No of subdomains in y
!       integer, parameter :: nx_gl = 6144 ! Number of grid points in X
!       integer, parameter :: ny_gl = 3072 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 96 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 64 ! No of subdomains in y
!       integer, parameter :: nx_gl = 2048 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1024 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 64 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 32 ! No of subdomains in y
!       integer, parameter :: nx_gl = 1536 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1536 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 32 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 48 ! No of subdomains in y
!       integer, parameter :: nx_gl = 1024 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1024 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 32 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 32 ! No of subdomains in y
!       integer, parameter :: nx_gl = 2048 ! Number of grid points in X
!       integer, parameter :: ny_gl = 2048 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 32 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 64 ! No of subdomains in y
!       integer, parameter :: nx_gl = 1440 ! Number of grid points in X
!       integer, parameter :: ny_gl = 720 ! Number of grid points in X
!       integer, parameter :: nz_gl = 72 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 10 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 12 ! No of subdomains in y


!------------------------------------------------------------------------------------- 

       ! define # of points in x and y direction to average for 
       !   output relating to statistical moments.
       ! For example, navgmom_x = 8 means the output will be   
       !  8 times coarser grid than the original.
       ! If don't wanna such output, just set them to -1 in both directions. 
       ! See Changes_log/README.UUmods for more details.

       integer, parameter :: navgmom_x = -1 
       integer, parameter :: navgmom_y = -1 
       integer, parameter :: ncoars_x = 8 
       integer, parameter :: ncoars_y = 8

       integer, parameter :: ntracers = 0 ! number of transported tracers (dotracers=.true.)

!       integer, parameter :: nz_ocn = 100 ! Number of ocean scalar levels
!       integer, parameter :: nz_ocn = 48 ! Number of ocean scalar levels
!       integer, parameter :: nz_ocn = 4 ! Number of ocean scalar levels

!       integer, parameter :: ocn_ntracers = 0 ! number of transported tracers (dotracers=.true.) in ocean



! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5. (see User's Guide). Note 2 is required (only even numbers!)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;
!  * if ntracers is > 0, don't forget to set dotracers to .true. in namelist 


! Acceptable values for nx_gl and ny_gl 
!	8	10	12	16	18
!	20	24	30	32	36	40	48	50
!	54	60	64	72	80	90	96	100
!	108	120	128	144	150	160	162	180
!	192	200	216	240	250	256	270	288
!	300	320	324	360	384	400	432	450
!	480	486	500	512	540	576	600	640
!	648	720	750	768	800	810	864	900
!	960	972	1000	1024	1080	1152	1200	1250
!	1280	1296	1350	1440	1458	1500	1536	1600
!	1620	1728	1800	1920	1944	2000	2048	2160
!	2250	2304	2400	2430	2500	2560	2592	2700
!	2880	2916	3000	3072	3200	3240	3456	3600
!	3750	3840	3888	4000	4050	4096	4320	4374
!	4500	4608	4800	4860	5000	5120	5184	5400
!	5760	5832	6000	6144	6250	6400	6480	6750
!	6912	7200	7290	7500	7680	7776	8000	8100
!	8192	8640	8748	9000	9216	9600	9720	10000
!	10240	10368	10800	11250	11520	11664	12000	12150
!	12288	12500	12800	12960	13122	13500	13824	14400
!	14580	15000	15360	15552	16000	16200	16384	17280
!	17496	18000	18432	18750	19200	19440		

end module domain
