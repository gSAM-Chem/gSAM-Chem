# Makefile to compile conversion utilities. You need to set include and library paths for NETCDF
#

PLATFORM := $(shell uname -s)

# uncomment line for your platform or add your own:

#---------------------------------------------------------
# AIX 

ifeq ($(PLATFORM),AIX)

LIB_NETCDF = /usr/local/xlf/lib
INC_NETCDF = /usr/local/xlf/include

FF = xlf90_r -qfixed -I$(INC_NETCDF) 
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

endif
#---------------------------------------------------------
# SGI IRIX64


ifeq ($(PLATFORM),IRIX64)


LIB_NETCDF = /usr/local/lib

FF = f90 -fixedform -extend_source -I$(INC_NETCDF) 
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

endif
#---------------------------------------------------------
# Mac OS X 


ifeq ($(PLATFORM),Darwin)


INC_NETCDF      := /usr/local/netcdf/include
LIB_NETCDF       := /usr/local/netcdf/lib

FF = mpif90 -fixed -g -O3 -traceback -pad -extend_source 
#FF = mpif90 -fixed -g -O3 -traceback -pad -extend_source -I$(INC_NETCDF)
#FF = mpif90 -fixed -extend_source -g -O0 -nowarn -ftrapuv -check all -debug full -traceback -check noarg_temp_created  -heap-arrays 1000  -I$(INC_NETCDF)
#FF = xlf90 -qfixed -I$(INC_NETCDF)
#FF = f95 -qfixed -kind=byte -gline -I$(INC_NETCDF)
#FF = gfortran  -O3 -ffixed-form -ffixed-line-length-0 -I$(INC_NETCDF)

LDFLAGS = -L${LIB_NETCDF} -lnetcdf


endif
#---------------------------------------------------------
#
ifeq ($(PLATFORM),Linux)

FF = ifx -O3 -m64 -traceback -fixed -extend-source -mcmodel=large -I${NETCDF}/include
LD = ifx  
LDFLAGS = -L${NETCDF}/lib -lnetcdff

endif

#---------------------------------------------------------
#---------------------------------------------------------

VPATH = ./SRC

all: 3Dsep2nc 3Dsep2nc_sepfields 2D2nc 2DZ2nc stat2nc isccp2nc  3D2nc 3D2nc_sepfields glue_movie_raw

.f:   
	$(FF) -o $@ -I./SRC $< ./SRC/hbuf_lib.f ./SRC/cape.f ./SRC/cin.f $(LDFLAGS) 

clean: 
	rm 2D* 3D* stat* isccp* *.o glue*
