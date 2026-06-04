# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
# NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

#SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)
SAM = gSAM

# Determine platform 
PLATFORM := $(shell uname -s)

#----------------------------------------------------------------------
# UNIVERSITY OF WASHINGTON
# Linux, Intel Compiler
#----------------------------------------------------------------------

ifeq ($(PLATFORM), Linux)

	FF77 = mpif90 -c -fixed -extend-source -r8
	FF90 = mpif90 -c -r8 
	CC = mpicc -c -DLINUX 

	# Trim HOSTNAME to strip off .atmos.washington.edu
	HOSTTRIM = $(basename $(basename $(basename $(basename $(HOSTNAME)))))

	ifeq ($(HOSTTRIM),olympus)
		# Compiling on olympus or challenger
		# module load intel/19.0.2 netcdf/4.6.1.i19 openmpi/3.1.3
		FFLAGS = -O2 -axSSE4.2 -fp-model source -debug -traceback 
		FFLAGS_O1 = -O1 -fp-model source 
		NCPATH = ${NETCDF}

	else ifeq ($(HOSTTRIM),hermes)
		FFLAGS = -O2 -axSSE4.2 -fp-model source -debug -traceback -lnetcdf -check all -fpp -check bounds -fpe0 -init=snan
		FFLAGS_O1 = -O1 -fp-model source
		NCPATH = ${NETCDF}
		
	else
	# Compiling on challenger or one of the olympus nodes.  Use older compiler setup.
	FFLAGS = -O2 -xHOST -fp-model source -mcmodel=large 
	FFLAGS_O1 = -O1 -fp-model source #	
	NCPATH = ${NETCDF}
	endif

	FFLAGS_NOOPT = -O0 -g -ftrapuv -fpe0 -check all -traceback -debug -gen-interfaces -warn interfaces -fp-model source
	#FFLAGS = ${FFLAGS_NOOPT}
	#FFLAGS_O1 = ${FFLAGS_NOOPT}

	LD = mpif90

	FFLAGS += -I$(NCPATH)/include
	FFLAGS_O1 += -I$(NCPATH)/include
	FFLAGS_NOOPT +=  -I$(NCPATH)/include
	LDFLAGS = -L$(NCPATH)/lib -Wl,-rpath $(NCPATH)/lib -lnetcdff -lnetcdf

	ifeq ($(MICRO_DIR),MICRO_M2005_PA)
	FFLAGS_MICRO = ${FFLAGS_O1}
	else
	FFLAGS_MICRO = ${FFLAGS}
	endif
endif

#----------------------------------
# AIX (tested only on IBM SP)
#

ifeq ($(PLATFORM),AIX)

#INC_MPI      := /usr/local/include
#LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib


FF77 = mpxlf90_r -c -qsuffix=f=f -qfixed=132
FF90 = mpxlf90_r -c -qsuffix=f=f90
CC = cc -c -DAIX
FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q -I${INC_NETCDF}
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en -I${INC_NETCDF}
LD = mpxlf90_r
LDFLAGS = -bmaxdata:512000000 -bmaxstack:256000000 -L${LIB_NETCDF} -lnetcdf

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(PLATFORM),IRIX64)

INC_MPI      := /usr/local/include
LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib

FF77 = f90 -c -fixedform  -extend_source
FF90 = f90 -c -freeform
CC = cc -c -DIRIX64
FFLAGS = -O3 
#FFLAGS = -g -DEBUG:subscript_check=ON:trap_uninitialized=ON 
FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
LD = f90 
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

endif
#----------------------------------------------------------------------
# Linux, Intel Compiler
#

#ifeq ($(PLATFORM),Linux)
#
#LIB_MPI = /usr/local/pkg/iopenmpi/lib
#INC_MPI = /usr/local/pkg/iopenmpi/include
#INC_NETCDF = /nfs/user08/marat/local/include
#LIB_NETCDF = /nfs/user08/marat/local/lib
#
#
#FF77 = /usr/local/pkg/iopenmpi/bin/mpif90 -c -fixed -extend_source
#FF90 = /usr/local/pkg/iopenmpi/bin/mpif90 -c
#CC = mpicc -c -DLINUX
#
#
#FFLAGS = -O3 
##FFLAGS = -g -ftrapuv -check all
#
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LD = /usr/local/pkg/iopenmpi/bin/mpif90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#

#endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at San Diego SC
#
#ifeq ($(PLATFORM),Linux)

#INC_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/include
#LIB_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/lib

#FF77 = mpxlf90  -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpxlf90  -qarch=440 -qsuffix=f=f90
#CC = mpcc -c -DLinux
#FFLAGS = -c -O3 -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpxlf90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf

#endif
#----------------------------------------------------------------------
# Linux,  handy computer at IACS
#

#ifeq ($(PLATFORM),Linux)

# ifort

#FF77 = mpiifort -c -fixed -extend_source
#FF90 = mpiifort -c
#CC = mpiicc -c -DLINUX


#FFLAGS = -O3 -pad
##FFLAGS = -g -ftrapuv -check all

#FFLAGS += -I${INC_NETCDF}
#LD = mpiifort
#LDFLAGS = -L${LIB_NETCDF} -lnetcdff 


#endif
#----------------------------------------------------------------------
# Linux, Derecho
# #

#  ifeq ($(PLATFORM),Linux)

#FF77 = mpipf90 -c -Mextend
#FF90 = mpipf90 -c -Mfreeform
#CC = mpipcc -c  -DLINUX
#
#FFLAGS = -Mnoframe -Mvect -Munroll -O2 -Mbyteswapio  
#
#FFLAGS += -I${INC_NETCDF}
#LD = mpipf90
#LDFLAGS =  -L${LIB_NETCDF} -lnetcdf


# IFORT
#  FF77 = mpif90 -c -fpp -fixed  -extend-source
#  FF90 = mpif90 -c -fpp -free 
#  FF95 = mpif90 -c -fpp -r8 -free 
#  CC = mpicc -c -O1 -DLINUX

# CRAY
# FF77 = mpif90 -c -f fixed -N 132
# FF90 = mpif90 -c -e Z -f free
# FF95 = mpif90 -c -s real64 -f free
# CC = mpicc -c -O1 -DLINUX
# FFLAGS = -O3  -g -hfp3 -O ipa5 -O aggress -O cache2 -eo


# #FFLAGS = -O3 -g 
# FFLAGS = -O3  -g -traceback -mcmodel=large -march=core-avx2 -assume buffered_io
# #FFLAGS = -O3  -g -traceback -fpe0 -mcmodel=large -march=core-avx2 -assume buffered_io
# #FFLAGS = -O3  -mcmodel=large -heap-arrays 1000
# #FFLAGS = -O3 -g -traceback -heap-arrays 1000 -pad -march=core-avx2 -assume buffered_io  -mcmodel=large
# #FFLAGS =  -Os -g -traceback -fpe0 -pad -mcmodel=large  
# #FFLAGS =  -Os -heap-arrays=1000 -g -traceback -pad -mcmodel=large  
# #FFLAGS = -Os -pad -no-prec-div -fp-model fast=2 -ipo -mcmodel=large
# #FFLAGS = -O3 -g -traceback -init=snan,arrays -pad -assume buffered_io  -mcmodel=large
# # debugging mode. Produces extreamly slow code.
# #FFLAGS = -g -O0 -fpe0 -heap-arrays 1000 -nowarn -ftrapuv -check all -init=snan,arrays -debug full -traceback -check noarg_temp_created -gen-interfaces -warn interfaces -W1 -mcmodel=large

# # LD = mpif90
#  LD = mpif90  -march=core-avx2 -mcmodel=large
# # LDFLAGS = -L${NETCDF}/lib -L${PNETCDF}/lib -lnetcdff -qmkl=sequential
# # LDFLAGS = -L${NETCDF}/lib -L${PNETCDF}/lib -lnetcdff -lpnetcdf -qmkl=sequential
#  LDFLAGS = -L${NETCDF}/lib -L${PNETCDF}/lib -lnetcdff  -lpnetcdf # -qmkl=sequential

# # NVFORTRAN
# #FF77 = mpif90 -c -Mpreprocess -Mfixed  
# #FF90 = mpif90 -c -Mpreprocess -Mfree 
# #FF95 = mpif90 -c -Mpreproc -r8 -Mfree 
# #CC = mpicc -c -O1 -DLINUX
# #FFLAGS = -O3 -mcmodel=medium
# #
# # LD = mpif90  -mcmodel=medium
# # LDFLAGS = -L${LIB_NETCDF} -lnetcdff -lblas -llapack

# FFLAGS += -I${NETCDF}/include -I${PNETCDF}/include

#  endif

#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#
#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#

# bloss:
#
#ifeq ($(PLATFORM),Linux)
#
#  # Default compiler flags
#  FFLAGS =  -Kieee -fastsse #-g -C -Ktrap=fp #
#
#  ifeq ($(HOSTNAME),olympus)
#    # special setup on olympus
#    MPIF90 = /usr/local/openmpi/bin/mpif90
#    FFLAGS += -tp k8-64
#
#  else
#    # Determine platform
#    PROCESSOR := $(shell uname -p)
#
#    # DEFAULT MPI LOCATION AND COMPILER
#    LIB_MPI = /usr/local/mpich/lib
#    INC_MPI = /usr/local/mpich/include
#    MPIF90 = pgf90
#
#    # CHANGE COMPILER IF RUNNING ON 64BIT MACHINE.
#    ifeq ($(PROCESSOR),x86_64)
#      MPIF90 = /usr/local/pgi/linux86-64/7.1-1/bin/pgf90
#      FFLAGS += -tp k8-64
#    else
#      # ADD -tp k8-32 IF RUNNING ON REX.
#      ifeq ($(HOSTNAME),rex)
#        FFLAGS += -tp k8-32
#      else
#        FFLAGS += -tp k7
#      endif
#    endif
#
#    ifeq ($(HOSTNAME),rex)
#      # UNCOMMENT TO USE LAHEY COMPILER -- USEFUL FOR DEBUGGING
#      LIB_MPI = /usr/local/mpich-lf/lib
#      INC_MPI = /usr/local/mpich-lf/include
#      MPIF90 = lf95
#      FFLAGS = -g --trap --chk aesu #--trap # --o2 #
#    endif
#
#    # UNCOMMENT FOR MYRINET RUN
#    #LIB_MPI = /usr/local/mpich-gm/lib
#    #INC_MPI = /usr/local/mpich-gm/include
#    #MPIF90 = /usr/local/mpich-gm/bin/mpif90
#
#    # ADD MPI FLAGS HERE IF NOT USING mpif90
#    FFLAGS += -I${INC_MPI}
#    LDFLAGS += -L${LIB_MPI} -lmpich
#  endif
#
#  FF77 = ${MPIF90} -c
#  FF90 = ${MPIF90} -c
#  CC = gcc -c -DLINUX -g
#
#  LD = ${MPIF90} ${FFLAGS}
#  FFLAGS += -I${INC_NETCDF}
#  LDFLAGS += -L${LIB_NETCDF} -lnetcdf
#
## end bloss:

# older options:

#LIB_MPI = /usr/pgi/linux86/5.1/lib
#INC_MPI = /usr/pgi/linux86/5.1/include
#LIB_MPI = /usr/local/lam-mpi/lib
#INC_MPI = /usr/local/lam-mpi/include

#FF77 = pgf90 -c
#FF90 = pgf90 -c -Mfreeform
#CC = cc -c  -DLINUX

#FFLAGS = -Mnoframe -Mvect -Munroll -O2 -Mbyteswapio  

#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LD = pgf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

#endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (Absoft Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/absoft/include
#LIB_NETCDF   := /usr/local/absoft/lib
#INC_MPI      := /usr/local/absoft/include
#LIB_MPI       := /usr/local/absoft/lib

#FF77 = f90 -c -f fixed
#FF90 = f90 -c -f free
#CC = cc -c -DMACOSX

#FFLAGS = -O3 -noconsole -nowdir -YEXT_NAMES=LCS -s -YEXT_SFX=_  -z4
#LD = f90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

#endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (NAG Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/nag/include
#LIB_NETCDF   := /usr/local/nag/lib
#INC_MPI      := /usr/local/nag/include
#LIB_MPI       := /usr/local/nag/lib

#FF77 = f95 -c -fixed -kind=byte
#FF90 = f95 -c -free -kind=byte
#CC = cc -c -DMACOSX

#FFLAGS =      # don't use any optimization -O* option! Will crash!
##FFLAGS = -gline -C=all -C=undefined  # use for debugging
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = f95 
##LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf 
#LDFLAGS =  

#endif


#----------------------------------
# Apple Mac OS X (Darwin) (XLF compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/xlf/include
#LIB_NETCDF   := /usr/local/xlf/lib
#INC_MPI      := /usr/local/xlf/include
#LIB_MPI       := /usr/local/xlf/lib

#FF77 = xlf90 -c -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -c -qsuffix=f=f90
#CC = cc -c -DMACOSX 
#FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = xlf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

#endif

#----------------------------------

#----------------------------------
# Apple Mac OS X (Darwin) (Intel compiler)
#

# ifeq ($(PLATFORM),Darwin)

# INC_NETCDF      := /usr/local/netcdf/include
# LIB_NETCDF       := /usr/local/netcdf/lib


# FF77 = mpif90 -c -fixed -extend_source
# FF90 = mpif90 -c 
# CC = mpicc -c -DLINUX


# FFLAGS = -Os -g -pad -traceback
# #FFLAGS = -g -ftrapuv -check all -traceback

# FFLAGS += -I${INC_NETCDF}
# LD = mpif90 
# LDFLAGS = -L${LIB_NETCDF} -lnetcdf -mkl=parallel

# endif
#
#----------------------------------
# Apple Mac OS X (Darwin) (GNU compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF := /usr/local/include
#LIB_NETCDF := /usr/local/lib
#
#FF77 = gfortran -c -ffixed-form -ffixed-line-length-0
#FF90 = gfortran -c -ffree-form -ffree-line-length-0
#CC = gcc -c -DLINUX
#
#
#FFLAGS = -O3
##FFLAGS = -g -fcheck=all
#
#FFLAGS += -I${INC_NETCDF}
#LD = gfortran
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#
#endif


#----------------------------------
#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .F90 .f .f90 .c .o 



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.F90.o:
	${FF90}  ${FFLAGS} $<
.F.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<


include Depends



clean: 
	rm ./OBJ/*


