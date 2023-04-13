# Makefile for REGCOIL
#
# For NERSC Edison and Cori you must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

LIBSTELL_DIR ?= mini_libstell
LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/mini_libstell.a
#LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/libstell.a

ifdef NERSC_HOST
  HOSTNAME = $(NERSC_HOST)
else
  HOSTNAME?="laptop"
endif

ifdef MACHINE
  HOSTNAME = $(MACHINE)
endif

ifeq ($(HOSTNAME),cori)
  REGCOIL_HOST = cori
  FC = ftn
  ## NERSC documentation recommends against specifying -O3
  ## -mkl MUST APPEAR AT THE END!!
  EXTRA_COMPILE_FLAGS = -qopenmp -mkl
  EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
  # Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32

else ifeq ($(HOSTNAME),marconi)
  REGCOIL_HOST = marconi
  FC = mpiifort
  EXTRA_COMPILE_FLAGS = -qopenmp -mkl $(shell nc-config --fflags)
  EXTRA_LINK_FLAGS =  -qopenmp -mkl $(shell nc-config --flibs)
  # Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32
  LIBSTELL_DIR=~/bin/libstell_dir
  LIBSTELL_FOR_REGCOIL=~/bin/libstell.a

else ifeq ($(HOSTNAME),pppl_gcc)
  REGCOIL_HOST=pppl_gcc
  NETCDF_F = $(NETCDF_FORTRAN_HOME)
  NETCDF_C = $(NETCDF_C_HOME)
  FC = mpifort
  EXTRA_COMPILE_FLAGS = -O2 -ffree-line-length-none -fexternal-blas -fopenmp -I$(NETCDF_F)/include -I$(NETCDF_C)/include
  EXTRA_LINK_FLAGS = -fopenmp -lopenblas -L$(NETCDF_C)/lib -lnetcdf -L$(NETCDF_F)/lib -lnetcdff
  LIBSTELL_DIR=$(STELLOPT_PATH)/LIBSTELL/Release
  LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/libstell.a
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -N 1 -n 1 -c 8 -q debug --mem 8G

else ifeq ($(HOSTNAME),pppl_intel)
  REGCOIL_HOST=pppl_intel
  FC = mpifort
  NETCDF_F = $(NETCDF_FORTRAN_HOME)
  NETCDF_C = $(NETCDF_C_HOME)
  EXTRA_COMPILE_FLAGS =-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback \
    -Wl,--end-group \
    -qopenmp \
    -lpthread \
    -I$(NETCDF_F)/include -I$(NETCDF_C)/include \
    -mkl	
  EXTRA_LINK_FLAGS = -qopenmp -mkl -L$(NETCDF_C)/lib -lnetcdf -L$(NETCDF_F)/lib -lnetcdff
  LIBSTELL_DIR=$(STELLOPT_PATH)/LIBSTELL/Release
  LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/libstell.a
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -N 1 -n 1 -c 8 -q debug --mem 8G
else ifeq ($(HOSTNAME),stellar-intel.princeton.edu)
  REGCOIL_HOST=stellar
  FC = ifort
  NETCDF_F = $(NETCDFDIR)
  NETCDF_C = $(NETCDFDIR)
  EXTRA_COMPILE_FLAGS =-mcmodel=large -O3 -m64 -unroll0 \
    -Wl,--end-group \
    -lpthread \
    -I$(NETCDF_F)/include -I$(NETCDF_C)/include
  EXTRA_LINK_FLAGS = -fopenmp -L$(NETCDF_C)/lib -lnetcdf -L$(NETCDF_F)/lib -lnetcdff -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -N 1 -n 1 -c 8 -q debug --mem 8G --time 00:10:00
  LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/mini_libstell.a
else ifeq ($(HOSTNAME),osx_brew)
  REGCOIL_HOST=osx_brew
  NETCDF_HOME ?= /usr/local/lib/netcdf
  LAPACK_HOME ?= /usr/local/lib/lapack
  FC = mpifort 
  EXTRA_COMPILE_FLAGS = -O2 -ffree-line-length-none -fopenmp -I$(NETCDF_HOME)/include
  EXTRA_LINK_FLAGS =  -fopenmp -L$(NETCDF_HOME)/lib -lnetcdf -lnetcdff -L$(LAPACK_HOME)/lib -lblas -llapack
  LIBSTELL_DIR = $(STELLOPT_PATH)/LIBSTELL/Release
  LIBSTELL_FOR_REGCOIL = $(LIBSTELL_DIR)/libstell.a
  REGCOIL_COMMAND_TO_SUBMIT_JOB = 

else ifeq ($(CLUSTER),DRACO)
  REGCOIL_HOST=draco
  FC = mpiifort
  #EXTRA_COMPILE_FLAGS = -qopenmp -mkl -I${NETCDF_HOME}/include
  #EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_ -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff
  EXTRA_COMPILE_FLAGS = -qopenmp  -I${MKLROOT}/include -I${NETCDF_HOME}/include
  EXTRA_LINK_FLAGS =  -qopenmp  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -Wl,-ydgemm_ -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff
  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun

else ifeq ($(CLUSTER),RAVEN)
  REGCOIL_HOST=raven
  FC = mpiifort
  EXTRA_COMPILE_FLAGS = -qopenmp  -I${MKLROOT}/include $(shell nc-config --fflags)
  EXTRA_LINK_FLAGS =  -qopenmp  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -Wl,-ydgemm_ $(shell nc-config --flibs)
  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
  REGCOIL_COMMAND_TO_SUBMIT_JOB = srun

else
  REGCOIL_HOST=macports
  FC = mpif90
  #EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
  #EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -O0 -g
  #EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf -framework Accelerate
  EXTRA_COMPILE_FLAGS = -fopenmp -I/usr/local/include -ffree-line-length-none -O0 -g -fallow-argument-mismatch
  EXTRA_LINK_FLAGS =  -fopenmp -L/usr/local/lib -lnetcdff  -lnetcdf -framework Accelerate

  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
  REGCOIL_COMMAND_TO_SUBMIT_JOB =
endif


# End of system-dependent variable assignments

TARGET = regcoil

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 ${LIBSTELL_FOR_REGCOIL}
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f ${LIBSTELL_FOR_REGCOIL}
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

lib$(TARGET).a: $(OBJ_FILES)
	ar rcs lib$(TARGET).a $(OBJ_FILES)

$(TARGET): lib$(TARGET).a $(TARGET).o $(LIBSTELL_FOR_REGCOIL)
	$(FC) -o $(TARGET) $(TARGET).o lib$(TARGET).a $(LIBSTELL_FOR_REGCOIL) $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

mini_libstell/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET) *.a
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export REGCOIL_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export REGCOIL_RETEST=yes && ./runExamples.py

test_make:
	@echo REGCOIL_HOST is $(REGCOIL_HOST)
	@echo HOSTNAME is $(HOSTNAME)
	@echo MACHINE is $(MACHINE)
	@echo FC is $(FC)
	@echo LIBSTELL_DIR is $(LIBSTELL_DIR)
	@echo LIBSTELL_FOR_REGCOIL is $(LIBSTELL_FOR_REGCOIL)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
