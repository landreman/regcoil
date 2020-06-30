# Makefile for REGCOIL
#
# For NERSC Edison and Cori you must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

LIBSTELL_DIR = mini_libstell
LIBSTELL_FOR_REGCOIL=$(LIBSTELL_DIR)/mini_libstell.a

ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else
        HOSTNAME?="laptop"
endif

HOSTNAME=draco

ifeq ($(HOSTNAME),edison)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -openmp -mkl
	EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 24

else ifeq ($(HOSTNAME),cori)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -qopenmp -mkl
	EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32

else ifeq ($(HOSTNAME),pppl)
	NETCDF_F = /usr/pppl/gcc/6.1-pkgs/netcdf-fortran-4.4.4
	NETCDF_C = /usr/pppl/gcc/6.1-pkgs/netcdf-c-4.4.1
	FC = mpifort
	EXTRA_COMPILE_FLAGS = -O2 -ffree-line-length-none -fexternal-blas -fopenmp -I$(NETCDF_F)/include -I$(NETCDF_C)/include
	EXTRA_LINK_FLAGS =  -fopenmp -L$(ACML_HOME)/lib -lacml -L$(NETCDF_C)/lib -lnetcdf -L$(NETCDF_F)/lib -lnetcdff
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32
        LIBSTELL_DIR=/u/slazerso/bin/libstell_dir
        LIBSTELL_FOR_REGCOIL=/u/slazerso/bin/libstell.a
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -N 1 -n 1 -c 8 -p dawson
else ifeq ($(HOSTNAME),draco)
	FC = mpiifort
	#EXTRA_COMPILE_FLAGS = -qopenmp -mkl -I${NETCDF_HOME}/include
	#EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_ -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff
	EXTRA_COMPILE_FLAGS = -qopenmp  -I${MKLROOT}/include -I${NETCDF_HOME}/include
	EXTRA_LINK_FLAGS =  -qopenmp  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -Wl,-ydgemm_ -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun

else
	FC = mpif90
	#EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
	EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -O0 -g
	EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf -framework Accelerate

	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB =
endif


# End of system-dependent variable assignments

TARGET = regcoil

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

lib$(TARGET).a: $(OBJ_FILES)
	ar rcs lib$(TARGET).a $(OBJ_FILES)

$(TARGET): lib$(TARGET).a $(TARGET).o $(LIBSTELL_FOR_REGCOIL)
	$(FC) -o $(TARGET) $(TARGET).o lib$(TARGET).a $(LIBSTELL_FOR_REGCOIL) $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET) *.a
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export REGCOIL_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export REGCOIL_RETEST=yes && ./runExamples.py

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo LIBSTELL_DIR is $(LIBSTELL_DIR)
	@echo LIBSTELL_FOR_REGCOIL is $(LIBSTELL_FOR_REGCOIL)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
