# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module:
#   module load cray-netcdf
# For Cori is is also necessary to run
#   module swap intel/16.0.0.109 intel/15.0.1.133
# to avoid a bug in the Intel MKL!!!
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

# For batch systems, set the following variable to the command used to run jobs.
# This variable is used by 'make test'.
REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 24
#REGCOIL_COMMAND_TO_SUBMIT_JOB =

ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else
        HOSTNAME="laptop"
endif

ifeq ($(HOSTNAME),edison)
	FC = ftn
	EXTRA_COMPILE_FLAGS = -openmp -mkl
	EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_
else ifeq ($(HOSTNAME),cori)
	FC = ftn
	EXTRA_COMPILE_FLAGS = -qopenmp -mkl
	EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
else
	#FC = gfortran
	FC = mpif90
	#EXTRA_COMPILE_FLAGS = -openmp -I/opt/local/include -ffree-form -ffree-line-length-none -ffixed-line-length-none -traditional
	EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
	EXTRA_LINK_FLAGS =  -fopenmp
endif
##LIBSTELL_DIR = /global/homes/l/landrema/20150410-02-stellinstall_245_edison/LIBSTELL/Release

#FC = ftn

## NERSC documentation recommends against specifying -O3
#EXTRA_COMPILE_FLAGS = -openmp -mkl
##EXTRA_COMPILE_FLAGS = -O3 -openmp -mkl
##EXTRA_COMPILE_FLAGS = -O0 -g -openmp
## -mkl MUST APPEAR AT THE END!!
#EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_

# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.

# End of system-dependent variable assignments

LIBSTELL_DIR = mini_libstell
TARGET = regcoil

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(TARGET): $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET)
	cd $(LIBSTELL_DIR); rm -f *.o *.mod *.MOD *.a

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export REGCOIL_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export REGCOIL_RETEST=yes && ./runExamples.py

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
