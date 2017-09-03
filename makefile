# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.


ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else
        # HOSTNAME="laptop"
        HOSTNAME=PPPL
endif

$(info HOSTNAME is set to: $(HOSTNAME))

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
else ifeq ($(HOSTNAME),PPPL)
	FC = $(MPIHOME)/bin/mpif90
	#EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
	EXTRA_COMPILE_FLAGS = -fopenmp -I$(NETCDF_HOME)/include -I/usr/pppl/gcc/4.6-pkgs/openblas-48f06dd/include -ffree-line-length-none
	EXTRA_LINK_FLAGS =  -fopenmp -L$(NETCDF_HOME)/lib -lnetcdff  -lnetcdf -L/usr/pppl/gcc/4.6-pkgs/openblas-48f06dd/lib -lopenblas

	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB =
else
	FC = mpif90
	#EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
	EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none
	EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf -framework Accelerate

	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB =
endif


# End of system-dependent variable assignments

LIBSTELL_DIR = mini_libstell
TARGET = regcoil

# JCS Modified to use regcoil_input_mod and renamed global_variables to regcoil_variables and more. see the diff file
# REGCOILLIB_OBJ_FILES = auto_regularization_solve.o compute_offset_surface_mod.o init_Fourier_modes_mod.o  \
#	read_bnorm.o regcoil.o validate_input.o build_matrices.o expand_plasma_surface.o  \
#	init_coil_surface.o read_efit_mod.o solve.o write_output.o  \
#	compute_diagnostics_for_nescout_potential.o  fzero.o  \
#	init_plasma_mod.o read_input.o splines.o compute_lambda.o global_variables.o  \
#	init_surface_mod.o read_nescin.o svd_scan.o
REGCOILLIB_OBJ_FILES = regcoil_auto_regularization_solve.o compute_offset_surface_mod.o init_Fourier_modes_mod.o  \
	read_regcoil_bnorm.o regcoil.o validate_regcoil_input.o build_regcoil_matrices.o expand_plasma_surface.o  \
	init_regcoil_coil_surface.o read_efit_mod.o solve.o write_regcoil_output.o  \
	compute_diagnostics_for_nescout_potential.o  fzero.o  \
	init_regcoil_plasma.o regcoil_input_mod.o splines.o compute_regcoil_lambda.o regcoil_variables.o  \
	init_surface_mod.o read_nescin.o svd_scan.o
REGCOILLIB_TARGET = libregcoil.a

export

.PHONY: all clean

all: $(TARGET) $(REGCOILLIB_TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(REGCOILLIB_TARGET): $(REGCOILLIB_OBJ_FILES)
	ar rcs $(REGCOILLIB_TARGET) $(REGCOILLIB_OBJ_FILES)

$(TARGET): $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET) $(REGCOIL_TARGET)
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
