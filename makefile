# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module:
#   module load cray-netcdf
# For Cori is is also necessary to run
#   module swap intel/16.0.0.109 intel/15.0.1.133
# to avoid a bug in the Intel MKL!!!
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

FC = ftn

#LIBSTELL_DIR = /global/homes/l/landrema/20150410-02-stellinstall_245_edison/LIBSTELL/Release

# NERSC documentation recommends against specifying -O3
EXTRA_COMPILE_FLAGS = -openmp -mkl
#EXTRA_COMPILE_FLAGS = -O3 -openmp -mkl
#EXTRA_COMPILE_FLAGS = -O0 -g -openmp
# -mkl MUST APPEAR AT THE END!!
EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_

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
