# REGCOIL Fortran sources linked into the Python extension.
#
# The extension boundary is ~4 stateless kernels (Phase 7 / ADR-020); the
# legacy regcoil_variables-coupled chain (namelist I/O, geometry init,
# build_matrices/prepare_solve/solve/diagnostics, lambda scan) and
# mini_libstell (ADR-005) are gone (Phase 13).

regcoil_lib_src = files(
  'regcoil_kinds_mod.f90',
  'regcoil_kernels_mod.f90',
  'regcoil_c_api.f90',
)

# Static library keeps Fortran objects out of a separate installable .so while
# still linking into the Python extension.
regcoil_fortran = static_library(
  'regcoil_fortran',
  regcoil_lib_src,
  dependencies: regcoil_deps,
  install: false,
)

regcoil_fortran_dep = declare_dependency(
  link_with: regcoil_fortran,
  include_directories: include_directories('.'),
  dependencies: regcoil_deps,
)
