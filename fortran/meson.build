# mini_libstell + REGCOIL Fortran sources linked into the Python extension.
#
# The extension boundary is now ~4 stateless kernels, so
# most of the legacy regcoil_variables-coupled chain (namelist I/O, geometry
# init, build_matrices/prepare_solve/solve/diagnostics, lambda scan) is no
# longer compiled in here.

mini_libstell_src = files(
  'mini_libstell/stel_kinds.f',
  'mini_libstell/stel_constants.f',
)

regcoil_lib_src = files(
  'regcoil_fzero.f',
  'regcoil_kernels_mod.f90',
  'regcoil_uniform_offset_surface_mod.f90',
  'regcoil_c_api.f90',
)

# Static library keeps Fortran objects out of a separate installable .so while
# still linking into the Python extension.
regcoil_fortran = static_library(
  'regcoil_fortran',
  mini_libstell_src + regcoil_lib_src,
  dependencies: regcoil_deps,
  install: false,
)

regcoil_fortran_dep = declare_dependency(
  link_with: regcoil_fortran,
  include_directories: include_directories('.'),
  dependencies: regcoil_deps,
)
