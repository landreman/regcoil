# mini_libstell + REGCOIL Fortran sources (linked into the Python extension).

mini_libstell_src = files(
  'mini_libstell/stel_kinds.f',
  'mini_libstell/stel_constants.f',
  'mini_libstell/vsvd0.f',
  'mini_libstell/vparams.f',
  'mini_libstell/system_mod.f',
  'mini_libstell/safe_open_mod.f',
  'mini_libstell/parse_extension.f',
  'mini_libstell/getcarg.f',
  'mini_libstell/vmec_getenv.f',
  'mini_libstell/handle_err.f90',
  'mini_libstell/v3_utilities.f',
  'mini_libstell/ezcdf_inqvar.f90',
  'mini_libstell/ezcdf_opncls.f90',
  'mini_libstell/ezcdf_attrib.f90',
  'mini_libstell/ezcdf_GenGet.f90',
  'mini_libstell/ezcdf_GenPut.f90',
  'mini_libstell/ezcdf.f90',
  'mini_libstell/vmec_input.f',
  'mini_libstell/mgrid_mod.f',
  'mini_libstell/read_wout_mod.F',
)

regcoil_lib_src = files(
  'regcoil_variables.f90',
  'regcoil_init_Fourier_modes_mod.f90',
  'regcoil_splines.f90',
  'regcoil_fzero.f',
  'regcoil_compute_offset_surface_mod.f90',
  'regcoil_read_efit_mod.f90',
  'regcoil_read_input.f90',
  'regcoil_validate_input.f90',
  'regcoil_compute_lambda.f90',
  'regcoil_init_plasma_mod.f90',
  'regcoil_evaluate_coil_surface.f90',
  'regcoil_expand_plasma_surface.f90',
  'regcoil_read_nescin.f90',
  'regcoil_write_nescin.f90',
  'regcoil_init_coil_surface.f90',
  'regcoil_read_bnorm.f90',
  'regcoil_build_matrices.f90',
  'regcoil_prepare_solve.f90',
  'regcoil_diagnostics.f90',
  'regcoil_solve.f90',
  'regcoil_lambda_scan.f90',
  'regcoil_compute_diagnostics_for_nescout_potential.f90',
  'regcoil_auto_regularization_solve.f90',
  'regcoil_write_input.f90',
  'regcoil_write_output.f90',
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
