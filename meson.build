# REGCOIL package build (Phase 4): Fortran library + Python iso_c_binding extension.
project(
  'regcoil',
  ['c', 'fortran'],
  version: '0.1.0',
  meson_version: '>=1.2.0',
  default_options: [
    'warning_level=1',
    'buildtype=release',
  ],
)

py = import('python').find_installation(pure: false)
fc = meson.get_compiler('fortran')

# Match flags used by the legacy makefile / CI hosts.
add_project_arguments(
  '-ffree-line-length-none',
  '-fallow-argument-mismatch',
  '-DNETCDF',
  language: 'fortran',
)

openmp_dep = dependency('openmp')
netcdf_dep = dependency('netcdf')
netcdff_dep = dependency('netcdf-fortran', required: false)
if not netcdff_dep.found()
  netcdff_dep = dependency('netcdff', required: false)
endif
if not netcdff_dep.found()
  # Fallback when pkg-config names differ (still need the Fortran library).
  netcdff_lib = fc.find_library('netcdff', required: true)
  netcdff_dep = declare_dependency(dependencies: netcdff_lib)
endif

lapack_dep = dependency('lapack', required: false)
blas_dep = dependency('blas', required: false)
if not lapack_dep.found() or not blas_dep.found()
  if host_machine.system() == 'darwin'
    accelerate_dep = declare_dependency(link_args: ['-framework', 'Accelerate'])
    lapack_dep = accelerate_dep
    blas_dep = declare_dependency()
  else
    lapack_dep = fc.find_library('lapack', required: true)
    blas_dep = fc.find_library('blas', required: true)
  endif
endif

regcoil_deps = [openmp_dep, netcdf_dep, netcdff_dep, blas_dep, lapack_dep]

subdir('fortran')
subdir('src/regcoil')
