# REGCOIL package build: Fortran library + Python iso_c_binding extension.
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
  language: 'fortran',
)

# OpenMP: dependency('openmp') detects via the C compiler, which fails on
# macOS with Apple Clang.  _core.c does not use OpenMP; only the Fortran
# library does.  On macOS, inject the flag per-language so that the C
# compiler is not involved.
if host_machine.system() == 'darwin'
  add_project_arguments('-fopenmp', language: 'fortran')
  add_project_link_arguments('-fopenmp', language: 'fortran')
  openmp_dep = declare_dependency()
else
  openmp_dep = dependency('openmp')
endif

blas_dep = dependency('blas', required: false)
if not blas_dep.found()
  if host_machine.system() == 'darwin'
    blas_dep = declare_dependency(link_args: ['-framework', 'Accelerate'])
  else
    blas_dep = fc.find_library('blas', required: true)
  endif
endif

regcoil_deps = [openmp_dep, blas_dep]

subdir('fortran')
subdir('src/regcoil')
