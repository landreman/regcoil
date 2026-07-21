# REGCOIL package build: Fortran library + Python iso_c_binding extension.
project(
  'regcoil',
  ['c', 'fortran'],
  # Placeholder only.  The package version is declared statically in
  # pyproject.toml ([project] version), which meson-python uses for all wheel
  # and sdist metadata; this field is ignored.  Do not bump it -- keeping it at
  # 0.0.0 avoids a stale number that looks authoritative but is not.
  version: '0.0.0',
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

# BLAS lookup, most specific first.  Debian/Ubuntu ship blas.pc (libblas-dev);
# the manylinux images are AlmaLinux, where openblas-devel installs openblas.pc
# but no blas.pc, and RHEL-family systems may present FlexiBLAS instead.  Trying
# all three by pkg-config keeps the wheel builds working without special-casing
# the CI images here.
blas_dep = dependency('blas', 'openblas', 'flexiblas', required: false)
if not blas_dep.found()
  if host_machine.system() == 'darwin'
    blas_dep = declare_dependency(link_args: ['-framework', 'Accelerate'])
  else
    # No .pc file: fall back to plain library search, again openblas-then-blas.
    blas_dep = fc.find_library('openblas', required: false)
    if not blas_dep.found()
      blas_dep = fc.find_library('blas', required: true)
    endif
  endif
endif

regcoil_deps = [openmp_dep, blas_dep]

subdir('fortran')
subdir('src/regcoil')
