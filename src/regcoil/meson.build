# Python package sources + compiled _core extension (iso_c_binding).

py.install_sources(
  '__init__.py',
  '_constants.py',
  '_io.py',
  'surface.py',
  'fourier_surface.py',
  'plasma_surface.py',
  'coil_surface.py',
  'regcoil.py',
  subdir: 'regcoil',
)

# _core.c crosses the boundary with NumPy arrays, so it needs
# NumPy's C headers. There is no meson "numpy" dependency; ask the
# interpreter that is doing the build for its include directory instead
# (standard meson-python idiom). Passed via c_args (a raw -I flag) rather
# than include_directories(), since the latter rejects absolute paths that
# happen to resolve inside the source tree (e.g. an in-repo venv).
numpy_incdir = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()

py.extension_module(
  '_core',
  '_core.c',
  include_directories: include_directories('.'),
  c_args: ['-I' + numpy_incdir],
  dependencies: [py.dependency(), regcoil_fortran_dep],
  link_language: 'fortran',
  install: true,
  subdir: 'regcoil',
)
