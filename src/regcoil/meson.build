# Python package sources + compiled _core extension (iso_c_binding).

py.install_sources(
  '__init__.py',
  '_constants.py',
  '_io.py',
  '_serialize.py',
  'log.py',
  'surface.py',
  'fourier_surface.py',
  'plasma_surface.py',
  'coil_surface.py',
  'regcoil.py',
  'plot.py',
  'cut.py',
  '_cli.py',
  'examples.py',
  subdir: 'regcoil',
)

# Bundled example datasets. The source files stay in the top-level equilibria/
# directory; meson copies them into regcoil/equilibria/ inside the installed
# package so regcoil.examples(...) can reach them via importlib.resources,
# both in editable installs and in built wheels. See src/regcoil/examples.py.
py.install_sources(
  '../../equilibria/wout_li383_1.4m.nc',
  '../../equilibria/bnorm.li383_1.4m',
  '../../equilibria/nescin.li383_realWindingSurface',
  '../../equilibria/wout_d23p4_tm.nc',
  '../../equilibria/bnorm.d23p4_tm',
  '../../equilibria/nescin.w7x_winding_surface_from_Drevlak',
  subdir: 'regcoil/equilibria',
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
