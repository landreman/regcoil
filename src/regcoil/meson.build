# Python package sources + compiled _core extension (iso_c_binding).

py.install_sources(
  '__init__.py',
  subdir: 'regcoil',
)

py.extension_module(
  '_core',
  '_core.c',
  include_directories: include_directories('.'),
  dependencies: [py.dependency(), regcoil_fortran_dep],
  link_language: 'fortran',
  install: true,
  subdir: 'regcoil',
)
