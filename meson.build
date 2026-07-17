# Packaging scaffold (Phase 2). Fortran extension wiring is Phase 4 (ADR-002).
project(
  'regcoil',
  version: '0.1.0',
  meson_version: '>=1.2.0',
  default_options: ['warning_level=1'],
)

py = import('python').find_installation(pure: true)

subdir('src/regcoil')

# Fortran sources live under fortran/ for the legacy makefile and for the
# future meson extension build (Phase 4). Not compiled by Meson yet.
