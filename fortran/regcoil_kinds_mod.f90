! Minimal replacement for mini_libstell's stel_kinds & stel_constants:
! the retained kernels only ever needed `dp`/`rprec` and a few constants.
module regcoil_kinds_mod

  implicit none

  integer, parameter :: rprec = selected_real_kind(12, 100)
  integer, parameter :: dp = rprec

  real(dp), parameter :: pi = 3.14159265358979323846264338328_dp
  real(dp), parameter :: twopi = 2 * pi
  real(dp), parameter :: mu0 = 2 * twopi * 1.0e-7_dp

end module regcoil_kinds_mod
