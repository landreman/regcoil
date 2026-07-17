module regcoil_compute_offset_surface_mod

  use stel_kinds

  implicit none

  private

  public :: regcoil_compute_offset_surface_xyz_of_thetazeta

contains

  subroutine regcoil_compute_offset_surface_xyz_of_thetazeta(prob, theta_rootSolve_in, &
       zeta_rootSolve_target_in, x_offsetSurface, y_offsetSurface, z_offsetSurface, separation_in)

    use regcoil_variables, only: regcoil_t
    use stel_kinds

    implicit none

    type(regcoil_t), intent(inout) :: prob
    real(dp), intent(in) :: theta_rootSolve_in, zeta_rootSolve_target_in, separation_in
    real(dp), intent(out) :: x_offsetSurface, y_offsetSurface, z_offsetSurface

    integer :: fzeroFlag
    real(dp) :: rootSolve_abserr, rootSolve_relerr, zeta_rootSolve_min, zeta_rootSolve_max
    real(dp) :: zeta_plasma_rootSolveSolution
    ! Locals for fzero nested residual (no module scratch — instance-safe / OpenMP-safe).
    real(dp) :: theta_rootSolve, zeta_rootSolve_target, separation

    theta_rootSolve = theta_rootSolve_in
    zeta_rootSolve_target = zeta_rootSolve_target_in
    separation = separation_in

    rootSolve_abserr = 1.0e-10_dp
    rootSolve_relerr = 1.0e-10_dp

    zeta_rootSolve_min = zeta_rootSolve_target - 1.0
    zeta_rootSolve_max = zeta_rootSolve_target + 1.0
    call regcoil_fzero(regcoil_fzero_residual, zeta_rootSolve_min, zeta_rootSolve_max, zeta_rootSolve_target, &
         rootSolve_relerr, rootSolve_abserr, fzeroFlag)
    ! Note: fzero returns its answer in zeta_rootSolve_min
    zeta_plasma_rootSolveSolution = zeta_rootSolve_min
    if (fzeroFlag == 4) then
       stop "ERROR: fzero returned error 4: no sign change in residual"
    else if (fzeroFlag > 2) then
       print *,"WARNING in cosm: fzero returned an error code:",fzeroFlag
    end if

    call regcoil_expand_plasma_surface(prob, theta_rootSolve, zeta_plasma_rootSolveSolution, separation, &
         x_offsetSurface, y_offsetSurface, z_offsetSurface)

  contains

    function regcoil_fzero_residual(zeta_plasma_test)

      use stel_constants

      implicit none

      real(dp) :: zeta_plasma_test, regcoil_fzero_residual
      real(dp) :: x_outer, y_outer, z_outer, zeta_outer_new, zeta_error

      call regcoil_expand_plasma_surface(prob, theta_rootSolve, zeta_plasma_test, separation, &
           x_outer, y_outer, z_outer)
      zeta_outer_new = atan2(y_outer,x_outer)
      zeta_error = zeta_outer_new - zeta_rootSolve_target
      if (zeta_error < -pi) then
         zeta_error = zeta_error + twopi
      end if
      if (zeta_error > pi) then
         zeta_error = zeta_error - twopi
      end if
      regcoil_fzero_residual = zeta_error

    end function regcoil_fzero_residual

  end subroutine regcoil_compute_offset_surface_xyz_of_thetazeta

end module regcoil_compute_offset_surface_mod
