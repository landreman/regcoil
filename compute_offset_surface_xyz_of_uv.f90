subroutine compute_offset_surface_xyz_of_uv(u_rootSolve,v_rootSolve_target,x_offsetSurface,y_offsetSurface,z_offsetSurface,separation)

  use stel_kinds

  implicit none

  real(dp), intent(in) :: u_rootSolve, v_rootSolve_target, separation
  real(dp), intent(out) :: x_offsetSurface, y_offsetSurface, z_offsetSurface
  
  integer :: fzeroFlag
  real(dp) :: rootSolve_abserr, rootSolve_relerr, v_rootSolve_min, v_rootSolve_max
  real(dp) :: v_plasma_rootSolveSolution
  
  !rootSolve_abserr = 0
  !rootSolve_relerr = 0
  rootSolve_abserr = 1.0e-10_dp
  rootSolve_relerr = 1.0e-10_dp
  
  v_rootSolve_min = v_rootSolve_target - 0.3
  v_rootSolve_max = v_rootSolve_target + 0.3
  call fzero(fzero_residual, v_rootSolve_min, v_rootSolve_max, v_rootSolve_target, &
       rootSolve_relerr, rootSolve_abserr, fzeroFlag)
  ! Note: fzero returns its answer in v_rootSolve_min
  v_plasma_rootSolveSolution = v_rootSolve_min
  if (fzeroFlag == 4) then
     stop "ERROR: fzero returned error 4: no sign change in residual"
  else if (fzeroFlag > 2) then
     print *,"WARNING: fzero returned an error code:",fzeroFlag
  end if
  
  call expand_plasma_surface(u_rootSolve, v_plasma_rootSolveSolution, separation, x_offsetSurface, y_offsetSurface, z_offsetSurface)
  
contains
  
  !------------------------------------------------------------------------------------
  
  function fzero_residual(v_plasma_test)
    
    use global_variables, only: nfp
    use stel_constants
    
    implicit none
    
    real(dp) :: v_plasma_test, fzero_residual
    real(dp) :: x_outer, y_outer, z_outer, v_outer_new, v_error
    
    call expand_plasma_surface(u_rootSolve, v_plasma_test, separation, x_outer, y_outer, z_outer)
    v_outer_new = atan2(y_outer,x_outer)*nfp/twopi
    v_error = v_outer_new - v_rootSolve_target
    if (v_error < -nfp/2.0_dp) then
       v_error = v_error + nfp
    end if
    if (v_error > nfp/2.0_dp) then
       v_error = v_error - nfp
    end if
    fzero_residual = v_error
    
  end function fzero_residual
  
end subroutine compute_offset_surface_xyz_of_uv
