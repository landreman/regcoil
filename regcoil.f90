! Main program

program regcoil

  use regcoil_variables, only: total_time, output_filename, general_option, &
    sensitivity_option, exit_code, fixed_norm_sensitivity_option
  use regcoil_init_plasma_mod

  implicit none

  integer :: tic, toc, countrate

  print *,"This is REGCOIL,"
  print *,"a regularized least-squares method for computing stellarator coil shapes."
  call system_clock(tic,countrate)

  call regcoil_read_input()
  call regcoil_validate_input()
  call regcoil_compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call regcoil_init_plasma()
  call regcoil_init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call regcoil_read_bnorm()

  !Initialize sensitivity arrays
  if (sensitivity_option > 1) then
    call regcoil_init_sensitivity()
  endif

  call regcoil_build_matrices()
  call regcoil_prepare_solve()

  select case (general_option)
  case (1)
     call regcoil_lambda_scan()
     if (sensitivity_option > 1) then
      call regcoil_adjoint_solve()
     end if
  case (2)
     call regcoil_compute_diagnostics_for_nescout_potential()
  case (3)
     call regcoil_svd_scan()
  case (4,5)
     call regcoil_auto_regularization_solve()
     if (sensitivity_option > 1 .and. exit_code == 0) then
       call regcoil_adjoint_solve()
     end if
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  if (fixed_norm_sensitivity_option .and. sensitivity_option > 1 .and. exit_code == 0) then
    call regcoil_fixed_norm_sensitivity()
  end if

  call system_clock(toc)
  total_time = real(toc-tic)/countrate

  call regcoil_write_output()
 
  print *,"REGCOIL complete. Total time=",total_time,"sec."
  print *,"You can run regcoilPlot ",trim(output_filename)," to plot results."

end program regcoil
