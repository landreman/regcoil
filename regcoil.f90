! Main program

program regcoil

  use regcoil_variables, only: total_time, output_filename, general_option
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
  call regcoil_build_matrices()
  call regcoil_prepare_solve()

  select case (general_option)
  case (1)
     call regcoil_lambda_scan()
  case (2)
     call regcoil_compute_diagnostics_for_nescout_potential()
  case (3)
     call regcoil_svd_scan()
  case (4,5)
     call regcoil_auto_regularization_solve()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  call system_clock(toc)
  total_time = real(toc-tic)/countrate

  call regcoil_write_output()
 
  print *,"REGCOIL complete. Total time=",total_time,"sec."
  print *,"You can run regcoilPlot ",trim(output_filename)," to plot results."

end program regcoil
