! Main program

program regcoil_main

  use regcoil_variables, only: total_time, outputFilename, general_option
  use regcoil_init_plasma
  use regcoil_input_mod
  use regcoil_validate_input
  use regcoil_compute_lambda
  use regcoil_init_coil_surface
  use regcoil_read_bnorm
  use regcoil_build_matrices
  use regcoil_auto_regularization_solve
  use regcoil_write_output

  implicit none

  integer :: tic, toc, countrate

  print *,"This is REGCOIL,"
  print *,"a regularized least-squares method for computing stellarator coil shapes."
  call system_clock(tic,countrate)

  call read_regcoil_cmd_input()
  call validate_input()
  call compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call init_plasma()
  call init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call read_bnorm()
  call build_matrices()

  select case (general_option)
  case (1)
     call solve()
  case (2)
     call compute_diagnostics_for_nescout_potential()
  case (3)
     call svd_scan()
  case (4,5)
     call auto_regularization_solve()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  call system_clock(toc)
  total_time = real(toc-tic)/countrate

  call write_output()
 
  print *,"REGCOIL complete. Total time=",total_time,"sec."
  print *,"You can run regcoilPlot ",trim(outputFilename)," to plot results."

end program regcoil_main
