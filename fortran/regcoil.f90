! Main program

program regcoil

  use regcoil_variables, only: regcoil_t
  use regcoil_init_plasma_mod

  implicit none

  type(regcoil_t) :: prob
  integer :: tic, toc, countrate

  print *,"This is REGCOIL,"
  print *,"a regularized least-squares method for computing stellarator coil shapes."
  call system_clock(tic,countrate)

  call regcoil_read_input(prob)
  call regcoil_validate_input(prob)
  call regcoil_compute_lambda(prob)

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call regcoil_init_plasma(prob)
  call regcoil_init_coil_surface(prob)

  ! Initialize some of the vectors and matrices needed:
  call regcoil_read_bnorm(prob)

  call regcoil_build_matrices(prob)
  call regcoil_prepare_solve(prob)

  select case (prob%input%general_option)
  case (1)
     call regcoil_lambda_scan(prob)
  case (2)
     call regcoil_compute_diagnostics_for_nescout_potential(prob)
  case (4,5)
     call regcoil_auto_regularization_solve(prob)
  case default
     print *,"Invalid general_option:",prob%input%general_option
     print *,"general_option=3 (SVD scan) has been removed."
     stop
  end select

  call system_clock(toc)
  prob%output%total_time = real(toc-tic)/countrate

  call regcoil_write_output(prob)

  print *,"REGCOIL complete. Total time=",prob%output%total_time,"sec."
  print *,"You can run regcoilPlot ",trim(prob%input%output_filename)," to plot results."

end program regcoil
