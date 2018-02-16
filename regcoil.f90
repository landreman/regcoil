! Main program

program regcoil_main

!<<<<<<< HEAD
!  use regcoil_variables, only: totalTime, outputFilename, general_option
  use regcoil_variables, only: outputFilename, general_option
  use init_regcoil_plasma
  use regcoil_input_mod
  use validate_regcoil_input
  use compute_regcoil_lambda
  use init_regcoil_coil_surface
  use read_regcoil_bnorm
  use build_regcoil_matrices
  use regcoil_auto_regularization_solve
  use write_regcoil_output
!=======
!  use global_variables, only: total_time, outputFilename, general_option
!  use init_plasma_mod
!>>>>>>> 5dc77ee6d32e1e0a0b643dd63edb6c91b610f3e4

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
