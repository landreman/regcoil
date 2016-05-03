! Main program

program regcoil

  use global_variables, only: totalTime, outputFilename

  implicit none

  integer :: tic, toc, countrate

  print *,"This is REGCOIL,"
  print *,"a regularized least-squares method for computing stellarator coils."
  call system_clock(tic,countrate)

  call read_input()
  call validate_input()
  call compute_alpha()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call init_plasma()
  call init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call read_bnorm()
  call build_matrices()

  ! Assemble transfer matrix and compute its SVD:
  call solve()

  call system_clock(toc)
  totalTime = real(toc-tic)/countrate

  call write_output()

  print *,"REGCOIL complete. Total time=",totalTime,"sec."
  print *,"You can run regcoilPlot ",trim(outputFilename)," to plot results."

end program regcoil
