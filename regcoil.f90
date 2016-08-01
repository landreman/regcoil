! Main program

program regcoil

  use global_variables, only: totalTime, outputFilename, general_option

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

  if (general_option==1) then
     call solve()
  elseif (general_option==2) then
     call compute_diagnostics_for_nescout_potential()
  else
     print *,"Invalid general_option:",general_option
     stop
  end if

  call system_clock(toc)
  totalTime = real(toc-tic)/countrate

  if (general_option==1) then
     call write_output()
  end if

  print *,"REGCOIL complete. Total time=",totalTime,"sec."
  print *,"You can run regcoilPlot ",trim(outputFilename)," to plot results."

end program regcoil
