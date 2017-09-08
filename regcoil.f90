




! Main program

program regcoil

  use global_variables, only: totalTime, outputFilename, general_option, sensitivity_option, normal_displacement_option, fixed_norm_sensitivity_option, exit_code
  use init_plasma_mod

  implicit none

  integer :: tic, toc, countrate

  print *,"This is REGCOIL,"
  print *,"a regularized least-squares method for computing stellarator coil shapes."
  call system_clock(tic,countrate)
  call read_input()
  call validate_input()
  call compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call init_plasma()
  call init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call read_bnorm()

  !Initialize sensitivity arrays
  if (sensitivity_option > 1) then
    select case (general_option)
    case (1,4,5)
      print *,"Initializing sensitivity."
      call init_sensitivity()
    end select
  endif

  call build_matrices()

  select case (general_option)
  case (1)
     call solve()
     if (sensitivity_option > 1) then
      call adjoint_solve()
     end if
  case (2)
     call compute_diagnostics_for_nescout_potential()
  case (3)
     call svd_scan()
  case (4,5)
     call auto_regularization_solve()
     if (sensitivity_option > 1 .and. exit_code == 0) then
       call adjoint_solve()
     end if
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  if (sensitivity_option > 1 .and. exit_code == 0) then
    if (normal_displacement_option == 1) then
      call normal_displacement()
      print *,"Normal displacement calculations complete."
    else if (normal_displacement_option == 2) then
      call normal_displacement_svd()
      print *,"Normal displacement calculations complete (with SVD)."
    endif
    if (fixed_norm_sensitivity_option > 1) then
      call lse_sensitivity()
    end if
  endif

  call system_clock(toc)
  totalTime = real(toc-tic)/countrate

  call write_output()

!  if (sensitivity_option > 1) then
!    call free_sensitivity()
!  endif

  print *,"REGCOIL complete. Total time=",totalTime,"sec."
  print *,"You can run regcoilPlot ",trim(outputFilename)," to plot results."

end program regcoil
