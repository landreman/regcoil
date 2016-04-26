subroutine validate_input

  use global_variables

  implicit none

  if (transfer_matrix_option < 1 .or. transfer_matrix_option > 2) then
     stop "Error! transfer_matrix_option must be 1 or 2."
  end if


  if (nu_plasma < 1) then
     stop "Error! nu_plasma must be >= 1."
  end if

  if (nu_middle < 1) then
     stop "Error! nu_middle must be >= 1."
  end if

  if (nu_outer < 1) then
     stop "Error! nu_outer must be >= 1."
  end if


  if (nv_plasma < 1) then
     stop "Error! nv_plasma must be >= 1."
  end if

  if (nv_middle < 1) then
     stop "Error! nv_middle must be >= 1."
  end if

  if (nv_outer < 1) then
     stop "Error! nv_outer must be >= 1."
  end if


  if (mpol_plasma < 0) then
     stop "Error! mpol_plasma must be >= 0."
  end if

  if (mpol_middle < 0) then
     stop "Error! mpol_middle must be >= 0."
  end if

  if (mpol_outer < 0) then
     stop "Error! mpol_outer must be >= 0."
  end if


  if (ntor_plasma < 0) then
     stop "Error! ntor_plasma must be >= 0."
  end if

  if (ntor_middle < 0) then
     stop "Error! ntor_middle must be >= 0."
  end if

  if (ntor_outer < 0) then
     stop "Error! ntor_outer must be >= 0."
  end if



  if (save_level < 0) then
     stop "Error! save_level must be >= 0."
  end if

  if (save_level > 4) then
     stop "Error! save_level must be <= 4."
  end if


  if (symmetry_option < 1) then
     stop "Error! symmetry_option must be >= 1."
  end if

  if (symmetry_option > 3) then
     stop "Error! symmetry_option must be <= 3."
  end if


  if (basis_option_plasma < 1) then
     stop "Error! basis_option_plasma must be >= 1."
  end if

  if (basis_option_plasma > 3) then
     stop "Error! basis_option_plasma must be <= 3."
  end if

  if (basis_option_middle < 1) then
     stop "Error! basis_option_middle must be >= 1."
  end if

  if (basis_option_middle > 3) then
     stop "Error! basis_option_middle must be <= 3."
  end if

  if (basis_option_outer < 1) then
     stop "Error! basis_option_outer must be >= 1."
  end if

  if (basis_option_outer > 3) then
     stop "Error! basis_option_outer must be <= 3."
  end if


  if (mode_order < 1) then
     stop "Error! mode_order must be >= 1."
  end if

  if (mode_order > 2) then
     stop "Error! mode_order must be <= 2."
  end if


  if (geometry_option_plasma < 0) then
     stop "Error! geometry_option_plasma must be >= 0."
  end if

  if (geometry_option_plasma > 6) then
     stop "Error! geometry_option_plasma must be <= 6."
  end if


  if (geometry_option_middle < 0) then
     stop "Error! geometry_option_middle must be >= 0."
  end if

  if (geometry_option_middle > 3) then
     stop "Error! geometry_option_middle must be <= 3."
  end if


  if (geometry_option_outer < 0) then
     stop "Error! geometry_option_outer must be >= 0."
  end if

  if (geometry_option_outer > 4) then
     stop "Error! geometry_option_outer must be <= 4."
  end if

  if (geometry_option_outer == 4 .and. geometry_option_middle .ne. 3) then
     stop "Error! If geometry_option_outer == 4 (offset from the middle nescin file), then you must have geometry_option_middle = 3 (nescin)."
  end if


  if (separation_middle < 0) then
     stop "Error! separation_middle must be >= 0."
  end if

  if (separation_outer <= 0 .and. geometry_option_outer==2) then
     stop "Error! For geometry_option_outer==2, separation_outer must be > 0."
  end if

  if ((separation_outer < separation_middle) .and. (geometry_option_outer .ne. 4)) then
     stop "Error! separation_outer must be >= separation_middle."
  end if


  if (load_bnorm) then
     ! To use bnorm data, we must have a VMEC file to set the overall normalization
     select case (geometry_option_plasma)
     case (2,3,4)
        ! Yes, we have a VMEC file available.
     case (0,1,5)
        stop "Error! If load_bnorm=.t., the plasma surface must come from a vmec wout file."
     case default
        stop "Error! Invalid geometry_option_plasma"
     end select
  end if

end subroutine validate_input
