subroutine validate_input

  use global_variables
  use safe_open_mod

  implicit none

  integer :: iunit = 7, istat, j
  character(300) :: myline
  character(*), parameter :: matchString = "---- Phi(m,n) for"

  if (ntheta_plasma < 1) then
     stop "Error! ntheta_plasma must be >= 1."
  end if

  if (ntheta_coil < 1) then
     stop "Error! ntheta_coil must be >= 1."
  end if


  if (nzeta_plasma < 1) then
     stop "Error! nzeta_plasma must be >= 1."
  end if

  if (nzeta_coil < 1) then
     stop "Error! nzeta_coil must be >= 1."
  end if


  if (mpol_coil < 0) then
     stop "Error! mpol_coil must be >= 0."
  end if


  if (ntor_coil < 0) then
     stop "Error! ntor_coil must be >= 0."
  end if


  if (save_level < 0) then
     stop "Error! save_level must be >= 0."
  end if

  if (save_level > 2) then
     stop "Error! save_level must be <= 2."
  end if


  if (symmetry_option < 1) then
     stop "Error! symmetry_option must be >= 1."
  end if

  if (symmetry_option > 3) then
     stop "Error! symmetry_option must be <= 3."
  end if


  if (geometry_option_plasma < 0) then
     stop "Error! geometry_option_plasma must be >= 0."
  end if

  if (geometry_option_plasma > 6) then
     stop "Error! geometry_option_plasma must be <= 6."
  end if


  if (geometry_option_coil < 0) then
     stop "Error! geometry_option_coil must be >= 0."
  end if

  if (geometry_option_coil > 3) then
     stop "Error! geometry_option_coil must be <= 3."
  end if



  if (separation < 0) then
     stop "Error! separation must be >= 0."
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

  if (nalpha < 1) then
     stop "nalpha must be at least 1."
  end if

  if (alpha_min <= 0) then
     stop "alpha_min must be greater than 0."
  end if

  if (alpha_max < alpha_min) then
     stop "alpha_max must be >= alpha_min."
  end if

  if (general_option<1) then
     stop "general_option must be at least 1."
  end if
  if (general_option>2) then
     stop "general_option must be no more than 2."
  end if

  if (general_option==2) then
     ! Replace nalpha with the number of current potentials saved in the nescout file.
     print *,"Opening nescout file",nescout_filename
     call safe_open(iunit, istat, trim(nescout_filename), 'old', 'formatted')
     if (istat .ne. 0) then
        stop 'Error opening nescout file'
     endif
     j = 0
     do
        read (iunit,"(a)",iostat=istat) myline
        if (istat<0) exit
        if (myline(:len(matchString)) == matchString) then
           j = j + 1
        end if
     end do
     print *,"Detected",j,"current potentials in the nescout file."
     nalpha = j
  end if

end subroutine validate_input
