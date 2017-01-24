subroutine validate_input

  use global_variables
  use safe_open_mod

  implicit none

  integer :: iunit = 7, istat, j
  character(300) :: myline
  character(*), parameter :: matchString = "---- Phi(m,n) for"
  character(len=*), parameter :: line="******************************************************************"

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

  if (save_level > 3) then
     stop "Error! save_level must be <= 3."
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

  if (nlambda < 1) then
     stop "nlambda must be at least 1."
  end if

  if (lambda_min <= 0) then
     stop "lambda_min must be greater than 0."
  end if

  if (lambda_max < lambda_min) then
     stop "lambda_max must be >= lambda_min."
  end if

  if (general_option<1) then
     stop "general_option must be at least 1."
  end if
  if (general_option>5) then
     stop "general_option must be no more than 5."
  end if

  if (general_option==2) then
     ! Replace nlambda with the number of current potentials saved in the nescout file.
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
     nlambda = j
  end if

  if (current_density_target<=0) then
     stop "current_density_target must be positive."
  end if

  if (current_density_target < 1e5) then
     print *,line
     print *,"Warning! The value of current_density_target you have set"
     print *,"is surprisingly small."
     if (general_option .ne. 5) then
        print *,"It is recommended that you run with general_option=5 to verify that this"
        print *,"value of current_density_target is attainable."
     end if
     print *,line
  end if

  if (current_density_target > 3e8) then
     print *,line
     print *,"Warning! The value of current_density_target you have set"
     print *,"is surprisingly large."
     if (general_option .ne. 5) then
        print *,"It is recommended that you run with general_option=5 to verify that this"
        print *,"value of current_density_target is attainable."
     end if
     print *,line
  end if

end subroutine validate_input
