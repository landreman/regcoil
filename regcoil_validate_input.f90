subroutine regcoil_validate_input

  use regcoil_variables
  use safe_open_mod

  implicit none

  integer :: iunit = 7, istat, j
  character(300) :: myline
  character(*), parameter :: matchString = "---- Phi(m,n) for"
  character(len=*), parameter :: line="******************************************************************"
  real(dp) :: typical_target_min, typical_target_max

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


  if (mpol_potential < 0) then
     stop "Error! mpol_potential must be >= 0."
  end if


  if (ntor_potential < 0) then
     stop "Error! ntor_potential must be >= 0."
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

  if (geometry_option_plasma > 7) then                  ! czhu changed this to option 7 for FOCUS
     stop "Error! geometry_option_plasma must be <= 7." 
  end if


  if (geometry_option_coil < 0) then
     stop "Error! geometry_option_coil must be >= 0."
  end if

  if (geometry_option_coil > 4) then
     stop "Error! geometry_option_coil must be <= 4."
  end if

  if (separation < 0) then
     stop "Error! separation must be >= 0."
  end if



  if (load_bnorm) then
     ! To use bnorm data, we must have a VMEC file to set the overall normalization
     select case (geometry_option_plasma)
     case (2,3,4,7)
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

  if ((general_option==2 .or. general_option==3) .and. sensitivity_option>1) then
    stop "sensitivity_option>1 must be used with general_option = 1, 4, or 5."
  end if

  if (general_option==2) then
     ! Replace nlambda with the number of current potentials saved in the nescout file.
     if (verbose) print *,"Opening nescout file",nescout_filename
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
     if (verbose) print *,"Detected",j,"current potentials in the nescout file."
     nlambda = j
  end if

  if (target_value<=0) then
     stop "target_value must be positive."
  end if

  if (general_option == 4 ) then
     print *,"It is recommended that you run with general_option=5 instead of 4"
     print *,"to verify that this value of target_value is attainable."
  end if

  if (general_option==4 .or. general_option==5) then
     select case (trim(target_option))
     case (target_option_max_K,target_option_rms_K,target_option_max_K_lse,target_option_lp_norm_K)
        typical_target_min = 1e5
        typical_target_max = 3e8
     case (target_option_chi2_K)
        typical_target_min = 1e14
        typical_target_max = 1e17
     case (target_option_max_Bnormal,target_option_rms_Bnormal)
        typical_target_min = 1e-5
        typical_target_max = 5
     case (target_option_chi2_B)
        typical_target_min = 1e-6
        typical_target_max = 100
     case default
        print *,"Invalid target_option: ",target_option
        stop
     end select

     if (target_value < typical_target_min) then
        print "(a)",line
        print "(a)","Warning! The value of target_value you have set is surprisingly small."
        print "(a)",line
     end if
     if (target_value > typical_target_max) then
        print "(a)",line
        print "(a)","Warning! The value of target_value you have set is surprisingly large."
        print "(a)",line
     end if
  end if

  select case (trim(regularization_term_option))
  case (regularization_term_option_chi2_K)
  case (regularization_term_option_Laplace_Beltrami)
    if (symmetry_option > 1) then
      stop "Error! sensitivity_option > 1 can only be used with regularization_term_option = chi2_K."
    end if
  case (regularization_term_option_K_xy)
    if (symmetry_option > 1) then
      stop "Error! sensitivity_option > 1 can only be used with regularization_term_option = chi2_K."
    end if
  case (regularization_term_option_K_zeta)
    if (symmetry_option > 1) then
      stop "Error! sensitivity_option > 1 can only be used with regularization_term_option = chi2_K."
    end if
  case default
     print *,"Error! Unrecognized regularization_term_option: ",trim(regularization_term_option)
     stop
  end select

  if (sensitivity_option>5 .or. sensitivity_option<1) then
    stop "sensitivity_option must be 1, 2, 3, 4, or 5."
  end if

  if (nmax_sensitivity<1) then
    stop "nmax_sensitivity must be >=1."
  end if

  if (mmax_sensitivity<1) then
    stop "mmax_sensitivity must be >=1."
  end if

  if (sensitivity_symmetry_option<1 .or. sensitivity_symmetry_option>2) then
    stop "sensitivity_symmetry_option must be 1 or 2."
  end if

  if (target_option_p <1) then
    stop "target_option_p must be >=1."
  end if

  if (coil_plasma_dist_lse_p <=1) then
    stop "coil_plasma_dist_lse_p must be > 1."
  end if

  if ((general_option==4 .or. general_option==5) .and. fixed_norm_sensitivity_option) then
     select case (trim(target_option))
       case (target_option_max_K_lse,target_option_lp_norm_K,target_option_chi2_B)
       case default
      print *,"fixed_norm_sensitivity_option must be used with target_option = 'max_K_lse', 'lp_norm_K', or 'chi2_B'"
        stop
     end select
  end if

end subroutine regcoil_validate_input
