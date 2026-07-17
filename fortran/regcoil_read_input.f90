subroutine regcoil_read_input(prob)

  use regcoil_variables, only: regcoil_t

  implicit none

  type(regcoil_t), intent(inout) :: prob
  integer :: numargs
  character(len=200) :: inputFilename
  integer :: ios

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named regcoil_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if

  call regcoil_read_input_file(prob, trim(inputFilename), ios)
  if (ios /= 0) stop

end subroutine regcoil_read_input


subroutine regcoil_read_input_file(prob, inputFilename, ios)

  use regcoil_variables, only: regcoil_t
  use stel_kinds

  implicit none

  type(regcoil_t), intent(inout) :: prob
  character(len=*), intent(in) :: inputFilename
  integer, intent(out) :: ios
  integer :: fileUnit, didFileAccessWork
  character(len=200) :: linebuf
  integer :: slashpos
  character(len=200) :: basename

  ! Locals for namelist I/O (copied into prob after read).
  integer :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil
  integer :: geometry_option_plasma, geometry_option_coil
  real(dp) :: R0_plasma, R0_coil, a_plasma, a_coil, separation
  character(len=200) :: wout_filename, nescin_filename, shape_filename_plasma
  character(len=200) :: nescout_filename, bnorm_filename, regularization_term_option, target_option
  integer :: save_level, nfp_imposed, symmetry_option
  integer :: mpol_potential, ntor_potential, mpol_coil_filter, ntor_coil_filter
  integer :: max_mpol_coil, max_ntor_coil, nlambda, general_option
  real(dp) :: mpol_transform_refinement, ntor_transform_refinement
  real(dp) :: net_poloidal_current_Amperes, net_toroidal_current_Amperes
  real(dp) :: lambda_min, lambda_max, target_value, lambda_search_tolerance, target_option_p
  logical :: load_bnorm, verbose

  namelist / regcoil_nml / ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, &
       geometry_option_plasma, geometry_option_coil, &
       R0_plasma, R0_coil, a_plasma, a_coil, &
       separation, wout_filename, &
       save_level, nfp_imposed, symmetry_option, &
       mpol_potential, ntor_potential, mpol_coil_filter, ntor_coil_filter, &
       nescin_filename, &
       mpol_transform_refinement, ntor_transform_refinement, max_mpol_coil, max_ntor_coil, &
       net_poloidal_current_Amperes, net_toroidal_current_Amperes, &
       load_bnorm, bnorm_filename, &
       shape_filename_plasma, nlambda, lambda_min, lambda_max, general_option, &
       regularization_term_option, verbose, nescout_filename, &
       target_option, target_value, lambda_search_tolerance, target_option_p

  ios = 0

  ! Seed locals from instance defaults.
  ntheta_plasma = prob%plasma%ntheta_plasma
  nzeta_plasma = prob%plasma%nzeta_plasma
  ntheta_coil = prob%coil%ntheta_coil
  nzeta_coil = prob%coil%nzeta_coil
  geometry_option_plasma = prob%plasma%geometry_option_plasma
  geometry_option_coil = prob%coil%geometry_option_coil
  R0_plasma = prob%plasma%R0_plasma
  R0_coil = prob%coil%R0_coil
  a_plasma = prob%plasma%a_plasma
  a_coil = prob%coil%a_coil
  separation = prob%coil%separation
  wout_filename = prob%plasma%wout_filename
  shape_filename_plasma = prob%plasma%shape_filename_plasma
  nescin_filename = prob%coil%nescin_filename
  nescout_filename = prob%coil%nescout_filename
  save_level = prob%input%save_level
  nfp_imposed = prob%input%nfp_imposed
  symmetry_option = prob%input%symmetry_option
  mpol_potential = prob%input%mpol_potential
  ntor_potential = prob%input%ntor_potential
  mpol_coil_filter = prob%coil%mpol_coil_filter
  ntor_coil_filter = prob%coil%ntor_coil_filter
  mpol_transform_refinement = prob%plasma%mpol_transform_refinement
  ntor_transform_refinement = prob%plasma%ntor_transform_refinement
  max_mpol_coil = prob%coil%max_mpol_coil
  max_ntor_coil = prob%coil%max_ntor_coil
  net_poloidal_current_Amperes = prob%input%net_poloidal_current_Amperes
  net_toroidal_current_Amperes = prob%input%net_toroidal_current_Amperes
  load_bnorm = prob%input%load_bnorm
  bnorm_filename = prob%input%bnorm_filename
  nlambda = prob%input%nlambda
  lambda_min = prob%input%lambda_min
  lambda_max = prob%input%lambda_max
  general_option = prob%input%general_option
  regularization_term_option = prob%input%regularization_term_option
  verbose = prob%input%verbose
  target_option = prob%input%target_option
  target_value = prob%input%target_value
  lambda_search_tolerance = prob%input%lambda_search_tolerance
  target_option_p = prob%input%target_option_p

  ! Prefer legacy naming when the basename starts with regcoil_in.
  slashpos = index(inputFilename, '/', back=.true.)
  if (slashpos == 0) then
     basename = trim(inputFilename)
  else
     basename = trim(inputFilename(slashpos+1:))
  end if
  if (basename(1:11) == "regcoil_in.") then
     prob%input%output_filename = "regcoil_out" // trim(basename(11:)) // ".nc"
  else
     prob%input%output_filename = "regcoil_out.nc"
  end if

  fileUnit=11
  open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     ios = 1
     return
  end if

  read(fileUnit, nml=regcoil_nml, iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error!  I was able to open the file ", trim(inputFilename), &
            " but not read data from the regcoil_nml namelist in it."
     backspace(fileUnit)
     read(fileUnit, fmt="(A)") linebuf
     print *, "Invalid line in namelist: ", trim(linebuf)
     print *, "(Error may be earlier)"
     if (didFileAccessWork==-1) then
        print *,"Make sure there is a carriage return after the / at the end of the namelist!"
     end if
     close(unit = fileUnit)
     ios = 2
     return
  end if
  if (verbose) print *,"Successfully read parameters from regcoil_nml namelist in ", trim(inputFilename), "."
  close(unit = fileUnit)

  ! Copy namelist locals into the instance.
  prob%plasma%ntheta_plasma = ntheta_plasma
  prob%plasma%nzeta_plasma = nzeta_plasma
  prob%coil%ntheta_coil = ntheta_coil
  prob%coil%nzeta_coil = nzeta_coil
  prob%plasma%geometry_option_plasma = geometry_option_plasma
  prob%coil%geometry_option_coil = geometry_option_coil
  prob%plasma%R0_plasma = R0_plasma
  prob%coil%R0_coil = R0_coil
  prob%plasma%a_plasma = a_plasma
  prob%coil%a_coil = a_coil
  prob%coil%separation = separation
  prob%plasma%wout_filename = wout_filename
  prob%plasma%shape_filename_plasma = shape_filename_plasma
  prob%coil%nescin_filename = nescin_filename
  prob%coil%nescout_filename = nescout_filename
  prob%input%save_level = save_level
  prob%input%nfp_imposed = nfp_imposed
  prob%input%symmetry_option = symmetry_option
  prob%input%mpol_potential = mpol_potential
  prob%input%ntor_potential = ntor_potential
  prob%coil%mpol_coil_filter = mpol_coil_filter
  prob%coil%ntor_coil_filter = ntor_coil_filter
  prob%plasma%mpol_transform_refinement = mpol_transform_refinement
  prob%plasma%ntor_transform_refinement = ntor_transform_refinement
  prob%coil%max_mpol_coil = max_mpol_coil
  prob%coil%max_ntor_coil = max_ntor_coil
  prob%input%net_poloidal_current_Amperes = net_poloidal_current_Amperes
  prob%input%net_toroidal_current_Amperes = net_toroidal_current_Amperes
  prob%input%load_bnorm = load_bnorm
  prob%input%bnorm_filename = bnorm_filename
  prob%input%nlambda = nlambda
  prob%input%lambda_min = lambda_min
  prob%input%lambda_max = lambda_max
  prob%input%general_option = general_option
  prob%input%regularization_term_option = regularization_term_option
  prob%input%verbose = verbose
  prob%input%target_option = target_option
  prob%input%target_value = target_value
  prob%input%lambda_search_tolerance = lambda_search_tolerance
  prob%input%target_option_p = target_option_p

  if (verbose) then
     print *,"Resolution parameters:"
     print "(a,i5)","   ntheta_plasma  =",ntheta_plasma
     print "(a,i5)","   ntheta_coil    =",ntheta_coil
     print "(a,i5)","   nzeta_plasma   =",nzeta_plasma
     print "(a,i5)","   nzeta_coil     =",nzeta_coil
     print "(a,i5)","   mpol_potential =",mpol_potential
     print "(a,i5)","   ntor_potential =",ntor_potential

     select case (symmetry_option)
     case (1)
        print *,"Symmetry: sin(m*theta - n*zeta) modes only"
     case (2)
        print *,"Symmetry: cos(m*theta - n*zeta) modes only"
     case (3)
        print *,"Symmetry: both sin(m*theta - n*zeta) and cos(m*theta - n*zeta) modes"
     case default
        print *,"Error! Invalid setting for symmetry_option: ",symmetry_option
        ios = 3
        return
     end select
  end if

end subroutine regcoil_read_input_file
