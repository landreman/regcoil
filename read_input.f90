subroutine read_input

  use global_variables

  implicit none

  integer :: numargs
  character(len=200) :: inputFilename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999

  namelist / regcoil / ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, &
       geometry_option_plasma, geometry_option_coil, &
       R0_plasma, R0_coil, a_plasma, a_coil, &
       separation, wout_filename, &
       save_level, nfp_imposed, symmetry_option, &
       mpol_coil, ntor_coil, &
       nescin_filename, efit_filename, efit_psiN, efit_num_modes, &
       mpol_transform_refinement, ntor_transform_refinement, &
       net_poloidal_current_Amperes, net_toroidal_current_Amperes, &
       load_bnorm, bnorm_filename, &
       shape_filename_plasma, nlambda, lambda_min, lambda_max, general_option, nescout_filename, &
       target_option, current_density_target, lambda_search_tolerance

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named regcoil_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (inputFilename(1:11) .ne. "regcoil_in.") then
     stop "Input file must be named regcoil_in.XXX for some extension XXX"
  end if

  outputFilename = "regcoil_out" // trim(inputFilename(11:)) // ".nc"

  fileUnit=11
  open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     stop
  else
     read(fileUnit, nml=regcoil, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(inputFilename), &
               " but not read data from the regcoil namelist in it."
        if (didFileAccessWork==-1) then
           print *,"Make sure there is a carriage return after the / at the end of the namelist!"
        end if
        stop
     end if
     print *,"Successfully read parameters from regcoil namelist in ", trim(inputFilename), "."
  end if
  close(unit = fileUnit)


  print *,"Resolution parameters:"
  print "(a,i4)","   ntheta_plasma =",ntheta_plasma
  print "(a,i4)","   ntheta_coil   =",ntheta_coil
  print "(a,i4)","   nzeta_plasma  =",nzeta_plasma
  print "(a,i4)","   nzeta_coil    =",nzeta_coil
  print "(a,i4)","   mpol_coil     =",mpol_coil
  print "(a,i4)","   ntor_coil     =",ntor_coil

  select case (symmetry_option)
  case (1)
     print *,"Symmetry: sin(m*theta - n*zeta) modes only"
  case (2)
     print *,"Symmetry: cos(m*theta - n*zeta) modes only"
  case (3)
     print *,"Symmetry: both sin(m*theta - n*zeta) and cos(m*theta - n*zeta) modes"
  case default
     print *,"Error! Invalid setting for symmetry_option: ",symmetry_option
     stop
  end select

end subroutine read_input
