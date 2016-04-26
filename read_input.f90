subroutine read_input

  use global_variables

  implicit none

  integer :: numargs
  character(len=200) :: inputFilename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999

  namelist / regcoil / nu_plasma, nv_plasma, nu_middle, nv_middle, nu_outer, nv_outer, &
       geometry_option_plasma, geometry_option_middle, geometry_option_outer, &
       R0_plasma, R0_middle, R0_outer, a_plasma, a_middle, a_outer, &
       separation_middle, separation_outer, wout_filename, pseudoinverse_thresholds, &
       save_level, n_singular_vectors_to_save, nfp_imposed, symmetry_option, mode_order, &
       mpol_plasma, ntor_plasma, mpol_middle, ntor_middle, mpol_outer, ntor_outer, &
       nescin_filename_middle, nescin_filename_outer, efit_filename, efit_psiN, efit_num_modes, &
       mpol_transform_refinement, ntor_transform_refinement, &
       basis_option_plasma, basis_option_middle, basis_option_outer, check_orthogonality, transfer_matrix_option, &
       zero_first_transfer_vector_in_overlap, net_poloidal_current_Amperes, load_bnorm, bnorm_filename, &
       shape_filename_plasma, N_alphas, alpha_min, alpha_max

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named regcoil_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (inputFilename(1:12) .ne. "regcoil_in.") then
     stop "Input file must be named regcoil_in.XXX for some extension XXX"
  end if

  outputFilename = "regcoil_out" // trim(inputFilename(12:)) // ".nc"

  pseudoinverse_thresholds = uninitialized
  ! Default: a single threshold
  pseudoinverse_thresholds(1) = 1e-12

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
        stop
     end if
     print *,"Successfully read parameters from regcoil namelist in ", trim(inputFilename), "."
  end if
  close(unit = fileUnit)

  do i=1,nmax_pseudoinverse_thresholds
     if (pseudoinverse_thresholds(i) == uninitialized) then
        n_pseudoinverse_thresholds = i-1
        exit
     end if
  end do

  save_vectors_in_uv_format = (save_level<4)

  if (transfer_matrix_option==2) then
     print *,"Overriding settings for outer surface. It will be the same as the middle surface."
     ! Override all settings for the outer surface. It will be the same as the middle surface.
     nu_outer = nu_middle
     nv_outer = nv_middle
     ! While nu and nv must be the same on the outer and middle surfaces, in principle ntor and mpol coud be different.
     ! For now, though, for simplicity I will force them to be the same.
     mpol_outer = mpol_middle
     ntor_outer = ntor_middle
     geometry_option_outer = geometry_option_middle
     nescin_filename_outer = nescin_filename_middle
     separation_outer = separation_middle
     R0_outer = R0_middle
     a_outer = a_middle
  end if


  print *,"Resolution parameters:"
  print *,"  nu_plasma =",nu_plasma
  print *,"  nu_middle =",nu_middle
  print *,"  nu_outer  =",nu_outer
  print *,"  nv_plasma =",nv_plasma
  print *,"  nv_middle =",nv_middle
  print *,"  nv_outer  =",nv_outer
  print *,"  mpol_plasma =",mpol_plasma
  print *,"  mpol_middle =",mpol_middle
  print *,"  mpol_outer  =",mpol_outer
  print *,"  ntor_plasma =",ntor_plasma
  print *,"  ntor_middle =",ntor_middle
  print *,"  ntor_outer  =",ntor_outer

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
