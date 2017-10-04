!-----------------------------------------------------------------------
!     Module:        regcoil_input_mod
!     Authors:       J.C. Schmitt (jcschmitt@auburn.edu)
!     Date:          08/31/2017
!     Description:   This module provides reading and writing functionality
!                    for the RECOIL input namelist.
!                    Used by REGCOIL and STELLOPT.
!-----------------------------------------------------------------------

module regcoil_input_mod

  use regcoil_variables

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

contains

subroutine read_regcoil_cmd_input
  ! This routine is called by REGCOIL. This is NOT called by external
  ! programs, such as STELLOPT.

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

end subroutine read_regcoil_cmd_input

! This routine is called from external programs (i.e. STELLOPT) to read the
! REGCOIL namelist from a file named "input.filename" where filename is 
! passed in
!subroutine read_regcoilin_input(filename, istat)
!  character(*), intent(in) :: filename
!  integer, intent(out) :: istat
!  logical :: lexist
!  integer :: i, k, iunit
! 
!  inquire(file='input.'//trim(filename), exist=lexist)
!  if (.not. lexist) then
!    istat=-3; return
!  end if
!  call safe_open(iunit,istat,'input.'//trim(filename),'old','formatted')
!  if (istat /= 0) return
!  read(iunit, nml=regcoil, iostat=istat)
!  if (istat /= 0) return
!  close(iunit)
!
!end subroutine read_regcoil_input

! This routine is called from external programs (i.e. STELLOPT) to read the
! REGCOIL namelist from a file associated with iunit.
subroutine read_regcoil_input(iunit, istat)
  integer, intent(out) :: istat
  integer, intent(in) :: iunit

  ntheta_plasma = 0
 
  print "(a,i4)","   iunit =",iunit
  print "(a,i4)","   istat =",istat
  read(iunit, nml=regcoil, iostat=istat)
  print "(a,i4)","   iunit =",iunit
  print "(a,i4)","   istat =",istat
!  if (istat /= 0) then
!    print *,"Error!  I was able to open the file ", trim(inputFilename), &
!            " but not read data from the regcoil namelist in it."
!  end if
!  if (istat==-1) then
!    print *,"Make sure there is a carriage return after the / at the end of the namelist!"
!  end if
!  print *,"Successfully read parameters from regcoil namelist"

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


end subroutine read_regcoil_input


subroutine write_regcoil_input(proc_string,iunit,istat)
  !-----------------------------------------------------------------------
  !     Local Variables
  !        ier         Error flag
  !        iunit       File unit number
  !----------------------------------------------------------------------
  ! proc_string is normally located in stellopt_runtime.f90
  integer, intent(in) :: istat
  character(256), intent(in) ::  proc_string
  character(256) ::  wout_filename_new
  integer ::  ier, iunit
  character(LEN=*), parameter :: OUTBOO  = "(2X,A,1X,'=',1X,L1)"
  character(LEN=*), parameter :: OUTINT  = "(2X,A,1X,'=',1X,I0)"
  character(LEN=*), parameter :: OUTFLT  = "(2X,A,1X,'=',1X,E22.14)"
  character(LEN=*), parameter :: OUTSTR  = "(2X,A,1X,'=',1X,'''',A,'''')"

  ! print "(a,a)", "proc_string = ", proc_string
  wout_filename_new = 'wout_'//TRIM(proc_string)//'.nc'
  ! print "(a,a)", "wout_filename_new = ", wout_filename_new
  ! print "(a,i12)", "iunit=", iunit


  write(iunit, '(A)') '&regcoil'
  write(iunit, OUTSTR) 'wout_filename', trim(wout_filename_new)
  write(iunit, OUTINT) 'ntheta_plasma', ntheta_plasma
  write(iunit, OUTINT) 'nzeta_plasma', nzeta_plasma
  write(iunit, OUTINT) 'ntheta_coil', ntheta_coil
  write(iunit, OUTINT) 'nzeta_coil', nzeta_coil
  write(iunit, OUTINT) 'mpol_coil', mpol_coil
  write(iunit, OUTINT) 'ntor_coil', ntor_coil
  write(iunit, OUTINT) 'general_option', general_option
  write(iunit, OUTINT) 'target_option', target_option
  !write(iunit, OUTFLT) 'current_density_target', current_density_target
  write(iunit, OUTINT) 'geometry_option_plasma', geometry_option_plasma
  write(iunit, OUTINT) 'geometry_option_coil', geometry_option_coil
  write(iunit, OUTFLT) 'separation', separation
  write(iunit, OUTINT) 'nlambda', nlambda
  write(iunit, OUTINT) 'symmetry_option', symmetry_option
  write(iunit, OUTFLT) 'net_poloidal_current_Amperes', net_poloidal_current_Amperes
  write(iunit, OUTFLT) 'net_toroidal_current_Amperes', net_toroidal_current_Amperes
  write(iunit, '(A)') '/'
!  write(iunit, '(A)') ''
  return
end subroutine write_regcoil_input

end module regcoil_input_mod
