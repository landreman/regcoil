! This subroutine is used only when REGCOIL is called from STELLOPT,
! and it not used when REGCOIL is run as a standalone code.

subroutine regcoil_write_input(proc_string,iunit,istat)
  ! proc_string is normally located in stellopt_runtime.f90
  
  ! istat is not used for anything. Why include it?

  integer, intent(in) :: istat
  character(256), intent(in) ::  proc_string
  character(256) ::  wout_filename_new
  integer ::  iunit
  character(LEN=*), parameter :: OUTBOO  = "(2X,A,1X,'=',1X,L1)"
  character(LEN=*), parameter :: OUTINT  = "(2X,A,1X,'=',1X,I0)"
  character(LEN=*), parameter :: OUTFLT  = "(2X,A,1X,'=',1X,E22.14)"
  character(LEN=*), parameter :: OUTSTR  = "(2X,A,1X,'=',1X,'''',A,'''')"

  wout_filename_new = 'wout_'//TRIM(proc_string)//'.nc'

  write(iunit, '(A)') '&regcoil_nml'
  write(iunit, OUTINT) 'general_option', general_option
  write(iunit, OUTINT) 'nlambda', nlambda
  write(iunit, OUTINT) 'ntheta_plasma', ntheta_plasma
  write(iunit, OUTINT) 'nzeta_plasma', nzeta_plasma
  write(iunit, OUTINT) 'ntheta_coil', ntheta_coil
  write(iunit, OUTINT) 'nzeta_coil', nzeta_coil
  write(iunit, OUTINT) 'mpol_potential', mpol_potential
  write(iunit, OUTINT) 'ntor_potential', ntor_potential
  write(iunit, OUTINT) 'target_option', target_option
  !write(iunit, OUTFLT) 'current_density_target', current_density_target
  write(iunit, OUTINT) 'geometry_option_plasma', geometry_option_plasma
  write(iunit, OUTSTR) 'wout_filename', trim(wout_filename_new)
  write(iunit, OUTINT) 'geometry_option_coil', geometry_option_coil
  write(iunit, OUTFLT) 'separation', separation
  write(iunit, OUTINT) 'symmetry_option', symmetry_option
  write(iunit, OUTFLT) 'net_poloidal_current_Amperes', net_poloidal_current_Amperes
  write(iunit, OUTFLT) 'net_toroidal_current_Amperes', net_toroidal_current_Amperes
  write(iunit, '(A)') '/'
!  write(iunit, '(A)') ''
  return

end subroutine regcoil_write_input
