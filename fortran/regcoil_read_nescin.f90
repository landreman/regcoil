subroutine regcoil_read_nescin()

  use regcoil_variables
  use safe_open_mod
  use stel_kinds
  
  implicit none
  
  integer :: iunit = 7, k, istat
  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"


  call safe_open(iunit, istat, trim(nescin_filename), 'old', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file'
  endif


  ! Skip down to the line in the file that begins with matchString
  do
     read (iunit,"(a)") myline
     if (myline(:len(matchString)) == matchString) then
        exit
     end if
  end do
  read (iunit, *)

  read (iunit, *) mnmax_coil
  if (verbose) print *,"  Reading",mnmax_coil,"Fourier modes from nescin"

  allocate(xn_coil(mnmax_coil))
  allocate(xm_coil(mnmax_coil))
  allocate(rmnc_coil(mnmax_coil))
  allocate(rmns_coil(mnmax_coil))
  allocate(zmnc_coil(mnmax_coil))
  allocate(zmns_coil(mnmax_coil))

  read (iunit, *)
  read (iunit, *)
  do k = 1, mnmax_coil
     read (iunit, *) xm_coil(k), xn_coil(k), rmnc_coil(k), zmns_coil(k), rmns_coil(k), zmnc_coil(k)
  end do

  ! nescin format uses cos(m * theta + n * nfp * zeta) whereas we use the vmec convention cos(m * theta - n * zeta), so convert:
  xn_coil = -nfp * xn_coil

  close(iunit)

end subroutine  regcoil_read_nescin
