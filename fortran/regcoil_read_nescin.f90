subroutine regcoil_read_nescin(prob)

  use regcoil_variables, only: regcoil_t
  use safe_open_mod
  use stel_kinds
  
  implicit none
  

  type(regcoil_t), intent(inout) :: prob
  integer :: iunit = 7, k, istat
  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"


    associate ( &
       nfp => prob%plasma%nfp, &
       nescin_filename => prob%coil%nescin_filename, &
       mnmax_coil => prob%coil%mnmax_coil, &
       verbose => prob%input%verbose &
       )
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

  allocate(prob%coil%xn_coil(mnmax_coil))
  allocate(prob%coil%xm_coil(mnmax_coil))
  allocate(prob%coil%rmnc_coil(mnmax_coil))
  allocate(prob%coil%rmns_coil(mnmax_coil))
  allocate(prob%coil%zmnc_coil(mnmax_coil))
  allocate(prob%coil%zmns_coil(mnmax_coil))

  read (iunit, *)
  read (iunit, *)
  do k = 1, mnmax_coil
     read (iunit, *) prob%coil%xm_coil(k), prob%coil%xn_coil(k), prob%coil%rmnc_coil(k), prob%coil%zmns_coil(k), prob%coil%rmns_coil(k), prob%coil%zmnc_coil(k)
  end do

  ! nescin format uses cos(m * theta + n * nfp * zeta) whereas we use the vmec convention cos(m * theta - n * zeta), so convert:
  prob%coil%xn_coil = -nfp * prob%coil%xn_coil

  close(iunit)


  end associate
end subroutine  regcoil_read_nescin
