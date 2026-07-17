subroutine regcoil_compute_diagnostics_for_nescout_potential

  use regcoil_variables
  use safe_open_mod
  use stel_constants
  use stel_kinds

  implicit none

  integer :: iflag, tic, toc, countrate
  integer :: ilambda

  integer :: iunit = 7, temp1, status, m, n, mm, nn, index, istat
  character(300) :: myline
  real(dp) :: curpol_nescout, temp2, temp3, amplitude
  character(*), parameter :: matchString_curpol = "np, iota_edge, phip_edge, curpol"
  character(*), parameter :: matchString_phi = "---- Phi(m,n) for"
  character(*), parameter :: matchString_phi2 = "---- end Phi(m,n)"


  matrix = matrix_B
  RHS    =    RHS_B
  ilambda=0

  print *,"Opening nescout file",nescout_filename
  call safe_open(iunit, istat, trim(nescout_filename), 'old', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescout file'
  endif


  ! Skip down to the line in the file that begins with matchString
  do
     read (iunit,"(a)") myline
     if (myline(:len(matchString_curpol)) == matchString_curpol) then
        exit
     end if
  end do
  ! Read in curpol
  read (iunit, *) temp1, temp2, temp3, curpol_nescout
  print *,"Read curpol=",curpol_nescout

  outerLoop: do
     ! Skip down to the line in the file that begins with matchString
     innerLoop: do
        read (iunit,"(a)",iostat=status) myline
        if (status<0) exit outerLoop
        if (myline(:len(matchString_phi)) == matchString_phi) then
           exit innerLoop
        end if
     end do innerLoop
     ! If we make it here, then we have found a current potential.
     print *,"Found a current potential in the nescout file."
     ilambda = ilambda + 1
     solution = 0
     do mm = 0,mpol_potential
        do nn = -ntor_potential, ntor_potential
           read (iunit, *) m, n, amplitude
           if ((m .ne. mm) .or. (n .ne. nn)) then
              print *,"Something weird happened:"
              print *,"regcoil m=",mm
              print *,"nescoil m=",m
              print *,"regcoil n=",nn
              print *,"nescoil n=",n
              print *,"Are you sure the mpol and ntor for regcoil match those from nescoil?"
              stop
           end if
           if (mm==0 .and. nn>0) then
              solution(nn) = amplitude*(-2) ! Need *2 here since nescoil prints n<0 values for m=0.
           elseif (mm>0) then
              ! Nescoil convention is m*u+n*v
              ! Regcoil convention is m*u-n*v
              ! so need to swap sign of n.
              index = ntor_potential + (mm-1)*(ntor_potential*2+1) + ntor_potential - nn + 1
              if ((xm_potential(index) .ne. m) .or. (xn_potential(index) .ne. -nfp*n)) then
                 print *,"Indexing error:"
                 print *,"index=",index
                 print *,"m =",m
                 print *,"xm_potential=",xm_potential(index)
                 print *,"n =",n
                 print *,"xn_potential=",xn_potential(index)
                 print *,"Are you sure the mpol and ntor for regcoil match those from nescoil?"
                 stop
              end if
              solution(index) = amplitude
           end if
        end do
     end do
     ! Convert from nescoil normalization to regcoil normalization:
     solution = solution * curpol_nescout / mu0
     ! Verify next line in the nescout file is what we expect:
     read (iunit,"(a)",iostat=status) myline
     if (myline(:len(matchString_phi2)) .ne. matchString_phi2) then
        print *,"Error! In the nescout file I expected the line"
        print *,matchString_phi2
        print *,"but instead I found the line"
        print *,myline
        stop
     end if
     print *,"Successfully read current potential"

     call regcoil_diagnostics(ilambda)

  end do outerLoop


end subroutine regcoil_compute_diagnostics_for_nescout_potential
