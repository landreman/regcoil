subroutine compute_diagnostics_for_nescout_potential

  use global_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer :: iflag, tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: matrix, this_current_potential
  real(dp), dimension(:), allocatable :: RHS, solution
  real(dp), dimension(:), allocatable :: JDifference_x, JDifference_y, Jdifference_z
  real(dp), dimension(:,:), allocatable :: this_J2_times_N
  real(dp) :: factor_theta, factor_zeta
  integer :: ialpha, itheta, izeta

  integer :: iunit = 7, temp1, status, m, n, mm, nn, index, istat
  character(300) :: myline
  real(dp) :: curpol_nescout, temp2, temp3, amplitude
  character(*), parameter :: matchString_curpol = "np, iota_edge, phip_edge, curpol"
  character(*), parameter :: matchString_phi = "---- Phi(m,n) for"
  character(*), parameter :: matchString_phi2 = "---- end Phi(m,n)"



  allocate(matrix(num_basis_functions, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(chi2_B(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(chi2_J(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(current_potential(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_mn(num_basis_functions,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_total(ntheta_plasma,nzeta_plasma,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(J2(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_J2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  matrix = matrix_B
  RHS    =    RHS_B
  ialpha=1

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
     solution = 0
     do mm = 0,mpol_coil
        do nn = -ntor_coil, ntor_coil
           read (iunit, *) m, n, amplitude
           if ((m .ne. mm) .or. (n .ne. nn)) then
              print *,"Something weird happened:"
              print *,"mm=",mm
              print *,"m =",m
              print *,"nn=",nn
              print *,"n =",n
              stop
           end if
           if (mm==0 .and. nn>0) then
              solution(nn) = amplitude*2 ! Need *2 here since nescoil prints n<0 values for m=0.
           elseif (mm>0) then
              ! Nescoil convention is m*u+n*v
              ! Regcoil convention is m*u-n*v
              ! so need to swap sign of n.
              index = ntor_coil + mm*(ntor_coil*2+1) + ntor_coil - nn + 1
              if ((xm(index) .ne. m) .or. (xn(index) .ne. nfp*n)) then
                 print *,"Indexing error:"
                 print *,"index=",index
                 print *,"m =",m
                 print *,"xm=",xm(index)
                 print *,"n =",n
                 print *,"xn=",xn(index)
                 print *,"Are you sure the mpol and ntor for regcoil match those from nescoil?"
                 stop
              end if
              solution(index) = amplitude
           end if
        end do
     end do
     ! Convert from nescoil normalization to regcoil normalization:
     solution = solution * mu0/curpol_nescout
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

     call system_clock(tic,countrate)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Now compute diagnostics.
     ! This code should be verbatim copied from solve.f90
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     single_valued_current_potential_mn(:,ialpha) = solution
     this_current_potential = reshape(matmul(basis_functions, solution), (/ ntheta_coil, nzeta_coil /)) ! Could use BLAS2 for this line for extra speed.
     single_valued_current_potential_thetazeta(:,:,ialpha) = this_current_potential
     do izeta = 1,nzeta_coil
        do itheta = 1,ntheta_coil
           this_current_potential(itheta,izeta) = this_current_potential(itheta,izeta) &
                + factor_zeta*zeta_coil(izeta) + factor_theta*theta_coil(itheta)
        end do
     end do
     current_potential(:,:,ialpha) = this_current_potential

     JDifference_x = d_x - matmul(f_x, solution)
     JDifference_y = d_y - matmul(f_y, solution)
     JDifference_z = d_z - matmul(f_z, solution)
     this_J2_times_N = reshape(JDifference_x*Jdifference_x + JDifference_y*JDifference_y + JDifference_z*JDifference_z, (/ ntheta_coil, nzeta_coil /)) &
          / norm_normal_coil
     chi2_J(ialpha) = nfp * dtheta_coil * dzeta_coil * sum(this_J2_times_N)
     J2(:,:,ialpha) = this_J2_times_N / norm_normal_coil

     Bnormal_total(:,:,ialpha) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
          + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

     chi2_B(ialpha) = nfp * dtheta_plasma * dzeta_plasma &
          * sum(Bnormal_total(:,:,ialpha) * Bnormal_total(:,:,ialpha) * norm_normal_plasma)

     call system_clock(toc)
     print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
     print *,"  chi2_B:",chi2_B(ialpha),"chi2_J:",chi2_J(ialpha)

  end do outerLoop


end subroutine compute_diagnostics_for_nescout_potential
