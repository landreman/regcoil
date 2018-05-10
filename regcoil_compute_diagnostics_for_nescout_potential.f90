subroutine regcoil_compute_diagnostics_for_nescout_potential

  use regcoil_variables
  use safe_open_mod
  use stel_constants
  use stel_kinds

  implicit none

  integer :: iflag, tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: matrix, this_current_potential
  real(dp), dimension(:), allocatable :: RHS, solution
  real(dp), dimension(:), allocatable :: KDifference_x, KDifference_y, KDifference_z
  real(dp), dimension(:,:), allocatable :: this_K2_times_N
  real(dp) :: factor_theta, factor_zeta
  integer :: ilambda, itheta, izeta

  integer :: iunit = 7, temp1, status, m, n, mm, nn, index, istat
  character(300) :: myline
  real(dp) :: curpol_nescout, temp2, temp3, amplitude
  character(*), parameter :: matchString_curpol = "np, iota_edge, phip_edge, curpol"
  character(*), parameter :: matchString_phi = "---- Phi(m,n) for"
  character(*), parameter :: matchString_phi2 = "---- end Phi(m,n)"

  if (allocated(matrix)) deallocate(matrix)
  allocate(matrix(num_basis_functions, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(RHS)) deallocate(RHS)
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(solution)) deallocate(solution)
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(chi2_B)) deallocate(chi2_B)
  allocate(chi2_B(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(chi2_K)) deallocate(chi2_K)
  allocate(chi2_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(max_Bnormal)) deallocate(max_Bnormal)
  allocate(max_Bnormal(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(max_K)) deallocate(max_K)
  allocate(max_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(current_potential)) deallocate(current_potential)
  allocate(current_potential(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(single_valued_current_potential_thetazeta)) &
        deallocate(single_valued_current_potential_thetazeta)
  allocate(single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(this_current_potential)) deallocate(this_current_potential)
  allocate(this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(single_valued_current_potential_mn)) deallocate(single_valued_current_potential_mn)
  allocate(single_valued_current_potential_mn(num_basis_functions,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(Bnormal_total)) deallocate(Bnormal_total)
  allocate(Bnormal_total(ntheta_plasma,nzeta_plasma,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(K2)) deallocate(K2)
  allocate(K2(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(KDifference_x)) deallocate(KDifference_x)
  allocate(KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(KDifference_y)) deallocate(KDifference_y)
  allocate(KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(KDifference_z)) deallocate(KDifference_z)
  allocate(KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (allocated(this_K2_times_N)) deallocate(this_K2_times_N)
  allocate(this_K2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

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
     do mm = 0,mpol_coil
        do nn = -ntor_coil, ntor_coil
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
              index = ntor_coil + (mm-1)*(ntor_coil*2+1) + ntor_coil - nn + 1
              if ((xm_coil(index) .ne. m) .or. (xn_coil(index) .ne. -nfp*n)) then
                 print *,"Indexing error:"
                 print *,"index=",index
                 print *,"m =",m
                 print *,"xm_coil=",xm_coil(index)
                 print *,"n =",n
                 print *,"xn_coil=",xn_coil(index)
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

     call system_clock(tic,countrate)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Now compute diagnostics.
     ! This code should be verbatim copied from regcoil_solve.f90
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     single_valued_current_potential_mn(:,ilambda) = solution
     this_current_potential = reshape(matmul(basis_functions, solution), (/ ntheta_coil, nzeta_coil /)) ! Could use BLAS2 for this line for extra speed.
     single_valued_current_potential_thetazeta(:,:,ilambda) = this_current_potential
     do izeta = 1,nzeta_coil
        do itheta = 1,ntheta_coil
           this_current_potential(itheta,izeta) = this_current_potential(itheta,izeta) &
                + factor_zeta*zeta_coil(izeta) + factor_theta*theta_coil(itheta)
        end do
     end do
     current_potential(:,:,ilambda) = this_current_potential

     KDifference_x = d_x - matmul(f_x, solution)
     KDifference_y = d_y - matmul(f_y, solution)
     KDifference_z = d_z - matmul(f_z, solution)
     this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
          / norm_normal_coil
     chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)
     K2(:,:,ilambda) = this_K2_times_N / norm_normal_coil

     Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
          + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

     max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
     max_K(ilambda) = sqrt(maxval(K2(:,:,ilambda)))

     chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
          * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)

     call system_clock(toc)
     print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
     print *,"  chi2_B:",chi2_B(ilambda),"chi2_K:",chi2_K(ilambda)

  end do outerLoop


end subroutine regcoil_compute_diagnostics_for_nescout_potential
