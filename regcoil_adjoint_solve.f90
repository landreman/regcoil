subroutine regcoil_adjoint_solve

  use regcoil_variables
  use stel_constants
  use omp_lib

  implicit none

  integer :: ilambda
  integer :: iflag, tic, toc, countrate

  real(dp), dimension(:,:), allocatable :: term1, term2, dBnormaldomega, &
    adjoint_Ax, adjoint_Ay, adjoint_Az, dKDifferencedomega, dFKdomega, dFBdomega, &
    dchi2Kdphi, dRHSdomega
  integer :: iomega
  real(dp), dimension(:), allocatable :: adjoint_bx, adjoint_by, adjoint_bz, &
    adjoint_c
  integer :: minLambda, maxLambda
  real(dp), dimension(:,:,:), allocatable :: dmatrixdomega
  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  ! In case of a lambda search, sensitivity only computed for Nlambda
  if (general_option == 1) then
    minLambda = 1
    maxLambda = NLambda
  else if (general_option > 3) then
    minLambda = NLambda
    maxLambda = NLambda
  end if

  allocate(dRMSKdomega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2domega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Kdomega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Kdomega_withoutadjoint(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Bdomega(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Bdomega_withoutadjoint(nomega_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dKDifferencedomega(3,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(term1(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(term2(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dBnormaldomega(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (sensitivity_option == 3 .or. sensitivity_option == 4 .or. fixed_norm_sensitivity_option) then
    allocate(dchi2Kdphi(nlambda, num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_bx(ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_by(ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_bz(ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Ax(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Ay(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Az(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  end if
  if (sensitivity_option == 3 .or. sensitivity_option == 4) then
    allocate(q_K(num_basis_functions,nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 5 .or. fixed_norm_sensitivity_option) then
    allocate(dchi2Bdphi(nlambda,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_c(ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  end if
  if (sensitivity_option == 3 .or. sensitivity_option == 5) then
    allocate(q_B(num_basis_functions,nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif
  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    allocate(dmatrixdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dRHSdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  end if
  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    allocate(dFdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif

  allocate(WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  if (verbose) then
    print *,"Optimal LWORK:",LWORK
  end if
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do ilambda=minLambda,maxLambda

    matrix = matrix_B + lambda(ilambda) * matrix_regularization
    RHS    =    RHS_B + lambda(ilambda) *    RHS_regularization
    solution = single_valued_current_potential_mn(:,ilambda)

    KDifference_x = d_x - matmul(f_x, solution)
    KDifference_y = d_y - matmul(f_y, solution)
    KDifference_z = d_z - matmul(f_z, solution)

    ! dmatrixdomega and dRHSdomega needed for adjoint solve
    if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
      dmatrixdomega = dmatrix_Bdomega + lambda(ilambda)*dmatrix_Kdomega
      dRHSdomega = dRHS_Bdomega + lambda(ilambda)*dRHS_Kdomega
    endif

    this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
        / norm_normal_coil

    call system_clock(tic,countrate)
    !$OMP PARALLEL
    !$OMP MASTER
    if (verbose) then
      print *,"  Number of OpenMP threads:",omp_get_num_threads()
    end if
    !$OMP END MASTER
    !$OMP DO PRIVATE(dKDifferencedomega,term1,term2)
    do iomega = 1, nomega_coil
      dKDifferencedomega(1,:) = dddomega(1,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfxdomega(iomega,:,:), solution)
      dKDifferencedomega(2,:) = dddomega(2,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfydomega(iomega,:,:), solution)
      dKDifferencedomega(3,:) = dddomega(3,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfzdomega(iomega,:,:), solution)
      term1 = -dnorm_normaldomega(iomega,:,:)*this_K2_times_N/norm_normal_coil
      term2 = reshape(KDifference_x*dKDifferencedomega(1,:) + KDifference_y*dKDifferencedomega(2,:) &
        + KDifference_z*dKDifferencedomega(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
      dchi2Kdomega_withoutadjoint(iomega,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(term1) + sum(term2))
      dchi2Kdomega(iomega,ilambda) = dchi2Kdomega_withoutadjoint(iomega,ilambda)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call system_clock(toc)
    if (verbose) then
      print *,"chi2_K sensitivity in regcoil_adjoint_solve :",real(toc-tic)/countrate," sec."
    end if

    if (sensitivity_option == 3 .or. sensitivity_option == 4 .or. fixed_norm_sensitivity_option) then
      ! Adjoint chi2_K calculation - compute dchi2Kdphi
      adjoint_bx = Kdifference_x/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_by = Kdifference_y/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_bz = Kdifference_z/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_Ax = transpose(f_x)
      adjoint_Ay = transpose(f_y)
      adjoint_Az = transpose(f_z)
      dchi2Kdphi(ilambda,:) = -2*nfp*dtheta_coil*dzeta_coil*(matmul(adjoint_Ax,adjoint_bx) + matmul(adjoint_Ay,adjoint_by) + matmul(adjoint_Az,adjoint_bz))
    end if

    if (sensitivity_option == 3 .or. sensitivity_option == 4) then
      call system_clock(tic,countrate)
      ! Solve Adjoint Equation
      ! matrix was overwritten previously
      matrix = matrix_B + lambda(ilambda) * matrix_regularization
      call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, dchi2Kdphi(ilambda,:), num_basis_functions, WORK, LWORK, INFO)
      if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
      end if
      q_K(:,ilambda) = dchi2Kdphi(ilambda,:)
      call system_clock(toc)
      if (verbose) then
        print *,"Adjoint chi2_K calculation in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
      end if
    endif

    call system_clock(tic,countrate)
    !$OMP PARALLEL
    !$OMP MASTER
    if (verbose) then
      print *,"  Number of OpenMP threads:",omp_get_num_threads()
    end if
    !$OMP END MASTER
    !$OMP DO PRIVATE(dBnormaldomega)
    do iomega = 1, nomega_coil
     dBnormaldomega = reshape(matmul(dgdomega(:,:,iomega),solution),(/ ntheta_plasma, nzeta_plasma /))/norm_normal_plasma + reshape(dhdomega(iomega,:),(/ ntheta_plasma, nzeta_plasma /))/norm_normal_plasma
     dchi2Bdomega_withoutadjoint(iomega,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldomega*norm_normal_plasma)
     dchi2Bdomega(iomega,ilambda) = dchi2Bdomega_withoutadjoint(iomega,ilambda)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call system_clock(toc)
    if (verbose) then
      print *,"chi2_B sensitivity in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
    end if
    dchi2domega(:,ilambda) = dchi2Bdomega(:,ilambda) + lambda(ilambda)*dchi2Kdomega(:,ilambda)

    ! Compute dFdomega for adjoint solution
    if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
      call system_clock(tic,countrate)
      do iomega = 1, nomega_coil
        dFdomega(iomega,:) = matmul(dmatrixdomega(iomega,:,:),solution) - dRHSdomega(iomega,:)
      enddo
      call system_clock(toc)
      if (verbose) then
        print *,"dFdomega in regcoil_adjoint_solve:", real(toc-tic)/countrate," sec."
      end if
    endif
    if (sensitivity_option == 3 .or. sensitivity_option == 5 .or. fixed_norm_sensitivity_option) then
      adjoint_c = reshape(Bnormal_total(:,:,ilambda), (/ ntheta_plasma * nzeta_plasma /))
      dchi2Bdphi(ilambda,:) = 2*nfp*dtheta_plasma*dzeta_plasma*matmul(transpose(g), adjoint_c)
    end if
    if (sensitivity_option == 3 .or. sensitivity_option == 5) then
      ! Compute dchi2Bdphi and q_B
      call system_clock(tic,countrate)
      ! Solve Adjoint Equation
      ! matrix was overwritten previously
      matrix = matrix_B + lambda(ilambda) * matrix_regularization
      call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, dchi2Bdphi(ilambda,:), num_basis_functions, WORK, LWORK, INFO)
      if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
      end if
      q_B(:,ilambda) = dchi2Bdphi(ilambda,:)
      call system_clock(toc)
      if (verbose) then
        print *,"Adjoint chi2_B calculation in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
      end if
    endif
    if (sensitivity_option == 3) then
      call system_clock(tic,countrate)
      dchi2Kdomega(:,ilambda) = dchi2Kdomega(:,ilambda) - matmul(dFdomega, q_K(:,ilambda))
      dchi2Bdomega(:,ilambda) = dchi2Bdomega(:,ilambda) - matmul(dFdomega, q_B(:,ilambda))
      call system_clock(toc)
      if (verbose) then
        print *,"dchi2K and dchi2B in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
      end if
    endif
    if (sensitivity_option == 4) then
      call system_clock(tic,countrate)
      dchi2Kdomega(:,ilambda) = dchi2Kdomega(:,ilambda) - matmul(dFdomega, q_K(:,ilambda))
      dchi2Bdomega(:,ilambda) = dchi2domega(:,ilambda) - lambda(ilambda)*dchi2Kdomega(:,ilambda)
      call system_clock(toc)
      if (verbose) then
        print *,"dchi2K and dchi2B in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
      end if
    endif
    if (sensitivity_option == 5) then
      call system_clock(tic,countrate)
      dchi2Bdomega(:,ilambda) = dchi2Bdomega(:,ilambda) - matmul(dFdomega, q_B(:,ilambda))
      if (lambda(ilambda) /= 0) then
        dchi2Kdomega(:,ilambda) = (dchi2domega(:,ilambda) - dchi2Bdomega(:,ilambda))/lambda(ilambda)
      else
        dchi2Kdomega(:,ilambda) = 0
      endif
      call system_clock(toc)
      if (verbose) then
        print *,"dchi2K and dchi2B in regcoil_adjoint_solve:",real(toc-tic)/countrate," sec."
      end if
    endif
    if (sensitivity_option > 2) then
      dRMSKdomega(:,ilambda)= 0.5*(chi2_K(ilambda)/area_coil)**(-0.5)*(dchi2Kdomega(:,ilambda)/area_coil - chi2_K(ilambda)*darea_coildomega/area_coil**2)
    end if
  end do
  
  deallocate(term1)
  deallocate(term2)
  deallocate(dBnormaldomega)
  if (sensitivity_option == 3 .or. sensitivity_option == 4) then
    deallocate(adjoint_Ax)
    deallocate(adjoint_bx)
    deallocate(adjoint_Ay)
    deallocate(adjoint_by)
    deallocate(adjoint_Az)
    deallocate(adjoint_bz)
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 5) then
    deallocate(adjoint_c)
  endif


end subroutine regcoil_adjoint_solve
