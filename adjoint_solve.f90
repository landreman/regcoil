subroutine adjoint_solve

  use global_variables
  use stel_constants
  use omp_lib

  implicit none

  integer :: ilambda
  integer :: iflag, tic, toc, countrate
  integer :: minLambda, MaxLambda

  real(dp), dimension(:,:), allocatable :: dKDifferencedomega, term1, term2, dBnormaldomega
  integer :: iomega
  real(dp), dimension(:,:), allocatable :: dchi2Kdphi, dchi2Bdphi
  real(dp), dimension(:,:), allocatable :: adjoint_bx, adjoint_by, adjoint_bz
  real(dp), dimension(:,:), allocatable :: adjoint_Ax, adjoint_Ay, adjoint_Az
  real(dp), dimension(:,:), allocatable :: adjoint_c
  real(dp), dimension(:,:), allocatable :: dFKdomega, dFBdomega, dFdomega
  real(dp), dimension(:,:), allocatable :: matrix, this_K2_times_N
  real(dp), dimension(:), allocatable :: RHS, KDifference_x, KDifference_y, KDifference_z, solution

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  if (general_option == 1) then
    minLambda = 1
    maxLambda = NLambda
    nlambda_sensitivity = Nlambda
  else if (general_option > 3) then
    minLambda = NLambda
    maxLambda = NLambda
    nlambda_sensitivity = 1
  end if

  allocate(dchi2domega(nomega_coil,nlambda_sensitivity),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Kdomega(nomega_coil,nlambda_sensitivity),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2Bdomega(nomega_coil,nlambda_sensitivity),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dKDifferencedomega(3,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(term1(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(term2(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dBnormaldomega(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (sensitivity_option == 3 .or. sensitivity_option == 4) then
    allocate(dchi2Kdphi(num_basis_functions,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(q_K(num_basis_functions,nlambda_sensitivity),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_bx(ntheta_coil*nzeta_coil,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_by(ntheta_coil*nzeta_coil,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_bz(ntheta_coil*nzeta_coil,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Ax(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Ay(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_Az(num_basis_functions,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 5) then
    allocate(dchi2Bdphi(num_basis_functions,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(adjoint_c(ntheta_coil*nzeta_coil,1),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(q_B(num_basis_functions,nlambda_sensitivity),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif
  if (sensitivity_option > 2) then
    allocate(dmatrixdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dRHSdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dFdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif

  allocate(matrix(num_basis_functions,num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_K2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  print *,"Optimal LWORK:",LWORK
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do ilambda=minLambda,maxLambda

    matrix = matrix_B + lambda(ilambda) * matrix_K
    RHS    =    RHS_B + lambda(ilambda) *    RHS_K
    solution = single_valued_current_potential_mn(:,ilambda)

    ! dmatrixdomega and dRHSdomega needed for adjoint solve
    if (sensitivity_option > 2) then
      dmatrixdomega = dmatrix_Bdomega + lambda(ilambda)*dmatrix_Kdomega
      dRHSdomega = dRHS_Bdomega + lambda(ilambda)*dRHS_Kdomega
    endif

    KDifference_x = d_x - matmul(f_x, solution)
    KDifference_y = d_y - matmul(f_y, solution)
    KDifference_z = d_z - matmul(f_z, solution)
    this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
        / norm_normal_coil

    call system_clock(tic,countrate)
    !$OMP PARALLEL
    !$OMP MASTER
    print *,"  Number of OpenMP threads:",omp_get_num_threads()
    !$OMP END MASTER
    !$OMP DO PRIVATE(dKDifferencedomega,term1,term2)
    do iomega = 1, nomega_coil
      dKDifferencedomega(1,:) = dddomega(1,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfxdomega(iomega,:,:), solution)
      dKDifferencedomega(2,:) = dddomega(2,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfydomega(iomega,:,:), solution)
      dKDifferencedomega(3,:) = dddomega(3,iomega,1:ntheta_coil*nzeta_coil)-matmul(dfzdomega(iomega,:,:), solution)
      term1 = -dnorm_normaldomega(iomega,:,:)*this_K2_times_N/norm_normal_coil
      term2 = reshape(KDifference_x*dKDifferencedomega(1,:) + KDifference_y*dKDifferencedomega(2,:) &
        + KDifference_z*dKDifferencedomega(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
      dchi2Kdomega(iomega,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(term1) + sum(term2))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call system_clock(toc)
    print *,"chi2_K sensitivity in solve:",real(toc-tic)/countrate," sec."

    if (sensitivity_option == 3 .or. sensitivity_option == 4) then
      call system_clock(tic,countrate)
      ! Adjoint chi2_K calculation - compute dchi2Kdphi
      ! dchi2Kdphi(num_basis_functions,1)
      adjoint_bx(:,1) = Kdifference_x/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_by(:,1) = Kdifference_y/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_bz(:,1) = Kdifference_z/reshape(norm_normal_coil, (/ ntheta_coil*nzeta_coil /))
      adjoint_Ax = transpose(f_x)
      adjoint_Ay = transpose(f_y)
      adjoint_Az = transpose(f_z)
      dchi2Kdphi = -2*nfp*dtheta_coil*dzeta_coil*(matmul(adjoint_Ax,adjoint_bx) + matmul(adjoint_Ay,adjoint_by) + matmul(adjoint_Az,adjoint_bz))
      ! Solve Adjoint Equation
      ! matrix was overwritten previously
      matrix = matrix_B + lambda(ilambda) * matrix_K
      call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, dchi2Kdphi, num_basis_functions, WORK, LWORK, INFO)
      if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
      end if
      q_K(:,ilambda) = dchi2Kdphi(:,1)
      call system_clock(toc)
      print *,"Adjoint chi2_K calculation in solve:",real(toc-tic)/countrate," sec."
    endif

    call system_clock(tic,countrate)
    !$OMP PARALLEL
    !$OMP MASTER
    print *,"  Number of OpenMP threads:",omp_get_num_threads()
    !$OMP END MASTER
    !$OMP DO PRIVATE(dBnormaldomega)
    do iomega = 1, nomega_coil
     dBnormaldomega = reshape(matmul(dgdomega(iomega,:,:),solution),(/ ntheta_plasma, nzeta_plasma /))/norm_normal_plasma + reshape(dhdomega(iomega,:),(/ ntheta_plasma, nzeta_plasma /))/norm_normal_plasma
     dchi2Bdomega(iomega,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldomega*norm_normal_plasma)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call system_clock(toc)
    print *,"chi2_B sensitivity in solve:",real(toc-tic)/countrate," sec."
    dchi2domega(:,ilambda) = dchi2Bdomega(:,ilambda) + lambda(ilambda)*dchi2Kdomega(:,ilambda)

    ! Compute dFdomega for adjoint solution
    if (sensitivity_option > 2) then
      call system_clock(tic,countrate)
      do iomega = 1, nomega_coil
        dFdomega(iomega,:) = matmul(dmatrixdomega(iomega,:,:),solution) - dRHSdomega(iomega,:)
      enddo
      call system_clock(toc)
      print *,"dFdomega in solve:", real(toc-tic)/countrate," sec."
    endif

    if (sensitivity_option == 3 .or. sensitivity_option == 5) then
      ! Compute dchi2Bdphi and q_B
      call system_clock(tic,countrate)
      adjoint_c(:,1) = reshape(Bnormal_total(:,:,ilambda), (/ ntheta_plasma * nzeta_plasma /))
      dchi2Bdphi = 2*nfp*dtheta_plasma*dzeta_plasma*matmul(transpose(g), adjoint_c)
      ! Solve Adjoint Equation
      ! matrix was overwritten previously
      matrix = matrix_B + lambda(ilambda) * matrix_K
      call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, dchi2Bdphi(:,1), num_basis_functions, WORK, LWORK, INFO)
      if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
      end if
      q_B(:,ilambda) = dchi2Bdphi(:,1)
      call system_clock(toc)
      print *,"Adjoint chi2_B calculation in solve:",real(toc-tic)/countrate," sec."
    endif
    if (sensitivity_option == 3) then
      call system_clock(tic,countrate)
      dchi2Kdomega(:,ilambda) = dchi2Kdomega(:,ilambda) - matmul(dFdomega, q_K(:,ilambda))
      dchi2Bdomega(:,ilambda) = dchi2Bdomega(:,ilambda) - matmul(dFdomega, q_B(:,ilambda))
      call system_clock(toc)
      print *,"dchi2K and dchi2B in solve:",real(toc-tic)/countrate," sec."
    endif
    if (sensitivity_option == 4) then
      call system_clock(tic,countrate)
      dchi2Kdomega(:,ilambda) = dchi2Kdomega(:,ilambda) - matmul(dFdomega, q_K(:,ilambda))
      dchi2Bdomega(:,ilambda) = dchi2domega(:,ilambda) - lambda(ilambda)*dchi2Kdomega(:,ilambda)
      call system_clock(toc)
      print *,"dchi2K and dchi2B in solve:",real(toc-tic)/countrate," sec."
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
      print *,"dchi2K and dchi2B in solve:",real(toc-tic)/countrate," sec."
    endif

  end do

  deallocate(term1)
  deallocate(term2)
  deallocate(dKDifferencedomega)
  deallocate(dBnormaldomega)
  if (sensitivity_option == 3 .or. sensitivity_option == 4) then
    deallocate(dchi2Kdphi)
    deallocate(adjoint_Ax)
    deallocate(adjoint_bx)
    deallocate(adjoint_Ay)
    deallocate(adjoint_by)
    deallocate(adjoint_Az)
    deallocate(adjoint_bz)
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 5) then
    deallocate(dchi2Bdphi)
    deallocate(adjoint_c)
  endif
  if (sensitivity_option > 2) then
    deallocate(dFdomega)
  endif


end subroutine adjoint_solve
