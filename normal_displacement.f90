subroutine normal_displacement()

  use global_variables

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK, M, N, LDA, LDU, LDVT
  real(dp), dimension(:), allocatable :: WORK, singular_values
  real(dp), dimension(:,:), allocatable :: U, VT
  integer, dimension(:), allocatable :: IPIV
  integer, dimension(:), allocatable :: IWORK
  character :: JOBZ

  real(dp), dimension(:,:), allocatable :: dxdomega, dydomega, dzdomega
  real(dp), dimension(:,:), allocatable :: domegadx, domegady, domegadz, sigma_UT
  real(dp), dimension(:,:), allocatable :: dchi2dx, dchi2dy,dchi2dz
  integer :: izeta_coil, itheta_coil, index_coil, izetal_coil, indexl_coil

  allocate(dxdomega(ntheta_coil*nzetal_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dydomega(ntheta_coil*nzetal_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dzdomega(ntheta_coil*nzetal_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadx(nomega_coil,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegady(nomega_coil,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadz(nomega_coil,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do izeta_coil = 1,nzeta_coil
    do itheta_coil = 1,ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      do l_coil = 0, (nfp-1)
        izetal_coil = izeta_coil + l_coil*nzeta_coil
        indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil
        dxdomega(indexl_coil,:) = drdomega(1,index_coil,l_coil+1,:)
        dydomega(indexl_coil,:) = drdomega(2,index_coil,l_coil+1,:)
        dzdomega(indexl_coil,:) = drdomega(3,index_coil,l_coil+1,:)
      enddo
    enddo
  enddo

  ! A has dimension N x N (to obtain pseudoinverse we must have M >= N)
  ! n_singular values = min(M,N) = N
  ! Sigma must be square (has dimension N x N)
  ! VT must have dimension N x N
  ! U must have dimension M x N
  ! Define variables needed for LAPACK SVD decomposition
  JOBZ='S'
  M = ntheta_coil*nzetal_coil ! number of rows of A
  N = nomega_coil ! number of columns of A
  LDA = M ! leading dimension of A
  LDU = M ! leading dimension of U
  LDVT = N ! leading dimension of VT
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*N + max(M,7*N), 3*N + max(M,5*N*N + 4*N),N*(6 + 4*N)+M)

  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  n_singular_values = N
  allocate(singular_values(n_singular_values),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(U(M,N),stat=iflag) ! If all singular vectors were computed, U would be M*M. But here we only compute the first N singular vectors,
  ! so U is M*N.
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  ! Call LAPACK to do the SVD:
  ! Note that svd_matrix is destroyed by LAPACK!
  call DGESDD(JOBZ, M, N, dxdomega, LDA, singular_values, U, LDU, &
    VT, LDVT, WORK, LWORK, IWORK, INFO)

  if (INFO==0) then
    print *,"SVD (DGESDD) successful."
  else if (INFO>0) then
    print *,"Error in SVD (DGESDD): Did not converge."
  else
    print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
  end if

 !Compute pseudoinverse
 !PA = V*Sigma^{-1}*UT
 !UT(N,M)
  allocate(sigma_UT(N,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  do i = 1,M
    sigma_UT(:,i) = singular_values(i)*U(i,:)
  enddo
  domegadx = matmul(transpose(VT),sigma_UT)

  ! Call LAPACK to do the SVD:
  ! Note that svd_matrix is destroyed by LAPACK!
  call DGESDD(JOBZ, M, N, dydomega, LDA, singular_values, U, LDU, &
  VT, LDVT, WORK, LWORK, IWORK, INFO)

  if (INFO==0) then
    print *,"SVD (DGESDD) successful."
  else if (INFO>0) then
    print *,"Error in SVD (DGESDD): Did not converge."
  else
    print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
  end if

  !Compute pseudoinverse
  !PA = V*Sigma^{-1}*UT
  !UT(N,M)
  do i = 1,M
    sigma_UT(:,i) = singular_values(i)*U(i,:)
  enddo
  domegady = matmul(transpose(VT),sigma_UT)

  ! Call LAPACK to do the SVD:
  ! Note that svd_matrix is destroyed by LAPACK!
  call DGESDD(JOBZ, M, N, dzdomega, LDA, singular_values, U, LDU, &
  VT, LDVT, WORK, LWORK, IWORK, INFO)

  if (INFO==0) then
    print *,"SVD (DGESDD) successful."
  else if (INFO>0) then
    print *,"Error in SVD (DGESDD): Did not converge."
  else
    print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
  end if

  !Compute pseudoinverse
  !PA = V*Sigma^{-1}*UT
  !UT(N,M)
  do i = 1,M
    sigma_UT(:,i) = singular_values(i)*U(i,:)
  enddo
  domegadz = matmul(transpose(VT),sigma_UT)

  allocate(dchi2dx(nlambda,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dy(nlambda,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dz(nlambda,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dr_normal(nlambda,ntheta_coil*nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  dchi2dx = matmul(transpose(domegadx),dchi2domega)
  dchi2dy = matmul(transpose(domegady),dchi2domega)
  dchi2dz = matmul(transpose(domegadz),dchi2domega)

  do itheta = 1, ntheta_coil
      do izeta = 1, nzeta_coil
        indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil
        dchi2dr_normal(:,indexl_coil) = &
            (normal_coil(1,itheta_coil,izetal_coil)*dchi2dx(:,indexl_coil) &
          + normal_coil(2,itheta_coil,izetal_coil)*dchi2dy(:,indexl_coil) &
          + normal_coil(3,itheta_coil,izetal_coil)*dchi2dz(:,indexl_coil)) &
          /norm_normal_coil(itheta_coil,izetal_coil)
      enddo
  enddo

end subroutine normal_displacement
