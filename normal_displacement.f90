subroutine normal_displacement()

  use global_variables

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK, M, N, LDA, LDU, LDVT
  real(dp), dimension(:), allocatable :: WORK, singular_values
  real(dp), dimension(:,:), allocatable :: U, VT
  integer, dimension(:), allocatable :: IPIV
  integer, dimension(:), allocatable :: IWORK
  character :: JOBZ

  real(dp), dimension(:,:), allocatable :: dxdomega, dydomega, dzdomega, sigma_UT
  real(dp), dimension(:,:), allocatable :: dchi2dx, dchi2dy,dchi2dz
  integer :: izeta_coil, itheta_coil, index_coil, izetal_coil, indexl_coil, l_coil

  allocate(dxdomega(ntheta_coil*nzeta_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dydomega(ntheta_coil*nzeta_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dzdomega(ntheta_coil*nzeta_coil,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadx(nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegady(nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadz(nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do izeta_coil = 1,nzeta_coil
    do itheta_coil = 1,ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      dxdomega(index_coil,:) = drdomega(1,index_coil,1,:)
      dydomega(index_coil,:) = drdomega(2,index_coil,1,:)
      dzdomega(index_coil,:) = drdomega(3,index_coil,1,:)
    enddo
  enddo

  ! A has dimension N x N (to obtain pseudoinverse we must have M >= N)
  ! n_singular values = min(M,N) = N
  ! Sigma must be square (has dimension N x N)
  ! VT must have dimension N x N
  ! U must have dimension M x N
  ! Define variables needed for LAPACK SVD decomposition
  JOBZ='A'
  M = ntheta_coil*nzeta_coil ! number of rows of A
  N = nomega_coil ! number of columns of A
  LDA = M ! leading dimension of A
  LDU = M ! leading dimension of U
  LDVT = N ! leading dimension of VT
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
  3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
  min(M,N)*(6+4*min(M,N))+max(M,N))

  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  n_singular_values = min(N,M)
  allocate(singular_values(n_singular_values),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(U(M,M),stat=iflag)
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
  sigma_UT = 0
  do i = 1,n_singular_values
    if (singular_values(i) /= 0) then
      sigma_UT(i,:) = (1/singular_values(i))*U(:,i)
    endif
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
  !PA = V*sigma^{-1}*UT
  !UT(N,M)
  ! Below we are multiplying the rows of UT by diagonal elements of sigma^{-1}
  sigma_UT = 0
  do i = 1,n_singular_values
    if (singular_values(i) /= 0) then
      sigma_UT(i,:) = (1/singular_values(i))*U(:,i)
    endif
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
  sigma_UT = 0
  do i = 1,n_singular_values
    if (singular_values(i) /= 0) then
      sigma_UT(i,:) = (1/singular_values(i))*U(:,i)
    endif
  enddo
  domegadz = matmul(transpose(VT),sigma_UT)

  allocate(dchi2dx(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dy(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dz(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dr_normal(nlambda,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  dchi2dx = matmul(transpose(domegadx),dchi2domega)
  dchi2dy = matmul(transpose(domegady),dchi2domega)
  dchi2dz = matmul(transpose(domegadz),dchi2domega)

  do itheta_coil = 1, ntheta_coil
		do izeta_coil = 1, nzeta_coil
			index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
			dchi2dr_normal(:,index_coil) = &
				(normal_coil(1,itheta_coil,izeta_coil)*dchi2dx(index_coil,:) &
				+ normal_coil(2,itheta_coil,izeta_coil)*dchi2dy(index_coil,:) &
				+ normal_coil(3,itheta_coil,izeta_coil)*dchi2dz(index_coil,:)) &
				/norm_normal_coil(itheta_coil,izeta_coil)
      enddo
  enddo

end subroutine normal_displacement
