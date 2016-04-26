! Documentation of LAPACK's SVD subroutine DGESDD is copied at the end of this file for convenience.
subroutine svd_inductance_matrices()
  ! This subroutine finds the singular values of the two inductance matrices

  use global_variables, only: inductance_plasma_outer, inductance_middle_outer, inductance_plasma_middle, &
       num_basis_functions_plasma, num_basis_functions_middle, num_basis_functions_outer, &
       n_singular_values_inductance_plasma_outer,  svd_s_inductance_plasma_outer, &
       n_singular_values_inductance_middle_outer,  svd_s_inductance_middle_outer, &
       n_singular_values_inductance_plasma_middle, svd_s_inductance_plasma_middle, &
       svd_uT_inductance_middle_outer, svd_v_inductance_middle_outer, allSVDsSucceeded, &
       n_singular_vectors_to_save, nu_plasma, nv_plasma, nu_middle, nv_middle, &
       basis_functions_plasma, basis_functions_middle, &
       svd_u_inductance_plasma_middle, svd_v_inductance_plasma_middle, &
       svd_u_inductance_plasma_middle_uv, svd_v_inductance_plasma_middle_uv, &
       xm_plasma, xn_plasma, mnmax_plasma, save_vectors_in_uv_format, &
       svd_u_inductance_plasma_middle_dominant_m, svd_u_inductance_plasma_middle_dominant_n, &
       svd_u_inductance_plasma_middle_all, svd_v_inductance_plasma_middle_all, &
       Bnormal_from_1_over_R_field_inductance, Bnormal_from_1_over_R_field, &
       Bnormal_from_const_v_coils_inductance,  Bnormal_from_const_v_coils, &
       Bnormal_from_plasma_current_inductance, Bnormal_from_plasma_current
  
  use stel_kinds
  
  implicit none
  
  character :: JOBZ
  integer :: INFO, LDA, LDU, LDVT, LWORK, M, N, iflag
  integer :: tic, toc, countrate, tic1, toc1
  real(dp), dimension(:,:), allocatable :: A, U, VT
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IWORK
  integer :: i,j,index
  real(dp) :: val

  !*************************************************************************
  ! Beginning of section related to the plasma-to-outer inductance matrix.
  !*************************************************************************
  
  allSVDsSucceeded = .true.

  print *,"Beginning SVD of the inductance matrix between the plasma and outer surfaces."
  call system_clock(tic,countrate)
  
  JOBZ='N'  ! For now compute none of the singular vectors. We could change this.
  M = num_basis_functions_plasma
  N = num_basis_functions_outer
  LDA = M
  LDU = M
  LDVT = N
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
       3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
       min(M,N)*(6+4*min(M,N))+max(M,N))
  
  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*min(M,N)),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  n_singular_values_inductance_plasma_outer = min(M,N)
  allocate(svd_s_inductance_plasma_outer(n_singular_values_inductance_plasma_outer),stat=iflag)
  
  ! Matrix is destroyed by LAPACK, so make a copy:
  allocate(A(M,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  A = inductance_plasma_outer
  
  allocate(U(M,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  ! Call LAPACK to do the SVD:
  call DGESDD(JOBZ, M, N, A, LDA, svd_s_inductance_plasma_outer, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)
  
  if (INFO==0) then
     print *,"SVD (DGESDD) successful."
     if (n_singular_values_inductance_plasma_outer<5) then
        print *,"Singular values:",svd_s_inductance_plasma_outer
     else
        print *,"First 5 singular values:",svd_s_inductance_plasma_outer(1:5)
        print *,"Last 5 singular values:", &
             svd_s_inductance_plasma_outer(n_singular_values_inductance_plasma_outer-4:n_singular_values_inductance_plasma_outer)
     end if
  else if (INFO>0) then
     print *,"Error in SVD (DGESDD): Did not converge."
     allSVDsSucceeded = .false.
  else
     print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
     allSVDsSucceeded = .false.
  end if
  
  deallocate(A,U,VT,WORK,IWORK)
  
  call system_clock(toc)
  print *,"Done with SVD. Took ",real(toc-tic)/countrate," sec."
  
  !*************************************************************************
  ! End of section related to the plasma-to-outer inductance matrix.
  !*************************************************************************
  
  !*************************************************************************
  ! Beginning of section related to the middle-to-outer inductance matrix.
  !*************************************************************************
  
  print *,"Beginning SVD of the inductance matrix between the middle and outer surfaces."
  call system_clock(tic,countrate)
  
  JOBZ='A'  ! For the middle-outer inductance matrix, we need all the singular vectors.
  M = num_basis_functions_middle
  N = num_basis_functions_outer
  LDA = M
  LDU = M
  LDVT = N
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
       3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
       min(M,N)*(6+4*min(M,N))+max(M,N))
  
  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*min(M,N)),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  n_singular_values_inductance_middle_outer = min(M,N)
  allocate(svd_s_inductance_middle_outer(n_singular_values_inductance_middle_outer),stat=iflag)
  
  ! Matrix is destroyed by LAPACK, so make a copy:
  allocate(A(M,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  A = inductance_middle_outer
  
  allocate(U(M,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_uT_inductance_middle_outer(M,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_v_inductance_middle_outer(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  ! Call LAPACK to do the SVD:
  call DGESDD(JOBZ, M, N, A, LDA, svd_s_inductance_middle_outer, U, LDU, &
       VT, LDVT, WORK, LWORK, IWORK, INFO)
  
  if (INFO==0) then
     print *,"SVD (DGESDD) successful."
     if (n_singular_values_inductance_middle_outer<5) then
        print *,"Singular values:",svd_s_inductance_middle_outer
     else
        print *,"First 5 singular values:",svd_s_inductance_middle_outer(1:5)
        print *,"Last 5 singular values:", &
             svd_s_inductance_middle_outer(n_singular_values_inductance_middle_outer-4:n_singular_values_inductance_middle_outer)
     end if
  else if (INFO>0) then
     print *,"Error in SVD (DGESDD): Did not converge."
     allSVDsSucceeded = .false.
  else
     print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
     allSVDsSucceeded = .false.
  end if
  
  svd_uT_inductance_middle_outer = transpose(U)
  svd_v_inductance_middle_outer = transpose(VT)
  deallocate(A,U,VT,WORK,IWORK)
  
  call system_clock(toc)
  print *,"Done with SVD. Took ",real(toc-tic)/countrate," sec."
  
  !*************************************************************************
  ! End of section related to the middle-to-outer inductance matrix.
  !*************************************************************************
  
  !*************************************************************************
  ! Beginning of section related to the plasma-to-middle inductance matrix.
  !*************************************************************************
  
  print *,"Beginning SVD of the inductance matrix between the plasma and middle surfaces."
  call system_clock(tic,countrate)
  
  JOBZ='A'  ! For the middle-outer inductance matrix, we need all the singular vectors.
  M = num_basis_functions_plasma
  N = num_basis_functions_middle
  LDA = M
  LDU = M
  LDVT = N
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
       3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
       min(M,N)*(6+4*min(M,N))+max(M,N))
  
  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*min(M,N)),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  n_singular_values_inductance_plasma_middle = min(M,N)
  allocate(svd_s_inductance_plasma_middle(n_singular_values_inductance_plasma_middle),stat=iflag)
  
  ! Matrix is destroyed by LAPACK, so make a copy:
  allocate(A(M,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  A = inductance_plasma_middle
  
  allocate(U(M,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_u_inductance_plasma_middle(num_basis_functions_plasma,n_singular_vectors_to_save),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_v_inductance_plasma_middle(num_basis_functions_middle,n_singular_vectors_to_save),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_u_inductance_plasma_middle_all(num_basis_functions_plasma,num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_v_inductance_plasma_middle_all(num_basis_functions_middle,num_basis_functions_middle),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (save_vectors_in_uv_format) then
     allocate(svd_u_inductance_plasma_middle_uv(nu_plasma*nv_plasma,n_singular_vectors_to_save),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(svd_v_inductance_plasma_middle_uv(nu_middle*nv_middle,n_singular_vectors_to_save),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
  end if
  
  ! Call LAPACK to do the SVD:
  call DGESDD(JOBZ, M, N, A, LDA, svd_s_inductance_plasma_middle, U, LDU, &
       VT, LDVT, WORK, LWORK, IWORK, INFO)
  
  if (INFO==0) then
     print *,"SVD (DGESDD) successful."
     if (n_singular_values_inductance_plasma_middle<5) then
        print *,"Singular values:",svd_s_inductance_plasma_middle
     else
        print *,"First 5 singular values:",svd_s_inductance_plasma_middle(1:5)
        print *,"Last 5 singular values:", &
             svd_s_inductance_plasma_middle(n_singular_values_inductance_plasma_middle-4:n_singular_values_inductance_plasma_middle)
     end if
  else if (INFO>0) then
     print *,"Error in SVD (DGESDD): Did not converge."
     allSVDsSucceeded = .false.
  else
     print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
     allSVDsSucceeded = .false.
  end if
  
  call system_clock(tic1)
  ! Convert singular vectors from basis functions to functions of (u,v):
  svd_u_inductance_plasma_middle = U(:,1:n_singular_vectors_to_save)
  svd_v_inductance_plasma_middle = transpose(VT(1:n_singular_vectors_to_save,:))
  svd_u_inductance_plasma_middle_all = U
  svd_v_inductance_plasma_middle_all = transpose(VT)
  if (save_vectors_in_uv_format) then
     svd_u_inductance_plasma_middle_uv = matmul(basis_functions_plasma, svd_u_inductance_plasma_middle)
     svd_v_inductance_plasma_middle_uv = matmul(basis_functions_middle, svd_v_inductance_plasma_middle)
  end if
  Bnormal_from_1_over_R_field_inductance = matmul(Bnormal_from_1_over_R_field,U)
  Bnormal_from_const_v_coils_inductance  = matmul(Bnormal_from_const_v_coils, U)
  Bnormal_from_plasma_current_inductance = matmul(Bnormal_from_plasma_current,U)
  call system_clock(toc1)
  print *,"  Final matmuls: ",real(toc1-tic1)/countrate," sec."

  ! Compute dominant (m,n) for each left singular vector.
  allocate(svd_u_inductance_plasma_middle_dominant_m(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_u_inductance_plasma_middle_dominant_n(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  do i = 1,num_basis_functions_plasma
     index = 1
     val = abs(U(1,i))
     do j = 2,num_basis_functions_plasma
        if (abs(U(j,i))>val) then
           index = j
           val = abs(U(j,i))
        end if
     end do
     if (index > mnmax_plasma) then
        svd_u_inductance_plasma_middle_dominant_m(i) = xm_plasma(index-mnmax_plasma)
        svd_u_inductance_plasma_middle_dominant_n(i) = xn_plasma(index-mnmax_plasma)
     else
        svd_u_inductance_plasma_middle_dominant_m(i) = xm_plasma(index)
        svd_u_inductance_plasma_middle_dominant_n(i) = xn_plasma(index)
     end if
  end do

  deallocate(A,U,VT,WORK,IWORK)
  
  call system_clock(toc)
  print *,"Done with SVD. Took ",real(toc-tic)/countrate," sec."
  
  !*************************************************************************
  ! End of section related to the middle-to-outer inductance matrix.
  !*************************************************************************
  
  
end subroutine svd_inductance_matrices



    ! Here is the LAPACK documentation for the relevant SVD subroutine:

!!$*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
!!$*                          LWORK, IWORK, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          JOBZ
!!$*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IWORK( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
!!$*      $                   VT( LDVT, * ), WORK( * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGESDD computes the singular value decomposition (SVD) of a real
!!$*> M-by-N matrix A, optionally computing the left and right singular
!!$*> vectors.  If singular vectors are desired, it uses a
!!$*> divide-and-conquer algorithm.
!!$*>
!!$*> The SVD is written
!!$*>
!!$*>      A = U * SIGMA * transpose(V)
!!$*>
!!$*> where SIGMA is an M-by-N matrix which is zero except for its
!!$*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
!!$*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
!!$*> are the singular values of A; they are real and non-negative, and
!!$*> are returned in descending order.  The first min(m,n) columns of
!!$*> U and V are the left and right singular vectors of A.
!!$*>
!!$*> Note that the routine returns VT = V**T, not V.
!!$*>
!!$*> The divide and conquer algorithm makes very mild assumptions about
!!$*> floating point arithmetic. It will work on machines with a guard
!!$*> digit in add/subtract, or on those binary machines without guard
!!$*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!!$*> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!!$*> without guard digits, but we know of none.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] JOBZ
!!$*> \verbatim
!!$*>          JOBZ is CHARACTER*1
!!$*>          Specifies options for computing all or part of the matrix U:
!!$*>          = 'A':  all M columns of U and all N rows of V**T are
!!$*>                  returned in the arrays U and VT;
!!$*>          = 'S':  the first min(M,N) columns of U and the first
!!$*>                  min(M,N) rows of V**T are returned in the arrays U
!!$*>                  and VT;
!!$*>          = 'O':  If M >= N, the first N columns of U are overwritten
!!$*>                  on the array A and all rows of V**T are returned in
!!$*>                  the array VT;
!!$*>                  otherwise, all columns of U are returned in the
!!$*>                  array U and the first M rows of V**T are overwritten
!!$*>                  in the array A;
!!$*>          = 'N':  no columns of U or rows of V**T are computed.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] M
!!$*> \verbatim
!!$*>          M is INTEGER
!!$*>          The number of rows of the input matrix A.  M >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of columns of the input matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the M-by-N matrix A.
!!$*>          On exit,
!!$*>          if JOBZ = 'O',  A is overwritten with the first N columns
!!$*>                          of U (the left singular vectors, stored
!!$*>                          columnwise) if M >= N;
!!$*>                          A is overwritten with the first M rows
!!$*>                          of V**T (the right singular vectors, stored
!!$*>                          rowwise) otherwise.
!!$*>          if JOBZ .ne. 'O', the contents of A are destroyed.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,M).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] S
!!$*> \verbatim
!!$*>          S is DOUBLE PRECISION array, dimension (min(M,N))
!!$*>          The singular values of A, sorted so that S(i) >= S(i+1).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] U
!!$*> \verbatim
!!$*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
!!$*>          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
!!$*>          UCOL = min(M,N) if JOBZ = 'S'.
!!$*>          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
!!$*>          orthogonal matrix U;
!!$*>          if JOBZ = 'S', U contains the first min(M,N) columns of U
!!$*>          (the left singular vectors, stored columnwise);
!!$*>          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDU
!!$*> \verbatim
!!$*>          LDU is INTEGER
!!$*>          The leading dimension of the array U.  LDU >= 1; if
!!$*>          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] VT
!!$*> \verbatim
!!$*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!!$*>          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!!$*>          N-by-N orthogonal matrix V**T;
!!$*>          if JOBZ = 'S', VT contains the first min(M,N) rows of
!!$*>          V**T (the right singular vectors, stored rowwise);
!!$*>          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDVT
!!$*> \verbatim
!!$*>          LDVT is INTEGER
!!$*>          The leading dimension of the array VT.  LDVT >= 1; if
!!$*>          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
!!$*>          if JOBZ = 'S', LDVT >= min(M,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] WORK
!!$*> \verbatim
!!$*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LWORK
!!$*> \verbatim
!!$*>          LWORK is INTEGER
!!$*>          The dimension of the array WORK. LWORK >= 1.
!!$*>          If JOBZ = 'N',
!!$*>            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
!!$*>          If JOBZ = 'O',
!!$*>            LWORK >= 3*min(M,N) + 
!!$*>                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
!!$*>          If JOBZ = 'S' or 'A'
!!$*>            LWORK >= min(M,N)*(6+4*min(M,N))+max(M,N)
!!$*>          For good performance, LWORK should generally be larger.
!!$*>          If LWORK = -1 but other input arguments are legal, WORK(1)
!!$*>          returns the optimal LWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IWORK
!!$*> \verbatim
!!$*>          IWORK is INTEGER array, dimension (8*min(M,N))
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit.
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!!$*>          > 0:  DBDSDC did not converge, updating process failed.
!!$*> \endverbatim
