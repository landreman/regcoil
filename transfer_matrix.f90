! Documentation of LAPACK's SVD subroutine DGESDD is copied at the end of this file for convenience.

subroutine transfer_matrix

  use global_variables
  use stel_kinds

  implicit none

  integer :: whichThreshold, iflag, tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: Mpc_V, tempMatrix, transferMatrix
  real(dp) :: threshold, val
  integer :: max_singular_value_index, i, j, index

  ! Variables needed by LAPACK:
  character :: JOBZ
  integer :: INFO, LDA, LDU, LDVT, LWORK, M, N
  real(dp), dimension(:,:), allocatable :: U, VT
  real(dp), dimension(:), allocatable :: WORK, svd_s_transferMatrix_single
  integer, dimension(:), allocatable :: IWORK

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA='N', TRANSB='N'
  integer :: M_DGEMM, N_DGEMM, K_DGEMM, LDA_DGEMM, LDB_DGEMM, LDC_DGEMM
  real(dp) :: ALPHA=1, BETA=0

  allocate(Mpc_V(num_basis_functions_plasma, num_basis_functions_outer), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(tempMatrix(num_basis_functions_outer, num_basis_functions_middle), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(transferMatrix(num_basis_functions_plasma, num_basis_functions_middle), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(n_singular_values_retained(n_pseudoinverse_thresholds), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(overlap_plasma(num_basis_functions_plasma,num_basis_functions_plasma), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(overlap_middle(num_basis_functions_middle,num_basis_functions_middle), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  print *,"Forming matrix product M_po * V_mo"
  call system_clock(tic,countrate)
  ! Slow but simpler method using "matmul":
  !Mpc_V = matmul(inductance_plasma_outer, svd_v_inductance_middle_outer)

  ! A = inductance_plasma
  ! B = svd_v_inductance_middle_outer
  ! C = Mpc_V
  M_DGEMM = size(inductance_plasma_outer,1) ! # rows of A
  N_DGEMM = size(svd_v_inductance_middle_outer,2) ! # cols of B
  K_DGEMM = size(inductance_plasma_outer,2) ! Common dimension of A and B
  LDA_DGEMM = M_DGEMM
  LDB_DGEMM = K_DGEMM
  LDC_DGEMM = M_DGEMM
  TRANSA = 'N' ! No transpose
  TRANSB = 'N' ! No transpose
  call DGEMM(TRANSA,TRANSB,M_DGEMM,N_DGEMM,K_DGEMM,ALPHA,inductance_plasma_outer,LDA_DGEMM,&
       svd_v_inductance_middle_outer,LDB_DGEMM,BETA,Mpc_V,LDC_DGEMM)

  call system_clock(toc)
  print *,"Done. Took",real(toc-tic)/countrate," sec."


  ! Before beginning the loop over thresholds, do the initialization that will be needed for the SVD.

  JOBZ='A'  ! Compute all the singular vectors.  Consider changing this if performance is an issue.
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
  if (iflag .ne. 0) stop 'Allocation error A!'
  allocate(IWORK(8*min(M,N)),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error B!'
     
  n_singular_values_transferMatrix = min(M,N)
  allocate(svd_s_transferMatrix(n_singular_values_transferMatrix,n_pseudoinverse_thresholds),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error C!'
  allocate(svd_s_transferMatrix_single(n_singular_values_transferMatrix),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error D!'
  
  allocate(U(M,M),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error E!'
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error F!'
  allocate(svd_u_transferMatrix(num_basis_functions_plasma,n_singular_vectors_to_save,n_pseudoinverse_thresholds),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error G1!'
  allocate(svd_v_transferMatrix(num_basis_functions_middle,n_singular_vectors_to_save,n_pseudoinverse_thresholds),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error H1!'
  if (save_vectors_in_uv_format) then
     allocate(svd_u_transferMatrix_uv(nu_plasma*nv_plasma,n_singular_vectors_to_save,n_pseudoinverse_thresholds),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error G1!'
     allocate(svd_v_transferMatrix_uv(nu_middle*nv_middle,n_singular_vectors_to_save,n_pseudoinverse_thresholds),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error H1!'
  end if

  allocate(svd_u_transferMatrix_dominant_m(num_basis_functions_plasma,n_pseudoinverse_thresholds),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(svd_u_transferMatrix_dominant_n(num_basis_functions_plasma,n_pseudoinverse_thresholds),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  ! Done with initialization. Now begin the loop over thresholds:

  do whichThreshold = 1,n_pseudoinverse_thresholds
     threshold = pseudoinverse_thresholds(whichThreshold)
     print *,"Beginning to build transfer matrix",whichThreshold,"of",n_pseudoinverse_thresholds
     print *,"  Beginning pseudoinverse with relative threshold",threshold

     ! Determine how many singular values from M_mo to include
     max_singular_value_index = n_singular_values_inductance_middle_outer
     do i = 1,n_singular_values_inductance_middle_outer
        if (svd_s_inductance_middle_outer(i)/svd_s_inductance_middle_outer(1) < threshold) then
           max_singular_value_index = i-1
           exit
        end if
     end do
     n_singular_values_retained(whichThreshold) = max_singular_value_index
     print *,"  Retaining",max_singular_value_index," of ",n_singular_values_inductance_middle_outer,&
          "singular values in the pseudoinverse."

     ! Next, form the matrix product s_mo^-1 * u_mo^-1
     call system_clock(tic)
     tempMatrix = 0
     do i = 1,max_singular_value_index
        tempMatrix(i,:) = (1/svd_s_inductance_middle_outer(i))*svd_uT_inductance_middle_outer(i,:)
     end do
     call system_clock(toc)
     print *,"  Assembly of s_mo^-1 * u_mo^-1 took",real(toc-tic)/countrate," sec."

     print *,"  Beginning final assembly of transfer matrix."
     call system_clock(tic)
     
     ! Transparent but slow method using "matmul":
     !transferMatrix = matmul(Mpc_V, tempMatrix)

     ! A = inductance_plasma_outer
     ! B = svd_v_inductance_middle_outer
     ! C = Mpc_V
     M_DGEMM = size(Mpc_V,1) ! # rows of A
     N_DGEMM = size(tempMatrix,2) ! # cols of B
     K_DGEMM = size(Mpc_V,2) ! Common dimension of A and B
     LDA_DGEMM = M_DGEMM
     LDB_DGEMM = K_DGEMM
     LDC_DGEMM = M_DGEMM
     TRANSA = 'N' ! No transpose
     TRANSB = 'N' ! No transpose
     call DGEMM(TRANSA,TRANSB,M_DGEMM,N_DGEMM,K_DGEMM,ALPHA,Mpc_V,LDA_DGEMM, &
          tempMatrix,LDB_DGEMM,BETA,transferMatrix,LDC_DGEMM)


     call system_clock(toc)
     print *,"  Done. Took",real(toc-tic)/countrate," sec."


     ! **********************************************************
     ! Done assembling the transfer matrix. Now do the SVD.
     ! **********************************************************

     print *,"  Beginning SVD of the transfer matrix."
     call system_clock(tic,countrate)  
  
     ! Call LAPACK to do the SVD. This destroys transferMatrix!
     call DGESDD(JOBZ, M, N, transferMatrix, LDA, svd_s_transferMatrix_single, U, LDU, &
          VT, LDVT, WORK, LWORK, IWORK, INFO)
  
     if (INFO==0) then
        print *,"  SVD (DGESDD) successful."
        if (n_singular_values_transferMatrix<5) then
           print *,"  Singular values:",svd_s_transferMatrix_single
        else
           print *,"  First 5 singular values:",svd_s_transferMatrix_single(1:5)
           print *,"  Last 5 singular values:", &
                svd_s_transferMatrix_single(n_singular_values_transferMatrix-4:n_singular_values_transferMatrix)
        end if
     else if (INFO>0) then
        print *,"  Error in SVD (DGESDD): Did not converge."
        allSVDsSucceeded = .false.
     else
        print *,"  Error in SVD (DGESDD): Argument",INFO," was invalid."
        allSVDsSucceeded = .false.
     end if

     call system_clock(toc)
     print *,"  Done with SVD. Took ",real(toc-tic)/countrate," sec."
  
     svd_s_transferMatrix(:,whichThreshold) = svd_s_transferMatrix_single

     call system_clock(tic)
     svd_u_transferMatrix(:,:,whichThreshold) = U(:,1:n_singular_vectors_to_save)
     svd_v_transferMatrix(:,:,whichThreshold) = transpose(VT(1:n_singular_vectors_to_save,:))
     if (save_vectors_in_uv_format) then
        svd_u_transferMatrix_uv(:,:,whichThreshold) = matmul(basis_functions_plasma, svd_u_transferMatrix(:,:,whichThreshold))
        svd_v_transferMatrix_uv(:,:,whichThreshold) = matmul(basis_functions_middle, svd_v_transferMatrix(:,:,whichThreshold))
     end if
     call system_clock(toc)
     if (whichThreshold==1) then
        Bnormal_from_1_over_R_field_transfer = matmul(Bnormal_from_1_over_R_field,U)
        Bnormal_from_const_v_coils_transfer  = matmul(Bnormal_from_const_v_coils, U)
        Bnormal_from_plasma_current_transfer = matmul(Bnormal_from_plasma_current,U)
     end if
     print *,"  Final matmuls: ",real(toc-tic)/countrate," sec."

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
           svd_u_transferMatrix_dominant_m(i,whichThreshold) = xm_plasma(index-mnmax_plasma)
           svd_u_transferMatrix_dominant_n(i,whichThreshold) = xn_plasma(index-mnmax_plasma)
        else
           svd_u_transferMatrix_dominant_m(i,whichThreshold) = xm_plasma(index)
           svd_u_transferMatrix_dominant_n(i,whichThreshold) = xn_plasma(index)
        end if
     end do

     if (whichThreshold==1) then
        call system_clock(tic)
        if (zero_first_transfer_vector_in_overlap) then
           U(:,1) = 0
           VT(1,:) = 0
        end if

        ! Transparent but slow method using matmul:
!!$        overlap_middle = matmul(VT, svd_v_inductance_plasma_middle_all)
!!$        overlap_plasma = matmul(transpose(U), svd_u_inductance_plasma_middle_all)

        ! Carry out the matrix-matrix multiplication C = A * B
        ! A = VT
        ! B = svd_v_inductance_plasma_middle_all
        ! C = overlap_middle
        M_DGEMM = num_basis_functions_middle ! # rows of A
        N_DGEMM = num_basis_functions_middle ! # cols of B
        K_DGEMM = num_basis_functions_middle ! Common dimension of A and B
        LDA_DGEMM = M_DGEMM
        LDB_DGEMM = K_DGEMM
        LDC_DGEMM = M_DGEMM
        TRANSA = 'N' ! No transpose
        TRANSB = 'N' ! No transpose
        call DGEMM(TRANSA,TRANSB,M_DGEMM,N_DGEMM,K_DGEMM,ALPHA,VT,LDA_DGEMM, &
             svd_v_inductance_plasma_middle_all,LDB_DGEMM,BETA,overlap_middle,LDC_DGEMM)

        ! Carry out the matrix-matrix multiplication C = A * B
        ! A = U
        ! B = svd_u_inductance_plasma_middle_all
        ! C = overlap_plasma
        M_DGEMM = num_basis_functions_plasma ! # rows of A
        N_DGEMM = num_basis_functions_plasma ! # cols of B
        K_DGEMM = num_basis_functions_plasma ! Common dimension of A and B
        LDA_DGEMM = M_DGEMM
        LDB_DGEMM = K_DGEMM
        LDC_DGEMM = M_DGEMM
        TRANSA = 'T' ! DO take transpose
        TRANSB = 'N' ! No transpose
        call DGEMM(TRANSA,TRANSB,M_DGEMM,N_DGEMM,K_DGEMM,ALPHA,U,LDA_DGEMM, &
             svd_u_inductance_plasma_middle_all,LDB_DGEMM,BETA,overlap_plasma,LDC_DGEMM)

        call system_clock(toc)
        print *,"  Compute overlap: ",real(toc-tic)/countrate," sec."
     end if

  end do

  deallocate(U,VT,WORK,IWORK)

end subroutine transfer_matrix

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
