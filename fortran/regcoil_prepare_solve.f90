subroutine regcoil_prepare_solve(prob)

  use regcoil_variables, only: regcoil_t, target_option_max_K_lse, target_option_lp_norm_K

  implicit none


  type(regcoil_t), intent(inout) :: prob
  integer :: iflag

    associate ( &
       ntheta_plasma => prob%plasma%ntheta_plasma, &
       nzeta_plasma => prob%plasma%nzeta_plasma, &
       ntheta_coil => prob%coil%ntheta_coil, &
       nzeta_coil => prob%coil%nzeta_coil, &
       verbose => prob%input%verbose, &
       nlambda => prob%input%nlambda, &
       target_option => prob%input%target_option, &
       num_basis_functions => prob%work%num_basis_functions, &
       LAPACK_INFO => prob%work%LAPACK_INFO, &
       LAPACK_LWORK => prob%work%LAPACK_LWORK &
       )
  if (allocated(prob%work%matrix)) deallocate(prob%work%matrix)
  allocate(prob%work%matrix(num_basis_functions, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 1!'

  if (allocated(prob%work%RHS)) deallocate(prob%work%RHS)
  allocate(prob%work%RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 2!'

  if (allocated(prob%work%solution)) deallocate(prob%work%solution)
  allocate(prob%work%solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 3!'

  if (allocated(prob%work%LAPACK_WORK)) deallocate(prob%work%LAPACK_WORK)
  allocate(prob%work%LAPACK_WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 4!'

  if (allocated(prob%work%LAPACK_IPIV)) deallocate(prob%work%LAPACK_IPIV)
  allocate(prob%work%LAPACK_IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 5!'

  if (allocated(prob%output%chi2_B)) deallocate(prob%output%chi2_B)
  allocate(prob%output%chi2_B(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 6!'

  if (allocated(prob%output%chi2_K)) deallocate(prob%output%chi2_K)
  allocate(prob%output%chi2_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 7!'

  if (allocated(prob%output%chi2_Laplace_Beltrami)) deallocate(prob%output%chi2_Laplace_Beltrami)
  allocate(prob%output%chi2_Laplace_Beltrami(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 7!'

  if (allocated(prob%output%max_Bnormal)) deallocate(prob%output%max_Bnormal)
  allocate(prob%output%max_Bnormal(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 8!'

  if (allocated(prob%output%max_K)) deallocate(prob%output%max_K)
  allocate(prob%output%max_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 9!'

  if (allocated(prob%output%current_potential)) deallocate(prob%output%current_potential)
  allocate(prob%output%current_potential(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 10!'

  if (allocated(prob%output%single_valued_current_potential_thetazeta)) &
        deallocate(prob%output%single_valued_current_potential_thetazeta)
  allocate(prob%output%single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 11!'

  if (allocated(prob%work%this_current_potential)) deallocate(prob%work%this_current_potential)
  allocate(prob%work%this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 12!'

  if (allocated(prob%output%single_valued_current_potential_mn)) &
        deallocate(prob%output%single_valued_current_potential_mn)
  allocate(prob%output%single_valued_current_potential_mn(num_basis_functions,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 13!'

  if (allocated(prob%output%Bnormal_total)) deallocate(prob%output%Bnormal_total)
  allocate(prob%output%Bnormal_total(ntheta_plasma,nzeta_plasma,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 14!'

  if (allocated(prob%output%K2)) deallocate(prob%output%K2)
  allocate(prob%output%K2(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 15!'

  if (allocated(prob%output%Laplace_Beltrami2)) deallocate(prob%output%Laplace_Beltrami2)
  allocate(prob%output%Laplace_Beltrami2(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 15!'

  if (allocated(prob%work%KDifference_x)) deallocate(prob%work%KDifference_x)
  allocate(prob%work%KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 16!'

  if (allocated(prob%work%KDifference_y)) deallocate(prob%work%KDifference_y)
  allocate(prob%work%KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 17!'

  if (allocated(prob%work%KDifference_z)) deallocate(prob%work%KDifference_z)
  allocate(prob%work%KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 18!'

  if (allocated(prob%work%KDifference_Laplace_Beltrami)) deallocate(prob%work%KDifference_Laplace_Beltrami)
  allocate(prob%work%KDifference_Laplace_Beltrami(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 18!'

  if (allocated(prob%work%this_K2_times_N)) deallocate(prob%work%this_K2_times_N)
  allocate(prob%work%this_K2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 19!'

  if (allocated(prob%work%this_Laplace_Beltrami2_times_N)) deallocate(prob%work%this_Laplace_Beltrami2_times_N)
  allocate(prob%work%this_Laplace_Beltrami2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 19!'

  if (trim(target_option)==target_option_max_K_lse) then
    if (allocated(prob%output%max_K_lse)) deallocate(prob%output%max_K_lse)
    allocate(prob%output%max_K_lse(nlambda), stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 20!'
  end if

  if (trim(target_option)==target_option_lp_norm_K) then
    if (allocated(prob%output%lp_norm_K)) deallocate(prob%output%lp_norm_K)
    allocate(prob%output%lp_norm_K(nlambda), stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_prepare_solve Allocation error 20!'
  end if

  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, prob%work%matrix, num_basis_functions, prob%work%LAPACK_IPIV, prob%work%RHS, num_basis_functions, prob%work%LAPACK_WORK, -1, LAPACK_INFO)
  LAPACK_LWORK = int(prob%work%LAPACK_WORK(1))
  if (verbose) print *,"Optimal LWORK:",LAPACK_LWORK
  deallocate(prob%work%LAPACK_WORK)
  allocate(prob%work%LAPACK_WORK(LAPACK_LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_prepare_solve LAPACK error!'


  end associate
end subroutine regcoil_prepare_solve

    ! Here is the LAPACK documentation for solving a symmetric linear system:

!!$
!!$
!!$*> \brief <b> DSYSV computes the solution to system of linear equations A * X = B for SY matrices</b>
!!$*
!!$*  =========== DOCUMENTATION ===========
!!$*
!!$* Online html documentation available at 
!!$*            http://www.netlib.org/lapack/explore-html/ 
!!$*
!!$*> \htmlonly
!!$*> Download DSYSV + dependencies 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [TGZ]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [ZIP]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [TXT]</a>
!!$*> \endhtmlonly 
!!$*
!!$*  Definition:
!!$*  ===========
!!$*
!!$*       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
!!$*                         LWORK, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          UPLO
!!$*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IPIV( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DSYSV computes the solution to a real system of linear equations
!!$*>    A * X = B,
!!$*> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!!$*> matrices.
!!$*>
!!$*> The diagonal pivoting method is used to factor A as
!!$*>    A = U * D * U**T,  if UPLO = 'U', or
!!$*>    A = L * D * L**T,  if UPLO = 'L',
!!$*> where U (or L) is a product of permutation and unit upper (lower)
!!$*> triangular matrices, and D is symmetric and block diagonal with
!!$*> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!!$*> used to solve the system of equations A * X = B.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] UPLO
!!$*> \verbatim
!!$*>          UPLO is CHARACTER*1
!!$*>          = 'U':  Upper triangle of A is stored;
!!$*>          = 'L':  Lower triangle of A is stored.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of linear equations, i.e., the order of the
!!$*>          matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] NRHS
!!$*> \verbatim
!!$*>          NRHS is INTEGER
!!$*>          The number of right hand sides, i.e., the number of columns
!!$*>          of the matrix B.  NRHS >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$*>          N-by-N upper triangular part of A contains the upper
!!$*>          triangular part of the matrix A, and the strictly lower
!!$*>          triangular part of A is not referenced.  If UPLO = 'L', the
!!$*>          leading N-by-N lower triangular part of A contains the lower
!!$*>          triangular part of the matrix A, and the strictly upper
!!$*>          triangular part of A is not referenced.
!!$*>
!!$*>          On exit, if INFO = 0, the block diagonal matrix D and the
!!$*>          multipliers used to obtain the factor U or L from the
!!$*>          factorization A = U*D*U**T or A = L*D*L**T as computed by
!!$*>          DSYTRF.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IPIV
!!$*> \verbatim
!!$*>          IPIV is INTEGER array, dimension (N)
!!$*>          Details of the interchanges and the block structure of D, as
!!$*>          determined by DSYTRF.  If IPIV(k) > 0, then rows and columns
!!$*>          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!!$*>          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!!$*>          then rows and columns k-1 and -IPIV(k) were interchanged and
!!$*>          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!!$*>          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!!$*>          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!!$*>          diagonal block.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*>          On entry, the N-by-NRHS right hand side matrix B.
!!$*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>          The leading dimension of the array B.  LDB >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] WORK
!!$*> \verbatim
!!$*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LWORK
!!$*> \verbatim
!!$*>          LWORK is INTEGER
!!$*>          The length of WORK.  LWORK >= 1, and for best performance
!!$*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!!$*>          DSYTRF.
!!$*>          for LWORK < N, TRS will be done with Level BLAS 2
!!$*>          for LWORK >= N, TRS will be done with Level BLAS 3
!!$*>
!!$*>          If LWORK = -1, then a workspace query is assumed; the routine
!!$*>          only calculates the optimal size of the WORK array, returns
!!$*>          this value as the first entry of the WORK array, and no error
!!$*>          message related to LWORK is issued by XERBLA.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0: successful exit
!!$*>          < 0: if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!!$*>               has been completed, but the block diagonal matrix D is
!!$*>               exactly singular, so the solution could not be computed.
!!$*> \endverbatim
