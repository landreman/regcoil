subroutine auto_regularization_solve

  use global_variables
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
  integer :: stage, next_stage
  logical :: initial_above_target, last_above_target
  real(dp) :: Brendt_a, Brendt_b, Brendt_c, Brendt_fa, Brendt_fb, Brendt_fc, Brendt_d, Brendt_e
  real(dp) :: Brendt_p, Brendt_q, Brendt_r, Brendt_s, Brendt_tol1, Brendt_xm, Brendt_EPS

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  allocate(matrix(num_basis_functions, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(chi2_B(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(chi2_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_Bnormal(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(current_potential(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_mn(num_basis_functions,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_total(ntheta_plasma,nzeta_plasma,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(K2(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_K2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'


  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  print *,"Optimal LWORK:",LWORK
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'  

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  if (general_option==4) then
     stage = 1
  else
     stage = 10
  end if
  exit_code = -1
  do ilambda = 1,nlambda

     ! First, pick the next value of lambda, depending on what stage we are in in the search algorithm:
     ! ------------------------------------------------------------------------------------------------
     select case (stage)
     case (1)
        ! Initial guess for lambda:
        ! guess lambda = chi^2_B / chi^2_K, where the right hand side is evaluated taking the
        ! single-valued part of the current potential to be 0.

        KDifference_x = d_x !- matmul(f_x, solution)
        KDifference_y = d_y !- matmul(f_y, solution)
        KDifference_z = d_z !- matmul(f_z, solution)
        this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
             / norm_normal_coil
        chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)

        Bnormal_total(:,:,ilambda) = & ! (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) +
             Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

        chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
             * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)
        ! chi2_B, chi2_K, and Bnormal_total for this ilambda will be over-written with the real values below.

        lambda(ilambda) = chi2_B(ilambda) / chi2_K(ilambda) / 1000
        next_stage = 2

     case (2)
        ! Vary lambda by factors of 100 until we bracket the target current density.
        next_stage = 2

        if (initial_above_target) then
           ! If the current density from the first solve was too high, then increase lambda
           lambda(ilambda) = lambda(ilambda-1) * 100
        else
           ! If the current density from the first solve was too low, then decrease lambda
           lambda(ilambda) = lambda(ilambda-1) * 0.01
        end if

     case (3)
        ! Now that we've bracketed the target current density, search for the target using Brendt's algorithm.
        next_stage = 3

        if (abs(Brendt_e) >= Brendt_tol1 .and. abs(Brendt_fa)>abs(Brendt_fb)) then
           ! Attempt inverse quadratic interpolation
           Brendt_s = Brendt_fb / Brendt_fa
           if (Brendt_a == Brendt_c) then
              Brendt_p = 2.0*Brendt_xm*Brendt_s
              Brendt_q = 1.0-Brendt_s
           else
              Brendt_q = Brendt_fa / Brendt_fc
              Brendt_r = Brendt_fb / Brendt_fc
              Brendt_p = Brendt_s*(2.0*Brendt_xm*Brendt_q*(Brendt_q-Brendt_r)-(Brendt_b-Brendt_a)*(Brendt_r-1.0))
              Brendt_q = (Brendt_q-1.0)*(Brendt_r-1.0)*(Brendt_s-1.0)
           end if
           if (Brendt_p > 0) Brendt_q = -Brendt_q ! Check whether in bounds
           Brendt_p = abs(Brendt_p)
           if (2.0*Brendt_p < min(3.0*Brendt_xm*Brendt_q-abs(Brendt_tol1*Brendt_q),abs(Brendt_e*Brendt_q))) then
              ! Accept interpolation
              Brendt_e = Brendt_d
              Brendt_d = Brendt_p / Brendt_q
           else
              ! Interpolation failed, so use bisection
              Brendt_d = Brendt_xm
              Brendt_e = Brendt_d
           end if
        else
           ! Bounds are decreasing too slowly, so use bisection
           Brendt_d = Brendt_xm
           Brendt_e = Brendt_d
        end if
        ! Move last best guess to a.
        Brendt_a = Brendt_b
        Brendt_fa = Brendt_fb
        if (abs(Brendt_d) > Brendt_tol1) then
           ! Evaluate new trial root
           Brendt_b = Brendt_b + Brendt_d
        else
           Brendt_b = Brendt_b + sign(Brendt_tol1,Brendt_xm)
        end if

        lambda(ilambda) = exp(Brendt_b)

     case (10)
        ! Try lambda = infinity
        lambda(ilambda) = 1.0d200
        next_stage = 11

     case (11)
        lambda(ilambda) = 0
        next_stage = 1

     case default
        print *,"Invalid stage in auto_regularization_solve:",stage
        stop
     end select

     print "(a,es10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of at most ",nlambda,")"
     call system_clock(tic,countrate)

     ! Done choosing the next lambda. Now comes the main solve.
     ! ------------------------------------------------------------------------------------------------

     if (stage==10) then
        matrix = matrix_K
        RHS    =    RHS_K
     else
        matrix = matrix_B + lambda(ilambda) * matrix_K
        RHS    =    RHS_B + lambda(ilambda) *    RHS_K
     end if

     call system_clock(toc)
     print *,"  Additions: ",real(toc-tic)/countrate," sec."
     call system_clock(tic)

     ! Compute solution = matrix \ RHS.
     ! Use LAPACK's DSYSV since matrix is symmetric.
     ! Note: RHS will be over-written with the solution.
     call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, LWORK, INFO)
     if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
        !stop
     end if
     solution = RHS

     call system_clock(toc)
     print *,"  DSYSV: ",real(toc-tic)/countrate," sec."
     call system_clock(tic)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Now compute diagnostics
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
     print "(a,es10.3,a,es10.3)","   chi2_B:",chi2_B(ilambda),",  chi2_K:",chi2_K(ilambda)
     print "(a,es10.3,a,es10.3,a,es10.3)","   max(B_n):",max_Bnormal(ilambda),",  max(K):",max_K(ilambda),",  rms K:",sqrt(chi2_K(ilambda)/area_coil)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Done computing diagnostics.
     ! Now analyze the results to determine what to do next.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     last_above_target = (target_function(ilambda) > current_density_target)
     if (stage==1) initial_above_target = last_above_target
     if (stage==2 .and. (last_above_target .neqv. initial_above_target)) then
        ! If we've bracketed the target, move on to stage 3.
        next_stage = 3  
        print *,"Target current density has been bracketed."
     end if
     if (stage==10 .and. last_above_target) then
        print *,"*******************************************************************************"
        print *,"*******************************************************************************"
        print *,"Error! The current_density_target you have set is not achievable because"
        print *,"it is too low."
        print *,"*******************************************************************************"
        print *,"*******************************************************************************"
        Nlambda = ilambda
        exit_code = -2
        exit
     end if
     if (stage==11 .and. (.not. last_above_target)) then
        print *,"*******************************************************************************"
        print *,"*******************************************************************************"
        print *,"Error! The current_density_target you have set is not achievable because"
        print *,"it is too high."
        print *,"*******************************************************************************"
        print *,"*******************************************************************************"
        Nlambda = ilambda
        exit_code = -3
        exit
     end if
     if (stage==2 .and. next_stage == 3) then
        ! Initialize Brendt's algorithm for root-finding
        Brendt_a = log(lambda(ilambda-1))
        Brendt_b = log(lambda(ilambda))
        Brendt_fa = log(target_function(ilambda-1) / current_density_target)
        Brendt_fb = log(target_function(ilambda) / current_density_target)
        Brendt_c = Brendt_b
        Brendt_fc = Brendt_fb
        Brendt_d = Brendt_b - Brendt_a
        Brendt_e = Brendt_d
     end if
     if (stage==3) Brendt_fb = log(target_function(ilambda) / current_density_target)
     if (next_stage==3) then
        ! Analyze the most recent diagnostics for Brendt's algorithm.
        if ((Brendt_fb > 0 .and. Brendt_fc > 0) .or. (Brendt_fb < 0 .and. Brendt_fc < 0)) then
           Brendt_c = Brendt_a
           Brendt_fc = Brendt_fa
           Brendt_d = Brendt_b - Brendt_a
           Brendt_e = Brendt_d
        end if
        if (abs(Brendt_fc) < abs(Brendt_fb)) then
           Brendt_a = Brendt_b
           Brendt_b = Brendt_c
           Brendt_c = Brendt_a
           Brendt_fa = Brendt_fb
           Brendt_fb = Brendt_fc
           Brendt_fc = Brendt_fa
        end if
        Brendt_EPS = 1d-15
        Brendt_tol1 = 2.0*Brendt_EPS*abs(Brendt_b) + 0.5*lambda_search_tolerance
        Brendt_xm = 0.5*(Brendt_c - Brendt_b)
        if (abs(Brendt_xm) <= Brendt_tol1 .or. (Brendt_fb==0)) then
           ! We met the requested tolerance
           print *,"Requested tolerance has been met."
           exit_code=0
           Nlambda = ilambda
           exit
        end if
     end if
     stage = next_stage
  end do

  chi2_B_target = chi2_B(Nlambda)

  if (exit_code == -1) then
     print *,"*******************************************************************************"
     print *,"*******************************************************************************"
     print *,"The lambda search did not converge within Nlambda iterations!"
     print *,"*******************************************************************************"
     print *,"*******************************************************************************"
  end if

contains
 
  function target_function(jlambda)

    implicit none

    integer, intent(in) :: jlambda
    real(dp) :: target_function

    target_function = 0
    select case (target_option)
    case (1)
       ! maximum current density
       target_function =  max_K(jlambda)
    case (2)
       ! root-mean-squared current density
       target_function = sqrt(chi2_K(jlambda) / area_coil)
    case default
       print *,"Invalid target_option: ",target_option
       stop
    end select

  end function target_function

end subroutine auto_regularization_solve

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
