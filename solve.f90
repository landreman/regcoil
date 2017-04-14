subroutine solve

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

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  ! Needed for sensitivity calculation
  real(dp), dimension(:,:), allocatable :: dKDifferencedrmnc, dKDifferencedrmns, dKDifferencedzmnc, dKDifferencedzmns
  real(dp), dimension(:,:), allocatable :: rmncterm1, rmnsterm1, zmncterm1, zmnsterm1, rmncterm2, rmnsterm2, zmncterm2, zmnsterm2
  real(dp), dimension(:,:), allocatable :: dBnormaldrmnc, dBnormaldrmns, dBnormaldzmnc, dBnormaldzmns
  integer :: imn

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
! Sensitivity matrices
  if (sensitivity_option > 1) then
    print *,"Allocating sensitivity matrices."
    allocate(dchi2Kdrmnc(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Kdrmns(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Kdzmnc(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Kdzmns(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Bdrmnc(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Bdrmns(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Bdzmnc(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dchi2Bdzmns(mnmax_sensitivity, nlambda),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dKDifferencedrmnc(3,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dKDifferencedrmns(3,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dKDifferencedzmnc(3,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dKDifferencedzmns(3,ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(rmncterm1(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(rmnsterm1(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(zmncterm1(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(zmnsterm1(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(rmncterm2(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(rmnsterm2(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(zmncterm2(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(zmnsterm2(ntheta_coil,nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dBnormaldrmnc(ntheta_plasma,nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dBnormaldrmns(ntheta_plasma,nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dBnormaldzmnc(ntheta_plasma,nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dBnormaldzmns(ntheta_plasma,nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    print *,"Allocation complete."
  endif

  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  print *,"Optimal LWORK:",LWORK
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  do ilambda = 1,nlambda
     print "(a,es10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of ",nlambda,")"
     call system_clock(tic,countrate)

     matrix = matrix_B + lambda(ilambda) * matrix_K
     RHS    =    RHS_B + lambda(ilambda) *    RHS_K

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

     if (sensitivity_option > 1) then

       print *,"Computing chi2_K sensitivity."
       ! sensitivity computation
       ! dddrmnc(3, mnmax_sensitivity, ntheta_coil*nzeta_coil)
       ! dfdrmnc(3, mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions)
       ! f_x(ntheta_coil*nzeta_coil, num_basis_functions)
       ! norm_normal_coil(ntheta_coil,nzeta_coil)
       do imn=1,mnmax_sensitivity
         dKDifferencedrmnc(1,:) = dddrmnc(1,imn,:)-matmul(dfxdrmnc(imn,:,:), solution)
         dKDifferencedrmnc(2,:) = dddrmnc(2,imn,:)-matmul(dfydrmnc(imn,:,:), solution)
         dKDifferencedrmnc(3,:) = dddrmnc(3,imn,:)-matmul(dfzdrmnc(imn,:,:), solution)
         dKDifferencedrmns(1,:) = dddrmns(1,imn,:)-matmul(dfxdrmns(imn,:,:), solution)
         dKDifferencedrmns(2,:) = dddrmns(2,imn,:)-matmul(dfydrmns(imn,:,:), solution)
         dKDifferencedrmns(3,:) = dddrmns(3,imn,:)-matmul(dfzdrmns(imn,:,:), solution)
         dKDifferencedzmnc(1,:) = dddrmnc(1,imn,:)-matmul(dfxdzmnc(imn,:,:), solution)
         dKDifferencedzmnc(2,:) = dddrmnc(2,imn,:)-matmul(dfydzmnc(imn,:,:), solution)
         dKDifferencedzmnc(3,:) = dddrmnc(3,imn,:)-matmul(dfzdzmnc(imn,:,:), solution)
         dKDifferencedzmns(1,:) = dddzmns(1,imn,:)-matmul(dfxdzmns(imn,:,:), solution)
         dKDifferencedzmns(2,:) = dddzmns(2,imn,:)-matmul(dfydzmns(imn,:,:), solution)
         dKDifferencedzmns(3,:) = dddzmns(3,imn,:)-matmul(dfzdzmns(imn,:,:), solution)
         rmncterm1(:,:) = -dnorm_normaldrmnc(imn,:,:)*this_K2_times_N/norm_normal_coil
         rmnsterm1(:,:) = -dnorm_normaldrmns(imn,:,:)*this_K2_times_N/norm_normal_coil
         zmncterm1(:,:) = -dnorm_normaldzmnc(imn,:,:)*this_K2_times_N/norm_normal_coil
         zmnsterm1(:,:) = -dnorm_normaldzmns(imn,:,:)*this_K2_times_N/norm_normal_coil
         rmncterm2(:,:) = reshape(KDifference_x*dKDifferencedrmnc(1,:) + KDifference_y*dKDifferencedrmnc(2,:) &
           + KDifference_z*dKDifferencedrmnc(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
         rmnsterm2(:,:) = reshape(KDifference_x*dKDifferencedrmns(1,:) + KDifference_y*dKDifferencedrmns(2,:) &
           + KDifference_z*dKDifferencedrmns(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
         zmncterm2(:,:) = reshape(KDifference_x*dKDifferencedzmnc(1,:) + KDifference_y*dKDifferencedzmnc(2,:) &
           + KDifference_z*dKDifferencedzmnc(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
         zmnsterm2(:,:) = reshape(KDifference_x*dKDifferencedzmns(1,:) + KDifference_y*dKDifferencedzmns(2,:) &
           + KDifference_z*dKDifferencedzmns(3,:),(/ ntheta_coil, nzeta_coil/))*(2/norm_normal_coil)
         dchi2Kdrmnc(imn,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(rmncterm1) + sum(rmncterm2))
         dchi2Kdrmns(imn,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(rmnsterm1) + sum(rmnsterm2))
         dchi2Kdzmnc(imn,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(zmncterm1) + sum(zmncterm2))
         dchi2Kdzmns(imn,ilambda) = nfp*dtheta_coil*dzeta_coil*(sum(zmnsterm1) + sum(zmnsterm2))
       enddo
     endif

     Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
          + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

     max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
     max_K(ilambda) = sqrt(maxval(K2    (:,:,ilambda)))

     chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
          * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)

     if (sensitivity_option > 1) then
       print *,"Computing chi2_B sensitivity."
       ! Compute chi2_B sensitivity
       ! dgdrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions)
       ! solution(num_basis_functions)
       ! dBnormaldrmnc(ntheta_plasma,nzeta_plasma)
       do imn=1,mnmax_sensitivity
         dBnormaldrmnc(:,:) = reshape(matmul(dgdrmnc(imn,:,:),solution),(/ ntheta_plasma, nzeta_plasma /))*norm_normal_plasma
         dBnormaldrmns(:,:) = reshape(matmul(dgdrmns(imn,:,:),solution),(/ ntheta_plasma, nzeta_plasma /))*norm_normal_plasma
         dBnormaldzmnc(:,:) = reshape(matmul(dgdzmnc(imn,:,:),solution),(/ ntheta_plasma, nzeta_plasma /))*norm_normal_plasma
         dBnormaldzmns(:,:) = reshape(matmul(dgdzmns(imn,:,:),solution),(/ ntheta_plasma, nzeta_plasma /))*norm_normal_plasma
         dchi2Bdrmnc(imn,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldrmnc*norm_normal_plasma)
         dchi2Bdrmns(imn,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldrmns*norm_normal_plasma)
         dchi2Bdzmnc(imn,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldzmnc*norm_normal_plasma)
         dchi2Bdzmns(imn,ilambda) = 2*nfp*dtheta_plasma*dzeta_plasma*sum(Bnormal_total(:,:,ilambda)*dBnormaldzmns*norm_normal_plasma)
       enddo
       print *,"Sensitivity computation complete."
     endif

     call system_clock(toc)
     print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
     print "(a,es10.3,a,es10.3)","   chi2_B:",chi2_B(ilambda),",  chi2_K:",chi2_K(ilambda)
     print "(a,es10.3,a,es10.3,a,es10.3)","   max(B_n):",max_Bnormal(ilambda),",  max(K):",max_K(ilambda),",  rms K:",sqrt(chi2_K(ilambda)/area_coil)
  end do

  print *,"Deallocating dynamic memory in solve."
  if (sensitivity_option > 1) then
    deallocate(rmncterm1)
    deallocate(rmnsterm1)
    deallocate(zmncterm1)
    deallocate(zmnsterm1)
    deallocate(rmncterm2)
    deallocate(rmnsterm2)
    deallocate(zmncterm2)
    deallocate(zmnsterm2)
    deallocate(dKDifferencedrmnc)
    deallocate(dKDifferencedrmns)
    deallocate(dKDifferencedzmnc)
    deallocate(dKDifferencedzmns)
    deallocate(dBnormaldrmnc)
    deallocate(dBnormaldrmns)
    deallocate(dBnormaldzmnc)
    deallocate(dBnormaldzmns)
  endif

  deallocate(matrix)
  deallocate(RHS)
  deallocate(solution)
  deallocate(WORK)
  deallocate(IPIV)
  deallocate(this_current_potential)
  deallocate(KDifference_x)
  deallocate(KDifference_y)
  deallocate(KDifference_z)
  deallocate(this_K2_times_N)

end subroutine solve

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
