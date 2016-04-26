module init_basis_functions_mod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine init_basis_functions(basis_option, mpol, ntor, mnmax, xm, xn, u, v, &
         num_basis_functions, basis_functions, area, norm_normal, should_be_identity)

      use global_variables, only: symmetry_option, nfp, check_orthogonality
      use init_Fourier_modes_mod
      use stel_kinds
      use stel_constants

      implicit none

      integer, intent(in) :: basis_option, mpol, ntor
      integer, intent(out) :: mnmax, num_basis_functions
      integer, dimension(:), allocatable :: xm, xn
      real(dp), dimension(:,:), allocatable :: basis_functions
      real(dp), intent(in) :: area
      real(dp), dimension(:,:), allocatable, intent(in) :: norm_normal
      real(dp), dimension(:), allocatable, intent(in) :: u, v
      real(dp), dimension(:,:), allocatable :: should_be_identity

      integer :: nu, nv, iu, iv, index, imn, imn2, iflag
      integer :: tic, toc, countrate
      integer :: whichSymmetry, minSymmetry, maxSymmetry, offset
      real(dp), dimension(:,:), allocatable :: tempMatrix, basis_to_Fourier, basis_functions_transpose
      real(dp) :: du, dv, factor, temp
      real(dp), dimension(:), allocatable :: basis_function_times_N_weight, norm_normal_1D
      logical :: any_problems

      ! Variables needed by LAPACK:
      character :: UPLO
      integer :: INFO
      integer, dimension(:), allocatable :: IPIV

      nu = size(u)
      nv = size(v)
      du = u(2)-u(1)
      dv = v(2)-v(1)

      select case (basis_option)
      case (1)
         print *,"  Weight w = 1 / (nfp * |N|), so the basis functions are sqrt(2) times sin/cos."
      case (2)
         print *,"  Weight w = 1 / area. Basis functions are sin/cos times sqrt(2*area/[nfp*|N|])."
      case (3)
         print *,"  Weight w = 1 / area. Basis functions are linear combinations of sin/cos determined by Cholesky decomposition."
      case default
         print *,"  Error! Invalid setting for basis_option: ",basis_option
         stop
      end select

      ! Initialize Fourier arrays
      call init_Fourier_modes(mpol, ntor, mnmax, xm, xn)

      select case (symmetry_option)
      case (1,2)
         num_basis_functions = mnmax
      case (3)
         num_basis_functions = mnmax * 2
      case default
         print *,"Error! Invalid setting for symmetry_option:",symmetry_option
         stop
      end select

      allocate(basis_functions(nu*nv, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      select case (symmetry_option)
      case (1)
         minSymmetry = 1
         maxSymmetry = 1
      case (2)
         minSymmetry = 2
         maxSymmetry = 2
      case (3)
         minSymmetry = 1
         maxSymmetry = 2
      end select

      call system_clock(tic,countrate)

      ! This loop could be made faster
      ! by using the sum-angle trig identities and pretabulating the trig functions.
      ! But these loops are not the rate-limiting step, so I'll use the more transparent direct method here.
      do whichSymmetry = minSymmetry, maxSymmetry
         
         if (whichSymmetry==2 .and. symmetry_option==3) then
            offset = mnmax
         else
            offset = 0
         end if
         
         do iu = 1, nu
            do iv = 1, nv
               index = (iv-1)*nu + iu
               ! This next "select case" block is where the normalization of the basis functions is set
               select case (basis_option)
               case (1,3)
                  factor = sqrt2
               case (2)
                  factor = sqrt(2.0_dp*area/(nfp*norm_normal(iu,iv)))
               case default
                  stop "Invalid basis_option!"
               end select
               do imn = 1, mnmax
                  if (whichSymmetry==1) then
                     basis_functions(index, imn) = factor * sin(twopi*(xm(imn)*u(iu)+xn(imn)*v(iv)))
                  else
                     basis_functions(index, imn+offset) = factor * cos(twopi*(xm(imn)*u(iu)+xn(imn)*v(iv)))
                  end if
               end do
            end do
         end do
      end do

      call system_clock(toc)
      print *,"  main loop:",real(toc-tic)/countrate,"sec."

      ! If needed for the Cholesky approach, compute the transformation between the Fourier functions and the 'real' basis functions.
      if (basis_option == 3) then
         call system_clock(tic)
         allocate(basis_to_Fourier(num_basis_functions, num_basis_functions),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'
         allocate(tempMatrix(num_basis_functions, nu*nv),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'

         index = 0
         do iv = 1,nv
            do iu = 1,nu
               index = index + 1
               ! The 2nd dimension of norm_normal is nvl rather than nv, but we can just ignore all periods after the 1st.
               tempMatrix(:,index) = basis_functions(index,:) * norm_normal(iu,iv)
            end do
         end do
         basis_to_Fourier = matmul(tempMatrix,basis_functions)*du*dv*nfp/area
         ! The 'basis_to_Fourier' array is not yet what the name suggests, since we will shortly do an in-place Cholesky decomposition
         ! Rather, at this point it is the matrix $\int d^2a f_i f_j$ where $f_i$ are the sqrt(2)*sin(..) functions.

         call system_clock(toc)
         print *,"  Assemble C:",real(toc-tic)/countrate,"sec."
         call system_clock(tic)

         ! Compute Cholesky factorization:
         UPLO = 'L'
         call DPOTRF(UPLO, num_basis_functions, basis_to_Fourier, num_basis_functions, INFO)
         if (INFO < 0) then
            print *,"Error in Cholesky decomposition DPOTRF. The i-th argument had an illegal value. INFO=",INFO
         elseif (INFO > 0) then
            print *,"Error in Cholesky decomposition DPOTRF. The leading minor of order i is not positive definite, and the factorization could not be completed. INFO=",INFO
         end if
         
         call system_clock(toc)
         print *,"  Cholesky decomp:",real(toc-tic)/countrate,"sec."
         deallocate(tempMatrix)

         ! LAPACK's DPOTRF leaves the upper-triangular part nonzero, so clean it up now.
         do imn = 2,num_basis_functions
            basis_to_Fourier(1:imn-1, imn) = 0
         end do

         call system_clock(tic)
         allocate(IPIV(num_basis_functions),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'
         allocate(tempMatrix(num_basis_functions, num_basis_functions),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'
         allocate(basis_functions_transpose(num_basis_functions, nu*nv),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'

         tempMatrix = basis_to_Fourier
         basis_functions_transpose = transpose(basis_functions)
         ! DGESV will overwrite tempMatrix with the L & U factors,
         ! and overwrite 'basis_functions_transpose' with the solution.
!!$*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         ! We probably could call a faster subroutine since basis_to_Fourier is already lower-triangular.
         call DGESV(num_basis_functions, nu*nv, tempMatrix, num_basis_functions, &
              IPIV, basis_functions_transpose, num_basis_functions, INFO)

         basis_functions = transpose(basis_functions_transpose)

         if (INFO < 0) then
            print *,"Error in DGESV. The i-th argument had an illegal value. INFO=",INFO
         elseif (INFO > 0) then
            print *,"Error in DGESV. Matrix is singular. INFO=",INFO
         end if
         deallocate(IPIV,tempMatrix, basis_functions_transpose)
         call system_clock(toc)
         print *,"  Convert basis:",real(toc-tic)/countrate,"sec."
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      ! Done assembling basis functions.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


      if (check_orthogonality) then
         call system_clock(tic)
         allocate(should_be_identity(num_basis_functions, num_basis_functions),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'
         allocate(basis_function_times_N_weight(nu*nv),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'
         allocate(norm_normal_1D(nu*nv),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error!'

         norm_normal_1D = reshape(norm_normal(:,1:nv), (/ nu*nv /))

         any_problems = .false.
         factor = du * dv * nfp
         ! This loop is not at all optimized, but it's only for testing, so probably no need.
         do imn = 1,num_basis_functions
            select case (basis_option)
            case (1)
               ! w = 1 / (nfp * |N|)
               ! Factors of |N| cancel
               basis_function_times_N_weight = basis_functions(:,imn) / nfp
            case (2,3)
               ! w = 1 / area
               basis_function_times_N_weight = basis_functions(:,imn) * norm_normal_1D / area
            case default
               stop "Invalid basis_option!"
            end select

            do imn2 = 1,num_basis_functions
               temp = dot_product(basis_function_times_N_weight, basis_functions(:,imn2)) * factor
               should_be_identity(imn,imn2) = temp
               if (imn==imn2) then
                  ! Diagonal elements should be 1
                  if (abs(temp-1) > 1e-12) then
                     any_problems = .true.
                  end if
               else
                  ! Off-diagonal elements should be 0
                  if (abs(temp) > 1e-12) then
                     any_problems = .true.
                  end if
               end if
            end do
         end do

         if (any_problems) then
            print *,"  WARNING!!!  Orthogonality test failed!!"
         else
            print *,"  Orthogonality test passed."
         end if

         deallocate(basis_function_times_N_weight, norm_normal_1D)
         call system_clock(toc)
         print *,"  Check orthogonality:",real(toc-tic)/countrate,"sec."
      end if


  end subroutine init_basis_functions
  
end module init_basis_functions_mod
  

! Documentation for LAPACK subroutine for Cholesky decomposition:

!!$*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          UPLO
!!$*       INTEGER            INFO, LDA, N
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       DOUBLE PRECISION   A( LDA, * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DPOTRF computes the Cholesky factorization of a real symmetric
!!$*> positive definite matrix A.
!!$*>
!!$*> The factorization has the form
!!$*>    A = U**T * U,  if UPLO = 'U', or
!!$*>    A = L  * L**T,  if UPLO = 'L',
!!$*> where U is an upper triangular matrix and L is lower triangular.
!!$*>
!!$*> This is the block version of the algorithm, calling Level 3 BLAS.
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
!!$*>          The order of the matrix A.  N >= 0.
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
!!$*>          On exit, if INFO = 0, the factor U or L from the Cholesky
!!$*>          factorization A = U**T*U or A = L*L**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0:  if INFO = i, the leading minor of order i is not
!!$*>                positive definite, and the factorization could not be
!!$*>                completed.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Documentation for LAPACK's subroutine for solving a linear system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!$
!!$
!!$*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       INTEGER            INFO, LDA, LDB, N, NRHS
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IPIV( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGESV computes the solution to a real system of linear equations
!!$*>    A * X = B,
!!$*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!!$*>
!!$*> The LU decomposition with partial pivoting and row interchanges is
!!$*> used to factor A as
!!$*>    A = P * L * U,
!!$*> where P is a permutation matrix, L is unit lower triangular, and U is
!!$*> upper triangular.  The factored form of A is then used to solve the
!!$*> system of equations A * X = B.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
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
!!$*>          On entry, the N-by-N coefficient matrix A.
!!$*>          On exit, the factors L and U from the factorization
!!$*>          A = P*L*U; the unit diagonal elements of L are not stored.
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
!!$*>          The pivot indices that define the permutation matrix P;
!!$*>          row i of the matrix was interchanged with row IPIV(i).
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*>          On entry, the N-by-NRHS matrix of right hand side matrix B.
!!$*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>          The leading dimension of the array B.  LDB >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!!$*>                has been completed, but the factor U is exactly
!!$*>                singular, so the solution could not be computed.
