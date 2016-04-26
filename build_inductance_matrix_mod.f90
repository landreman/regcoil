module build_inductance_matrix_mod

  implicit none

contains

  subroutine build_inductance_matrix(inductance, u_1, v_1, r_1, normal_1, basis_functions_1,  &
       u_2, v_2, r_2, normal_2, basis_functions_2)

    use global_variables, only: nfp
    use stel_constants
    use stel_kinds
    use omp_lib

    implicit none

    real(dp), dimension(:,:), allocatable, intent(out) :: inductance
    real(dp), dimension(:,:,:), allocatable, intent(in) :: r_1, r_2, normal_1, normal_2
    real(dp), dimension(:), allocatable, intent(in) :: u_1, u_2, v_1, v_2
    real(dp), dimension(:,:), allocatable, intent(in) :: basis_functions_1, basis_functions_2

    integer :: nu_1, nu_2, nv_1, nv_2
    real(dp) :: du_1, du_2, dv_1, dv_2
    integer :: num_basis_functions_1, num_basis_functions_2

    integer :: l_2, iu_1, iv_1, iu_2, iv_2, ivl_2, iflag
    real(dp) :: x, y, z, dx, dy, dz, dr2, dr32
    integer :: imn_1, imn_2, index_1, index_2
    real(dp), dimension(:,:), allocatable :: inductance_xbasis
    integer :: tic, toc, countrate

    ! Variables needed by BLAS DGEMM:
    character :: TRANSA, TRANSB
    integer :: M, N, K, LDA, LDB, LDC
    real(dp) :: ALPHA=1, BETA=0
    real(dp), dimension(:,:), allocatable :: tempMatrix

    nu_1 = size(u_1)
    nu_2 = size(u_2)
    nv_1 = size(v_1)
    nv_2 = size(v_2)
    du_1 = u_1(2)-u_1(1)
    du_2 = u_2(2)-u_2(1)
    dv_1 = v_1(2)-v_1(1)
    dv_2 = v_2(2)-v_2(1)
    num_basis_functions_1 = size(basis_functions_1,2)
    num_basis_functions_2 = size(basis_functions_2,2)

    allocate(inductance(num_basis_functions_1, num_basis_functions_2),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(inductance_xbasis(nu_1*nv_1, nu_2*nv_2),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now compute S from Boozer's eq (39)-(40).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    inductance_xbasis = 0
    call system_clock(tic,countrate)

    !$OMP PARALLEL

    !$OMP MASTER
    print *,"  Number of OpenMP threads:",omp_get_num_threads()
    !$OMP END MASTER

    !$OMP DO PRIVATE(index_1,index_2,x,y,z,ivl_2,dx,dy,dz,dr2,dr32)
    do iv_2 = 1, nv_2
       do iu_2 = 1, nu_2
          index_2 = (iv_2-1)*nu_2 + iu_2
          do iv_1 = 1, nv_1
             do iu_1 = 1, nu_1
                index_1 = (iv_1-1)*nu_1 + iu_1
                x = r_1(1,iu_1,iv_1)
                y = r_1(2,iu_1,iv_1)
                z = r_1(3,iu_1,iv_1)
                do l_2 = 0, (nfp-1)
                   ivl_2 = iv_2 + l_2*nv_2
                   dx = x - r_2(1,iu_2,ivl_2)
                   dy = y - r_2(2,iu_2,ivl_2)
                   dz = z - r_2(3,iu_2,ivl_2)

                   dr2 = dx*dx + dy*dy + dz*dz
                   dr32 = dr2*sqrt(dr2)

                   inductance_xbasis(index_1,index_2) = inductance_xbasis(index_1,index_2) + &
                        (normal_1(1,iu_1,iv_1)*normal_2(1,iu_2,ivl_2) &
                        +normal_1(2,iu_1,iv_1)*normal_2(2,iu_2,ivl_2) &
                        +normal_1(3,iu_1,iv_1)*normal_2(3,iu_2,ivl_2) &
                        - (3/dr2) * &
                        (normal_1(1,iu_1,iv_1)*dx + normal_1(2,iu_1,iv_1)*dy + normal_1(3,iu_1,iv_1)*dz) * &
                        (normal_2(1,iu_2,ivl_2)*dx &
                        +normal_2(2,iu_2,ivl_2)*dy &
                        +normal_2(3,iu_2,ivl_2)*dz)) / dr32

!!$                   inductance_xbasis(index_1,index_2) = inductance_xbasis(index_1,index_2) + &
!!$                        (normal_1(1,iu_1,iv_1)*dx + normal_1(2,iu_1,iv_1)*dy + normal_1(3,iu_1,iv_1)*dz) / dr32
                end do
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call system_clock(toc)
    print *,"  inductance_xbasis:",real(toc-tic)/countrate,"sec."

    call system_clock(tic)

    ! Next, convert from the integrand as a function of (u,v) to the modal matrix
    ! by sandwiching with the basis_functions matrices.

    ! For some reason, the BLAS matrix-matrix multiplication function DGEMM sometimes causes the
    ! program to crash on Edison unless you are careful to use the Intel MKL instead of Cray LibSci.
    ! If you like, you can use the following method which is slower but more reliable:
!    inductance = matmul(transpose(basis_functions_1), matmul(inductance_xbasis, basis_functions_2))

    !*******************************************************
    ! Call BLAS3 subroutine DGEMM for matrix multiplications:
    !*******************************************************

    ! Here we carry out tempMatrix = inductance_xbasis * basis_functions_2
    ! A = inductance_xbasis
    ! B = basis_functions_2
    ! C = tempMatrix
    M = nu_1*nv_1 ! # rows of A
    N = num_basis_functions_2 ! # cols of B
    K = nu_2*nv_2 ! Common dimension of A and B
    LDA = M
    LDB = K
    LDC = M
    TRANSA = 'N' ! No transposes
    TRANSB = 'N'
    allocate(tempMatrix(M,N),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    tempMatrix = 0
    ALPHA=1
    BETA=0
    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,inductance_xbasis,LDA,basis_functions_2,LDB,BETA,tempMatrix,LDC)

    ! Here we carry out inductance = (basis_functions_1 ^ T) * tempMatrix
    ! A = basis_functions_1
    ! B = tempMatrix
    ! C = inductance
    M = num_basis_functions_1 ! # rows of A^T
    N = num_basis_functions_2 ! # cols of B
    K = nu_1*nv_1 ! Common dimension of A^T and B
    LDA = K ! Would be M if not taking the transpose.
    LDB = K
    LDC = M
    TRANSA = 'T' ! DO take a transpose!
    TRANSB = 'N'
    ALPHA=1
    BETA=0
    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,basis_functions_1,LDA,tempMatrix,LDB,BETA,inductance,LDC)

    deallocate(tempMatrix)
    
    call system_clock(toc)
    print *,"  matmul:",real(toc-tic)/countrate,"sec."

    ! Multiply by some overall constants:
    inductance = inductance * (nfp * du_1 * dv_1 * du_2 * dv_2 * mu0 / (4*pi))

    deallocate(inductance_xbasis)

  end subroutine build_inductance_matrix


!**********************************************************************************************
!**********************************************************************************************

  subroutine build_NESTOR_like_matrix(inductance, u_1, v_1, r_1, normal_1, basis_functions_1,  &
       u_2, v_2, r_2, normal_2, basis_functions_2, &
       drdu_1, drdv_1, d2rdu2_1, d2rdudv_1, d2rdv2_1, subtract_singularity)

    use global_variables, only: nfp, Merkel_Kmn, nu_middle, nv_middle, mnmax_outer, xm_outer, xn_outer, symmetry_option
    use stel_constants
    use stel_kinds
    use omp_lib

    implicit none

    real(dp), dimension(:,:), allocatable, intent(out) :: inductance
    real(dp), dimension(:,:,:), allocatable, intent(in) :: r_1, r_2, normal_1, normal_2
    real(dp), dimension(:), allocatable, intent(in) :: u_1, u_2, v_1, v_2
    real(dp), dimension(:,:), allocatable, intent(in) :: basis_functions_1, basis_functions_2
    real(dp), dimension(:,:,:), allocatable, intent(in) :: drdu_1, drdv_1, d2rdu2_1, d2rdudv_1, d2rdv2_1
    logical :: subtract_singularity

    integer :: nu_1, nu_2, nv_1, nv_2
    real(dp) :: du_1, du_2, dv_1, dv_2
    integer :: num_basis_functions_1, num_basis_functions_2

    integer :: l_2, iu_1, iv_1, iu_2, iv_2, ivl_2, iflag, min_l_2
    real(dp) :: x, y, z, dx, dy, dz, dr2, dr32
    integer :: imn_1, imn_2, index_1, index_2
    real(dp), dimension(:,:), allocatable :: inductance_xbasis
    integer :: tic, toc, countrate
    real(dp) :: Merkel_little_a, Merkel_little_b, Merkel_little_c
    real(dp) :: Merkel_big_A, Merkel_big_B, Merkel_big_C
    real(dp) :: temp, Ys
    real(dp), dimension(:,:), allocatable :: tanU, tanV, tan2U, tan2V
    integer :: whichSymmetry, minSymmetry, maxSymmetry, offset
    real(dp) :: epsilon, big_number

    ! Variables needed by BLAS DGEMM:
    character :: TRANSA, TRANSB
    integer :: M, N, K, LDA, LDB, LDC
    real(dp) :: ALPHA=1, BETA=0
    real(dp), dimension(:,:), allocatable :: tempMatrix, sincos_on_middle

    nu_1 = size(u_1)
    nu_2 = size(u_2)
    nv_1 = size(v_1)
    nv_2 = size(v_2)
    du_1 = u_1(2)-u_1(1)
    du_2 = u_2(2)-u_2(1)
    dv_1 = v_1(2)-v_1(1)
    dv_2 = v_2(2)-v_2(1)
    num_basis_functions_1 = size(basis_functions_1,2)
    num_basis_functions_2 = size(basis_functions_2,2)

    allocate(inductance(num_basis_functions_1, num_basis_functions_2),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(inductance_xbasis(nu_1*nv_1, nu_2*nv_2),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    inductance_xbasis = 0
    call system_clock(tic,countrate)

    if (subtract_singularity) then

       allocate(Merkel_Kmn(nu_middle*nv_middle, mnmax_outer),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       call compute_Kmn()
       print *,"Back in build_inductance_matrix_mod."

       allocate(tanU(nu_1, nu_2),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(tanV(nv_1, nv_2),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(tan2U(nu_1, nu_2),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(tan2V(nv_1, nv_2),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(sincos_on_middle(nu_1*nv_1,num_basis_functions_2),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'

       ! Instead of forming a 2D array for tanU and tanV, we could just form 1D arrays.
       ! But the 2D method is not too slow, and it is less prone to indexing errors.
       ! If the argument of tan() is close to (within epsilon of) +/- pi/2, where tan() is singular, replace
       ! tan() with big_number. 
       epsilon = 1e-12
       big_number = 1e12
       do iu_1 = 1,nu_1
          do iu_2 = 1,nu_2
             if (abs(abs(u_2(iu_2) - u_1(iu_1))-0.5) < epsilon) then
                temp = big_number
             else
                temp = tan(pi*(u_2(iu_2) - u_1(iu_1)))
             end if
             tanU(iu_1,iu_2) = temp
             tan2U(iu_1,iu_2) = temp*temp
          end do
       end do
       do iv_1 = 1,nv_1
          do iv_2 = 1,nv_2
             if (abs(abs(v_2(iv_2) - v_1(iv_1))-0.5) < epsilon) then
                temp = big_number
             else
                temp = tan(pi*(v_2(iv_2) - v_1(iv_1)))
             end if
             tanV(iv_1,iv_2) = temp
             tan2V(iv_1,iv_2) = temp*temp
          end do
       end do

       ! We need sin/cos functions evaluated on the middle (1) (u,v) grid
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

       ! This loop could be made faster
       ! by using the sum-angle trig identities and pretabulating the trig functions.
       ! But these loops are not the rate-limiting step, so I'll use the more transparent direct method here.
       do whichSymmetry = minSymmetry, maxSymmetry
         
          if (whichSymmetry==2 .and. symmetry_option==3) then
             offset = mnmax_outer
          else
             offset = 0
          end if
         
          do iu_1 = 1, nu_1
             do iv_1 = 1, nv_1
                index_1 = (iv_1-1)*nu_1 + iu_1
                do imn_2 = 1, mnmax_outer
                   if (whichSymmetry==1) then
                      sincos_on_middle(index_1, imn_2)        = sqrt2 * sin(twopi*(xm_outer(imn_2)*u_1(iu_1)+xn_outer(imn_2)*v_1(iv_1)))
                   else
                      sincos_on_middle(index_1, imn_2+offset) = sqrt2 * cos(twopi*(xm_outer(imn_2)*u_1(iu_1)+xn_outer(imn_2)*v_1(iv_1)))
                   end if
                end do
             end do
          end do
       end do

       call system_clock(toc)
       print *,"  Build tan and sincos matrices:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Now compute (r - r') dot N / |r - r'|^(3/2)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !$OMP PARALLEL

       !$OMP MASTER
       print *,"  Number of OpenMP threads:",omp_get_num_threads()
       !$OMP END MASTER

       !$OMP DO PRIVATE(index_1,index_2,x,y,z,ivl_2,dx,dy,dz,dr2,dr32,Merkel_little_a,Merkel_little_b,Merkel_little_c,Merkel_big_A,Merkel_big_B,Merkel_big_C,temp,Ys,min_l_2)
       do iv_1 = 1, nv_1
          do iu_1 = 1, nu_1
             index_1 = (iv_1-1)*nu_1 + iu_1
             x = r_1(1,iu_1,iv_1)
             y = r_1(2,iu_1,iv_1)
             z = r_1(3,iu_1,iv_1)

             Merkel_little_a = drdu_1(1,iu_1,iv_1)*drdu_1(1,iu_1,iv_1) + drdu_1(2,iu_1,iv_1)*drdu_1(2,iu_1,iv_1) + drdu_1(3,iu_1,iv_1)*drdu_1(3,iu_1,iv_1)
             Merkel_little_b = drdu_1(1,iu_1,iv_1)*drdv_1(1,iu_1,iv_1) + drdu_1(2,iu_1,iv_1)*drdv_1(2,iu_1,iv_1) + drdu_1(3,iu_1,iv_1)*drdv_1(3,iu_1,iv_1)
             Merkel_little_c = drdv_1(1,iu_1,iv_1)*drdv_1(1,iu_1,iv_1) + drdv_1(2,iu_1,iv_1)*drdv_1(2,iu_1,iv_1) + drdv_1(3,iu_1,iv_1)*drdv_1(3,iu_1,iv_1)

             Merkel_big_A = 0.5*( d2rdu2_1(1,iu_1,iv_1)*normal_1(1,iu_1,iv_1) &
                                + d2rdu2_1(2,iu_1,iv_1)*normal_1(2,iu_1,iv_1) &
                                + d2rdu2_1(3,iu_1,iv_1)*normal_1(3,iu_1,iv_1) )

             Merkel_big_B = 0.5*( d2rdudv_1(1,iu_1,iv_1)*normal_1(1,iu_1,iv_1) &
                                + d2rdudv_1(2,iu_1,iv_1)*normal_1(2,iu_1,iv_1) &
                                + d2rdudv_1(3,iu_1,iv_1)*normal_1(3,iu_1,iv_1) )

             Merkel_big_C = 0.5*( d2rdv2_1(1,iu_1,iv_1)*normal_1(1,iu_1,iv_1) &
                                + d2rdv2_1(2,iu_1,iv_1)*normal_1(2,iu_1,iv_1) &
                                + d2rdv2_1(3,iu_1,iv_1)*normal_1(3,iu_1,iv_1) )

             do iv_2 = 1, nv_2
                do iu_2 = 1, nu_2
                   index_2 = (iv_2-1)*nu_2 + iu_2
                   if ((iv_1==iv_2) .and. (iu_1==iu_2)) then
                      ! Avoid divide-by-0 singularity
                      min_l_2 = 1
                   else
                      min_l_2 = 0
                      temp = Merkel_little_a*tan2U(iu_1,iu_2) + 2*Merkel_little_b*tanU(iu_1,iu_2)*tanV(iv_1,iv_2) + Merkel_little_c*tan2V(iv_1,iv_2)
                      Ys = -pi*(Merkel_big_A*tan2U(iu_1,iu_2) + 2*Merkel_big_B*tanU(iu_1,iu_2)*tanV(iv_1,iv_2) + Merkel_big_C*tan2V(iv_1,iv_2)) &
                           / (temp*sqrt(temp))

                      inductance_xbasis(index_1,index_2) = inductance_xbasis(index_1,index_2) - Ys
                   end if
                   do l_2 = min_l_2, (nfp-1)
                      ivl_2 = iv_2 + l_2*nv_2
                      dx = x - r_2(1,iu_2,ivl_2)
                      dy = y - r_2(2,iu_2,ivl_2)
                      dz = z - r_2(3,iu_2,ivl_2)
                      
                      dr2 = dx*dx + dy*dy + dz*dz
                      dr32 = dr2*sqrt(dr2)
                      
                      inductance_xbasis(index_1,index_2) = inductance_xbasis(index_1,index_2) + &
                           (dx*normal_1(1,iu_1,iv_1) + dy*normal_1(2,iu_1,iv_1) + dz*normal_1(3,iu_1,iv_1)) / dr32
                      
                   end do
                end do
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

       deallocate(tanU,tanV,tan2U,tan2V)

       call system_clock(toc)
       print *,"  inductance_xbasis:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

       !*******************************************************
       ! Call BLAS3 subroutine DGEMM for matrix multiplications:
       !*******************************************************

       ! Here we carry out tempMatrix = inductance_xbasis * basis_functions_2
       ! A = inductance_xbasis
       ! B = basis_functions_2
       ! C = tempMatrix
       M = nu_1*nv_1 ! # rows of A
       N = num_basis_functions_2 ! # cols of B
       K = nu_2*nv_2 ! Common dimension of A and B
       LDA = M
       LDB = K
       LDC = M
       TRANSA = 'N' ! No transposes
       TRANSB = 'N'
       allocate(tempMatrix(M,N),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       tempMatrix = 0
       call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,inductance_xbasis,LDA,basis_functions_2,LDB,BETA,tempMatrix,LDC)

       tempMatrix = tempMatrix * (du_2*dv_2)

       call system_clock(toc)
       print *,"  1st matmul:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

       !*******************************************************
       ! Now add back the singularity, using Merkel's analytic Kmn integrals.
       !*******************************************************

!!$       do imn_2 = 1,mnmax_outer
!!$          do index_1 = 1,(nu_1*nv_1)
!!$             tempMatrix(index_1,imn_2) = tempMatrix(index_1,imn_2) &
!!$                  - Merkel_Kmn(index_1,imn_2)*sincos_on_middle(index_1,imn_2)
!!$          end do
!!$       end do
!!$       if (symmetry_option==3) then
!!$          do imn_2 = 1,mnmax_outer
!!$             do index_1 = 1,(nu_1*nv_1)
!!$                tempMatrix(index_1,imn_2+mnmax_outer) = tempMatrix(index_1,imn_2+mnmax_outer) &
!!$                     - Merkel_Kmn(index_1,imn_2)*sincos_on_middle(index_1,imn_2+mnmax_outer)
!!$             end do
!!$          end do
!!$       end if

       do imn_2 = 1,mnmax_outer
          index_1 = 0
          do iv_1 = 1,nv_1
             do iu_1 = 1,nu_1
                index_1 = index_1+1
                tempMatrix(index_1,imn_2) = tempMatrix(index_1,imn_2) &
                     - (Merkel_Kmn(index_1,imn_2) + twopi) * sincos_on_middle(index_1,imn_2)
!                tempMatrix(index_1,imn_2) = twopi*norm_normal_middle(iu_1,iv_1) * sincos_on_middle(index_1,imn_2)
             end do
          end do
       end do
       if (symmetry_option==3) then
          do imn_2 = 1,mnmax_outer
             index_1 = 0
             do iv_1 = 1,nv_1
                do iu_1 = 1,nu_1
                   index_1 = index_1 + 1
                   tempMatrix(index_1,imn_2+mnmax_outer) = tempMatrix(index_1,imn_2+mnmax_outer) &
                        - (Merkel_Kmn(index_1,imn_2) + twopi) * sincos_on_middle(index_1,imn_2+mnmax_outer)
!                   tempMatrix(index_1,imn_2+mnmax_outer) = twopi*norm_normal_middle(iu_1,iv_1) * sincos_on_middle(index_1,imn_2+mnmax_outer)
                end do
             end do
          end do
       end if

       deallocate(sincos_on_middle)

       call system_clock(toc)
       print *,"  Add Kmn:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

    else

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Now compute (r - r') dot N / |r - r'|^(3/2)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !$OMP PARALLEL

       !$OMP MASTER
       print *,"  Number of OpenMP threads:",omp_get_num_threads()
       !$OMP END MASTER

       !$OMP DO PRIVATE(index_1,index_2,x,y,z,ivl_2,dx,dy,dz,dr2,dr32)
       do iv_1 = 1, nv_1
          do iu_1 = 1, nu_1
             index_1 = (iv_1-1)*nu_1 + iu_1
             x = r_1(1,iu_1,iv_1)
             y = r_1(2,iu_1,iv_1)
             z = r_1(3,iu_1,iv_1)
             do iv_2 = 1, nv_2
                do iu_2 = 1, nu_2
                   index_2 = (iv_2-1)*nu_2 + iu_2
                   do l_2 = 0, (nfp-1)
                      ivl_2 = iv_2 + l_2*nv_2
                      dx = x - r_2(1,iu_2,ivl_2)
                      dy = y - r_2(2,iu_2,ivl_2)
                      dz = z - r_2(3,iu_2,ivl_2)
                      
                      dr2 = dx*dx + dy*dy + dz*dz
                      dr32 = dr2*sqrt(dr2)
                      
                      inductance_xbasis(index_1,index_2) = inductance_xbasis(index_1,index_2) + &
                           (dx*normal_1(1,iu_1,iv_1) + dy*normal_1(2,iu_1,iv_1) + dz*normal_1(3,iu_1,iv_1)) / dr32
                      
                   end do
                end do
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

       call system_clock(toc)
       print *,"  inductance_xbasis:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

       !*******************************************************
       ! Call BLAS3 subroutine DGEMM for matrix multiplications:
       !*******************************************************

       ! Here we carry out tempMatrix = inductance_xbasis * basis_functions_2
       ! A = inductance_xbasis
       ! B = basis_functions_2
       ! C = tempMatrix
       M = nu_1*nv_1 ! # rows of A
       N = num_basis_functions_2 ! # cols of B
       K = nu_2*nv_2 ! Common dimension of A and B
       LDA = M
       LDB = K
       LDC = M
       TRANSA = 'N' ! No transposes
       TRANSB = 'N'
       allocate(tempMatrix(M,N),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       tempMatrix = 0
       call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,inductance_xbasis,LDA,basis_functions_2,LDB,BETA,tempMatrix,LDC)

       tempMatrix = tempMatrix * (du_2*dv_2)

       call system_clock(toc)
       print *,"  1st matmul:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

    end if


    ! Next, convert from the integrand as a function of (u,v) to the modal matrix
    ! by sandwiching with the basis_functions matrices.

    ! For some reason, the BLAS matrix-matrix multiplication function DGEMM sometimes causes the
    ! program to crash on Edison unless you are careful to use the Intel MKL instead of Cray LibSci.
    ! If you like, you can use the following method which is slower but more reliable:
!    inductance = matmul(transpose(basis_functions_1), matmul(inductance_xbasis, basis_functions_2))

    !*******************************************************
    ! Call BLAS3 subroutine DGEMM for matrix multiplications:
    !*******************************************************
    ! Here we carry out inductance = (basis_functions_1 ^ T) * tempMatrix
    ! A = basis_functions_1
    ! B = tempMatrix
    ! C = inductance
    M = num_basis_functions_1 ! # rows of A^T
    N = num_basis_functions_2 ! # cols of B
    K = nu_1*nv_1 ! Common dimension of A^T and B
    LDA = K ! Would be M if not taking the transpose.
    LDB = K
    LDC = M
    TRANSA = 'T' ! DO take a transpose!
    TRANSB = 'N'
    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,basis_functions_1,LDA,tempMatrix,LDB,BETA,inductance,LDC)

    deallocate(tempMatrix)
    
    call system_clock(toc)
    print *,"  2nd matmul:",real(toc-tic)/countrate,"sec."

    ! Multiply by some overall constants:
    inductance = inductance * (-nfp * du_1 * dv_1)

    deallocate(inductance_xbasis)

  end subroutine build_NESTOR_like_matrix

end module build_inductance_matrix_mod


! Documentation of BLAS3 DGEMM subroutine for matrix-matrix multiplication:

!!$*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       DOUBLE PRECISION ALPHA,BETA
!!$*       INTEGER K,LDA,LDB,LDC,M,N
!!$*       CHARACTER TRANSA,TRANSB
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGEMM  performs one of the matrix-matrix operations
!!$*>
!!$*>    C := alpha*op( A )*op( B ) + beta*C,
!!$*>
!!$*> where  op( X ) is one of
!!$*>
!!$*>    op( X ) = X   or   op( X ) = X**T,
!!$*>
!!$*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!!$*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] TRANSA
!!$*> \verbatim
!!$*>          TRANSA is CHARACTER*1
!!$*>           On entry, TRANSA specifies the form of op( A ) to be used in
!!$*>           the matrix multiplication as follows:
!!$*>
!!$*>              TRANSA = 'N' or 'n',  op( A ) = A.
!!$*>
!!$*>              TRANSA = 'T' or 't',  op( A ) = A**T.
!!$*>
!!$*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] TRANSB
!!$*> \verbatim
!!$*>          TRANSB is CHARACTER*1
!!$*>           On entry, TRANSB specifies the form of op( B ) to be used in
!!$*>           the matrix multiplication as follows:
!!$*>
!!$*>              TRANSB = 'N' or 'n',  op( B ) = B.
!!$*>
!!$*>              TRANSB = 'T' or 't',  op( B ) = B**T.
!!$*>
!!$*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] M
!!$*> \verbatim
!!$*>          M is INTEGER
!!$*>           On entry,  M  specifies  the number  of rows  of the  matrix
!!$*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>           On entry,  N  specifies the number  of columns of the matrix
!!$*>           op( B ) and the number of columns of the matrix C. N must be
!!$*>           at least zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] K
!!$*> \verbatim
!!$*>          K is INTEGER
!!$*>           On entry,  K  specifies  the number of columns of the matrix
!!$*>           op( A ) and the number of rows of the matrix op( B ). K must
!!$*>           be at least  zero.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] ALPHA
!!$*> \verbatim
!!$*>          ALPHA is DOUBLE PRECISION.
!!$*>           On entry, ALPHA specifies the scalar alpha.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!!$*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!!$*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!!$*>           part of the array  A  must contain the matrix  A,  otherwise
!!$*>           the leading  k by m  part of the array  A  must contain  the
!!$*>           matrix A.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>           On entry, LDA specifies the first dimension of A as declared
!!$*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!!$*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!!$*>           least  max( 1, k ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!!$*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!!$*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!!$*>           part of the array  B  must contain the matrix  B,  otherwise
!!$*>           the leading  n by k  part of the array  B  must contain  the
!!$*>           matrix B.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>           On entry, LDB specifies the first dimension of B as declared
!!$*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!!$*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!!$*>           least  max( 1, n ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] BETA
!!$*> \verbatim
!!$*>          BETA is DOUBLE PRECISION.
!!$*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!!$*>           supplied as zero then C need not be set on input.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] C
!!$*> \verbatim
!!$*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!!$*>           Before entry, the leading  m by n  part of the array  C must
!!$*>           contain the matrix  C,  except when  beta  is zero, in which
!!$*>           case C need not be set on entry.
!!$*>           On exit, the array  C  is overwritten by the  m by n  matrix
!!$*>           ( alpha*op( A )*op( B ) + beta*C ).
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDC
!!$*> \verbatim
!!$*>          LDC is INTEGER
!!$*>           On entry, LDC specifies the first dimension of C as declared
!!$*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!!$*>           max( 1, m ).
!!$*> \endverbatim
