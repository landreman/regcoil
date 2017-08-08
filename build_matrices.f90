subroutine build_matrices()

  use global_variables
  use stel_constants
  use stel_kinds
  use omp_lib
  use init_Fourier_modes_mod
  
  implicit none

  integer :: l_coil, itheta_plasma, izeta_plasma, itheta_coil, izeta_coil, izetal_coil,iomega
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, dr52inv
  real(dp) :: cosangle_xn, cosangle_xm, sinangle_xn, sinangle_xm
  real(dp) :: dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
  real(dp) :: dy_norm3, dy_norm1, dx_norm2, dx_norm3, dz_norm1, dz_norm2, this_h
  integer :: index_plasma, index_coil, j, imn
  integer :: tic, toc, countrate, iflag, indexl_coil
  integer :: minSymmetry, maxSymmetry, whichSymmetry, offset
  real(dp) :: angle, sinangle, cosangle
  real(dp), dimension(:,:,:), allocatable :: factor_for_h
  real(dp), dimension(:,:,:), allocatable :: f_xdNdomega_over_N_coil2,f_ydNdomega_over_N_coil2,f_zdNdomega_over_N_coil2
  real(dp), dimension(:), allocatable :: norm_normal_plasma_inv1D, norm_normal_coil_inv1D
  real(dp), dimension(:,:), allocatable :: g_over_N_plasma, f_x_over_N_coil, f_y_over_N_coil, f_z_over_N_coil
  real(dp), dimension(:), allocatable :: dinductancednorm, dinductancedr
  
  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0
  real(dp), dimension(:,:), allocatable :: tempMatrix

  call system_clock(tic,countrate)
  print *,"Initializing basis functions and f"
  
  ! Initialize Fourier arrays
  call init_Fourier_modes(mpol_coil, ntor_coil, mnmax_coil, xm_coil, xn_coil)
  xn_coil = xn_coil * nfp
  
  select case (symmetry_option)
  case (1,2)
     num_basis_functions = mnmax_coil
  case (3)
     num_basis_functions = mnmax_coil * 2
  case default
     print *,"Error! Invalid setting for symmetry_option:",symmetry_option
     stop
  end select
  print *,"num_basis_function = ", num_basis_functions
  
  allocate(basis_functions(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_x(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_y(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_z(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(d_x(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(d_y(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(d_z(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (sensitivity_option > 1) then
    allocate(dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dfydomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dfzdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif
  if (sensitivity_option > 2) then
    allocate(dmatrix_Kdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dmatrix_Bdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(f_xdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(f_ydNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(f_zdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dRHS_Kdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dRHS_Bdomega(nomega_coil,num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif

  d_x = reshape((net_poloidal_current_Amperes * drdtheta_coil(1,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(1,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  d_y = reshape((net_poloidal_current_Amperes * drdtheta_coil(2,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(2,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  d_z = reshape((net_poloidal_current_Amperes * drdtheta_coil(3,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(3,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))

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
        offset = mnmax_coil
     else
        offset = 0
     end if

    !$OMP PARALLEL
    !$OMP MASTER
    print *,"  Number of OpenMP threads:",omp_get_num_threads()
    !$OMP END MASTER
    !$OMP DO PRIVATE(itheta_coil,index_coil,imn,angle,sinangle,cosangle,cosangle_xn,cosangle_xm,sinangle_xn,sinangle_xm)
     do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
           index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
           do imn = 1, mnmax_coil
              angle = xm_coil(imn)*theta_coil(itheta_coil)-xn_coil(imn)*zeta_coil(izeta_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)
              cosangle_xn = cosangle*xn_coil(imn)
              sinangle_xn = sinangle*xn_coil(imn)
              cosangle_xm = cosangle*xm_coil(imn)
              sinangle_xm = sinangle*xm_coil(imn)
              if (whichSymmetry==1) then
                 basis_functions(index_coil, imn) = sinangle
                 f_x(index_coil, imn) = cosangle_xn*drdtheta_coil(1,itheta_coil,izeta_coil) &
                    + cosangle_xm*drdzeta_coil(1,itheta_coil,izeta_coil)
                 f_y(index_coil, imn) = cosangle_xn*drdtheta_coil(2,itheta_coil,izeta_coil) &
                    + cosangle_xm*drdzeta_coil(2,itheta_coil,izeta_coil)
                 f_z(index_coil, imn) = cosangle_xn*drdtheta_coil(3,itheta_coil,izeta_coil) &
                    + cosangle_xm*drdzeta_coil(3,itheta_coil,izeta_coil)
                 if (sensitivity_option > 1) then
                    dfxdomega(:, index_coil,imn) = &
                      cosangle_xn*domegadxdtheta(:,itheta_coil,izeta_coil) &
                      + cosangle_xm*domegadxdzeta(:,itheta_coil,izeta_coil)
                    dfydomega(:, index_coil,imn) = &
                      cosangle_xn*domegadydtheta(:,itheta_coil,izeta_coil) &
                      + cosangle_xm*domegadydzeta(:,itheta_coil,izeta_coil)
                    dfzdomega(:, index_coil,imn) = &
                      cosangle_xn*domegadzdtheta(:,itheta_coil,izeta_coil) &
                      + cosangle_xm*domegadzdzeta(:,itheta_coil,izeta_coil)
                 endif
              else
                 basis_functions(index_coil, imn+offset) = cosangle
                 f_x(index_coil, imn+offset) = -sinangle_xn*drdtheta_coil(1,itheta_coil,izeta_coil) &
                   -sinangle_xm*drdzeta_coil(1,itheta_coil,izeta_coil)
                 f_y(index_coil, imn+offset) = -sinangle_xn*drdtheta_coil(2,itheta_coil,izeta_coil) &
                   -sinangle_xm*drdzeta_coil(2,itheta_coil,izeta_coil)
                 f_z(index_coil, imn+offset) = -sinangle_xn*drdtheta_coil(3,itheta_coil,izeta_coil) &
                   -sinangle_xm*drdzeta_coil(3,itheta_coil,izeta_coil)
                if (sensitivity_option > 1) then
                  dfxdomega(:, index_coil,imn+offset) = &
                    -sinangle_xn*domegadxdtheta(:,itheta_coil,izeta_coil) &
                    -sinangle_xm*domegadxdzeta(:,itheta_coil,izeta_coil)
                  dfydomega(:, index_coil,imn+offset) = &
                    -sinangle_xn*domegadydtheta(:,itheta_coil,izeta_coil) &
                    -sinangle_xm*domegadydzeta(:,itheta_coil,izeta_coil)
                  dfzdomega(:, index_coil,imn+offset) = &
                    -sinangle_xn*domegadzdtheta(:,itheta_coil,izeta_coil) &
                    -sinangle_xm*domegadzdzeta(:,itheta_coil,izeta_coil)
                endif
              endif
           end do
        end do
     end do
    !$OMP END DO
    !$OMP END PARALLEL
  end do
  
  call system_clock(toc)
  print *,"Done. Took",real(toc-tic)/countrate,"sec."

  allocate(g(ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(inductance(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(h(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(factor_for_h(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_net_coil_currents(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (sensitivity_option > 1) then
    allocate(dgdomega(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancednorm(3),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancedr(3),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancedomega(nomega_coil, ntheta_plasma*nzeta_plasma, &
      ntheta_coil*nzeta_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dhdomega(nomega_coil,ntheta_plasma*nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now compute g and h
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  inductance = 0
  h=0
  if (sensitivity_option > 1) then
    dinductancedomega = 0
    dhdomega = 0
  endif
  factor_for_h = net_poloidal_current_Amperes * drdtheta_coil - net_toroidal_current_Amperes * drdzeta_coil

  call system_clock(tic,countrate)
  print *,"Building inductance matrix and h."
 !$OMP PARALLEL
 !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
 !$OMP END MASTER

  ! Note: the outermost loop below must be over the plasma variables rather than over the coil variables.
  ! This ensures the multiple threads write to different indices in h() rather than to the same indices in h(),
  ! in which case the h(index+plasma)=h(index_plasma)+... update does not work properly.
 !$OMP DO PRIVATE(index_plasma,index_coil,x,y,z,izetal_coil,dx,dy,dz,dr2inv,dr32inv,indexl_coil,dr52inv,dr_dot_norm_coil,dr_dot_norm_plasma,norm_plasma_dot_norm_coil,dx_norm2,dx_norm3,dy_norm1,dy_norm3,dz_norm1,dz_norm2,this_h)
  do izeta_plasma = 1, nzeta_plasma
     do itheta_plasma = 1, ntheta_plasma
        index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
        x = r_plasma(1,itheta_plasma,izeta_plasma)
        y = r_plasma(2,itheta_plasma,izeta_plasma)
        z = r_plasma(3,itheta_plasma,izeta_plasma)
        do izeta_coil = 1, nzeta_coil
           do itheta_coil = 1, ntheta_coil
              index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
              do l_coil = 0, (nfp-1)
                 izetal_coil = izeta_coil + l_coil*nzeta_coil
                 dx = x - r_coil(1,itheta_coil,izetal_coil)
                 dy = y - r_coil(2,itheta_coil,izetal_coil)
                 dz = z - r_coil(3,itheta_coil,izetal_coil)
                 dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                 dr32inv = dr2inv*sqrt(dr2inv)
                 ! Pre-multiplying these factors for optimization purposes
                 dy_norm3 = dy * normal_plasma(3,itheta_plasma,izeta_plasma)
                 dz_norm1 = dz * normal_plasma(1,itheta_plasma,izeta_plasma)
                 dx_norm2 = dx * normal_plasma(2,itheta_plasma,izeta_plasma)
                 dy_norm1 = dy * normal_plasma(1,itheta_plasma,izeta_plasma)
                 dz_norm2 = dz * normal_plasma(2,itheta_plasma,izeta_plasma)
                 dx_norm3 = dx * normal_plasma(3,itheta_plasma,izeta_plasma)
                 norm_plasma_dot_norm_coil = normal_plasma(1,itheta_plasma,izeta_plasma) &
                   *normal_coil(1,itheta_coil,izetal_coil) + normal_plasma(2,itheta_plasma,izeta_plasma) &
                   *normal_coil(2,itheta_coil,izetal_coil) + normal_plasma(3,itheta_plasma,izeta_plasma) &
                   *normal_coil(3,itheta_coil,izetal_coil)
                 dr_dot_norm_coil = dx*normal_coil(1,itheta_coil,izetal_coil) &
                  + dy*normal_coil(2,itheta_coil,izetal_coil) + dz*normal_coil(3,itheta_coil,izetal_coil)
                 dr_dot_norm_plasma = dx*normal_plasma(1,itheta_plasma,izeta_plasma) &
                  + dy*normal_plasma(2,itheta_plasma,izeta_plasma) + dz*normal_plasma(3, itheta_plasma,izeta_plasma)
                 this_h = (factor_for_h(1,itheta_coil,izetal_coil) * dy_norm3 + &
                   factor_for_h(2,itheta_coil,izetal_coil) * dz_norm1 + &
                   factor_for_h(3,itheta_coil,izetal_coil) * dx_norm2  &
                   - factor_for_h(3,itheta_coil,izetal_coil) * dy_norm1 &
                   - factor_for_h(1,itheta_coil,izetal_coil) * dz_norm2 &
                   - factor_for_h(2,itheta_coil,izetal_coil) * dx_norm3 ) * dr32inv

                 inductance(index_plasma,index_coil) = inductance(index_plasma,index_coil) + &
                      (norm_plasma_dot_norm_coil - (3*dr2inv) * &
                      (dr_dot_norm_plasma * dr_dot_norm_coil)) * dr32inv
                 
                 h(index_plasma) = h(index_plasma) + this_h

                if (sensitivity_option > 1) then
                  indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil
                  dr52inv = dr2inv*dr32inv

                 dinductancedomega(:,index_plasma,index_coil) = dinductancedomega(:,index_plasma,index_coil) &
                    + (normal_plasma(1,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dx)*(dr32inv*mu0/(4*pi))*dnormxdomega(:,index_coil,l_coil+1) &
                    + (normal_plasma(2,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dy)*(dr32inv*mu0/(4*pi))*dnormydomega(:,index_coil,l_coil+1) &
                    + (normal_plasma(3,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dz)*(dr32inv*mu0/(4*pi))*dnormzdomega(:,index_coil,l_coil+1) &
                    + (3*dx*norm_plasma_dot_norm_coil &
                    - 15*dx*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                    + 3*(normal_plasma(1,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                    + normal_coil(1,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))&
                    *drdomega(1,index_coil,l_coil+1,:) &
                    + (3*dy*norm_plasma_dot_norm_coil &
                    - 15*dy*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                    + 3*(normal_plasma(2,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                    + normal_coil(2,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))&
                    *drdomega(2,index_coil,l_coil+1,:) &
                    + (3*dz*norm_plasma_dot_norm_coil &
                    - 15*dz*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                    + 3*(normal_plasma(3,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                    + normal_coil(3,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))&
                    *drdomega(3,index_coil,l_coil+1,:)
                 dhdomega(:,index_plasma) = dhdomega(:,index_plasma) + &
                     (dddomega(1,:,indexl_coil)*dy_norm3 &
                    +  dddomega(2,:,indexl_coil)*dz_norm1 &
                    +  dddomega(3,:,indexl_coil)*dx_norm2 &
                    -  dddomega(3,:,indexl_coil)*dy_norm1 &
                    -  dddomega(1,:,indexl_coil)*dz_norm2 &
                    -  dddomega(2,:,indexl_coil)*dx_norm3)*2*pi*dr32inv
                 dhdomega(:,index_plasma) = dhdomega(:,index_plasma) &
                    - (factor_for_h(1,itheta_coil,izetal_coil)*drdomega(2,index_coil,l_coil+1,:) &
                    *normal_plasma(3,itheta_plasma,izeta_plasma) &
                    +  factor_for_h(2,itheta_coil,izetal_coil)*drdomega(3,index_coil,l_coil+1,:) &
                    *normal_plasma(1,itheta_plasma,izeta_plasma) &
                    +  factor_for_h(3,itheta_coil,izetal_coil)*drdomega(1,index_coil,l_coil+1,:) &
                    *normal_plasma(2,itheta_plasma,izeta_plasma) &
                    -  factor_for_h(3,itheta_coil,izetal_coil)*drdomega(2,index_coil,l_coil+1,:) &
                    *normal_plasma(1,itheta_plasma,izeta_plasma) &
                    -  factor_for_h(1,itheta_coil,izetal_coil)*drdomega(3,index_coil,l_coil+1,:) &
                    *normal_plasma(2,itheta_plasma,izeta_plasma) &
                    -  factor_for_h(2,itheta_coil,izetal_coil)*drdomega(1,index_coil,l_coil+1,:) &
                    *normal_plasma(3,itheta_plasma,izeta_plasma))*dr32inv
                  dhdomega(:,index_plasma) = dhdomega(:,index_plasma) &
                    + (drdomega(1,index_coil,l_coil+1,:)*dx + drdomega(2,index_coil,l_coil+1,:)*dy &
                    + drdomega(3,index_coil,l_coil+1,:)*dz)*this_h*3*dr2inv
                endif
              end do
           end do
        end do
     end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL

  call system_clock(toc)
  print *,"Done. Took",real(toc-tic)/countrate,"sec."
  
  h = h * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi))
  inductance = inductance * (mu0/(4*pi))
  if (sensitivity_option > 1) then
    dhdomega = dhdomega * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi))
  endif
  deallocate(factor_for_h)
  Bnormal_from_net_coil_currents = reshape(h, (/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma
  !Bnormal_from_net_coil_currents = transpose(reshape(h, (/ nzeta_plasma, ntheta_plasma /))) / norm_normal_plasma
  
  call system_clock(tic)

  ! For some reason, the BLAS matrix-matrix multiplication function DGEMM sometimes causes the
  ! program to crash on Edison unless you are careful to use the Intel MKL instead of Cray LibSci.
  ! If you like, you can use the following method which is slower but more reliable:
  !    inductance = matmul(transpose(basis_functions_1), matmul(inductance_xbasis, basis_functions_2))

  !*******************************************************
  ! Call BLAS3 subroutine DGEMM for matrix multiplications:
  !*******************************************************

  ! Here we carry out g = inductance * basis_functions
  ! A = inductance
  ! B = basis_functions
  ! C = g
  M = ntheta_plasma*nzeta_plasma ! # rows of A
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A and B
  LDA = M
  LDB = K
  LDC = M
  TRANSA = 'N' ! No transposes
  TRANSB = 'N'
  g = 0
  BLAS_ALPHA=dtheta_coil*dzeta_coil
  BLAS_BETA=0
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,inductance,LDA,basis_functions,LDB,BLAS_BETA,g,LDC)

  !$OMP PARALLEL
  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER
  !$OMP DO 
  do iomega = 1, nomega_coil
    ! With matmuls is about 10x slower so commenting out
    !dgdomega(iomega,:,:) = dtheta_coil*dzeta_coil*matmul(dinductancedomega(iomega,:,:), basis_functions)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedomega(iomega,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdomega(iomega,:,:),LDC)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call system_clock(toc)
  print *,"inductance*basis_functions:",real(toc-tic)/countrate,"sec."

  allocate(matrix_B(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(matrix_K(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS_B(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS_K(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_normal_plasma_inv1D(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_normal_coil_inv1D(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(g_over_N_plasma(ntheta_plasma*nzeta_plasma,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_x_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_y_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(f_z_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  call system_clock(tic)

  RHS_B = -dtheta_plasma*dzeta_plasma*matmul( &
       reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)), g)

  call system_clock(toc)
  print *,"Form RHS_B:",real(toc-tic)/countrate,"sec."

  norm_normal_plasma_inv1D = reshape(1/norm_normal_plasma, (/ ntheta_plasma*nzeta_plasma /))
  norm_normal_coil_inv1D   = reshape(1/norm_normal_coil,   (/ ntheta_coil  *nzeta_coil /))
  do j = 1,num_basis_functions
     g_over_N_plasma(:,j) = g(:,j) * norm_normal_plasma_inv1D
     f_x_over_N_coil(:,j) = f_x(:,j) * norm_normal_coil_inv1D
     f_y_over_N_coil(:,j) = f_y(:,j) * norm_normal_coil_inv1D
     f_z_over_N_coil(:,j) = f_z(:,j) * norm_normal_coil_inv1D
     if (sensitivity_option > 2) then
        ! I'm premultiplying this
        do iomega = 1, nomega_coil
          f_xdNdomega_over_N_coil2(iomega,:,j) = f_x(:,j)*norm_normal_coil_inv1D*norm_normal_coil_inv1D &
          * reshape(dnorm_normaldomega(iomega,:,:), (/ntheta_coil*nzeta_coil/))
          f_ydNdomega_over_N_coil2(iomega,:,j) = f_y(:,j)*norm_normal_coil_inv1D*norm_normal_coil_inv1D &
          * reshape(dnorm_normaldomega(iomega,:,:), (/ntheta_coil*nzeta_coil/))
          f_zdNdomega_over_N_coil2(iomega,:,j) = f_z(:,j)*norm_normal_coil_inv1D*norm_normal_coil_inv1D &
          * reshape(dnorm_normaldomega(iomega,:,:), (/ntheta_coil*nzeta_coil/))
        enddo
     endif
  end do

  if (sensitivity_option > 2) then
    call system_clock(tic)

    do iomega = 1, nomega_coil
      dRHS_Bdomega(iomega,:) = -dtheta_plasma*dzeta_plasma*matmul( &
      reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)), dgdomega(iomega,:,:))
      dRHS_Bdomega(iomega,:) = dRHS_Bdomega(iomega,:) - &
        dtheta_plasma*dzeta_plasma*matmul(transpose(g_over_N_plasma),dhdomega(iomega,:))
    enddo

    call system_clock(toc)
    print *,"Form dRHS_Bdomega: ",real(toc-tic)/countrate,"sec."

  endif

  matrix_B = 0

  deallocate(norm_normal_plasma_inv1D)
  deallocate(norm_normal_coil_inv1D)

  call system_clock(toc)
  print *,"Prepare for matrix_B:",real(toc-tic)/countrate,"sec."
  call system_clock(tic)

  ! Here we carry out matrix_B = (dtheta*dzeta)*(g ^ T) * g_over_N_plasma
  ! A = g
  ! B = g_over_N_plasma
  ! C = inductance
  M = num_basis_functions ! # rows of A^T
  N = num_basis_functions ! # cols of B
  K = ntheta_plasma*nzeta_plasma ! Common dimension of A^T and B
  LDA = K ! Would be M if not taking the transpose.
  LDB = K
  LDC = M
  TRANSA = 'T' ! DO take a transpose!
  TRANSB = 'N'
  BLAS_ALPHA = dtheta_plasma*dzeta_plasma
  BLAS_BETA=0
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,g,LDA,g_over_N_plasma,LDB,BLAS_BETA,matrix_B,LDC)

  call system_clock(toc)
  print *,"matmul for matrix_B:",real(toc-tic)/countrate,"sec."

  call system_clock(tic)

  if (sensitivity_option > 2) then
    do iomega = 1, nomega_coil
      ! g_over_N_plasma(ntheta_plasma*nzeta_plasma,num_basis_functions)
      ! dgdomega(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions)
      !dmatrix_Bdomega(iomega,:,:) = 2*dtheta_plasma*dzeta_plasma*matmul(transpose(g_over_N_plasma),dgdomega(iomega,:,:))
      dmatrix_Bdomega(iomega,:,:) = dtheta_plasma*dzeta_plasma*(matmul(transpose(g_over_N_plasma),dgdomega(iomega,:,:)) &
        + matmul(transpose(dgdomega(iomega,:,:)),g_over_N_plasma))
      !call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dgdomega(iomega,:,:),LDA,g_over_N_plasma,LDB,2*BLAS_BETA,dmatrix_Bdomega(iomega,:,:),LDC)
    enddo
  endif

  call system_clock(toc)
  print *,"matmul for dmatrix_Bdomega:",real(tic-toc)/countrate,"sec."

  deallocate(g_over_N_plasma)

  matrix_K = 0
  dmatrix_Kdomega = 0

  call system_clock(tic)
  ! Here we carry out matrix_K += dtheta*dzeta*(f_x ^ T) * f_x_over_N_coil
  ! A = f_x
  ! B = f_x_over_N_plasma
  ! C = matrix_K
  M = num_basis_functions ! # rows of A^T
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A^T and B
  LDA = K ! Would be M if not taking the transpose.
  LDB = K
  LDC = M
  TRANSA = 'T' ! DO take a transpose!
  TRANSB = 'N'
  BLAS_ALPHA = dtheta_coil*dzeta_coil
  BLAS_BETA=1
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_x,LDA,f_x_over_N_coil,LDB,BLAS_BETA,matrix_K,LDC)

  call system_clock(toc)
  print *,"matmul for matrix_K:",real(tic-toc)/countrate,"sec."


  ! Construct dmatrix_Kdomega
  if (sensitivity_option > 2) then
    call system_clock(tic)
    do iomega = 1, nomega_coil
      ! f_xdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions)
      ! f_x(ntheta_coil*nzeta_coil, num_basis_functions)
      ! dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions)
      ! f_x_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions)
      ! This should be changed to DGEMM eventually
      dmatrix_Kdomega(iomega,:,:) = dtheta_coil*dzeta_coil*(matmul(transpose(dfxdomega(iomega,:,:)),f_x_over_N_coil) &
        + matmul(transpose(dfydomega(iomega,:,:)),f_y_over_N_coil) &
        + matmul(transpose(dfzdomega(iomega,:,:)),f_z_over_N_coil) &
        + matmul(transpose(f_x_over_N_coil),dfxdomega(iomega,:,:)) &
        + matmul(transpose(f_y_over_N_coil),dfydomega(iomega,:,:)) &
        + matmul(transpose(f_z_over_N_coil),dfzdomega(iomega,:,:)) &
        - matmul(transpose(f_x),f_xdNdomega_over_N_coil2(iomega,:,:)) &
        - matmul(transpose(f_y),f_ydNdomega_over_N_coil2(iomega,:,:)) &
        - matmul(transpose(f_z),f_zdNdomega_over_N_coil2(iomega,:,:)))
      !call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dfxdomega(iomega,:,:),LDA,f_x_over_N_coil,LDB,2*BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
      !call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_x,LDA,f_xdNdomega_over_N_coil2(iomega,:,:),LDB,-BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
    enddo
    call system_clock(toc)
    print *,"matmul for dmatrix_Kdomega:",real(tic-toc)/countrate,"sec."
  endif



  call system_clock(tic)
  ! Here we carry out matrix_K += dtheta*dzeta*(f_y ^ T) * f_y_over_N_coil
  ! A = f_y
  ! B = f_y_over_N_plasma
  ! C = matrix_K
  M = num_basis_functions ! # rows of A^T
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A^T and B
  LDA = K ! Would be M if not taking the transpose.
  LDB = K
  LDC = M
  TRANSA = 'T' ! DO take a transpose!
  TRANSB = 'N'
  BLAS_ALPHA = dtheta_coil*dzeta_coil
  BLAS_BETA=1
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_y,LDA,f_y_over_N_coil,LDB,BLAS_BETA,matrix_K,LDC)

  call system_clock(toc)
  print *,"matmul 2 for matrix_K:",real(toc-tic)/countrate,"sec."

!  ! Construct dmatrix_Kdomega
!  if (sensitivity_option > 2) then
!    do iomega = 1, nomega_coil
!      call DGEMM(TRANSA,TRANSB,M,N,K,2*BLAS_ALPHA,dfydomega(iomega,:,:),LDA,f_y_over_N_coil,LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!      call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_y,LDA,f_ydNdomega_over_N_coil2(iomega,:,:),LDB,-BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!    enddo
!  endif


  call system_clock(tic)
  ! Here we carry out matrix_K += dtheta*dzeta*(f_z ^ T) * f_z_over_N_coil
  ! A = f_z
  ! B = f_z_over_N_plasma
  ! C = matrix_K
  M = num_basis_functions ! # rows of A^T
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A^T and B
  LDA = K ! Would be M if not taking the transpose.
  LDB = K
  LDC = M
  TRANSA = 'T' ! DO take a transpose!
  TRANSB = 'N'
  BLAS_ALPHA = dtheta_coil*dzeta_coil
  BLAS_BETA=1
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_z,LDA,f_z_over_N_coil,LDB,BLAS_BETA,matrix_K,LDC)

  call system_clock(toc)
  print *,"matmul 3 for matrix_K:",real(toc-tic)/countrate,"sec."

!  ! Construct dmatrix_Kdomega
!  if (sensitivity_option > 2) then
!    do iomega = 1, nomega_coil
!      call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dfzdomega(iomega,:,:),LDA,f_z_over_N_coil,LDB,2*BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!      call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_z,LDA,f_zdNdomega_over_N_coil2(iomega,:,:),LDB,-BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!    enddo
!  endif

  call system_clock(tic)

  RHS_K = (matmul(d_x, f_x_over_N_coil) + matmul(d_y, f_y_over_N_coil) + matmul(d_z, f_z_over_N_coil)) &
       * (dtheta_coil*dzeta_coil)

  call system_clock(toc)
  print *,"Matmuls for RHS_K:",real(toc-tic)/countrate,"sec."

  ! Compute dRHS_Kdomega
  ! f_x_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions)
  ! d_x(ntheta_coil*nzeta_coil)
  ! f_xdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions)
  ! dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions
  if (sensitivity_option > 2) then
    call system_clock(tic)
    do iomega = 1,nomega_coil
      dRHS_Kdomega(iomega,:) = dtheta_coil*dzeta_coil*(matmul(dddomega(1,iomega,1:ntheta_coil*nzeta_coil),f_x_over_N_coil) &
        + matmul(dddomega(2,iomega,1:ntheta_coil*nzeta_coil),f_y_over_N_coil) &
        + matmul(dddomega(3,iomega,1:ntheta_coil*nzeta_coil),f_z_over_N_coil) &
        - matmul(transpose(f_xdNdomega_over_N_coil2(iomega,:,:)),d_x) &
        - matmul(transpose(f_ydNdomega_over_N_coil2(iomega,:,:)),d_y) &
        - matmul(transpose(f_zdNdomega_over_N_coil2(iomega,:,:)),d_z) &
        + matmul(transpose(dfxdomega(iomega,:,:)),d_x/reshape(norm_normal_coil, (/ntheta_coil*nzeta_coil/))) &
        + matmul(transpose(dfydomega(iomega,:,:)),d_y/reshape(norm_normal_coil, (/ntheta_coil*nzeta_coil/))) &
        + matmul(transpose(dfzdomega(iomega,:,:)),d_z/reshape(norm_normal_coil, (/ntheta_coil*nzeta_coil/))))
    enddo
    call system_clock(toc)
    print *,"Matmuls for dRHS_Kdomega:",real(toc-tic)/countrate,"sec."
  endif

  deallocate(f_x_over_N_coil)
  deallocate(f_y_over_N_coil)
  deallocate(f_z_over_N_coil)
  if (sensitivity_option > 1) then
    deallocate(dinductancednorm)
    deallocate(dinductancedr)
  endif
  if (sensitivity_option > 2) then
    deallocate(f_xdNdomega_over_N_coil2)
    deallocate(f_ydNdomega_over_N_coil2)
    deallocate(f_zdNdomega_over_N_coil2)
  endif

end subroutine build_matrices


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
