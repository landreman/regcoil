subroutine init_sensitivity()

  use global_variables
  use stel_constants
  use init_Fourier_modes_mod

  implicit none

  integer :: imn, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry, tic, toc, countrate
  !real(dp), dimension(:,:,:,:), allocatable :: dnormdrmnc, dnormdrmns, dnormdzmnc, dnormdzmns
  real(dp), dimension(:), allocatable :: drdrmnc, drdrmns, drdzmnc, drdzmns
  real(dp), dimension(:), allocatable :: dinductancednorm, dinductancedr
  !real(dp), dimension(:,:,:), allocatable :: dinductancedrmnc, dinductancedrmns, dinductancedzmnc, dinductancedzmns

  real(dp) :: angle, angle2, sinangle, cosangle
  real(dp) :: sinangle2, cosangle2, dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
  real(dp) :: drmncdxdtheta, drmncdxdzeta, drmncdydtheta, drmncdydzeta
  real(dp) :: drmncdzdtheta, drmncdzdzeta, drmnsdxdtheta, drmnsdxdzeta
  real(dp) :: drmnsdydtheta, drmnsdydzeta, drmnsdzdtheta, drmnsdzdzeta
  real(dp) :: dzmncdxdtheta, dzmncdxdzeta, dzmncdydtheta, dzmncdydzeta
  real(dp) :: dzmncdzdtheta, dzmncdzdzeta, dzmnsdydzeta, dzmnsdzdtheta
  real(dp) :: dzmnsdzdzeta, dzmnsdxdtheta, dzmnsdxdzeta, dzmnsdydtheta
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
  integer :: izeta_plasma, itheta_plasma, izeta_coil, itheta_coil
  integer :: index_coil, l_coil, izetal_coil, index_plasma, imn_coil
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, dr52inv
  real(dp) :: angle_coil, sinangle_coil, cosangle_coil

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0

  ! Initialize Fourier arrays - need to include m = 0, n = 0 mode
  call init_Fourier_modes_sensitivity(mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, &
xm_sensitivity, xn_sensitivity)

  print *, "mnmax_sensitivity: ", mnmax_sensitivity

  allocate(drdrmnc(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdrmns(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdzmnc(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdzmns(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dinductancednorm(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedr(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dinductancedrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
  ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedrmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
  ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedzmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
  ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedzmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
  ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dnorm_normaldrmnc(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldrmns(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldzmnc(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldzmns(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dddrmnc(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dddrmns(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dddzmnc(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dddzmns(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dgdrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dgdrmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dgdzmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dgdzmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dnormxdrmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormxdrmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormxdzmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormxdzmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dnormydrmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormydrmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormydzmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormydzmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dnormzdrmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormzdrmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormzdzmnc(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormzdzmns(mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dfxdrmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfxdrmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfxdzmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfxdzmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dfydrmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfydrmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfydzmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfydzmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(dfzdrmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfzdrmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfzdzmnc(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfzdzmns(mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  print *,"Allocation complete."

  ! Computing chi2_B sensitivity
  print *,"Computing chi2_B sensitivity."
  call system_clock(tic,countrate)

  ! These quantities used to compute dnormd(...)
  drmncdzdtheta = 0
  drmncdzdzeta = 0
  drmnsdzdtheta = 0
  drmnsdzdzeta = 0

  dzmncdxdtheta = 0
  dzmncdxdzeta = 0
  dzmncdydtheta = 0
  dzmncdydzeta = 0

  dzmnsdxdtheta = 0
  dzmnsdxdzeta = 0
  dzmnsdydtheta = 0
  dzmnsdydzeta = 0

  do izeta_coil = 1,nzeta_coil

    do itheta_coil = 1,ntheta_coil

      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil

      do l_coil = 0, (nfp-1)

        izetal_coil = izeta_coil + l_coil*nzeta_coil

        angle2 = zetal_coil(izetal_coil)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)

        dxdtheta = drdtheta_coil(1,itheta_coil,izetal_coil)
        dxdzeta = drdzeta_coil(1,itheta_coil,izetal_coil)
        dydzeta = drdzeta_coil(2,itheta_coil,izetal_coil)
        dydtheta = drdtheta_coil(2,itheta_coil,izetal_coil)
        dzdtheta = drdtheta_coil(3,itheta_coil,izetal_coil)
        dzdzeta = drdzeta_coil(3,itheta_coil,izetal_coil)

        do imn = 1,mnmax_sensitivity

          ! For Fourier decomposition of surface, need to index using izetal_coil
          angle = xm_sensitivity(imn)*theta_coil(itheta_coil) - xn_sensitivity(imn)*zetal_coil(izetal_coil)
          sinangle = sin(angle)
          cosangle = cos(angle)

          drmncdxdtheta = -xm_sensitivity(imn)*sinangle*cosangle2
          drmncdxdzeta = xn_sensitivity(imn)*sinangle*cosangle2 - cosangle*sinangle2
          drmncdydtheta = -xm_sensitivity(imn)*sinangle*sinangle2
          drmncdydzeta = xn_sensitivity(imn)*sinangle*sinangle2 + cosangle*cosangle2

          drmnsdxdtheta = xm_sensitivity(imn)*cosangle*cosangle2
          drmnsdxdzeta = -xn_sensitivity(imn)*cosangle*cosangle2 - sinangle*sinangle2!
          drmnsdydtheta = xm_sensitivity(imn)*cosangle*sinangle2
          drmnsdydzeta = -xn_sensitivity(imn)*cosangle*sinangle2 + sinangle*cosangle2

          dzmncdzdtheta = -xm_sensitivity(imn)*sinangle
          dzmncdzdzeta = xn_sensitivity(imn)*sinangle
          dzmnsdzdtheta = xm_sensitivity(imn)*cosangle
          dzmnsdzdzeta = -xn_sensitivity(imn)*cosangle

          ! dnormd(...) stored for use in chi_k^2 sensitivity
          ! dnormd(...)(3, mnmax_sensitivity,ntheta_coil,nzetal_coil)
          dnormxdrmnc(imn, itheta_coil, izetal_coil) = drmncdydzeta*dzdtheta + drmncdzdtheta*dydzeta &
            - drmncdydtheta*dzdzeta - drmncdzdzeta*dydtheta
          dnormydrmnc(imn, itheta_coil, izetal_coil) = drmncdzdzeta*dxdtheta + drmncdxdtheta*dzdzeta &
            - drmncdzdtheta*dxdzeta - drmncdxdzeta*dzdtheta
          dnormzdrmnc(imn, itheta_coil, izetal_coil) = drmncdxdzeta*dydtheta + drmncdydtheta*dxdzeta &
            - drmncdxdtheta*dydzeta - drmncdydzeta*dxdtheta

          dnormxdrmns(imn, itheta_coil, izetal_coil) = drmnsdydzeta*dzdtheta + drmnsdzdtheta*dydzeta &
            - drmnsdydtheta*dzdzeta - drmnsdzdzeta*dydtheta
          dnormydrmns(imn, itheta_coil, izetal_coil) = drmnsdzdzeta*dxdtheta + drmnsdxdtheta*dzdzeta &
            - drmnsdzdtheta*dxdzeta - drmnsdxdzeta*dzdtheta
          dnormzdrmns(imn, itheta_coil, izetal_coil) = drmnsdxdzeta*dydtheta + drmnsdydtheta*dxdzeta &
            - drmnsdxdtheta*dydzeta - drmnsdydzeta*dxdtheta
          
          dnormxdzmnc(imn, itheta_coil, izetal_coil) = dzmncdydzeta*dzdtheta + dzmncdzdtheta*dydzeta &
            - dzmncdydtheta*dzdzeta - dzmncdzdzeta*dydtheta
          dnormydzmnc(imn, itheta_coil, izetal_coil) = dzmncdzdzeta*dxdtheta + dzmncdxdtheta*dzdzeta &
            - dzmncdzdtheta*dxdzeta - dzmncdxdzeta*dzdtheta
          dnormzdzmnc(imn, itheta_coil, izetal_coil) = dzmncdxdzeta*dydtheta + dzmncdydtheta*dxdzeta &
            - dzmncdxdtheta*dydzeta - dzmncdydzeta*dxdtheta

          dnormxdzmns(imn, itheta_coil, izetal_coil) = dzmnsdydzeta*dzdtheta + dzmnsdzdtheta*dydzeta &
            - dzmnsdydtheta*dzdzeta - dzmnsdzdzeta*dydtheta
          dnormydzmns(imn, itheta_coil, izetal_coil) = dzmnsdzdzeta*dxdtheta + dzmnsdxdtheta*dzdzeta &
            - dzmnsdzdtheta*dxdzeta - dzmnsdxdzeta*dzdtheta
          dnormzdzmns(imn, itheta_coil, izetal_coil) = dzmnsdxdzeta*dydtheta + dzmnsdydtheta*dxdzeta &
            - dzmnsdxdtheta*dydzeta - dzmnsdydzeta*dxdtheta

          ! dxdOmega
          drdrmnc(1) = cosangle*cosangle2
          drdrmns(1) = sinangle*cosangle2
          drdzmnc(1) = 0
          drdzmns(1) = 0
          ! dydOmega
          drdrmnc(2) = cosangle*sinangle2
          drdrmns(2) = sinangle*sinangle2
          drdzmnc(2) = 0
          drdzmns(2) = 0
          ! dzdOmega
          drdrmnc(3) = 0
          drdrmns(3) = 0
          drdzmnc(3) = cosangle2
          drdzmns(3) = sinangle2

          do izeta_plasma = 1, nzeta_plasma
            do itheta_plasma = 1, ntheta_plasma
              index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma

              x = r_plasma(1,itheta_plasma,izeta_plasma)
              y = r_plasma(2,itheta_plasma,izeta_plasma)
              z = r_plasma(3,itheta_plasma,izeta_plasma)

              dx = x - r_coil(1,itheta_coil,izetal_coil)
              dy = y - r_coil(2,itheta_coil,izetal_coil)
              dz = z - r_coil(3,itheta_coil,izetal_coil)

              dr2inv = 1/(dx*dx + dy*dy + dz*dz)
              dr32inv = dr2inv*sqrt(dr2inv)
              dr52inv = dr2inv*dr32inv

              dr_dot_norm_coil = dx*normal_coil(1,itheta_coil,izetal_coil) &
                + dy*normal_coil(2,itheta_coil,izetal_coil) + dz*normal_coil(3,itheta_coil,izetal_coil)
              dr_dot_norm_plasma = dx*normal_plasma(1,itheta_plasma,izeta_plasma) &
                + dy*normal_plasma(2,itheta_plasma,izeta_plasma) + dz*normal_plasma(3, itheta_coil,izeta_coil)
              norm_plasma_dot_norm_coil = normal_plasma(1,itheta_plasma,izeta_plasma)*normal_coil(1,itheta_coil, izetal_coil) &
                + normal_plasma(2,itheta_plasma,izeta_plasma)*normal_coil(2,itheta_coil,izetal_coil) &
                + normal_plasma(3,itheta_plasma,izeta_plasma)*normal_coil(3,itheta_coil,izetal_coil)

              ! dgdr before sum over l_coil
              dinductancedr(1) = (3*dx*norm_plasma_dot_norm_coil &
                - 15*dx*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                + 3*(normal_plasma(1,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                + normal_coil(1,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
              dinductancedr(2) = (3*dy*norm_plasma_dot_norm_coil &
                - 15*dy*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                + 3*(normal_plasma(2,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                + normal_coil(2,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
              dinductancedr(3) = (3*dz*norm_plasma_dot_norm_coil &
                - 15*dz*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                + 3*(normal_plasma(3,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                + normal_coil(3,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))

              ! dgdN' before sum over l_coil
              dinductancednorm(1) = (normal_plasma(1,itheta_plasma,izeta_plasma) &
                - 3*dr2inv*dr_dot_norm_plasma*dx)*(dr32inv*mu0/(4*pi))
              dinductancednorm(2) = (normal_plasma(2,itheta_plasma,izeta_plasma) &
                - 3*dr2inv*dr_dot_norm_plasma*dy)*(dr32inv*mu0/(4*pi))
              dinductancednorm(3) = (normal_plasma(3,itheta_plasma,izeta_plasma) &
                - 3*dr2inv*dr_dot_norm_plasma*dz)*(dr32inv*mu0/(4*pi))

              ! These quantity are summed over l_coil
              dinductancedrmnc(imn, index_plasma,index_coil) = dinductancedrmnc(imn, index_plasma, index_coil) &
                + dinductancednorm(1)*dnormxdrmnc(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(2)*dnormydrmnc(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(3)*dnormzdrmnc(imn,itheta_coil,izetal_coil) &
                + dinductancedr(1)*drdrmnc(1) + dinductancedr(2)*drdrmnc(2) &
                + dinductancedr(3)*drdrmnc(3)
              dinductancedrmns(imn, index_plasma,index_coil) = dinductancedrmns(imn, index_plasma,index_coil) &
                + dinductancednorm(1)*dnormxdrmns(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(2)*dnormydrmns(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(3)*dnormzdrmns(imn,itheta_coil,izetal_coil) &
                + dinductancedr(1)*drdrmns(1) + dinductancedr(2)*drdrmns(2) &
                + dinductancedr(3)*drdrmns(3)
              dinductancedzmnc(imn, index_plasma,index_coil) = dinductancedzmnc(imn, index_plasma,index_coil) &
                + dinductancednorm(1)*dnormxdzmnc(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(2)*dnormydzmnc(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(3)*dnormzdzmnc(imn,itheta_coil,izetal_coil) &
                + dinductancedr(1)*drdzmnc(1) + dinductancedr(2)*drdzmnc(2) &
                + dinductancedr(3)*drdzmnc(3)
              dinductancedzmns(imn, index_plasma,index_coil) = dinductancedzmns(imn, index_plasma,index_coil) &
                + dinductancednorm(1)*dnormxdzmns(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(2)*dnormydzmns(imn,itheta_coil,izetal_coil) &
                + dinductancednorm(3)*dnormzdzmns(imn,itheta_coil,izetal_coil) &
                + dinductancedr(1)*drdzmns(1) + dinductancedr(2)*drdzmns(2) &
                + dinductancedr(3)*drdzmns(3)
            end do
          end do
        end do
      end do
    end do
  end do

  call system_clock(toc)
  print *,"chi2_B sensitivity in init_sensitivity: ",real(toc-tic)/countrate," sec."

  print *,"Loop with matrix multiplications."
  call system_clock(tic,countrate)

  ! Now need to multiply dinductanced(...) by basis_functions to compute sensitivity of g_j
  ! dincutanced(...)(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil)
  ! basis_functions(ntheta_coil*nzeta_coil, num_basis_functions)
  ! dgdrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions)

  ! Now we need to compute dgd(...) = dinductanced(...)*basis_functions
   M = ntheta_plasma*nzeta_plasma ! # rows of A
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A and B
  LDA = M
  LDB = K
  LDC = M
  TRANSA = 'N' ! No transposes
  TRANSB = 'N'
  BLAS_ALPHA=dtheta_coil*dzeta_coil
  BLAS_BETA=0
  dgdrmnc = 0
  dgdrmns = 0
  dgdzmnc = 0
  dgdzmns = 0

  do imn = 1, mnmax_sensitivity
    !dgdrmnc(imn,:,:) = matmul(dinductancedrmnc(imn,:,:),basis_functions)
    !dgdrmns(imn,:,:) = matmul(dinductancedrmns(imn,:,:),basis_functions)
    !dgdzmnc(imn,:,:) = matmul(dinductancedzmnc(imn,:,:),basis_functions)
    !dgdzmns(imn,:,:) = matmul(dinductancedzmns(imn,:,:),basis_functions)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedrmnc(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdrmnc(imn,:,:),LDC)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedrmns(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdrmns(imn,:,:),LDC)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedzmnc(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdzmnc(imn,:,:),LDC)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedzmns(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdzmns(imn,:,:),LDC)
  enddo

  call system_clock(toc)
  print *,"sensitivity MM in init_sensitivity: ",real(toc-tic)/countrate," sec."

  ! Computing quantities needed chi2_K sensitivity
  ! Partial derivatives of d, f, and N' computed

  print *,"Computing chi2_K sensitivity."
  call system_clock(tic,countrate)

  ! Needed for computing sensitivity of f
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

  do izeta_coil = 1, nzeta_coil

    ! For geometric components (drdtheta, etc.) need to use izetal_coil
    angle2 = zetal_coil(izeta_coil)
    sinangle2 = sin(angle2)
    cosangle2 = cos(angle2)

    do itheta_coil = 1, ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      do imn = 1, mnmax_sensitivity

        ! For geometric components (drdtheta, etc.) need to use izetal_coil
        angle = xm_sensitivity(imn)*theta_coil(itheta_coil) - xn_sensitivity(imn)*zetal_coil(izeta_coil)
        sinangle = sin(angle)
        cosangle = cos(angle)

        drmncdxdtheta = -xm_sensitivity(imn)*sinangle*cosangle2
        drmncdxdzeta = xn_sensitivity(imn)*sinangle*cosangle2 - cosangle*sinangle2
        drmncdydtheta = -xm_sensitivity(imn)*sinangle*sinangle2
        drmncdydzeta = xn_sensitivity(imn)*sinangle*sinangle2 + cosangle*cosangle2

        drmnsdxdtheta = xm_sensitivity(imn)*cosangle*cosangle2
        drmnsdxdzeta = -xn_sensitivity(imn)*cosangle*cosangle2 - sinangle*sinangle2
        drmnsdydtheta = xm_sensitivity(imn)*cosangle*sinangle2
        drmnsdydzeta = -xn_sensitivity(imn)*cosangle*sinangle2 + sinangle*cosangle2

        dzmncdzdtheta = -xm_sensitivity(imn)*sinangle
        dzmncdzdzeta = xn_sensitivity(imn)*sinangle
        dzmnsdzdtheta = xm_sensitivity(imn)*cosangle
        dzmnsdzdzeta = -xn_sensitivity(imn)*cosangle

        ! Here the sensitivity is l_coil periodic - indexed by izeta_coil
        dnorm_normaldrmnc(imn,itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdrmnc(imn, itheta_coil, izeta_coil) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydrmnc(imn, itheta_coil, izeta_coil) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdrmnc(imn, itheta_coil, izeta_coil)) &
            /norm_normal_coil(itheta_coil,izeta_coil)
        dnorm_normaldrmns(imn, itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdrmns(imn, itheta_coil, izeta_coil) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydrmns(imn, itheta_coil, izeta_coil) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdrmns(imn, itheta_coil, izeta_coil)) &
            /norm_normal_coil(itheta_coil,izeta_coil)
        dnorm_normaldzmnc(imn, itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdzmnc(imn, itheta_coil, izeta_coil) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydzmnc(imn, itheta_coil, izeta_coil) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdzmnc(imn, itheta_coil, izeta_coil)) &
           /norm_normal_coil(itheta_coil,izeta_coil)
        dnorm_normaldzmns(imn, itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdzmns(imn, itheta_coil, izeta_coil) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydzmns(imn, itheta_coil, izeta_coil) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdzmns(imn, itheta_coil, izeta_coil)) &
           /norm_normal_coil(itheta_coil,izeta_coil)

        dddrmnc(1, imn, index_coil) = (net_poloidal_current_Amperes*drmncdxdtheta &
          - net_toroidal_current_Amperes*drmncdxdzeta)/twopi
        dddrmnc(2, imn, index_coil) = (net_poloidal_current_Amperes*drmncdydtheta &
          - net_toroidal_current_Amperes*drmncdydzeta)/twopi
        dddrmnc(3, imn, index_coil) = (net_poloidal_current_Amperes*drmncdzdtheta &
          - net_toroidal_current_Amperes*drmncdzdzeta)/twopi
        dddrmns(1, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdxdtheta &
          - net_toroidal_current_Amperes*drmnsdxdzeta)/twopi
        dddrmns(2, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdydtheta &
          - net_toroidal_current_Amperes*drmnsdydzeta)/twopi
        dddrmns(3, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdzdtheta &
          - net_toroidal_current_Amperes*drmnsdzdzeta)/twopi
        dddzmnc(1, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdxdtheta &
          - net_toroidal_current_Amperes*dzmncdxdzeta)/twopi
        dddzmnc(2, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdydtheta &
          - net_toroidal_current_Amperes*dzmncdydzeta)/twopi
        dddzmnc(3, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdzdtheta &
          - net_toroidal_current_Amperes*dzmncdzdzeta)/twopi
        dddzmns(1, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdxdtheta &
          - net_toroidal_current_Amperes*dzmnsdxdzeta)/twopi
        dddzmns(2, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdydtheta &
          - net_toroidal_current_Amperes*dzmnsdydzeta)/twopi
        dddzmns(3, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdzdtheta &
          - net_toroidal_current_Amperes*dzmnsdzdzeta)/twopi

        do whichSymmetry = minSymmetry, maxSymmetry

          if (whichSymmetry==2 .and. symmetry_option==3) then
            offset = mnmax_coil
          else
            offset = 0
          end if

          do imn_coil = 1, mnmax_coil

            angle_coil = xm_coil(imn_coil)*theta_coil(itheta_coil) - xn_coil(imn_coil)*zeta_coil(izeta_coil)
            cosangle_coil = cos(angle_coil)
            sinangle_coil = sin(angle_coil)
            if (whichSymmetry==1) then
              dfxdrmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdxdtheta &
                + xm_coil(imn_coil)*drmncdxdzeta)
              dfxdrmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdxdtheta &
                + xm_coil(imn_coil)*drmnsdxdzeta)
              dfxdzmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdxdtheta &
                + xm_coil(imn_coil)*dzmncdxdzeta)
              dfxdzmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdxdtheta &
                + xm_coil(imn_coil)*dzmnsdxdzeta)

              dfydrmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdydtheta &
                + xm_coil(imn_coil)*drmncdydzeta)
              dfydrmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdydtheta &
                + xm_coil(imn_coil)*drmnsdydzeta)
              dfydzmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdydtheta &
                + xm_coil(imn_coil)*dzmncdydzeta)
              dfydzmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdydtheta &
                + xm_coil(imn_coil)*dzmnsdydzeta)

              dfzdrmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdzdtheta &
                + xm_coil(imn_coil)*drmncdzdzeta)
              dfzdrmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdzdtheta &
                + xm_coil(imn_coil)*drmnsdzdzeta)
              dfzdzmnc(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdzdtheta &
                + xm_coil(imn_coil)*dzmncdzdzeta)
              dfzdzmns(imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdzdtheta &
                + xm_coil(imn_coil)*dzmnsdzdzeta)

            else
              dfxdrmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmncdxdtheta &
                + xm_coil(imn_coil)*drmncdxdzeta)
              dfxdrmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdxdtheta &
                + xm_coil(imn_coil)*drmnsdxdzeta)
              dfxdzmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdxdtheta &
                + xm_coil(imn_coil)*dzmncdxdzeta)
              dfxdzmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdxdtheta &
                + xm_coil(imn_coil)*dzmnsdxdzeta)

              dfydrmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmncdydtheta &
                + xm_coil(imn_coil)*drmncdydzeta)
              dfydrmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdydtheta &
                + xm_coil(imn_coil)*drmnsdydzeta)
              dfydzmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdydtheta &
                + xm_coil(imn_coil)*dzmncdydzeta)
              dfydzmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdydtheta &
                + xm_coil(imn_coil)*dzmnsdydzeta)

              dfzdrmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmncdzdtheta &
                + xm_coil(imn_coil)*drmncdzdzeta)
              dfzdrmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdzdtheta &
                + xm_coil(imn_coil)*drmnsdzdzeta)
              dfzdzmnc(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdzdtheta &
                + xm_coil(imn_coil)*dzmncdzdzeta)
              dfzdzmns(imn, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdzdtheta &
                + xm_coil(imn_coil)*dzmnsdzdzeta)
            end if
          enddo
        enddo
      enddo
    enddo
  enddo

  call system_clock(toc)
  print *,"chi2_K sensitivity in init_sensitivity: ",real(toc-tic)/countrate," sec."
  print *,"max(dddzmns): ", real(maxval(dddzmns))
  print *,"Deallocating."
  !deallocate(dnormdrmnc)
  !deallocate(dnormdrmns)
  !deallocate(dnormdzmnc)
  !deallocate(dnormdzmns)
  ! in global variables for testing
  !deallocate(dinductancednorm)
  !deallocate(dinductancedr)
  !deallocate(dinductancedrmnc)
  !deallocate(dinductancedrmns)
  !deallocate(dinductancedzmnc)
  !deallocate(dinductancedzmns)
  deallocate(drdrmnc)
  deallocate(drdrmns)
  deallocate(drdzmnc)
  deallocate(drdzmns)

  print *,"Init sensitivity complete."

end subroutine init_sensitivity





