subroutine init_sensitivity()

  use global_variables
  use stel_constants
  use init_Fourier_modes_mod

  implicit none

  integer :: iomega, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry, tic, toc, countrate, j_basis
  !real(dp), dimension(:,:,:,:), allocatable :: dnormdrmnc, dnormdrmns, dnormdzmnc, dnormdzmns
  real(dp), dimension(:), allocatable :: drdomega, dinductancednorm, dinductancedr
  !real(dp), dimension(:,:,:), allocatable :: dinductancedrmnc, dinductancedrmns, dinductancedzmnc, dinductancedzmns

  real(dp) :: angle, angle2, sinangle, cosangle
  real(dp) :: sinangle2, cosangle2, dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
  real(dp) :: domegadxdtheta, domegadxdzeta, domegadydtheta, domegadydzeta
  real(dp) :: domegadzdtheta, domegadzdzeta
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
  integer :: izeta_plasma, itheta_plasma, izeta_coil, itheta_coil
  integer :: index_coil, l_coil, izetal_coil, index_plasma, imn_coil
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, dr52inv
  real(dp) :: angle_coil, sinangle_coil, cosangle_coil
  real(dp), dimension(:,:), allocatable :: Afactorx, Afactory, Afactorz, Afactor
  real(dp), dimension(:,:), allocatable :: bfactorx, bfactory, bfactorz
  real(dp), dimension(:), allocatable :: norm_normal_coil_inv1D
  real(dp), dimension(:), allocatable :: norm_normal_plasma_inv1D
  real(dp), dimension(:,:), allocatable :: dnorm_normaldomega2D

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0

  ! Initialize Fourier arrays - need to include m = 0, n = 0 mode
  call init_Fourier_modes_sensitivity(mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, nomega_coil, &
    xm_sensitivity, xn_sensitivity, omega_coil)

  print *, "nomega_coil: ", nomega_coil

  allocate(drdomega(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancednorm(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedr(3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedomega(nomega_coil, ntheta_plasma*nzeta_plasma, &
  ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldomega(nomega_coil,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dddomega(3, nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dgdomega(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormxdomega(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormydomega(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormzdomega(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfydomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfzdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dAKdomega(num_basis_functions,num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dABdomega(num_basis_functions,num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dbKdomega(num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dbBdomega(num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactorx(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactory(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactorz(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactor(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactorx(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactory(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactorz(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_normal_coil_inv1D(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldomega2D(nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  print *,"Allocation complete."

  ! Computing chi2_B sensitivity
  print *,"Computing chi2_B sensitivity."
  call system_clock(tic,countrate)

  ! These quantities used to compute dnormdomega

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

        do iomega = 1,nomega_coil

          ! For Fourier decomposition of surface, need to index using izetal_coil
          angle = xm_sensitivity(iomega)*theta_coil(itheta_coil) - xn_sensitivity(iomega)*zetal_coil(izetal_coil)
          sinangle = sin(angle)
          cosangle = cos(angle)

          ! omega = rmnc
          if (omega_coil(iomega) == 1) then
            domegadxdtheta = -xm_sensitivity(iomega)*sinangle*cosangle2
            domegadxdzeta = xn_sensitivity(iomega)*sinangle*cosangle2 - cosangle*sinangle2
            domegadydtheta = -xm_sensitivity(iomega)*sinangle*sinangle2
            domegadydzeta = xn_sensitivity(iomega)*sinangle*sinangle2 + cosangle*cosangle2
            domegadzdtheta = 0
            domegadzdzeta = 0
            drdomega(1) = cosangle*cosangle2
            drdomega(2) = cosangle*sinangle2
            drdomega(3) = 0
          endif
          ! omega = zmns
          if (omega_coil(iomega) == 2) then
            domegadxdtheta = 0
            domegadxdzeta = 0
            domegadydtheta = 0
            domegadydzeta = 0
            domegadzdtheta = xm_sensitivity(iomega)*cosangle
            domegadzdzeta = -xn_sensitivity(iomega)*cosangle
            drdomega(1) = 0
            drdomega(2) = 0
            drdomega(3) = sinangle2
          endif
          ! omega = rmns
          if (omega_coil(iomega) == 3) then
            domegadxdtheta = xm_sensitivity(iomega)*cosangle*cosangle2
            domegadxdzeta = -xn_sensitivity(iomega)*cosangle*cosangle2 - sinangle*sinangle2!
            domegadydtheta = xm_sensitivity(iomega)*cosangle*sinangle2
            domegadydzeta = -xn_sensitivity(iomega)*cosangle*sinangle2 + sinangle*cosangle2
            domegadzdtheta = 0
            domegadzdzeta = 0
            drdomega(1) = sinangle*cosangle2
            drdomega(2) = sinangle*sinangle2
            drdomega(3) = 0
          endif
          ! omega = zmnc
          if (omega_coil(iomega) == 4) then
            domegadxdtheta = 0
            domegadxdzeta = 0
            domegadydtheta = 0
            domegadydzeta = 0
            domegadzdtheta = -xm_sensitivity(iomega)*sinangle
            domegadzdzeta = xn_sensitivity(iomega)*sinangle
            drdomega(1) = 0
            drdomega(2) = 0
            drdomega(3) = cosangle2
          endif

          ! dnormdomega stored for use in chi_k^2 sensitivity
          ! dnormdomega(3, nomega_coil,ntheta_coil,nzetal_coil)
          dnormxdomega(iomega, itheta_coil, izetal_coil) = domegadydzeta*dzdtheta + domegadzdtheta*dydzeta &
            - domegadydtheta*dzdzeta - domegadzdzeta*dydtheta
          dnormydomega(iomega, itheta_coil, izetal_coil) = domegadzdzeta*dxdtheta + domegadxdtheta*dzdzeta &
            - domegadzdtheta*dxdzeta - domegadxdzeta*dzdtheta
          dnormzdomega(iomega, itheta_coil, izetal_coil) = domegadxdzeta*dydtheta + domegadydtheta*dxdzeta &
            - domegadxdtheta*dydzeta - domegadydzeta*dxdtheta

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

              ! This quantity is summed over l_coil
              dinductancedomega(iomega, index_plasma,index_coil) = dinductancedomega(iomega, index_plasma, index_coil) &
                + dinductancednorm(1)*dnormxdomega(iomega,itheta_coil,izetal_coil) &
                + dinductancednorm(2)*dnormydomega(iomega,itheta_coil,izetal_coil) &
                + dinductancednorm(3)*dnormzdomega(iomega,itheta_coil,izetal_coil) &
                + dinductancedr(1)*drdomega(1) + dinductancedr(2)*drdomega(2) &
                + dinductancedr(3)*drdomega(3)
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

  ! Now need to multiply dinductancedomega by basis_functions to compute sensitivity of g_j
  ! dincutancedomega(nomega_coil, ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil)
  ! basis_functions(ntheta_coil*nzeta_coil, num_basis_functions)
  ! dgdomega(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions)

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
  dgdomega = 0

  do iomega = 1, nomega_coil
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedomega(iomega,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdomega(iomega,:,:),LDC)
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
      do iomega = 1,nomega_coil

        ! For geometric components (drdtheta, etc.) need to use izetal_coil
        angle = xm_sensitivity(iomega)*theta_coil(itheta_coil) - xn_sensitivity(iomega)*zetal_coil(izeta_coil)
        sinangle = sin(angle)
        cosangle = cos(angle)

        ! omega = rmnc
        if (omega_coil(iomega) == 1) then
          domegadxdtheta = -xm_sensitivity(iomega)*sinangle*cosangle2
          domegadxdzeta = xn_sensitivity(iomega)*sinangle*cosangle2 - cosangle*sinangle2
          domegadydtheta = -xm_sensitivity(iomega)*sinangle*sinangle2
          domegadydzeta = xn_sensitivity(iomega)*sinangle*sinangle2 + cosangle*cosangle2
          domegadzdtheta = 0
          domegadzdzeta = 0
        endif
        ! omega = zmns
        if (omega_coil(iomega) == 2) then
          domegadxdtheta = 0
          domegadxdzeta = 0
          domegadydtheta = 0
          domegadydzeta = 0
          domegadzdtheta = xm_sensitivity(iomega)*cosangle
          domegadzdzeta = -xn_sensitivity(iomega)*cosangle
        endif
        ! omega = rmns
        if (omega_coil(iomega) == 3) then
          domegadxdtheta = xm_sensitivity(iomega)*cosangle*cosangle2
          domegadxdzeta = -xn_sensitivity(iomega)*cosangle*cosangle2 - sinangle*sinangle2!
          domegadydtheta = xm_sensitivity(iomega)*cosangle*sinangle2
          domegadydzeta = -xn_sensitivity(iomega)*cosangle*sinangle2 + sinangle*cosangle2
          domegadzdtheta = 0
          domegadzdzeta = 0
        endif
        ! omega = zmnc
        if (omega_coil(iomega) == 4) then
          domegadxdtheta = 0
          domegadxdzeta = 0
          domegadydtheta = 0
          domegadydzeta = 0
          domegadzdtheta = -xm_sensitivity(iomega)*sinangle
          domegadzdzeta = xn_sensitivity(iomega)*sinangle
        endif

        ! Here the sensitivity is l_coil periodic - indexed by izeta_coil
        dnorm_normaldomega(iomega,itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdomega(iomega, itheta_coil, izeta_coil) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydomega(iomega, itheta_coil, izeta_coil) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdomega(iomega, itheta_coil, izeta_coil)) &
            /norm_normal_coil(itheta_coil,izeta_coil)

        dddomega(1, iomega, index_coil) = (net_poloidal_current_Amperes*domegadxdtheta &
          - net_toroidal_current_Amperes*domegadxdzeta)/twopi
        dddomega(2, iomega, index_coil) = (net_poloidal_current_Amperes*domegadydtheta &
          - net_toroidal_current_Amperes*domegadydzeta)/twopi
        dddomega(3, iomega, index_coil) = (net_poloidal_current_Amperes*domegadzdtheta &
          - net_toroidal_current_Amperes*domegadzdzeta)/twopi

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
              dfxdomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadxdtheta &
                + xm_coil(imn_coil)*domegadxdzeta)
              dfydomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadydtheta &
                + xm_coil(imn_coil)*domegadydzeta)
              dfzdomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadzdtheta &
                + xm_coil(imn_coil)*domegadzdzeta)
            else
              dfxdomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadxdtheta &
                + xm_coil(imn_coil)*domegadxdzeta)
              dfydomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadydtheta &
                + xm_coil(imn_coil)*domegadydzeta)
              dfzdomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadzdtheta &
                + xm_coil(imn_coil)*domegadzdzeta)
            end if
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Compute matrices needed for dFKdomega
  ! dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions)
  ! dnorm_normaldrmnc(nomega_coil,ntheta_coil,nzeta_coil)
  ! f_x(ntheta_coil*nzeta_coil, num_basis_functions)
  ! norm_normal(ntheta_coil, nzeta_coil)
  ! dAKdrmnc(num_basis_functions,num_basis_functions,nomega_coil)
  ! Afactorx(ntheta*nzeta, num_basis_functions)
  ! dddrmnc(3, nomega_coil,ntheta_coil*nzeta_coil)
  ! d_x(ntheta_coil*nzeta_coil)
  ! bfactorx(ntheta_coil*nzeta_coil, num_basis_functions)
  norm_normal_coil_inv1D   = reshape(1/norm_normal_coil,   (/ ntheta_coil*nzeta_coil /))
  do iomega = 1, nomega_coil
    dAKdomega(:,:,iomega) = dtheta_coil*dzeta_coil*(2*(matmul(transpose(f_x), dfxdomega(iomega,:,:)) &
      + matmul(transpose(f_y), dfydomega(iomega,:,:)) + matmul(transpose(f_z),dfzdomega(iomega,:,:))))
    dnorm_normaldomega2D(iomega,:) = reshape(dnorm_normaldomega(iomega,:,:),(/ ntheta_coil * nzeta_coil /))
    do j_basis = 1, num_basis_functions
      Afactorx(:, j_basis) = f_x(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      Afactory(:, j_basis) = f_y(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      Afactorz(:, j_basis) = f_z(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      bfactorx(:,j_basis) = f_x(:,j_basis)*norm_normal_coil_inv1D
      bfactory(:,j_basis) = f_y(:,j_basis)*norm_normal_coil_inv1D
      bfactorz(:,j_basis) = f_z(:,j_basis)*norm_normal_coil_inv1D
    enddo
    dAKdomega(:,:,iomega) = dAKdomega(:,:,iomega) - dtheta_coil*dzeta_coil*(matmul(transpose(f_x),Afactorx) &
     + matmul(transpose(f_y),Afactory) + matmul(transpose(f_z), Afactorz))
    dbKdomega(:,iomega) = dtheta_coil*dzeta_coil*(matmul(transpose(bfactorx),dddomega(1,iomega,:)) + &
      matmul(transpose(bfactory),dddomega(2,iomega,:)) + matmul(transpose(bfactorz),dddomega(3,iomega,:)))
    dbKdomega(:,iomega) = dbKdomega(:,iomega) + dtheta_coil*dzeta_coil*(matmul(transpose(dfxdomega(iomega,:,:)), &
      d_x*norm_normal_coil_inv1D) + matmul(transpose(dfydomega(iomega,:,:)),d_y*norm_normal_coil_inv1D) &
      + matmul(transpose(dfzdomega(iomega,:,:)),d_z*norm_normal_coil_inv1D))
    dbKdomega(:,iomega) = dbKdomega(:,iomega) - dtheta_coil*dzeta_coil*(matmul(transpose(f_x), &
      d_x*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D) + &
    matmul(transpose(f_y),d_y*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D) &
      + matmul(transpose(f_z),d_z*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D))
  enddo

  ! Compute matrices needed for dFBdomega
  ! g(ntheta_plasma*nzeta_plasma, num_basis_functions)
  ! dgdrmnc(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions)
  norm_normal_plasma_inv1D   = reshape(1/norm_normal_plasma,   (/ ntheta_plasma*nzeta_plasma /))
  do iomega = 1, nomega_coil
    do j_basis = 1, num_basis_functions
      Afactor(:,j_basis) = dgdomega(iomega,:,j_basis)*norm_normal_plasma_inv1D
    enddo
    dABdomega(:,:,iomega) = 2*dtheta_plasma*dzeta_plasma*(matmul(transpose(g),Afactor))
    dbBdomega(:,iomega) = -dtheta_plasma*dzeta_plasma*matmul(transpose(dgdomega(iomega,:,:)), &
      reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)))
  enddo

  call system_clock(toc)
  print *,"chi2_K sensitivity in init_sensitivity: ",real(toc-tic)/countrate," sec."
  print *,"Deallocating."

  deallocate(drdomega)
  deallocate(Afactorx)
  deallocate(Afactory)
  deallocate(Afactorz)
  deallocate(bfactorx)
  deallocate(bfactory)
  deallocate(bfactorz)
  deallocate(Afactor)

  print *,"Init sensitivity complete."

end subroutine init_sensitivity





