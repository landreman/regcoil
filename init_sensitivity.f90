subroutine init_sensitivity()

  use global_variables
  use stel_constants
  use init_Fourier_modes_mod
  use omp_lib

  implicit none

  integer :: iomega, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry, tic, toc, countrate, j_basis
  real(dp), dimension(:,:), allocatable :: dinductancednorm, dinductancedr

  real(dp) :: angle, angle2, sinangle, cosangle
  real(dp) :: sinangle2, cosangle2, dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
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

  allocate(drdomega(3,ntheta_coil*nzeta_coil,nfp,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancednorm(3,nfp),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dinductancedr(3,nfp),stat=iflag)
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
  allocate(dnormxdomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormydomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnormzdomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfydomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dfzdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadxdtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadxdzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadydtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadydzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadzdtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(domegadzdzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  if (sensitivity_option == 3) then
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
  endif

  call system_clock(tic,countrate)

  do izeta_coil = 1,nzeta_coil
    do itheta_coil = 1,ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      do l_coil = 0, (nfp-1)
        izetal_coil = izeta_coil + l_coil*nzeta_coil

        angle2 = zetal_coil(izetal_coil)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)

        do iomega = 1,nomega_coil

          ! For Fourier decomposition of surface, need to index using izetal_coil
          ! We are differentiating something that is not n_p periodic and summing over nfp,
          ! so need to index by izetal_coil
          ! Same convention for angle as used by nescin input file
          angle = xm_sensitivity(iomega)*theta_coil(itheta_coil) + nfp*xn_sensitivity(iomega)*zetal_coil(izetal_coil)
          sinangle = sin(angle)
          cosangle = cos(angle)

          ! omega = rmnc
          if (omega_coil(iomega) == 1) then
            domegadxdtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle*cosangle2
            domegadxdzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle*cosangle2 - cosangle*sinangle2
            domegadydtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle*sinangle2
            domegadydzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle*sinangle2 + cosangle*cosangle2
            domegadzdtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadzdzeta(iomega,itheta_coil,izetal_coil) = 0
            drdomega(1,index_coil,l_coil+1,iomega) = cosangle*cosangle2
            drdomega(2,index_coil,l_coil+1,iomega) = cosangle*sinangle2
            drdomega(3,index_coil,l_coil+1,iomega) = 0
          endif
          ! omega = zmns
          if (omega_coil(iomega) == 2) then
            domegadxdtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadxdzeta(iomega,itheta_coil,izetal_coil) = 0
            domegadydtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadydzeta(iomega,itheta_coil,izetal_coil) = 0
            domegadzdtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle
            domegadzdzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle
            drdomega(1,index_coil,l_coil+1,iomega) = 0
            drdomega(2,index_coil,l_coil+1,iomega) = 0
            drdomega(3,index_coil,l_coil+1,iomega) = sinangle
          endif
          ! omega = rmns
          if (omega_coil(iomega) == 3) then
            domegadxdtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle*cosangle2
            domegadxdzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle*cosangle2 - sinangle*sinangle2
            domegadydtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle*sinangle2
            domegadydzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle*sinangle2 + sinangle*cosangle2
            domegadzdtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadzdzeta(iomega,itheta_coil,izetal_coil) = 0
            drdomega(1,index_coil,l_coil+1,iomega) = sinangle*cosangle2
            drdomega(2,index_coil,l_coil+1,iomega) = sinangle*sinangle2
            drdomega(3,index_coil,l_coil+1,iomega) = 0
          endif
          ! omega = zmnc
          if (omega_coil(iomega) == 4) then
            domegadxdtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadxdzeta(iomega,itheta_coil,izetal_coil) = 0
            domegadydtheta(iomega,itheta_coil,izetal_coil) = 0
            domegadydzeta(iomega,itheta_coil,izetal_coil) = 0
            domegadzdtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle
            domegadzdzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle
            drdomega(1,index_coil,l_coil+1,iomega) = 0
            drdomega(2,index_coil,l_coil+1,iomega) = 0
            drdomega(3,index_coil,l_coil+1,iomega) = cosangle
          endif
          dxdtheta = drdtheta_coil(1,itheta_coil,izetal_coil)
          dydtheta = drdtheta_coil(2,itheta_coil,izetal_coil)
          dzdtheta = drdtheta_coil(3,itheta_coil,izetal_coil)
          dxdzeta = drdzeta_coil(1,itheta_coil,izetal_coil)
          dydzeta = drdzeta_coil(2,itheta_coil,izetal_coil)
          dzdzeta = drdzeta_coil(3,itheta_coil,izetal_coil)
          dnormxdomega(iomega,index_coil,l_coil+1) = &
              domegadydzeta(iomega,itheta_coil,izetal_coil)*dzdtheta &
            + domegadzdtheta(iomega,itheta_coil,izetal_coil)*dydzeta &
            - domegadydtheta(iomega,itheta_coil,izetal_coil)*dzdzeta &
            - domegadzdzeta(iomega,itheta_coil,izetal_coil)*dydtheta
          dnormydomega(iomega,index_coil,l_coil+1) = &
              domegadzdzeta(iomega,itheta_coil,izetal_coil)*dxdtheta &
            + domegadxdtheta(iomega,itheta_coil,izetal_coil)*dzdzeta &
            - domegadzdtheta(iomega,itheta_coil,izetal_coil)*dxdzeta &
            - domegadxdzeta(iomega,itheta_coil,izetal_coil)*dzdtheta
          dnormzdomega(iomega,index_coil,l_coil+1) = &
              domegadxdzeta(iomega,itheta_coil,izetal_coil)*dydtheta &
            + domegadydtheta(iomega,itheta_coil,izetal_coil)*dxdzeta &
            - domegadxdtheta(iomega,itheta_coil,izetal_coil)*dydzeta &
            - domegadydzeta(iomega,itheta_coil,izetal_coil)*dxdtheta
        enddo
      enddo
    enddo
  enddo

  call system_clock(toc)
  print *,"Loop over coil geom in init_sensitivity:",real(toc-tic)/countrate," sec."

  call system_clock(tic,countrate)

  !$OMP PARALLEL
  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER
  !$OMP DO PRIVATE(itheta_plasma,index_plasma,x,y,z,izeta_coil,itheta_coil,index_coil,l_coil,izetal_coil,angle2,sinangle2,cosangle2,dx,dy,dz,dr2inv,dr32inv,dr52inv,dr_dot_norm_coil,dr_dot_norm_plasma,norm_plasma_dot_norm_coil,dinductancedr,dinductancednorm)
  do izeta_plasma = 1, nzeta_plasma
    do itheta_plasma = 1, ntheta_plasma
      index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
      x = r_plasma(1,itheta_plasma,izeta_plasma)
      y = r_plasma(2,itheta_plasma,izeta_plasma)
      z = r_plasma(3,itheta_plasma,izeta_plasma)
      do izeta_coil = 1,nzeta_coil
        do itheta_coil = 1,ntheta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
          do l_coil = 0, (nfp-1)
            izetal_coil = izeta_coil + l_coil*nzeta_coil

            angle2 = zetal_coil(izetal_coil)
            sinangle2 = sin(angle2)
            cosangle2 = cos(angle2)

            dx = x - r_coil(1,itheta_coil,izetal_coil)
            dy = y - r_coil(2,itheta_coil,izetal_coil)
            dz = z - r_coil(3,itheta_coil,izetal_coil)

            dr2inv = 1/(dx*dx + dy*dy + dz*dz)
            dr32inv = dr2inv*sqrt(dr2inv)
            dr52inv = dr2inv*dr32inv

            dr_dot_norm_coil = dx*normal_coil(1,itheta_coil,izetal_coil) &
              + dy*normal_coil(2,itheta_coil,izetal_coil) + dz*normal_coil(3,itheta_coil,izetal_coil)
            dr_dot_norm_plasma = dx*normal_plasma(1,itheta_plasma,izeta_plasma) &
              + dy*normal_plasma(2,itheta_plasma,izeta_plasma) + dz*normal_plasma(3, itheta_plasma,izeta_plasma)
            norm_plasma_dot_norm_coil = normal_plasma(1,itheta_plasma,izeta_plasma) &
              *normal_coil(1,itheta_coil,izetal_coil) + normal_plasma(2,itheta_plasma,izeta_plasma) &
              *normal_coil(2,itheta_coil,izetal_coil) + normal_plasma(3,itheta_plasma,izeta_plasma) &
              *normal_coil(3,itheta_coil,izetal_coil)

            dinductancedr(1,l_coil+1) = (3*dx*norm_plasma_dot_norm_coil &
              - 15*dx*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
              + 3*(normal_plasma(1,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
              + normal_coil(1,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
            dinductancedr(2,l_coil+1) = (3*dy*norm_plasma_dot_norm_coil &
              - 15*dy*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
              + 3*(normal_plasma(2,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
              + normal_coil(2,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
            dinductancedr(3,l_coil+1) = (3*dz*norm_plasma_dot_norm_coil &
              - 15*dz*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
              + 3*(normal_plasma(3,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
              + normal_coil(3,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))

            dinductancednorm(1,l_coil+1) = (normal_plasma(1,itheta_plasma,izeta_plasma) &
              - 3*dr2inv*dr_dot_norm_plasma*dx)*(dr32inv*mu0/(4*pi))
            dinductancednorm(2,l_coil+1) = (normal_plasma(2,itheta_plasma,izeta_plasma) &
              - 3*dr2inv*dr_dot_norm_plasma*dy)*(dr32inv*mu0/(4*pi))
            dinductancednorm(3,l_coil+1) = (normal_plasma(3,itheta_plasma,izeta_plasma) &
              - 3*dr2inv*dr_dot_norm_plasma*dz)*(dr32inv*mu0/(4*pi))

            dinductancedomega(:,index_plasma,index_coil) = &
                dinductancedomega(:,index_plasma,index_coil) &
              + dinductancednorm(1,l_coil+1)*dnormxdomega(:,index_coil,l_coil+1) &
              + dinductancednorm(2,l_coil+1)*dnormydomega(:,index_coil,l_coil+1) &
              + dinductancednorm(3,l_coil+1)*dnormzdomega(:,index_coil,l_coil+1) &
              + dinductancedr(1,l_coil+1)*drdomega(1,index_coil,l_coil+1,:) &
              + dinductancedr(2,l_coil+1)*drdomega(2,index_coil,l_coil+1,:) &
              + dinductancedr(3,l_coil+1)*drdomega(3,index_coil,l_coil+1,:)
          end do
          ! Using matmuls seems to be about twice as slow
          !dinductancedomega(:,index_plasma,index_coil) = &
           !   matmul(dnormxdomega(:, index_coil,:),dinductancednorm(1,:)) &
           ! + matmul(dnormydomega(:, index_coil,:),dinductancednorm(2,:)) &
           ! + matmul(dnormzdomega(:, index_coil,:),dinductancednorm(3,:)) &
           ! + matmul(dinductancedr(1,:),drdomega(1,index_coil,:,:)) &
           ! + matmul(dinductancedr(2,:),drdomega(2,index_coil,:,:)) &
           ! + matmul(dinductancedr(3,:),drdomega(3,index_coil,:,:))
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call system_clock(toc)
  print *,"Inductance sensitivity loop in init_sensitivity:",real(toc-tic)/countrate," sec."

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
    !dgdomega(iomega,:,:) = dtheta_coil*dzeta_coil*matmul(dinductancedomega(iomega,:,:), basis_functions)
    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedomega(iomega,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdomega(iomega,:,:),LDC)
  enddo

  call system_clock(toc)
  print *,"loop with DGEMM for g sensitivity: ",real(toc-tic)/countrate," sec."

  ! Computing quantities needed chi2_K sensitivity
  ! Partial derivatives of d, f, and N' computed

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
    do itheta_coil = 1, ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      do iomega = 1,nomega_coil

        ! Here the sensitivity is l_coil periodic - indexed by izeta_coil
        dnorm_normaldomega(iomega,itheta_coil,izeta_coil) = &
           (normal_coil(1, itheta_coil, izeta_coil)*dnormxdomega(iomega,index_coil,1) &
          + normal_coil(2, itheta_coil, izeta_coil)*dnormydomega(iomega,index_coil,1) &
          + normal_coil(3, itheta_coil, izeta_coil)*dnormzdomega(iomega,index_coil,1)) &
            /norm_normal_coil(itheta_coil,izeta_coil)

        dddomega(1, iomega, index_coil) = (net_poloidal_current_Amperes*domegadxdtheta(iomega,itheta_coil,izeta_coil) &
          - net_toroidal_current_Amperes*domegadxdzeta(iomega,itheta_coil,izeta_coil))/twopi
        dddomega(2, iomega, index_coil) = (net_poloidal_current_Amperes*domegadydtheta(iomega,itheta_coil,izeta_coil) &
          - net_toroidal_current_Amperes*domegadydzeta(iomega,itheta_coil,izeta_coil))/twopi
        dddomega(3, iomega, index_coil) = (net_poloidal_current_Amperes*domegadzdtheta(iomega,itheta_coil,izeta_coil) &
          - net_toroidal_current_Amperes*domegadzdzeta(iomega,itheta_coil,izeta_coil))/twopi

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
              dfxdomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadxdtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadxdzeta(iomega,itheta_coil,izeta_coil))
              dfydomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadydtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadydzeta(iomega,itheta_coil,izeta_coil))
              dfzdomega(iomega, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*domegadzdtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadzdzeta(iomega,itheta_coil,izeta_coil))
            else
              dfxdomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadxdtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadxdzeta(iomega,itheta_coil,izeta_coil))
              dfydomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadydtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadydzeta(iomega,itheta_coil,izeta_coil))
              dfzdomega(iomega, index_coil,imn_coil+offset) = -sinangle_coil*(xn_coil(imn_coil)*domegadzdtheta(iomega,itheta_coil,izeta_coil) &
                + xm_coil(imn_coil)*domegadzdzeta(iomega,itheta_coil,izeta_coil))
            end if
          enddo
        enddo
      enddo
    enddo
  enddo

  call system_clock(toc)
  print *,"Main chi2_K sensitivity loop in init_sensitivity: ",real(toc-tic)/countrate," sec."

  if (sensitivity_option == 3) then
    call system_clock(tic,countrate)

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

    call system_clock(toc)
    print *,"Computing matrices for dFKdomega: ",real(toc-tic)/countrate," sec."

    call system_clock(tic,countrate)

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
    print *,"Computing matrices for dFBdomega: ",real(toc-tic)/countrate," sec."
  endif

  deallocate(dinductancednorm)
  deallocate(dinductancedr)
  if (sensitivity_option == 3) then
    deallocate(Afactorx)
    deallocate(Afactory)
    deallocate(Afactorz)
    deallocate(bfactorx)
    deallocate(bfactory)
    deallocate(bfactorz)
    deallocate(Afactor)
    deallocate(norm_normal_coil_inv1D)
    deallocate(dnorm_normaldomega2D)
  endif
  print *,"Init sensitivity complete."

end subroutine init_sensitivity





