subroutine regcoil_build_matrices()

  use regcoil_variables
  use stel_constants
  use stel_kinds
  use omp_lib
  use regcoil_init_Fourier_modes_mod
  
  implicit none

  integer :: l_coil, itheta_plasma, izeta_plasma, itheta_coil, izeta_coil, izetal_coil
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, dr52inv
  integer :: index_plasma, index_coil, j, imn
  integer :: tic, toc, countrate, iflag
  integer :: minSymmetry, maxSymmetry, whichSymmetry, offset
  real(dp) :: angle, sinangle, cosangle
  real(dp), dimension(:,:,:), allocatable :: factor_for_h
  real(dp), dimension(:), allocatable :: norm_normal_plasma_inv1D, norm_normal_coil_inv1D
  real(dp), dimension(:,:), allocatable :: g_over_N_plasma, f_x_over_N_coil, f_y_over_N_coil, f_z_over_N_coil, f_Laplace_Beltrami_over_N_coil
  real(dp), dimension(:,:), allocatable :: g_theta_theta, g_theta_zeta, g_zeta_zeta
  real(dp), dimension(:,:), allocatable :: d_g_theta_theta_d_theta, d_g_theta_theta_d_zeta
  real(dp), dimension(:,:), allocatable :: d_g_theta_zeta_d_theta, d_g_theta_zeta_d_zeta
  real(dp), dimension(:,:), allocatable :: d_g_zeta_zeta_d_theta, d_g_zeta_zeta_d_zeta
  real(dp), dimension(:,:), allocatable :: d_N_d_theta, d_N_d_zeta
  real(dp), dimension(:,:), allocatable :: Laplace_Beltrami_d_Phi_d_theta_coefficient, Laplace_Beltrami_d_Phi_d_zeta_coefficient
  real(dp), dimension(:,:,:), allocatable :: f_xdNdomega_over_N_coil2,f_ydNdomega_over_N_coil2,f_zdNdomega_over_N_coil2
  real(dp), dimension(:), allocatable :: dinductancednorm, dinductancedr
  real(dp) :: cosangle_xm, sinangle_xm, cosangle_xn, &
    sinangle_xn
  real(dp) :: dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
  real(dp) :: dy_norm3, dy_norm1, dx_norm2, dx_norm3, dz_norm1, dz_norm2, this_h
  integer :: iomega, indexl_coil
  real(dp), dimension(:,:,:), allocatable :: dinductancedomega

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0
  real(dp), dimension(:,:), allocatable :: tempMatrix

  call system_clock(tic,countrate)
  if (verbose) print *,"Initializing basis functions and f"
  
  ! Initialize Fourier arrays
  call regcoil_init_Fourier_modes(mpol_potential, ntor_potential, mnmax_potential, xm_potential, xn_potential, .false., helicity_ratio)
  xn_potential = xn_potential * nfp
  
  select case (symmetry_option)
  case (1,2)
     num_basis_functions = mnmax_potential
  case (3)
     num_basis_functions = mnmax_potential * 2
  case default
     print *,"Error! Invalid setting for symmetry_option:",symmetry_option
     stop
  end select

  if (allocated(basis_functions)) deallocate(basis_functions)
  allocate(basis_functions(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 1!'

  if (allocated(f_x)) deallocate(f_x)
  allocate(f_x(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 2!'

  if (allocated(f_y)) deallocate(f_y)
  allocate(f_y(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 3!'

  if (allocated(f_z)) deallocate(f_z)
  allocate(f_z(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 4!'

  if (allocated(f_Laplace_Beltrami)) deallocate(f_Laplace_Beltrami)
  allocate(f_Laplace_Beltrami(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 4!'

  if (allocated(d_x)) deallocate(d_x)
  allocate(d_x(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 5!'

  if (allocated(d_y)) deallocate(d_y)
  allocate(d_y(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 6!'

  if (allocated(d_z)) deallocate(d_z)
  allocate(d_z(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 7!'

  if (allocated(d_Laplace_Beltrami)) deallocate(d_Laplace_Beltrami)
  allocate(d_Laplace_Beltrami(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 7!'


  allocate(g_theta_theta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(g_theta_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(g_zeta_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_theta_theta_d_theta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_theta_theta_d_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_theta_zeta_d_theta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_theta_zeta_d_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_zeta_zeta_d_theta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_g_zeta_zeta_d_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_N_d_theta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(d_N_d_zeta(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(Laplace_Beltrami_d_Phi_d_theta_coefficient(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  allocate(Laplace_Beltrami_d_Phi_d_zeta_coefficient(ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error!'

  if (sensitivity_option > 1) then
      allocate(dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dfydomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dfzdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dmatrix_Kdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dmatrix_Bdomega(nomega_coil,num_basis_functions,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dRHS_Kdomega(nomega_coil,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dRHS_Bdomega(nomega_coil,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
    endif
    if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
      allocate(f_xdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(f_ydNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(f_zdNdomega_over_N_coil2(nomega_coil,ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
    endif

  ! For the Laplace-Beltrami operator, compute the metric coefficients and their derivatives on the coil surface.
  ! See my notes from 20180516 for the derivation of the necessary quantities.

  g_theta_theta = drdtheta_coil(1,:,1:nzeta_coil) * drdtheta_coil(1,:,1:nzeta_coil) &
       +          drdtheta_coil(2,:,1:nzeta_coil) * drdtheta_coil(2,:,1:nzeta_coil) &
       +          drdtheta_coil(3,:,1:nzeta_coil) * drdtheta_coil(3,:,1:nzeta_coil)

  g_theta_zeta  = drdtheta_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +          drdtheta_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +          drdtheta_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil)

  g_zeta_zeta   = drdzeta_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +          drdzeta_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +          drdzeta_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil)

  d_g_theta_theta_d_theta = (d2rdtheta2_coil(1,:,1:nzeta_coil) * drdtheta_coil(1,:,1:nzeta_coil) &
       +                     d2rdtheta2_coil(2,:,1:nzeta_coil) * drdtheta_coil(2,:,1:nzeta_coil) &
       +                     d2rdtheta2_coil(3,:,1:nzeta_coil) * drdtheta_coil(3,:,1:nzeta_coil)) * 2

  d_g_theta_theta_d_zeta  = (d2rdthetadzeta_coil(1,:,1:nzeta_coil) * drdtheta_coil(1,:,1:nzeta_coil) &
       +                     d2rdthetadzeta_coil(2,:,1:nzeta_coil) * drdtheta_coil(2,:,1:nzeta_coil) &
       +                     d2rdthetadzeta_coil(3,:,1:nzeta_coil) * drdtheta_coil(3,:,1:nzeta_coil)) * 2

  d_g_theta_zeta_d_theta = d2rdtheta2_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +                   d2rdtheta2_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +                   d2rdtheta2_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil) &
       +                   drdtheta_coil(1,:,1:nzeta_coil) * d2rdthetadzeta_coil(1,:,1:nzeta_coil) &
       +                   drdtheta_coil(2,:,1:nzeta_coil) * d2rdthetadzeta_coil(2,:,1:nzeta_coil) &
       +                   drdtheta_coil(3,:,1:nzeta_coil) * d2rdthetadzeta_coil(3,:,1:nzeta_coil)

  d_g_theta_zeta_d_zeta = d2rdthetadzeta_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +                  d2rdthetadzeta_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +                  d2rdthetadzeta_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil) &
       +                  drdtheta_coil(1,:,1:nzeta_coil) * d2rdzeta2_coil(1,:,1:nzeta_coil) &
       +                  drdtheta_coil(2,:,1:nzeta_coil) * d2rdzeta2_coil(2,:,1:nzeta_coil) &
       +                  drdtheta_coil(3,:,1:nzeta_coil) * d2rdzeta2_coil(3,:,1:nzeta_coil)

  d_g_zeta_zeta_d_theta = (d2rdthetadzeta_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +                   d2rdthetadzeta_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +                   d2rdthetadzeta_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil)) * 2

  d_g_zeta_zeta_d_zeta = (d2rdzeta2_coil(1,:,1:nzeta_coil) * drdzeta_coil(1,:,1:nzeta_coil) &
       +                  d2rdzeta2_coil(2,:,1:nzeta_coil) * drdzeta_coil(2,:,1:nzeta_coil) &
       +                  d2rdzeta2_coil(3,:,1:nzeta_coil) * drdzeta_coil(3,:,1:nzeta_coil)) * 2

  d_N_d_theta = (d_g_zeta_zeta_d_theta * g_theta_theta + d_g_theta_theta_d_theta * g_zeta_zeta - 2 * d_g_theta_zeta_d_theta * g_theta_zeta) / (2 * norm_normal_coil)

  d_N_d_zeta  = (d_g_zeta_zeta_d_zeta  * g_theta_theta + d_g_theta_theta_d_zeta  * g_zeta_zeta - 2 * d_g_theta_zeta_d_zeta  * g_theta_zeta) / (2 * norm_normal_coil)

  Laplace_Beltrami_d_Phi_d_theta_coefficient = (d_g_zeta_zeta_d_theta - d_g_theta_zeta_d_zeta &
       + (-g_zeta_zeta * d_N_d_theta + g_theta_zeta * d_N_d_zeta) / norm_normal_coil) / norm_normal_coil

  Laplace_Beltrami_d_Phi_d_zeta_coefficient  = (d_g_theta_theta_d_zeta - d_g_theta_zeta_d_theta &
       + (g_theta_zeta * d_N_d_theta - g_theta_theta * d_N_d_zeta) / norm_normal_coil) / norm_normal_coil

  d_x = reshape((net_poloidal_current_Amperes * drdtheta_coil(1,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(1,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  d_y = reshape((net_poloidal_current_Amperes * drdtheta_coil(2,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(2,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  d_z = reshape((net_poloidal_current_Amperes * drdtheta_coil(3,:,1:nzeta_coil) - net_toroidal_current_Amperes * drdzeta_coil(3,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))

  d_Laplace_Beltrami = -reshape((net_poloidal_current_Amperes / twopi) * Laplace_Beltrami_d_Phi_d_zeta_coefficient &
       + (net_toroidal_current_Amperes / twopi) * Laplace_Beltrami_d_Phi_d_theta_coefficient, (/ ntheta_coil*nzeta_coil /))

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
        offset = mnmax_potential
     else
        offset = 0
     end if

     do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
           index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
           do imn = 1, mnmax_potential

              angle = xm_potential(imn)*theta_coil(itheta_coil)-xn_potential(imn)*zeta_coil(izeta_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)
              cosangle_xn = cosangle*xn_potential(imn)
              sinangle_xn = sinangle*xn_potential(imn)
              cosangle_xm = cosangle*xm_potential(imn)
              sinangle_xm = sinangle*xm_potential(imn)
              if (whichSymmetry==1) then
                 basis_functions(index_coil, imn) = sinangle
                 f_x(index_coil, imn) = cosangle*(xn_potential(imn)*drdtheta_coil(1,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(1,itheta_coil,izeta_coil))
                 f_y(index_coil, imn) = cosangle*(xn_potential(imn)*drdtheta_coil(2,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(2,itheta_coil,izeta_coil))
                 f_z(index_coil, imn) = cosangle*(xn_potential(imn)*drdtheta_coil(3,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(3,itheta_coil,izeta_coil))
                 f_Laplace_Beltrami(index_coil, imn) = ( &
                      xm_potential(imn) * Laplace_Beltrami_d_Phi_d_theta_coefficient(itheta_coil, izeta_coil) &
                      - xn_potential(imn) * Laplace_Beltrami_d_Phi_d_zeta_coefficient(itheta_coil, izeta_coil)) * cosangle &
                      + (   xm_potential(imn) * xm_potential(imn) * g_zeta_zeta(  itheta_coil, izeta_coil) &
                      +     xn_potential(imn) * xn_potential(imn) * g_theta_theta(itheta_coil, izeta_coil) &
                      + 2 * xm_potential(imn) * xn_potential(imn) * g_theta_zeta( itheta_coil, izeta_coil) ) * (-sinangle) / norm_normal_coil(itheta_coil, izeta_coil)
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
                 f_x(index_coil, imn+offset) = -sinangle*(xn_potential(imn)*drdtheta_coil(1,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(1,itheta_coil,izeta_coil))
                 f_y(index_coil, imn+offset) = -sinangle*(xn_potential(imn)*drdtheta_coil(2,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(2,itheta_coil,izeta_coil))
                 f_z(index_coil, imn+offset) = -sinangle*(xn_potential(imn)*drdtheta_coil(3,itheta_coil,izeta_coil) + xm_potential(imn)*drdzeta_coil(3,itheta_coil,izeta_coil))
                 f_Laplace_Beltrami(index_coil, imn + offset) = ( &
                      xm_potential(imn) * Laplace_Beltrami_d_Phi_d_theta_coefficient(itheta_coil, izeta_coil) &
                      - xn_potential(imn) * Laplace_Beltrami_d_Phi_d_zeta_coefficient(itheta_coil, izeta_coil)) * (-sinangle) &
                      + (   xm_potential(imn) * xm_potential(imn) * g_zeta_zeta(  itheta_coil, izeta_coil) &
                      +     xn_potential(imn) * xn_potential(imn) * g_theta_theta(itheta_coil, izeta_coil) &
                      + 2 * xm_potential(imn) * xn_potential(imn) * g_theta_zeta( itheta_coil, izeta_coil) ) * (-cosangle) / norm_normal_coil(itheta_coil, izeta_coil)
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
              end if
           end do
        end do
     end do
  end do
  
  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."
  

  if (allocated(g)) deallocate(g)
  allocate(g(ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 8!'

  if (allocated(inductance)) deallocate(inductance)
  allocate(inductance(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 9!'

  if (allocated(h)) deallocate(h)
  allocate(h(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 10!'

  allocate(factor_for_h(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 11!'

  if (allocated(Bnormal_from_net_coil_currents)) deallocate(Bnormal_from_net_coil_currents)
  allocate(Bnormal_from_net_coil_currents(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 12!'

  if (sensitivity_option > 1) then
    allocate(dgdomega(ntheta_plasma*nzeta_plasma,num_basis_functions,nomega_coil),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancednorm(3),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancedr(3),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(dinductancedomega(ntheta_plasma*nzeta_plasma, &
      ntheta_coil*nzeta_coil,nomega_coil),stat=iflag)
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
  if (verbose) print *,"Building inductance matrix and h."
  !$OMP PARALLEL

  !$OMP MASTER
  if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
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

                ! Pre-multiplying these factors
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
                      (normal_plasma(1,itheta_plasma,izeta_plasma)*normal_coil(1,itheta_coil,izetal_coil) &
                      +normal_plasma(2,itheta_plasma,izeta_plasma)*normal_coil(2,itheta_coil,izetal_coil) &
                      +normal_plasma(3,itheta_plasma,izeta_plasma)*normal_coil(3,itheta_coil,izetal_coil) &
                      - (3*dr2inv) * &
                      (normal_plasma(1,itheta_plasma,izeta_plasma)*dx &
                      + normal_plasma(2,itheta_plasma,izeta_plasma)*dy &
                      + normal_plasma(3,itheta_plasma,izeta_plasma)*dz) * &
                      (normal_coil(1,itheta_coil,izetal_coil)*dx &
                      +normal_coil(2,itheta_coil,izetal_coil)*dy &
                      +normal_coil(3,itheta_coil,izetal_coil)*dz)) * dr32inv
                 
                 h(index_plasma) = h(index_plasma) + this_h

                if (sensitivity_option > 1) then
                  indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil
                  dr52inv = dr2inv*dr32inv

                 dinductancedomega(index_plasma,index_coil,:) = dinductancedomega(index_plasma,index_coil,:) &
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

  ! Deallocate quantities no longer needed
  if (allocated(drdomega)) deallocate(drdomega)

  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."
  
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
  ! If you like, you can use "matmul" instead which is slower but more reliable.

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

  if (sensitivity_option > 1) then
   !$OMP PARALLEL
   !$OMP MASTER
   if (verbose) then
     print *,"  Number of OpenMP threads:",omp_get_num_threads()
   end if
   !$OMP END MASTER
   !$OMP DO
     do iomega = 1, nomega_coil
      call &
        DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedomega(:,:,iomega),&
        LDA,basis_functions,LDB,BLAS_BETA,dgdomega(:,:,iomega),LDC)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end if

  call system_clock(toc)
  if (verbose) print *,"inductance*basis_functions:",real(toc-tic)/countrate,"sec."

  if (allocated(matrix_B)) deallocate(matrix_B)
  allocate(matrix_B(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 13!'

  if (allocated(matrix_regularization)) deallocate(matrix_regularization)
  allocate(matrix_regularization(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 14!'

  if (allocated(RHS_B)) deallocate(RHS_B)
  allocate(RHS_B(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 15!'

  if (allocated(RHS_regularization)) deallocate(RHS_regularization)
  allocate(RHS_regularization(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 16!'

  ! if (allocated(norm_normal_plasma_inv1D)) deallocate(norm_normal_plasma_inv1D)
  allocate(norm_normal_plasma_inv1D(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 17!'

  ! if (allocated(norm_normal_coil_inv1D)) deallocate(norm_normal_coil_inv1D)
  allocate(norm_normal_coil_inv1D(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 18!'

  ! if (allocated(g_over_N_plasma)) deallocate(g_over_N_plasma)
  allocate(g_over_N_plasma(ntheta_plasma*nzeta_plasma,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 19!'

  ! if (allocated(f_x_over_N_coil)) deallocate(f_x_over_N_coil)
  allocate(f_x_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 20!'

  ! if (allocated(f_y_over_N_coil)) deallocate(f_y_over_N_coil)
  allocate(f_y_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 21!'

  ! if (allocated(f_z_over_N_coil)) deallocate(f_z_over_N_coil)
  allocate(f_z_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 22!'

  allocate(f_Laplace_Beltrami_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 22!'

  call system_clock(tic)

  RHS_B = -dtheta_plasma*dzeta_plasma*matmul( &
       reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)), g)

  call system_clock(toc)
  if (verbose) print *,"Form RHS_B:",real(toc-tic)/countrate,"sec."
  call system_clock(tic)

  norm_normal_plasma_inv1D = reshape(1/norm_normal_plasma, (/ ntheta_plasma*nzeta_plasma /))
  norm_normal_coil_inv1D   = reshape(1/norm_normal_coil,   (/ ntheta_coil  *nzeta_coil /))
  do j = 1,num_basis_functions
     g_over_N_plasma(:,j) = g(:,j) * norm_normal_plasma_inv1D
     f_x_over_N_coil(:,j) = f_x(:,j) * norm_normal_coil_inv1D
     f_y_over_N_coil(:,j) = f_y(:,j) * norm_normal_coil_inv1D
     f_z_over_N_coil(:,j) = f_z(:,j) * norm_normal_coil_inv1D
     f_Laplace_Beltrami_over_N_coil(:,j) = f_Laplace_Beltrami(:,j) * norm_normal_coil_inv1D
      if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
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

  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    call system_clock(tic)

    do iomega = 1, nomega_coil
      dRHS_Bdomega(iomega,:) = -dtheta_plasma*dzeta_plasma*matmul( &
      reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)), dgdomega(:,:,iomega))
      dRHS_Bdomega(iomega,:) = dRHS_Bdomega(iomega,:) - &
        dtheta_plasma*dzeta_plasma*matmul(transpose(g_over_N_plasma),dhdomega(iomega,:))
    enddo

    call system_clock(toc)
    if (verbose) then
      print *,"Form dRHS_Bdomega: ",real(toc-tic)/countrate,"sec."
    end if
  endif


  matrix_B = 0
  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    dmatrix_Bdomega = 0
  end if

  deallocate(norm_normal_plasma_inv1D)
  deallocate(norm_normal_coil_inv1D)

  call system_clock(toc)
  if (verbose) print *,"Prepare for matrix_B:",real(toc-tic)/countrate,"sec."
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
  if (verbose) print *,"matmul for matrix_B:",real(toc-tic)/countrate,"sec."

  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    BLAS_BETA=1
    call system_clock(tic)
    !$OMP PARALLEL
    !$OMP MASTER
    if (verbose) then
      print *,"  Number of OpenMP threads:",omp_get_num_threads()
    end if
    !$OMP END MASTER
    !$OMP DO
    do iomega = 1, nomega_coil
      call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dgdomega(:,:,iomega),LDA,g_over_N_plasma,LDB,BLAS_BETA,dmatrix_Bdomega(iomega,:,:),LDC)
      call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,g_over_N_plasma,LDA,dgdomega(:,:,iomega),LDB,BLAS_BETA,dmatrix_Bdomega(iomega,:,:),LDC)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call system_clock(toc)
    if (verbose) then
      print *,"matmul for dmatrix_Bdomega:",real(toc-tic)/countrate,"sec."
    end if
  endif

  deallocate(g_over_N_plasma)
    

  matrix_regularization = 0
  if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
    dmatrix_Kdomega = 0
  end if
     
  select case (trim(regularization_term_option))
  case (regularization_term_option_chi2_K, regularization_term_option_K_xy)
  
     call system_clock(tic)
     ! Here we carry out matrix_regularization += dtheta*dzeta*(f_x ^ T) * f_x_over_N_coil
     ! A = f_x
     ! B = f_x_over_N_plasma
     ! C = matrix_regularization
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_x,LDA,f_x_over_N_coil,LDB,BLAS_BETA,matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul 1 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     call system_clock(tic)
     ! Here we carry out matrix_regularization += dtheta*dzeta*(f_y ^ T) * f_y_over_N_coil
     ! A = f_y
     ! B = f_y_over_N_plasma
     ! C = matrix_regularization
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_y,LDA,f_y_over_N_coil,LDB,BLAS_BETA,matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul 2 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     if (trim(regularization_term_option) == regularization_term_option_chi2_K) then
        call system_clock(tic)
        ! Here we carry out matrix_regularization += dtheta*dzeta*(f_z ^ T) * f_z_over_N_coil
        ! A = f_z
        ! B = f_z_over_N_plasma
        ! C = matrix_regularization
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
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_z,LDA,f_z_over_N_coil,LDB,BLAS_BETA,matrix_regularization,LDC)
        
        call system_clock(toc)
        if (verbose) print *,"matmul 3 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     end if
     
     call system_clock(tic)
     
     if (trim(regularization_term_option) == regularization_term_option_chi2_K) then     
        RHS_regularization = (matmul(d_x, f_x_over_N_coil) + matmul(d_y, f_y_over_N_coil) + matmul(d_z, f_z_over_N_coil)) &
             * (dtheta_coil*dzeta_coil)
     else
        RHS_regularization = (matmul(d_x, f_x_over_N_coil) + matmul(d_y, f_y_over_N_coil)) &
             * (dtheta_coil*dzeta_coil)
     end if
     
     call system_clock(toc)
     if (verbose) print *,"Matmuls for RHS_regularization:",real(toc-tic)/countrate,"sec."

    ! Compute dRHS_Kdomega
    if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
      call system_clock(tic)
      !$OMP PARALLEL
      !$OMP MASTER
      if (verbose) then
        print *,"  Number of OpenMP threads:",omp_get_num_threads()
      end if
      !$OMP END MASTER
      !$OMP DO
      ! This would be faster with LAPACK, but constructing dinductancematrixdomega is much more expensive
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
      !$OMP END DO
      !$OMP END PARALLEL
      call system_clock(toc)
      if (verbose) then
        print *,"Matmuls for dRHS_Kdomega:",real(toc-tic)/countrate,"sec."
      end if
    endif

    ! Construct dmatrix_Kdomega
    if (sensitivity_option > 2 .or. fixed_norm_sensitivity_option) then
      call system_clock(tic)
      !$OMP PARALLEL
      !$OMP MASTER
      if (verbose) then
        print *,"  Number of OpenMP threads:",omp_get_num_threads()
      end if
      !$OMP END MASTER
      !$OMP DO
      do iomega = 1, nomega_coil
        ! Converted matmul to dot_product to avoid segfault issues for large matrices
        do j = 1,num_basis_functions
      ! dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions
      ! f_x_over_N_coil(ntheta_coil*nzeta_coil,num_basis_functions)
      ! fx(ntheta_coil*nzeta_coil, num_basis_functions)

          dmatrix_Kdomega(iomega,:,j) = matmul(transpose(dfxdomega(iomega,:,:)),f_x_over_N_coil(:,j)) &
          + matmul(transpose(dfydomega(iomega,:,:)),f_y_over_N_coil(:,j)) &
          + matmul(transpose(dfzdomega(iomega,:,:)),f_z_over_N_coil(:,j)) &
          + matmul(transpose(f_x_over_n_coil(:,:)),dfxdomega(iomega,:,j)) &
          + matmul(transpose(f_y_over_n_coil(:,:)),dfydomega(iomega,:,j)) &
          + matmul(transpose(f_z_over_n_coil(:,:)),dfzdomega(iomega,:,j)) &
          + matmul(transpose(f_x(:,:)),f_xdNdomega_over_N_coil2(iomega,:,j)) &
          + matmul(transpose(f_y(:,:)),f_ydNdomega_over_N_coil2(iomega,:,j)) &
          + matmul(transpose(f_z(:,:)),f_zdNdomega_over_N_coil2(iomega,:,j))

!          dmatrix_Kdomega(iomega,:,j) = 2*(dot_product(dfxdomega(iomega,:,:),f_x_over_N_coil(:,j)) &
!          + dot_product(dfydomega(iomega,:,j),f_y_over_N_coil(:,j)) &
!          + dot_product(dfzdomega(iomega,:,j),f_z_over_N_coil(:,j))) &
!          + dot_product(f_x(:,j),f_xdNdomega_over_N_coil2(iomega,:,j)) &
!          + dot_product(f_y(:,j),f_ydNdomega_over_N_coil2(iomega,:,j)) &
!          + dot_product(f_z(:,j),f_zdNdomega_over_N_coil2(iomega,:,j))
        end do
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dfxdomega(iomega,:,:),LDA,&
!          f_x_over_N_coil,LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dfydomega(iomega,:,:),LDA,&
!          f_y_over_N_coil,LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dfzdomega(iomega,:,:),LDA,&
!          f_z_over_N_coil,LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_x_over_N_coil,LDA,&
!          dfxdomega(iomega,:,:),LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_y_over_N_coil,LDA,&
!          dfydomega(iomega,:,:),LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_z_over_N_coil,LDA,&
!          dfzdomega(iomega,:,:),LDB,BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,-BLAS_ALPHA,f_x,LDA,f_xdNdomega_over_N_coil2(iomega,:,:),LDB,&
!          BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,-BLAS_ALPHA,f_y,LDA,f_ydNdomega_over_N_coil2(iomega,:,:),LDB,&
!          BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
!        call DGEMM(TRANSA,TRANSB,M,N,K,-BLAS_ALPHA,f_z,LDA,f_zdNdomega_over_N_coil2(iomega,:,:),LDB,&
!          BLAS_BETA,dmatrix_Kdomega(iomega,:,:),LDC)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      call system_clock(toc)
      if (verbose) then
        print *,"matmul for dmatrix_Kdomega in regcoil_build_matrices :",real(toc-tic)/countrate,"sec."
      end if
      dmatrix_Kdomega = dmatrix_Kdomega*dtheta_coil*dzeta_coil
    endif

  case (regularization_term_option_Laplace_Beltrami)
     ! ------------------------------------------------------------------
     ! Laplace-Beltrami matrix and RHS:
     
     call system_clock(tic)
     ! Here we carry out matrix = dtheta*dzeta*(f ^ T) * f_over_N_coil
     ! A = f
     ! B = f_over_N_plasma
     ! C = matrix
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,f_Laplace_Beltrami,LDA,f_Laplace_Beltrami_over_N_coil,LDB,BLAS_BETA,matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     call system_clock(tic)
     
     RHS_regularization = matmul(d_Laplace_Beltrami, f_Laplace_Beltrami_over_N_coil) * (dtheta_coil*dzeta_coil)
     
     call system_clock(toc)
     if (verbose) print *,"Matmul for RHS_regularization:",real(toc-tic)/countrate,"sec."

  case default
     print *,"Error! Unrecognized regularization_term_option: ",regularization_term_option
  end select

  ! ------------------------------------------------------------------

  deallocate(f_x_over_N_coil)
  deallocate(f_y_over_N_coil)
  deallocate(f_z_over_N_coil)
  deallocate(f_Laplace_Beltrami_over_N_coil)

  deallocate(g_theta_theta, g_theta_zeta, g_zeta_zeta)
  deallocate(d_g_theta_theta_d_theta, d_g_theta_theta_d_zeta)
  deallocate(d_g_theta_zeta_d_theta, d_g_theta_zeta_d_zeta)
  deallocate(d_g_zeta_zeta_d_theta, d_g_zeta_zeta_d_zeta)
  deallocate(d_N_d_theta, d_N_d_zeta)
  deallocate(Laplace_Beltrami_d_Phi_d_theta_coefficient, Laplace_Beltrami_d_Phi_d_zeta_coefficient)

end subroutine regcoil_build_matrices

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



