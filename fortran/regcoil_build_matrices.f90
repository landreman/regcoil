subroutine regcoil_build_matrices(prob)

  use regcoil_variables, only: regcoil_t, regularization_term_option_chi2_K, regularization_term_option_Laplace_Beltrami, regularization_term_option_K_xy
  use stel_constants
  use stel_kinds
  use omp_lib
  use regcoil_init_Fourier_modes_mod
  
  implicit none


  type(regcoil_t), intent(inout) :: prob
  integer :: l_coil, itheta_plasma, izeta_plasma, itheta_coil, izeta_coil, izetal_coil
  real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv
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
  real(dp) :: dy_norm3, dy_norm1, dx_norm2, dx_norm3, dz_norm1, dz_norm2, this_h

  ! Variables needed by BLAS DGEMM:
  character :: TRANSA, TRANSB
  integer :: M, N, K, LDA, LDB, LDC
  real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0
  real(dp), dimension(:,:), allocatable :: tempMatrix

    associate ( &
       ntheta_plasma => prob%plasma%ntheta_plasma, &
       nzeta_plasma => prob%plasma%nzeta_plasma, &
       dtheta_plasma => prob%plasma%dtheta_plasma, &
       dzeta_plasma => prob%plasma%dzeta_plasma, &
       nfp => prob%plasma%nfp, &
       ntheta_coil => prob%coil%ntheta_coil, &
       nzeta_coil => prob%coil%nzeta_coil, &
       nzetal_coil => prob%coil%nzetal_coil, &
       dtheta_coil => prob%coil%dtheta_coil, &
       dzeta_coil => prob%coil%dzeta_coil, &
       verbose => prob%input%verbose, &
       regularization_term_option => prob%input%regularization_term_option, &
       symmetry_option => prob%input%symmetry_option, &
       mpol_potential => prob%input%mpol_potential, &
       ntor_potential => prob%input%ntor_potential, &
       net_poloidal_current_Amperes => prob%input%net_poloidal_current_Amperes, &
       net_toroidal_current_Amperes => prob%input%net_toroidal_current_Amperes, &
       mnmax_potential => prob%work%mnmax_potential, &
       num_basis_functions => prob%work%num_basis_functions &
       )
  call system_clock(tic,countrate)
  if (verbose) print *,"Initializing basis functions and f"
  
  ! Initialize Fourier arrays
  call regcoil_init_Fourier_modes(mpol_potential, ntor_potential, prob%work%mnmax_potential, prob%work%xm_potential, prob%work%xn_potential, .false.)
  prob%work%xn_potential = prob%work%xn_potential * nfp
  
  select case (symmetry_option)
  case (1,2)
     num_basis_functions = mnmax_potential
  case (3)
     num_basis_functions = mnmax_potential * 2
  case default
     print *,"Error! Invalid setting for symmetry_option:",symmetry_option
     stop
  end select

  if (allocated(prob%work%basis_functions)) deallocate(prob%work%basis_functions)
  allocate(prob%work%basis_functions(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 1!'

  if (allocated(prob%work%f_x)) deallocate(prob%work%f_x)
  allocate(prob%work%f_x(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 2!'

  if (allocated(prob%work%f_y)) deallocate(prob%work%f_y)
  allocate(prob%work%f_y(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 3!'

  if (allocated(prob%work%f_z)) deallocate(prob%work%f_z)
  allocate(prob%work%f_z(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 4!'

  if (allocated(prob%work%f_Laplace_Beltrami)) deallocate(prob%work%f_Laplace_Beltrami)
  allocate(prob%work%f_Laplace_Beltrami(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 4!'

  if (allocated(prob%work%d_x)) deallocate(prob%work%d_x)
  allocate(prob%work%d_x(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 5!'

  if (allocated(prob%work%d_y)) deallocate(prob%work%d_y)
  allocate(prob%work%d_y(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 6!'

  if (allocated(prob%work%d_z)) deallocate(prob%work%d_z)
  allocate(prob%work%d_z(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 7!'

  if (allocated(prob%work%d_Laplace_Beltrami)) deallocate(prob%work%d_Laplace_Beltrami)
  allocate(prob%work%d_Laplace_Beltrami(ntheta_coil*nzeta_coil),stat=iflag)
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


  ! For the Laplace-Beltrami operator, compute the metric coefficients and their derivatives on the coil surface.
  ! See my notes from 20180516 for the derivation of the necessary quantities.

  g_theta_theta = prob%coil%drdtheta_coil(1,:,1:nzeta_coil) * prob%coil%drdtheta_coil(1,:,1:nzeta_coil) &
       +          prob%coil%drdtheta_coil(2,:,1:nzeta_coil) * prob%coil%drdtheta_coil(2,:,1:nzeta_coil) &
       +          prob%coil%drdtheta_coil(3,:,1:nzeta_coil) * prob%coil%drdtheta_coil(3,:,1:nzeta_coil)

  g_theta_zeta  = prob%coil%drdtheta_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +          prob%coil%drdtheta_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +          prob%coil%drdtheta_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil)

  g_zeta_zeta   = prob%coil%drdzeta_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +          prob%coil%drdzeta_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +          prob%coil%drdzeta_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil)

  d_g_theta_theta_d_theta = (prob%coil%d2rdtheta2_coil(1,:,1:nzeta_coil) * prob%coil%drdtheta_coil(1,:,1:nzeta_coil) &
       +                     prob%coil%d2rdtheta2_coil(2,:,1:nzeta_coil) * prob%coil%drdtheta_coil(2,:,1:nzeta_coil) &
       +                     prob%coil%d2rdtheta2_coil(3,:,1:nzeta_coil) * prob%coil%drdtheta_coil(3,:,1:nzeta_coil)) * 2

  d_g_theta_theta_d_zeta  = (prob%coil%d2rdthetadzeta_coil(1,:,1:nzeta_coil) * prob%coil%drdtheta_coil(1,:,1:nzeta_coil) &
       +                     prob%coil%d2rdthetadzeta_coil(2,:,1:nzeta_coil) * prob%coil%drdtheta_coil(2,:,1:nzeta_coil) &
       +                     prob%coil%d2rdthetadzeta_coil(3,:,1:nzeta_coil) * prob%coil%drdtheta_coil(3,:,1:nzeta_coil)) * 2

  d_g_theta_zeta_d_theta = prob%coil%d2rdtheta2_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +                   prob%coil%d2rdtheta2_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +                   prob%coil%d2rdtheta2_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil) &
       +                   prob%coil%drdtheta_coil(1,:,1:nzeta_coil) * prob%coil%d2rdthetadzeta_coil(1,:,1:nzeta_coil) &
       +                   prob%coil%drdtheta_coil(2,:,1:nzeta_coil) * prob%coil%d2rdthetadzeta_coil(2,:,1:nzeta_coil) &
       +                   prob%coil%drdtheta_coil(3,:,1:nzeta_coil) * prob%coil%d2rdthetadzeta_coil(3,:,1:nzeta_coil)

  d_g_theta_zeta_d_zeta = prob%coil%d2rdthetadzeta_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +                  prob%coil%d2rdthetadzeta_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +                  prob%coil%d2rdthetadzeta_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil) &
       +                  prob%coil%drdtheta_coil(1,:,1:nzeta_coil) * prob%coil%d2rdzeta2_coil(1,:,1:nzeta_coil) &
       +                  prob%coil%drdtheta_coil(2,:,1:nzeta_coil) * prob%coil%d2rdzeta2_coil(2,:,1:nzeta_coil) &
       +                  prob%coil%drdtheta_coil(3,:,1:nzeta_coil) * prob%coil%d2rdzeta2_coil(3,:,1:nzeta_coil)

  d_g_zeta_zeta_d_theta = (prob%coil%d2rdthetadzeta_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +                   prob%coil%d2rdthetadzeta_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +                   prob%coil%d2rdthetadzeta_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil)) * 2

  d_g_zeta_zeta_d_zeta = (prob%coil%d2rdzeta2_coil(1,:,1:nzeta_coil) * prob%coil%drdzeta_coil(1,:,1:nzeta_coil) &
       +                  prob%coil%d2rdzeta2_coil(2,:,1:nzeta_coil) * prob%coil%drdzeta_coil(2,:,1:nzeta_coil) &
       +                  prob%coil%d2rdzeta2_coil(3,:,1:nzeta_coil) * prob%coil%drdzeta_coil(3,:,1:nzeta_coil)) * 2

  d_N_d_theta = (d_g_zeta_zeta_d_theta * g_theta_theta + d_g_theta_theta_d_theta * g_zeta_zeta - 2 * d_g_theta_zeta_d_theta * g_theta_zeta) / (2 * prob%coil%norm_normal_coil)

  d_N_d_zeta  = (d_g_zeta_zeta_d_zeta  * g_theta_theta + d_g_theta_theta_d_zeta  * g_zeta_zeta - 2 * d_g_theta_zeta_d_zeta  * g_theta_zeta) / (2 * prob%coil%norm_normal_coil)

  Laplace_Beltrami_d_Phi_d_theta_coefficient = (d_g_zeta_zeta_d_theta - d_g_theta_zeta_d_zeta &
       + (-g_zeta_zeta * d_N_d_theta + g_theta_zeta * d_N_d_zeta) / prob%coil%norm_normal_coil) / prob%coil%norm_normal_coil

  Laplace_Beltrami_d_Phi_d_zeta_coefficient  = (d_g_theta_theta_d_zeta - d_g_theta_zeta_d_theta &
       + (g_theta_zeta * d_N_d_theta - g_theta_theta * d_N_d_zeta) / prob%coil%norm_normal_coil) / prob%coil%norm_normal_coil

  prob%work%d_x = reshape((net_poloidal_current_Amperes * prob%coil%drdtheta_coil(1,:,1:nzeta_coil) - net_toroidal_current_Amperes * prob%coil%drdzeta_coil(1,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  prob%work%d_y = reshape((net_poloidal_current_Amperes * prob%coil%drdtheta_coil(2,:,1:nzeta_coil) - net_toroidal_current_Amperes * prob%coil%drdzeta_coil(2,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))
  prob%work%d_z = reshape((net_poloidal_current_Amperes * prob%coil%drdtheta_coil(3,:,1:nzeta_coil) - net_toroidal_current_Amperes * prob%coil%drdzeta_coil(3,:,1:nzeta_coil)) / twopi, &
       (/ ntheta_coil*nzeta_coil /))

  prob%work%d_Laplace_Beltrami = -reshape((net_poloidal_current_Amperes / twopi) * Laplace_Beltrami_d_Phi_d_zeta_coefficient &
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

              angle = prob%work%xm_potential(imn)*prob%coil%theta_coil(itheta_coil)-prob%work%xn_potential(imn)*prob%coil%zeta_coil(izeta_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)
              if (whichSymmetry==1) then
                 prob%work%basis_functions(index_coil, imn) = sinangle
                 prob%work%f_x(index_coil, imn) = cosangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(1,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(1,itheta_coil,izeta_coil))
                 prob%work%f_y(index_coil, imn) = cosangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(2,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(2,itheta_coil,izeta_coil))
                 prob%work%f_z(index_coil, imn) = cosangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(3,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(3,itheta_coil,izeta_coil))
                 prob%work%f_Laplace_Beltrami(index_coil, imn) = ( &
                      prob%work%xm_potential(imn) * Laplace_Beltrami_d_Phi_d_theta_coefficient(itheta_coil, izeta_coil) &
                      - prob%work%xn_potential(imn) * Laplace_Beltrami_d_Phi_d_zeta_coefficient(itheta_coil, izeta_coil)) * cosangle &
                      + (   prob%work%xm_potential(imn) * prob%work%xm_potential(imn) * g_zeta_zeta(  itheta_coil, izeta_coil) &
                      +     prob%work%xn_potential(imn) * prob%work%xn_potential(imn) * g_theta_theta(itheta_coil, izeta_coil) &
                      + 2 * prob%work%xm_potential(imn) * prob%work%xn_potential(imn) * g_theta_zeta( itheta_coil, izeta_coil) ) * (-sinangle) / prob%coil%norm_normal_coil(itheta_coil, izeta_coil)
              else
                 prob%work%basis_functions(index_coil, imn+offset) = cosangle
                 prob%work%f_x(index_coil, imn+offset) = -sinangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(1,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(1,itheta_coil,izeta_coil))
                 prob%work%f_y(index_coil, imn+offset) = -sinangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(2,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(2,itheta_coil,izeta_coil))
                 prob%work%f_z(index_coil, imn+offset) = -sinangle*(prob%work%xn_potential(imn)*prob%coil%drdtheta_coil(3,itheta_coil,izeta_coil) + prob%work%xm_potential(imn)*prob%coil%drdzeta_coil(3,itheta_coil,izeta_coil))
                 prob%work%f_Laplace_Beltrami(index_coil, imn + offset) = ( &
                      prob%work%xm_potential(imn) * Laplace_Beltrami_d_Phi_d_theta_coefficient(itheta_coil, izeta_coil) &
                      - prob%work%xn_potential(imn) * Laplace_Beltrami_d_Phi_d_zeta_coefficient(itheta_coil, izeta_coil)) * (-sinangle) &
                      + (   prob%work%xm_potential(imn) * prob%work%xm_potential(imn) * g_zeta_zeta(  itheta_coil, izeta_coil) &
                      +     prob%work%xn_potential(imn) * prob%work%xn_potential(imn) * g_theta_theta(itheta_coil, izeta_coil) &
                      + 2 * prob%work%xm_potential(imn) * prob%work%xn_potential(imn) * g_theta_zeta( itheta_coil, izeta_coil) ) * (-cosangle) / prob%coil%norm_normal_coil(itheta_coil, izeta_coil)
              end if
           end do
        end do
     end do
  end do
  
  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."
  

  if (allocated(prob%work%g)) deallocate(prob%work%g)
  allocate(prob%work%g(ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 8!'

  if (allocated(prob%work%inductance)) deallocate(prob%work%inductance)
  allocate(prob%work%inductance(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 9!'

  if (allocated(prob%work%h)) deallocate(prob%work%h)
  allocate(prob%work%h(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 10!'

  allocate(factor_for_h(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 11!'

  if (allocated(prob%coil%Bnormal_from_net_coil_currents)) deallocate(prob%coil%Bnormal_from_net_coil_currents)
  allocate(prob%coil%Bnormal_from_net_coil_currents(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 12!'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now compute prob%work%g and prob%work%h
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  prob%work%inductance = 0
  prob%work%h=0
  factor_for_h = net_poloidal_current_Amperes * prob%coil%drdtheta_coil - net_toroidal_current_Amperes * prob%coil%drdzeta_coil

  call system_clock(tic,countrate)
  if (verbose) print *,"Building inductance matrix and h."
  !$OMP PARALLEL

  !$OMP MASTER
  if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER

  ! Note: the outermost loop below must be over the plasma variables rather than over the coil variables.
  ! This ensures the multiple threads write to different indices in prob%work%h() rather than to the same indices in prob%work%h(),
  ! in which case the prob%work%h(index+plasma)=prob%work%h(index_plasma)+... update does not work properly.
  !$OMP DO PRIVATE(index_plasma,index_coil,x,y,z,izetal_coil,dx,dy,dz,dr2inv,dr32inv,dx_norm2,dx_norm3,dy_norm1,dy_norm3,dz_norm1,dz_norm2,this_h)
  do izeta_plasma = 1, nzeta_plasma
     do itheta_plasma = 1, ntheta_plasma
        index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
        x = prob%plasma%r_plasma(1,itheta_plasma,izeta_plasma)
        y = prob%plasma%r_plasma(2,itheta_plasma,izeta_plasma)
        z = prob%plasma%r_plasma(3,itheta_plasma,izeta_plasma)
        do izeta_coil = 1, nzeta_coil
           do itheta_coil = 1, ntheta_coil
              index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
              do l_coil = 0, (nfp-1)
                 izetal_coil = izeta_coil + l_coil*nzeta_coil
                 dx = x - prob%coil%r_coil(1,itheta_coil,izetal_coil)
                 dy = y - prob%coil%r_coil(2,itheta_coil,izetal_coil)
                 dz = z - prob%coil%r_coil(3,itheta_coil,izetal_coil)
                 
                 dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                 dr32inv = dr2inv*sqrt(dr2inv)

                ! Pre-multiplying these factors
                 dy_norm3 = dy * prob%plasma%normal_plasma(3,itheta_plasma,izeta_plasma)
                 dz_norm1 = dz * prob%plasma%normal_plasma(1,itheta_plasma,izeta_plasma)
                 dx_norm2 = dx * prob%plasma%normal_plasma(2,itheta_plasma,izeta_plasma)
                 dy_norm1 = dy * prob%plasma%normal_plasma(1,itheta_plasma,izeta_plasma)
                 dz_norm2 = dz * prob%plasma%normal_plasma(2,itheta_plasma,izeta_plasma)
                 dx_norm3 = dx * prob%plasma%normal_plasma(3,itheta_plasma,izeta_plasma)
                 this_h = (factor_for_h(1,itheta_coil,izetal_coil) * dy_norm3 + &
                   factor_for_h(2,itheta_coil,izetal_coil) * dz_norm1 + &
                   factor_for_h(3,itheta_coil,izetal_coil) * dx_norm2  &
                   - factor_for_h(3,itheta_coil,izetal_coil) * dy_norm1 &
                   - factor_for_h(1,itheta_coil,izetal_coil) * dz_norm2 &
                   - factor_for_h(2,itheta_coil,izetal_coil) * dx_norm3 ) * dr32inv
                 
                 prob%work%inductance(index_plasma,index_coil) = prob%work%inductance(index_plasma,index_coil) + &
                      (prob%plasma%normal_plasma(1,itheta_plasma,izeta_plasma)*prob%coil%normal_coil(1,itheta_coil,izetal_coil) &
                      +prob%plasma%normal_plasma(2,itheta_plasma,izeta_plasma)*prob%coil%normal_coil(2,itheta_coil,izetal_coil) &
                      +prob%plasma%normal_plasma(3,itheta_plasma,izeta_plasma)*prob%coil%normal_coil(3,itheta_coil,izetal_coil) &
                      - (3*dr2inv) * &
                      (prob%plasma%normal_plasma(1,itheta_plasma,izeta_plasma)*dx &
                      + prob%plasma%normal_plasma(2,itheta_plasma,izeta_plasma)*dy &
                      + prob%plasma%normal_plasma(3,itheta_plasma,izeta_plasma)*dz) * &
                      (prob%coil%normal_coil(1,itheta_coil,izetal_coil)*dx &
                      +prob%coil%normal_coil(2,itheta_coil,izetal_coil)*dy &
                      +prob%coil%normal_coil(3,itheta_coil,izetal_coil)*dz)) * dr32inv
                 
                 prob%work%h(index_plasma) = prob%work%h(index_plasma) + this_h

              end do
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL


  call system_clock(toc)
  if (verbose) print *,"Done. Took",real(toc-tic)/countrate,"sec."
  
  prob%work%h = prob%work%h * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi))
  prob%work%inductance = prob%work%inductance * (mu0/(4*pi))
  deallocate(factor_for_h)
  prob%coil%Bnormal_from_net_coil_currents = reshape(prob%work%h, (/ ntheta_plasma, nzeta_plasma /)) / prob%plasma%norm_normal_plasma
  !prob%coil%Bnormal_from_net_coil_currents = transpose(reshape(prob%work%h, (/ nzeta_plasma, ntheta_plasma /))) / prob%plasma%norm_normal_plasma
  
  call system_clock(tic)

  ! For some reason, the BLAS prob%work%matrix-prob%work%matrix multiplication function DGEMM sometimes causes the
  ! program to crash on Edison unless you are careful to use the Intel MKL instead of Cray LibSci.
  ! If you like, you can use "matmul" instead which is slower but more reliable.

  !*******************************************************
  ! Call BLAS3 subroutine DGEMM for prob%work%matrix multiplications:
  !*******************************************************

  ! Here we carry out prob%work%g = prob%work%inductance * prob%work%basis_functions
  ! A = prob%work%inductance
  ! B = prob%work%basis_functions
  ! C = prob%work%g
  M = ntheta_plasma*nzeta_plasma ! # rows of A
  N = num_basis_functions ! # cols of B
  K = ntheta_coil*nzeta_coil ! Common dimension of A and B
  LDA = M
  LDB = K
  LDC = M
  TRANSA = 'N' ! No transposes
  TRANSB = 'N'
  prob%work%g = 0
  BLAS_ALPHA=dtheta_coil*dzeta_coil
  BLAS_BETA=0
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%inductance,LDA,prob%work%basis_functions,LDB,BLAS_BETA,prob%work%g,LDC)


  call system_clock(toc)
  if (verbose) print *,"inductance*basis_functions:",real(toc-tic)/countrate,"sec."

  if (allocated(prob%work%matrix_B)) deallocate(prob%work%matrix_B)
  allocate(prob%work%matrix_B(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 13!'

  if (allocated(prob%work%matrix_regularization)) deallocate(prob%work%matrix_regularization)
  allocate(prob%work%matrix_regularization(num_basis_functions, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 14!'

  if (allocated(prob%work%RHS_B)) deallocate(prob%work%RHS_B)
  allocate(prob%work%RHS_B(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_build_matrices Allocation error 15!'

  if (allocated(prob%work%RHS_regularization)) deallocate(prob%work%RHS_regularization)
  allocate(prob%work%RHS_regularization(num_basis_functions),stat=iflag)
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

  prob%work%RHS_B = -dtheta_plasma*dzeta_plasma*matmul( &
       reshape(prob%plasma%Bnormal_from_plasma_current+prob%coil%Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)), prob%work%g)

  call system_clock(toc)
  if (verbose) print *,"Form RHS_B:",real(toc-tic)/countrate,"sec."
  call system_clock(tic)

  norm_normal_plasma_inv1D = reshape(1/prob%plasma%norm_normal_plasma, (/ ntheta_plasma*nzeta_plasma /))
  norm_normal_coil_inv1D   = reshape(1/prob%coil%norm_normal_coil,   (/ ntheta_coil  *nzeta_coil /))
  do j = 1,num_basis_functions
     g_over_N_plasma(:,j) = prob%work%g(:,j) * norm_normal_plasma_inv1D
     f_x_over_N_coil(:,j) = prob%work%f_x(:,j) * norm_normal_coil_inv1D
     f_y_over_N_coil(:,j) = prob%work%f_y(:,j) * norm_normal_coil_inv1D
     f_z_over_N_coil(:,j) = prob%work%f_z(:,j) * norm_normal_coil_inv1D
     f_Laplace_Beltrami_over_N_coil(:,j) = prob%work%f_Laplace_Beltrami(:,j) * norm_normal_coil_inv1D
  end do



  prob%work%matrix_B = 0

  deallocate(norm_normal_plasma_inv1D)
  deallocate(norm_normal_coil_inv1D)

  call system_clock(toc)
  if (verbose) print *,"Prepare for matrix_B:",real(toc-tic)/countrate,"sec."
  call system_clock(tic)


  ! Here we carry out prob%work%matrix_B = (dtheta*dzeta)*(prob%work%g ^ T) * g_over_N_plasma
  ! A = prob%work%g
  ! B = g_over_N_plasma
  ! C = prob%work%inductance
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
  call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%g,LDA,g_over_N_plasma,LDB,BLAS_BETA,prob%work%matrix_B,LDC)

  call system_clock(toc)
  if (verbose) print *,"matmul for matrix_B:",real(toc-tic)/countrate,"sec."


  deallocate(g_over_N_plasma)
    

  prob%work%matrix_regularization = 0
     
  select case (trim(regularization_term_option))
  case (regularization_term_option_chi2_K, regularization_term_option_K_xy)
  
     call system_clock(tic)
     ! Here we carry out prob%work%matrix_regularization += dtheta*dzeta*(prob%work%f_x ^ T) * f_x_over_N_coil
     ! A = prob%work%f_x
     ! B = f_x_over_N_plasma
     ! C = prob%work%matrix_regularization
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%f_x,LDA,f_x_over_N_coil,LDB,BLAS_BETA,prob%work%matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul 1 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     call system_clock(tic)
     ! Here we carry out prob%work%matrix_regularization += dtheta*dzeta*(prob%work%f_y ^ T) * f_y_over_N_coil
     ! A = prob%work%f_y
     ! B = f_y_over_N_plasma
     ! C = prob%work%matrix_regularization
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%f_y,LDA,f_y_over_N_coil,LDB,BLAS_BETA,prob%work%matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul 2 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     if (trim(regularization_term_option) == regularization_term_option_chi2_K) then
        call system_clock(tic)
        ! Here we carry out prob%work%matrix_regularization += dtheta*dzeta*(prob%work%f_z ^ T) * f_z_over_N_coil
        ! A = prob%work%f_z
        ! B = f_z_over_N_plasma
        ! C = prob%work%matrix_regularization
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
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%f_z,LDA,f_z_over_N_coil,LDB,BLAS_BETA,prob%work%matrix_regularization,LDC)
        
        call system_clock(toc)
        if (verbose) print *,"matmul 3 for matrix_regularization:",real(toc-tic)/countrate,"sec."
     end if
     
     call system_clock(tic)
     
     if (trim(regularization_term_option) == regularization_term_option_chi2_K) then     
        prob%work%RHS_regularization = (matmul(prob%work%d_x, f_x_over_N_coil) + matmul(prob%work%d_y, f_y_over_N_coil) + matmul(prob%work%d_z, f_z_over_N_coil)) &
             * (dtheta_coil*dzeta_coil)
     else
        prob%work%RHS_regularization = (matmul(prob%work%d_x, f_x_over_N_coil) + matmul(prob%work%d_y, f_y_over_N_coil)) &
             * (dtheta_coil*dzeta_coil)
     end if
     
     call system_clock(toc)
     if (verbose) print *,"Matmuls for RHS_regularization:",real(toc-tic)/countrate,"sec."

    ! Compute dRHS_Kdomega

    ! Construct dmatrix_Kdomega

  case (regularization_term_option_Laplace_Beltrami)
     ! ------------------------------------------------------------------
     ! Laplace-Beltrami prob%work%matrix and prob%work%RHS:
     
     call system_clock(tic)
     ! Here we carry out prob%work%matrix = dtheta*dzeta*(f ^ T) * f_over_N_coil
     ! A = f
     ! B = f_over_N_plasma
     ! C = prob%work%matrix
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
     call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,prob%work%f_Laplace_Beltrami,LDA,f_Laplace_Beltrami_over_N_coil,LDB,BLAS_BETA,prob%work%matrix_regularization,LDC)
     
     call system_clock(toc)
     if (verbose) print *,"matmul for matrix_regularization:",real(toc-tic)/countrate,"sec."
     
     call system_clock(tic)
     
     prob%work%RHS_regularization = matmul(prob%work%d_Laplace_Beltrami, f_Laplace_Beltrami_over_N_coil) * (dtheta_coil*dzeta_coil)
     
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


  end associate
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



