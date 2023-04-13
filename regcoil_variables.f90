module regcoil_variables

  use stel_kinds

  implicit none

  integer :: general_option=1
  logical :: verbose = .true.

  character(len=*), parameter :: &
       regularization_term_option_chi2_K = "chi2_K", &
       regularization_term_option_Laplace_Beltrami = "Laplace-Beltrami", &
       regularization_term_option_K_xy = "K_xy", &
       regularization_term_option_K_zeta = "K_zeta"
  character(len=200) :: regularization_term_option = regularization_term_option_chi2_K

  integer :: ntheta_plasma=64, nzeta_plasma=64, nzetal_plasma
  integer :: ntheta_coil=64, nzeta_coil=64, nzetal_coil

  integer :: geometry_option_plasma = 0
  integer :: geometry_option_coil = 0

  real(dp) :: R0_plasma = 10.0, R0_coil = 10.0
  real(dp) :: a_plasma = 0.5, a_coil = 1.0
  real(dp) :: separation=0.2

  character(len=200) :: wout_filename=""
  character(len=200) :: shape_filename_plasma=""
  character(len=200) :: nescin_filename="nescin.out"
  character(len=200) :: nescout_filename=""
  character(len=200) :: efit_filename=""
  character(len=200) :: output_filename

  real(dp), dimension(:), allocatable :: theta_plasma, zeta_plasma, zetal_plasma
  real(dp), dimension(:,:,:), allocatable :: r_plasma, drdtheta_plasma, drdzeta_plasma, normal_plasma

  real(dp), dimension(:,:), allocatable :: g, f_x, f_y, f_z, f_Laplace_Beltrami
  real(dp), dimension(:), allocatable :: h, d_x, d_y, d_z, d_Laplace_Beltrami

  real(dp), dimension(:,:), allocatable :: Bnormal_from_plasma_current
  real(dp), dimension(:,:), allocatable :: Bnormal_from_net_coil_currents
  real(dp), dimension(:,:), allocatable :: matrix_B, matrix_regularization, inductance
  real(dp), dimension(:,:), allocatable :: single_valued_current_potential_mn
  real(dp), dimension(:,:,:), allocatable :: single_valued_current_potential_thetazeta
  real(dp), dimension(:,:,:), allocatable :: current_potential
  real(dp), dimension(:), allocatable :: RHS_B, RHS_regularization
  real(dp), dimension(:,:,:), allocatable :: Bnormal_total
  real(dp), dimension(:,:,:), allocatable :: K2, Laplace_Beltrami2
  real(dp), dimension(:), allocatable :: chi2_B, chi2_K, max_Bnormal, max_K, chi2_Laplace_Beltrami
  real(dp), dimension(:), allocatable :: lp_norm_K, max_K_lse

  real(dp), dimension(:), allocatable :: theta_coil, zeta_coil, zetal_coil
  real(dp), dimension(:,:,:), allocatable :: r_coil, drdtheta_coil, drdzeta_coil, normal_coil
  real(dp), dimension(:,:,:), allocatable :: d2rdtheta2_coil, d2rdthetadzeta_coil, d2rdzeta2_coil

  real(dp), dimension(:,:), allocatable :: norm_normal_plasma, norm_normal_coil
  real(dp), dimension(:,:), allocatable :: basis_functions

  real(dp) :: dtheta_plasma, dzeta_plasma, dtheta_coil, dzeta_coil

  integer :: mpol_potential=12
  integer :: ntor_potential=12
  integer :: mnmax_plasma, mnmax_coil, mnmax_potential
  integer :: num_basis_functions
  integer, dimension(:), allocatable :: xm_plasma, xn_plasma, xm_coil, xn_coil, xm_potential, xn_potential
  real(dp), dimension(:), allocatable :: rmns_plasma, zmnc_plasma, rmnc_plasma, zmns_plasma
  real(dp), dimension(:), allocatable :: rmns_coil, zmnc_coil, rmnc_coil, zmns_coil
  integer :: nfp
  integer :: helicity_ratio ! to produce coils with this helicity i.e. contours of constant curr potential lie on curves of 
!theta = 2*pi - helicity_ratio*Nfp*zeta
  logical :: lasym
  integer :: max_mpol_coil = 24, max_ntor_coil = 24 ! These variables are upper limits on the # of Fourier modes used to describe a uniform-offset coil surface.
  integer :: mpol_coil_filter = 24, ntor_coil_filter = 24

  integer :: save_level = 3
  integer :: nfp_imposed = 1

  integer :: symmetry_option = 1
  real(dp) :: total_time

  integer :: efit_num_modes = 10
  real(dp) :: efit_psiN = 0.98
  real(dp) :: constant_arclength_tolerance = 1.0e-6

  real(dp) :: mpol_transform_refinement=5, ntor_transform_refinement=1
  real(dp) :: area_plasma, area_coil, volume_plasma, volume_coil

  real(dp) :: net_poloidal_current_Amperes = 1
  real(dp) :: net_toroidal_current_Amperes = 0
  logical :: load_bnorm = .false.
  character(len=200) :: bnorm_filename=""
  real(dp) :: curpol = 1  ! number which multiplies data in bnorm file.
  integer :: nbf ! number of Fourier harmonics in FOCUS format boundary.
  integer, dimension(:), allocatable :: bfn, bfm
  real(dp), dimension(:), allocatable :: bfs, bfc

  integer :: nlambda = 4
  real(dp) :: lambda_min = 1.0d-19, lambda_max = 1.0d-13
  real(dp), dimension(:), allocatable :: lambda

  real(dp), dimension(:,:), allocatable :: matrix, this_current_potential
  real(dp), dimension(:), allocatable :: RHS, solution
  real(dp), dimension(:), allocatable :: KDifference_x, KDifference_y, KDifference_z, KDifference_Laplace_Beltrami
  real(dp), dimension(:,:), allocatable :: this_K2_times_N, this_Laplace_Beltrami2_times_N

  ! Variables needed by LAPACK:
  integer :: LAPACK_INFO, LAPACK_LWORK
  real(dp), dimension(:), allocatable :: LAPACK_WORK
  integer, dimension(:), allocatable :: LAPACK_IPIV

  real(dp) :: target_value = 8.0d+6
  real(dp) :: lambda_search_tolerance = 1.0d-5
  integer :: exit_code = 0
  real(dp) :: chi2_B_target = 0

  ! Input parameters for adjoint
  logical :: fixed_norm_sensitivity_option = .false.
  integer :: sensitivity_option = 1
  integer :: sensitivity_symmetry_option = 1
  integer :: nmax_sensitivity = 1
  integer :: mmax_sensitivity = 1
  real(dp) :: target_option_p = 4.0
  real(dp) :: coil_plasma_dist_lse_p = 1.0d4

  ! Needed for sensitivity calculation
  real(dp), dimension(:), allocatable :: L_p_norm_with_area
  real(dp), dimension(:,:), allocatable :: dRMSKdomega
  real(dp), dimension(:), allocatable :: darea_coildomega
  integer :: mnmax_sensitivity
  integer :: nomega_coil
  integer, dimension(:), allocatable :: omega_coil
  integer, dimension(:), allocatable :: xm_sensitivity, xn_sensitivity
  real(dp), dimension(:,:), allocatable :: dchi2Bdphi
  real(dp), dimension(:,:,:), allocatable :: dddomega, dfxdomega, dfydomega, dfzdomega
  real(dp), dimension(:,:,:), allocatable :: dnorm_normaldomega, dgdomega
  real(dp), dimension(:,:,:), allocatable :: dnormxdomega, dnormydomega, dnormzdomega
  real(dp), dimension(:,:), allocatable :: dchi2domega, dchi2Kdomega, dchi2Bdomega, dchi2Kdomega_withoutadjoint, dchi2Bdomega_withoutadjoint
  real(dp), dimension(:,:,:,:), allocatable :: drdomega
  real(dp), dimension(:,:,:), allocatable :: domegadxdtheta, domegadxdzeta, domegadydtheta, domegadydzeta, domegadzdtheta, domegadzdzeta
  real(dp), dimension(:,:), allocatable :: dhdomega
  real(dp), dimension(:,:), allocatable :: domegadx, domegady, domegadz
  real(dp), dimension(:,:), allocatable :: q_K, q_B, dFdomega, dRHS_Bdomega, dRHS_Kdomega
  real(dp), dimension(:,:,:), allocatable :: dmatrix_Kdomega, dmatrix_Bdomega
  real(dp), dimension(:), allocatable :: dvolume_coildomega
  real(dp), dimension(:), allocatable :: dcoil_plasma_dist_mindomega
  real(dp) :: coil_plasma_dist_min, coil_plasma_dist_max
  real(dp) :: coil_plasma_dist_min_lse, coil_plasma_dist_max_lse

  character(len=*), parameter :: &
       target_option_max_K = "max_K", &
       target_option_rms_K = "rms_K", &
       target_option_chi2_K = "chi2_K", &
       target_option_max_Bnormal = "max_Bnormal", &
       target_option_rms_Bnormal = "rms_Bnormal", &
       target_option_chi2_B = "chi2_B", &
       target_option_max_K_lse = "max_K_lse", &
       target_option_lp_norm_K = "lp_norm_K"
  character(len=200) :: target_option = target_option_max_K

  namelist / regcoil_nml / ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, &
       geometry_option_plasma, geometry_option_coil, &
       R0_plasma, R0_coil, a_plasma, a_coil, &
       separation, wout_filename, &
       save_level, nfp_imposed, symmetry_option, &
       mpol_potential, ntor_potential, mpol_coil_filter, ntor_coil_filter, &
       nescin_filename, efit_filename, efit_psiN, efit_num_modes, &
       mpol_transform_refinement, ntor_transform_refinement, max_mpol_coil, max_ntor_coil, &
       net_poloidal_current_Amperes, net_toroidal_current_Amperes, &
       load_bnorm, bnorm_filename, &
       shape_filename_plasma, nlambda, lambda_min, lambda_max, general_option, regularization_term_option, verbose, nescout_filename, &
       target_option, target_value, lambda_search_tolerance, &
       sensitivity_option, nmax_sensitivity, mmax_sensitivity, &
       sensitivity_symmetry_option, target_option_p, &
       fixed_norm_sensitivity_option, coil_plasma_dist_lse_p, helicity_ratio

end module regcoil_variables

