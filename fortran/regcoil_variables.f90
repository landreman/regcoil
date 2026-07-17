! Instance state for REGCOIL (Phase 5).
! Pattern for new routines: pass type(regcoil_t) as the first argument (name: prob).
! Nested components: plasma, coil, input, output, work.
! Immutable string option constants remain module parameters (not problem state).

module regcoil_variables

  use stel_kinds

  implicit none

  character(len=*), parameter :: &
       regularization_term_option_chi2_K = "chi2_K", &
       regularization_term_option_Laplace_Beltrami = "Laplace-Beltrami", &
       regularization_term_option_K_xy = "K_xy", &
       regularization_term_option_K_zeta = "K_zeta"

  character(len=*), parameter :: &
       target_option_max_K = "max_K", &
       target_option_rms_K = "rms_K", &
       target_option_chi2_K = "chi2_K", &
       target_option_max_Bnormal = "max_Bnormal", &
       target_option_rms_Bnormal = "rms_Bnormal", &
       target_option_chi2_B = "chi2_B", &
       target_option_max_K_lse = "max_K_lse", &
       target_option_lp_norm_K = "lp_norm_K"

  type :: regcoil_plasma_surface_t
     integer :: ntheta_plasma = 64, nzeta_plasma = 64, nzetal_plasma = 0
     integer :: geometry_option_plasma = 0
     real(dp) :: R0_plasma = 10.0_dp, a_plasma = 0.5_dp
     character(len=200) :: wout_filename = ""
     character(len=200) :: shape_filename_plasma = ""
     character(len=200) :: efit_filename = ""
     real(dp), dimension(:), allocatable :: theta_plasma, zeta_plasma, zetal_plasma
     real(dp), dimension(:,:,:), allocatable :: r_plasma, drdtheta_plasma, drdzeta_plasma, normal_plasma
     real(dp), dimension(:,:), allocatable :: Bnormal_from_plasma_current
     real(dp), dimension(:,:), allocatable :: norm_normal_plasma
     real(dp) :: dtheta_plasma = 0, dzeta_plasma = 0
     integer :: mnmax_plasma = 0
     integer, dimension(:), allocatable :: xm_plasma, xn_plasma
     real(dp), dimension(:), allocatable :: rmns_plasma, zmnc_plasma, rmnc_plasma, zmns_plasma
     integer :: nfp = 0
     logical :: lasym = .false.
     integer :: efit_num_modes = 10
     real(dp) :: efit_psiN = 0.98_dp
     real(dp) :: area_plasma = 0, volume_plasma = 0
     integer :: nbf = 0
     integer, dimension(:), allocatable :: bfn, bfm
     real(dp), dimension(:), allocatable :: bfs, bfc
     real(dp) :: constant_arclength_tolerance = 1.0e-6_dp
     real(dp) :: mpol_transform_refinement = 5, ntor_transform_refinement = 1
  end type regcoil_plasma_surface_t

  type :: regcoil_coil_surface_t
     integer :: ntheta_coil = 64, nzeta_coil = 64, nzetal_coil = 0
     integer :: geometry_option_coil = 0
     real(dp) :: R0_coil = 10.0_dp, a_coil = 1.0_dp
     real(dp) :: separation = 0.2_dp
     character(len=200) :: nescin_filename = "nescin.out"
     character(len=200) :: nescout_filename = ""
     real(dp), dimension(:), allocatable :: theta_coil, zeta_coil, zetal_coil
     real(dp), dimension(:,:,:), allocatable :: r_coil, drdtheta_coil, drdzeta_coil, normal_coil
     real(dp), dimension(:,:,:), allocatable :: d2rdtheta2_coil, d2rdthetadzeta_coil, d2rdzeta2_coil
     real(dp), dimension(:,:), allocatable :: Bnormal_from_net_coil_currents
     real(dp), dimension(:,:), allocatable :: norm_normal_coil
     real(dp) :: dtheta_coil = 0, dzeta_coil = 0
     integer :: mnmax_coil = 0
     integer, dimension(:), allocatable :: xm_coil, xn_coil
     real(dp), dimension(:), allocatable :: rmns_coil, zmnc_coil, rmnc_coil, zmns_coil
     integer :: max_mpol_coil = 24, max_ntor_coil = 24
     integer :: mpol_coil_filter = 24, ntor_coil_filter = 24
     real(dp) :: area_coil = 0, volume_coil = 0
  end type regcoil_coil_surface_t

  type :: regcoil_solver_input_t
     integer :: general_option = 1
     logical :: verbose = .true.
     character(len=200) :: regularization_term_option = regularization_term_option_chi2_K
     integer :: save_level = 3
     integer :: nfp_imposed = 1
     integer :: symmetry_option = 1
     integer :: mpol_potential = 12
     integer :: ntor_potential = 12
     real(dp) :: net_poloidal_current_Amperes = 1
     real(dp) :: net_toroidal_current_Amperes = 0
     logical :: load_bnorm = .false.
     character(len=200) :: bnorm_filename = ""
     real(dp) :: curpol = 1
     integer :: nlambda = 4
     real(dp) :: lambda_min = 1.0d-19, lambda_max = 1.0d-13
     real(dp), dimension(:), allocatable :: lambda
     real(dp) :: target_value = 8.0d+6
     real(dp) :: lambda_search_tolerance = 1.0d-5
     character(len=200) :: target_option = target_option_max_K
     real(dp) :: target_option_p = 4.0_dp
     real(dp) :: chi2_B_target = 0
     character(len=200) :: output_filename = ""
  end type regcoil_solver_input_t

  type :: regcoil_solver_output_t
     real(dp), dimension(:,:), allocatable :: single_valued_current_potential_mn
     real(dp), dimension(:,:,:), allocatable :: single_valued_current_potential_thetazeta
     real(dp), dimension(:,:,:), allocatable :: current_potential
     real(dp), dimension(:,:,:), allocatable :: Bnormal_total
     real(dp), dimension(:,:,:), allocatable :: K2, Laplace_Beltrami2
     real(dp), dimension(:), allocatable :: chi2_B, chi2_K, max_Bnormal, max_K, chi2_Laplace_Beltrami
     real(dp), dimension(:), allocatable :: lp_norm_K, max_K_lse
     integer :: exit_code = 0
     real(dp) :: total_time = 0
  end type regcoil_solver_output_t

  ! Matrices, basis, and LAPACK scratch for one problem instance.
  type :: regcoil_work_t
     real(dp), dimension(:,:), allocatable :: g, f_x, f_y, f_z, f_Laplace_Beltrami
     real(dp), dimension(:), allocatable :: h, d_x, d_y, d_z, d_Laplace_Beltrami
     real(dp), dimension(:,:), allocatable :: matrix_B, matrix_regularization, inductance
     real(dp), dimension(:), allocatable :: RHS_B, RHS_regularization
     real(dp), dimension(:,:), allocatable :: basis_functions
     integer :: mnmax_potential = 0
     integer :: num_basis_functions = 0
     integer, dimension(:), allocatable :: xm_potential, xn_potential
     real(dp), dimension(:,:), allocatable :: matrix
     real(dp), dimension(:), allocatable :: RHS, solution
     real(dp), dimension(:,:), allocatable :: this_current_potential
     real(dp), dimension(:), allocatable :: KDifference_x, KDifference_y, KDifference_z
     real(dp), dimension(:), allocatable :: KDifference_Laplace_Beltrami
     real(dp), dimension(:,:), allocatable :: this_K2_times_N, this_Laplace_Beltrami2_times_N
     integer :: LAPACK_INFO = 0, LAPACK_LWORK = 0
     real(dp), dimension(:), allocatable :: LAPACK_WORK
     integer, dimension(:), allocatable :: LAPACK_IPIV
  end type regcoil_work_t

  type :: regcoil_t
     type(regcoil_plasma_surface_t) :: plasma
     type(regcoil_coil_surface_t) :: coil
     type(regcoil_solver_input_t) :: input
     type(regcoil_solver_output_t) :: output
     type(regcoil_work_t) :: work
  end type regcoil_t

end module regcoil_variables
