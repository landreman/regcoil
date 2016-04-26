module global_variables

  use stel_kinds

  implicit none

  integer :: transfer_matrix_option = 1

  integer :: nu_plasma=64, nv_plasma=64, nvl_plasma
  integer :: nu_middle=64, nv_middle=64, nvl_middle
  integer :: nu_outer =64, nv_outer =64, nvl_outer

  integer :: geometry_option_plasma = 0
  integer :: geometry_option_middle = 0
  integer :: geometry_option_outer  = 0

  real(dp) :: R0_plasma = 10.0, R0_middle = 10.0, R0_outer = 10.0
  real(dp) :: a_plasma = 0.5, a_middle = 1.0, a_outer = 1.5
  real(dp) :: separation_middle=0.2, separation_outer=0.4

  character(len=200) :: wout_filename=""
  character(len=200) :: shape_filename_plasma=""
  character(len=200) :: nescin_filename_middle=""
  character(len=200) :: nescin_filename_outer=""
  character(len=200) :: efit_filename=""
  character(len=200) :: outputFilename

  real(dp), dimension(:), allocatable :: u_plasma, v_plasma, vl_plasma
  real(dp), dimension(:,:,:), allocatable :: r_plasma, drdu_plasma, drdv_plasma, normal_plasma

  real(dp), dimension(:), allocatable :: Bnormal_from_1_over_R_field
  real(dp), dimension(:), allocatable :: Bnormal_from_1_over_R_field_inductance
  real(dp), dimension(:), allocatable :: Bnormal_from_1_over_R_field_transfer
  real(dp), dimension(:,:), allocatable :: Bnormal_from_1_over_R_field_uv

  real(dp), dimension(:), allocatable :: Bnormal_from_const_v_coils
  real(dp), dimension(:), allocatable :: Bnormal_from_const_v_coils_inductance
  real(dp), dimension(:), allocatable :: Bnormal_from_const_v_coils_transfer
  real(dp), dimension(:,:), allocatable :: Bnormal_from_const_v_coils_uv

  real(dp), dimension(:), allocatable :: Bnormal_from_plasma_current
  real(dp), dimension(:), allocatable :: Bnormal_from_plasma_current_inductance
  real(dp), dimension(:), allocatable :: Bnormal_from_plasma_current_transfer
  real(dp), dimension(:,:), allocatable :: Bnormal_from_plasma_current_uv

  real(dp), dimension(:), allocatable :: u_middle, v_middle, vl_middle
  real(dp), dimension(:,:,:), allocatable :: r_middle, drdu_middle, drdv_middle, normal_middle

  real(dp), dimension(:), allocatable :: u_outer, v_outer, vl_outer
  real(dp), dimension(:,:,:), allocatable :: r_outer, drdu_outer, drdv_outer, normal_outer
  real(dp), dimension(:,:,:), allocatable :: d2rdu2_middle, d2rdudv_middle, d2rdv2_middle
  real(dp), dimension(:,:), allocatable :: Merkel_Kmn

  real(dp), dimension(:,:), allocatable :: norm_normal_plasma, norm_normal_middle, norm_normal_outer
  real(dp), dimension(:,:), allocatable :: basis_functions_plasma, basis_functions_middle, basis_functions_outer

  real(dp) :: du_plasma, dv_plasma, du_middle, dv_middle, du_outer, dv_outer

  integer :: mpol_plasma=8, mpol_middle=8, mpol_outer=8
  integer :: ntor_plasma=8, ntor_middle=8, ntor_outer=8
  integer :: mnmax_plasma, mnmax_middle, mnmax_outer
  integer :: num_basis_functions_plasma, num_basis_functions_middle, num_basis_functions_outer
  integer, dimension(:), allocatable :: xm_plasma, xm_middle, xm_outer
  integer, dimension(:), allocatable :: xn_plasma, xn_middle, xn_outer

  real(dp), dimension(:), allocatable :: rmns, zmnc, rmnc, zmns
  integer :: mnmax, nfp
  integer, dimension(:), allocatable :: xm, xn
  logical :: lasym
  real(dp), dimension(:,:), allocatable :: inductance_plasma_outer
  real(dp), dimension(:,:), allocatable :: inductance_middle_outer
  real(dp), dimension(:,:), allocatable :: inductance_plasma_middle
  integer :: n_singular_values_inductance_plasma_outer
  integer :: n_singular_values_inductance_middle_outer
  integer :: n_singular_values_inductance_plasma_middle
  real(dp), dimension(:), allocatable :: svd_s_inductance_plasma_outer
  real(dp), dimension(:), allocatable :: svd_s_inductance_middle_outer
  real(dp), dimension(:), allocatable :: svd_s_inductance_plasma_middle

  real(dp), dimension(:,:), allocatable :: svd_uT_inductance_middle_outer, svd_v_inductance_middle_outer
  real(dp), dimension(:,:), allocatable :: svd_u_inductance_plasma_middle, svd_v_inductance_plasma_middle
  real(dp), dimension(:,:), allocatable :: svd_u_inductance_plasma_middle_uv, svd_v_inductance_plasma_middle_uv
  real(dp), dimension(:,:), allocatable :: svd_u_inductance_plasma_middle_all, svd_v_inductance_plasma_middle_all
  real(dp), dimension(:,:), allocatable :: overlap_plasma, overlap_middle

  integer, parameter :: nmax_pseudoinverse_thresholds = 1000
  integer :: n_pseudoinverse_thresholds
  real(dp) :: pseudoinverse_thresholds(nmax_pseudoinverse_thresholds)

  integer :: save_level = 3
  logical :: save_vectors_in_uv_format
  integer :: n_singular_vectors_to_save = 5, n_singular_values_transferMatrix
  integer, dimension(:), allocatable :: n_singular_values_retained
  real(dp), dimension(:,:), allocatable :: svd_s_transferMatrix

  real(dp), dimension(:,:), allocatable :: should_be_identity_plasma
  real(dp), dimension(:,:), allocatable :: should_be_identity_middle
  real(dp), dimension(:,:), allocatable :: should_be_identity_outer

  ! EZCDF doesn't allow 4D arrays, so I can't include sin/cos as a 4th dimension.
  real(dp), dimension(:,:,:), allocatable :: svd_u_transferMatrix, svd_v_transferMatrix
  real(dp), dimension(:,:,:), allocatable :: svd_u_transferMatrix_uv, svd_v_transferMatrix_uv

  logical :: allSVDsSucceeded
  integer :: nfp_imposed = 1

  integer :: symmetry_option = 1, mode_order = 1
  integer :: basis_option_plasma = 1, basis_option_middle = 1, basis_option_outer = 1
  real(dp) :: totalTime

  integer :: efit_num_modes = 10
  real(dp) :: efit_psiN = 0.98

  real(dp) :: mpol_transform_refinement=5, ntor_transform_refinement=1
  real(dp) :: area_plasma, area_middle, area_outer
  logical :: check_orthogonality = .false.

  integer, dimension(:), allocatable :: svd_u_inductance_plasma_middle_dominant_m
  integer, dimension(:), allocatable :: svd_u_inductance_plasma_middle_dominant_n
  integer, dimension(:,:), allocatable :: svd_u_transferMatrix_dominant_m
  integer, dimension(:,:), allocatable :: svd_u_transferMatrix_dominant_n

  logical :: zero_first_transfer_vector_in_overlap = .false.

  real(dp) :: net_poloidal_current_Amperes = 1
  logical :: load_bnorm = .false.
  character(len=200) :: bnorm_filename=""
  real(dp) :: curpol = 1  ! number which multiplies data in bnorm file.

end module global_variables

