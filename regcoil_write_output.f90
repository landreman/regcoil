subroutine regcoil_write_output

  use regcoil_variables
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_nfp = "nfp", &
       vn_geometry_option_plasma = "geometry_option_plasma", &
       vn_geometry_option_coil = "geometry_option_coil", &
       vn_ntheta_plasma = "ntheta_plasma", &
       vn_nzeta_plasma = "nzeta_plasma", &
       vn_nzetal_plasma = "nzetal_plasma", &
       vn_ntheta_coil = "ntheta_coil", &
       vn_nzeta_coil = "nzeta_coil", &
       vn_nzetal_coil = "nzetal_coil", &
       vn_a_plasma  = "a_plasma", &
       vn_a_coil  = "a_coil", &
       vn_R0_plasma  = "R0_plasma", &
       vn_R0_coil  = "R0_coil", &
       vn_mpol_potential = "mpol_potential", &
       vn_ntor_potential = "ntor_potential", &
       vn_mnmax_potential = "mnmax_potential", &
       vn_mnmax_plasma = "mnmax_plasma", &
       vn_mnmax_coil = "mnmax_coil", &
       vn_num_basis_functions = "num_basis_functions", &
       vn_symmetry_option = "symmetry_option", &
       vn_area_plasma = "area_plasma", &
       vn_area_coil = "area_coil", &
       vn_volume_plasma = "volume_plasma", &
       vn_volume_coil = "volume_coil", &
       vn_net_poloidal_current_Amperes = "net_poloidal_current_Amperes", &
       vn_net_toroidal_current_Amperes = "net_toroidal_current_Amperes", &
       vn_curpol = "curpol", &
       vn_nlambda = "nlambda", &
       vn_total_time = "total_time", &
       vn_exit_code = "exit_code", &
       vn_chi2_B_target = "chi2_B_target", &
       vn_sensitivity_option = "sensitivity_option", &
       vn_nmax_sensitivity = "nmax_sensitivity", &
       vn_mmax_sensitivity = "mmax_sensitivity", &
       vn_mnmax_sensitivity = "mnmax_sensitivity", &
       vn_nomega_coil = "nomega_coil", &
       vn_sensitivity_symmetry_option = "sensitivity_symmetry_option", &
       vn_fixed_norm_sensitivity_option = "fixed_norm_sensitivity_option", &
       vn_coil_plasma_dist_min = "coil_plasma_dist_min", &
       vn_coil_plasma_dist_max = "coil_plasma_dist_max", &
       vn_coil_plasma_dist_min_lse = "coil_plasma_dist_min_lse", &
       vn_coil_plasma_dist_max_lse = "coil_plasma_dist_max_lse", &
       vn_coil_plasma_dist_lse_p = "coil_plasma_dist_lse_p", &
       vn_helicity_ratio = "helicity_ratio"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_theta_plasma = "theta_plasma", &
       vn_zeta_plasma = "zeta_plasma", &
       vn_zetal_plasma = "zetal_plasma", &
       vn_theta_coil = "theta_coil", &
       vn_zeta_coil = "zeta_coil", &
       vn_zetal_coil = "zetal_coil", &
       vn_xm_potential = "xm_potential", &
       vn_xn_potential = "xn_potential", &
       vn_xm_plasma = "xm_plasma", &
       vn_xn_plasma = "xn_plasma", &
       vn_xm_coil = "xm_coil", &
       vn_xn_coil = "xn_coil", &
       vn_rmnc_plasma = "rmnc_plasma", &
       vn_rmns_plasma = "rmns_plasma", &
       vn_zmnc_plasma = "zmnc_plasma", &
       vn_zmns_plasma = "zmns_plasma", &
       vn_rmnc_coil = "rmnc_coil", &
       vn_rmns_coil = "rmns_coil", &
       vn_zmnc_coil = "zmnc_coil", &
       vn_zmns_coil = "zmns_coil", &
       vn_h = "h", &
       vn_RHS_B = "RHS_B", &
       vn_RHS_regularization = "RHS_regularization", &
       vn_lambda = "lambda", &
       vn_chi2_B = "chi2_B", &
       vn_chi2_K = "chi2_K", &
       vn_chi2_Laplace_Beltrami = "chi2_Laplace_Beltrami", &
       vn_max_Bnormal = "max_Bnormal", &
       vn_max_K = "max_K", &
       vn_xn_sensitivity = "xn_sensitivity", &
       vn_xm_sensitivity = "xm_sensitivity", &
       vn_omega_coil = "omega_coil", &
       vn_dvolume_coildomega = "dvolume_coildomega", &
       vn_darea_coildomega = "darea_coildomega", &
       vn_dcoil_plasma_dist_mindomega = "dcoil_plasma_dist_mindomega", &
       vn_max_K_lse = "max_K_lse", &
       vn_lp_norm_K = "lp_norm_K"

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_norm_normal_plasma  = "norm_normal_plasma", &
       vn_norm_normal_coil  = "norm_normal_coil", &
       vn_Bnormal_from_plasma_current = "Bnormal_from_plasma_current", &
       vn_Bnormal_from_net_coil_currents = "Bnormal_from_net_coil_currents", &
       vn_inductance = "inductance", &
       vn_g = "g", &
       vn_matrix_B = "matrix_B", &
       vn_matrix_regularization = "matrix_regularization", &
       vn_single_valued_current_potential_mn = "single_valued_current_potential_mn", &
       vn_dchi2Bdomega = "dchi2Bdomega", &
       vn_dchi2Kdomega = "dchi2Kdomega", &
       vn_dchi2domega = "dchi2domega", &
       vn_dRMSKdomega = "dRMSKdomega"

  ! Arrays with dimension 3
  character(len=*), parameter :: &
       vn_r_plasma  = "r_plasma", &
       vn_r_coil  = "r_coil", &
       vn_drdtheta_plasma  = "drdtheta_plasma", &
       vn_drdtheta_coil  = "drdtheta_coil", &
       vn_drdzeta_plasma  = "drdzeta_plasma", &
       vn_drdzeta_coil  = "drdzeta_coil", &
       vn_normal_plasma = "normal_plasma", &
       vn_normal_coil = "normal_coil", &
       vn_single_valued_current_potential_thetazeta = "single_valued_current_potential_thetazeta", &
       vn_current_potential = "current_potential", &
       vn_Bnormal_total = "Bnormal_total", &
       vn_K2 = "K2", &
       vn_Laplace_Beltrami2 = "Laplace_Beltrami2"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       ntheta_plasma_dim = (/'ntheta_plasma'/), &
       nzeta_plasma_dim = (/'nzeta_plasma'/), &
       nzetal_plasma_dim = (/'nzetal_plasma'/), &
       ntheta_coil_dim = (/'ntheta_coil'/), &
       nzeta_coil_dim = (/'nzeta_coil'/), &
       nzetal_coil_dim = (/'nzetal_coil'/), &
       mnmax_potential_dim = (/'mnmax_potential'/), &
       mnmax_plasma_dim = (/'mnmax_plasma'/), &
       mnmax_coil_dim = (/'mnmax_coil'/), &
       nthetanzeta_plasma_dim = (/'ntheta_nzeta_plasma'/), &
       num_basis_functions_dim = (/'num_basis_functions'/), &
       nlambda_dim = (/'nlambda'/), &
       nomega_coil_dim = (/'nomega_coil'/), &
       ntheta_times_nzeta_coil_dim = (/'ntheta_times_nzeta_coil'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       ntheta_nzeta_plasma_dim = (/ character(len=50) :: 'ntheta_plasma','nzeta_plasma'/), &
       ntheta_nzeta_coil_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil'/), &
       nthetanzeta_plasma_nthetanzeta_coil_dim = (/ character(len=50) :: 'ntheta_nzeta_plasma','ntheta_nzeta_coil'/), &
       nthetanzeta_plasma_basis_dim = (/ character(len=50) :: 'ntheta_nzeta_plasma','num_basis_functions'/), &
       basis_basis_dim = (/ character(len=50) :: 'num_basis_functions','num_basis_functions'/), &
       basis_nlambda_dim = (/ character(len=50) :: 'num_basis_functions','nlambda'/), &
       nomega_coil_nlambda_dim = (/ character(len=50) :: 'nomega_coil', 'nlambda'/), &
       nthetanzeta_coil_basis_dim = (/ character(len=50) :: &
         'ntheta_times_nzeta_coil','num_basis_functions'/), &
       nomega_coil_nthetanzeta_plasma_dim = (/character(len=50) :: &
         'nomega_coil','ntheta_nzeta_plasma'/), &
       nlambda_nthetanzeta_coil_dim = (/character(len=50) :: &
          'nlambda', 'ntheta_times_nzeta_coil'/), &
       nomega_coil_num_basis_functions_dim = (/character(len=50) :: &
          'nomega_coil', 'num_basis_functions'/), &
       nomega_nlambda_dim = (/ character(len=50) :: &
          'nomega_coil', 'nlambda' /), &
       ntheta_nzetal_coil_dim = (/ character(len=50) :: &
          'ntheta_coil', 'nzetal_coil' /), &
       nlambda_nomega_dim = (/ character(len=50) :: &
          'nlambda', 'nomega_coil' /)

  ! Arrays with dimension 3:
  character(len=*), parameter, dimension(3) :: &
       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &
       xyz_ntheta_nzetal_coil_dim = (/ character(len=50) :: 'xyz','ntheta_coil','nzetal_coil'/), &
       ntheta_nzeta_coil_nlambda_dim = (/ character(len=50) :: 'ntheta_coil','nzeta_coil','nlambda'/), &
       ntheta_nzeta_plasma_nlambda_dim = (/ character(len=50) :: 'ntheta_plasma','nzeta_plasma','nlambda'/), &
       nomega_coil_ntheta_nzeta_coil_dim = (/ character(len=50) :: 'nomega_coil', 'ntheta_coil', 'nzeta_coil'/), &
       xyz_nomega_coil_ntheta_nzetal_coil_dim = (/character(len=50) :: 'xyz', 'nomega_coil', 'ntheta_times_nzetal_coil'/), &
       nomega_coil_ntheta_times_nzeta_num_basis_functions_dim = (/character(len=50) :: 'nomega_coil', 'ntheta_times_nzeta_plasma', 'num_basis_functions' /), &
       nomega_coil_ntheta_nzetal_coil_dim = (/character(len=50) :: 'nomega_coil', 'ntheta_coil', &
         'nzetal_coil'/), &
       nomega_coil_ntheta_times_nzeta_coil_basis_dim = (/character(len=50):: 'nomega_coil', 'ntheta_times_nzeta_coil','num_basis_functions' /), &
       nomega_coil_num_basis_num_basis_dim = (/character(len=50) :: 'nomega_coil', 'num_basis_functions', &
         'num_basis_functions' /), &
       nlambda_num_basis_nomega_coil_dim = (/character(len=50):: 'nlambda',  &
        'num_basis_functions', 'nomega_coil' /), &
       ntheta_times_nzeta_plasma_coil_nomega_coil_dim = (/character(len=50) :: 'ntheta_times_nzeta_plasma', &
        'ntheta_times_nzeta_coil', 'nomega_coil' /), &
       ntheta_times_nzeta_num_basis_functions_nomega_dim = (/character(len=50) :: 'ntheta_times_nzeta_plasma', &
        'num_basis_functions', 'nomega_coil' /)


  character(len=*), parameter :: input_parameter_text = ' See the user manual documentation for the input parameter of the same name.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.' // input_parameter_text)

  call cdf_define(ncid, vn_geometry_option_plasma, geometry_option_plasma)
  call cdf_setatt(ncid, vn_geometry_option_plasma, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_geometry_option_coil, geometry_option_coil)
  call cdf_setatt(ncid, vn_geometry_option_coil, 'Method used to define the geometry of the coil winding surface.' // input_parameter_text)

  call cdf_define(ncid, vn_ntheta_plasma, ntheta_plasma)
  call cdf_setatt(ncid, vn_ntheta_plasma, 'Number of grid points used in the poloidal angle theta on the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzeta_plasma, nzeta_plasma)
  call cdf_setatt(ncid, vn_nzeta_plasma, 'Number of grid points used in the toridal angle zeta per identical toroidal period, on the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzetal_plasma, nzetal_plasma)
  call cdf_setatt(ncid, vn_nzetal_plasma, 'Number of grid points used in the toridal angle, including all the nfp identical toroidal periods, on the plasma surface.')

  call cdf_define(ncid, vn_ntheta_coil, ntheta_coil)
  call cdf_setatt(ncid, vn_ntheta_coil, 'Number of grid points used in the poloidal angle theta on the coil surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzeta_coil, nzeta_coil)
  call cdf_setatt(ncid, vn_nzeta_coil, 'Number of grid points used in the toridal angle zeta per identical toroidal period, on the coil surface.' // input_parameter_text)

  call cdf_define(ncid, vn_nzetal_coil, nzetal_coil)
  call cdf_setatt(ncid, vn_nzetal_coil, 'Number of grid points used in the toridal angle, including all the nfp identical toroidal periods, on the coil surface.')

  call cdf_define(ncid, vn_a_plasma, a_plasma)
  call cdf_define(ncid, vn_a_coil, a_coil)
  call cdf_define(ncid, vn_R0_plasma, R0_plasma)
  call cdf_define(ncid, vn_R0_coil, R0_coil)

  call cdf_define(ncid, vn_mpol_potential, mpol_potential)
  call cdf_setatt(ncid, vn_mpol_potential, 'The maximum poloidal mode number retained in the current potential.' // input_parameter_text)

  call cdf_define(ncid, vn_ntor_potential, ntor_potential)
  call cdf_setatt(ncid, vn_ntor_potential, 'ntor_potential * nfp is the maximum toroidal mode number retained in the current potential.' // input_parameter_text)

  call cdf_define(ncid, vn_mnmax_potential, mnmax_potential)
  call cdf_setatt(ncid, vn_mnmax_potential, 'Number of unique (m,n) pairs for the Fourier modes retained in the single-valued part of the current potential. ' // &
       'Equal to mpol_potential*(ntor_potential*2+1) + ntor_potential.')

  call cdf_define(ncid, vn_mnmax_plasma, mnmax_plasma)
  call cdf_setatt(ncid, vn_mnmax_plasma, 'Number of unique (m,n) pairs for the Fourier modes retained in the plasma surface.')

  call cdf_define(ncid, vn_mnmax_coil, mnmax_coil)
  call cdf_setatt(ncid, vn_mnmax_coil, 'Number of unique (m,n) pairs for the Fourier modes retained in the coil winding surface.')

  call cdf_define(ncid, vn_num_basis_functions, num_basis_functions)
  call cdf_setatt(ncid, vn_num_basis_functions, 'Number of cos(m*theta-n*zeta) and/or sin(m*theta-n*zeta) Fourier modes retained ' // &
       'in the single-valued part of the current potential. Equal to mnmax_potential * 1 or 2, depending on symmetry_option.')

  call cdf_define(ncid, vn_symmetry_option, symmetry_option)
  call cdf_setatt(ncid, vn_symmetry_option, '1 = The single-valued part of the current potential was forced to be stellarator-symmetric. ' // &
       '2 = The single-valued part of the current potential was forced to be stellarator-antisymmetric. ' // &
       '3 = The single-valued part of the current potential was allowed to have both stellarator-symmetric and antisymmetric components')

  call cdf_define(ncid, vn_area_plasma, area_plasma)
  call cdf_setatt(ncid, vn_area_plasma, 'Area of the plasma surface in meters^2')

  call cdf_define(ncid, vn_area_coil, area_coil)
  call cdf_setatt(ncid, vn_area_coil, 'Area of the coil winding surface in meters^2')

  call cdf_define(ncid, vn_volume_plasma, volume_plasma)
  call cdf_setatt(ncid, vn_volume_plasma, 'Volume of the plasma surface in meters^3')

  call cdf_define(ncid, vn_volume_coil, volume_coil)
  call cdf_setatt(ncid, vn_volume_coil, 'Volume of the coil winding surface in meters^3')

  call cdf_define(ncid, vn_net_poloidal_current_Amperes, net_poloidal_current_Amperes)
  call cdf_setatt(ncid, vn_net_poloidal_current_Amperes, 'Net current (in Amperes) that flows on the coil winding surface in the poloidal direction. ' // &
       'This quantity corresponds to G in the 2017 Nuclear Fusion paper. For modular coils (as opposed to helical or wavy-PF coils), ' // &
       'this quantity is then the sum of the currents in all the modular coils.')

  call cdf_define(ncid, vn_net_toroidal_current_Amperes, net_toroidal_current_Amperes)
  call cdf_setatt(ncid, vn_net_toroidal_current_Amperes, 'Net current (in Amperes) that flows on the coil winding surface in the toroidal direction. ' // &
       'This quantity corresponds to I in the 2017 Nuclear Fusion paper. For modular coils, this quantity is 0.')

  call cdf_define(ncid, vn_helicity_ratio, helicity_ratio)
  call cdf_setatt(ncid, vn_helicity_ratio, 'Helicity of the helical coils (i.e. number of poloidal turns per field period).')

  call cdf_define(ncid, vn_curpol, curpol)

  call cdf_define(ncid, vn_nlambda, nlambda)
  call cdf_setatt(ncid, vn_nlambda, 'Number of values of the regularization parameter lambda examined.')

  call cdf_define(ncid, vn_total_time, total_time)
  call cdf_setatt(ncid, vn_total_time, 'Total time it took regcoil to run, in seconds.')

  call cdf_define(ncid, vn_exit_code, exit_code)
  call cdf_setatt(ncid, vn_exit_code, "Only meaningful when general_option = 4 or 5 so a lambda search is performed. " // &
       "exit_code = 0 means the lambda search was successful. " // &
       "exit_code = -1 means the lambda search did not converge to the requested tolerance within nlambda iterations. " // &
       "exit_code = -2 means the current_density_target you have set is not achievable because it is too low. " // &
       "exit_code = -3 means the current_density_target you have set is not achievable because it is too high.")

  if (general_option==4 .or. general_option==5) then
     call cdf_define(ncid, vn_chi2_B_target, chi2_B_target)
     call cdf_setatt(ncid, vn_chi2_B_target, 'The value of chi^2_B at the final value of regularization parameter lambda resulting from the lambda search. ' // &
          'Units = Tesla^2 meters^2.')
  end if

  call cdf_define(ncid, vn_sensitivity_option, sensitivity_option)
  if (sensitivity_option > 1) then
    call cdf_define(ncid, vn_mmax_sensitivity, mmax_sensitivity)
    call cdf_define(ncid, vn_nmax_sensitivity, nmax_sensitivity)
    call cdf_define(ncid, vn_mnmax_sensitivity, mnmax_sensitivity)
    call cdf_define(ncid, vn_nomega_coil, nomega_coil)
    call cdf_define(ncid, vn_sensitivity_symmetry_option, sensitivity_symmetry_option)
    call cdf_define(ncid, vn_fixed_norm_sensitivity_option, fixed_norm_sensitivity_option)
  endif
  if (sensitivity_option > 1) then
    call cdf_define(ncid, vn_coil_plasma_dist_min, coil_plasma_dist_min)
    call cdf_define(ncid, vn_coil_plasma_dist_max, coil_plasma_dist_max)
    call cdf_define(ncid, vn_coil_plasma_dist_lse_p, coil_plasma_dist_lse_p)
    call cdf_define(ncid, vn_coil_plasma_dist_min_lse, coil_plasma_dist_min_lse)
    call cdf_define(ncid, vn_coil_plasma_dist_max_lse, coil_plasma_dist_max_lse)
  end if

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_theta_plasma, theta_plasma, dimname=ntheta_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the poloidal angle on the plasma surface.')

  call cdf_define(ncid, vn_zeta_plasma, zeta_plasma, dimname=nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the toroidal angle on the plasma surface. ' // &
       'Only the first of the nfp identical toroidal periods is included')

  call cdf_define(ncid, vn_zetal_plasma, zetal_plasma, dimname=nzetal_plasma_dim)
  call cdf_setatt(ncid, vn_theta_plasma, 'Grid points of the toroidal angle on the plasma surface, including all nfp toroidal periods.')

  call cdf_define(ncid, vn_theta_coil, theta_coil, dimname=ntheta_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the poloidal angle on the coil surface.')

  call cdf_define(ncid, vn_zeta_coil, zeta_coil, dimname=nzeta_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the toroidal angle on the coil surface. ' // &
       'Only the first of the nfp identical toroidal periods is included')

  call cdf_define(ncid, vn_zetal_coil, zetal_coil, dimname=nzetal_coil_dim)
  call cdf_setatt(ncid, vn_theta_coil, 'Grid points of the toroidal angle on the coil surface, including all nfp toroidal periods.')

  call cdf_define(ncid, vn_xm_potential, xm_potential, dimname=mnmax_potential_dim)
  call cdf_setatt(ncid, vn_xm_potential, 'Values of poloidal mode number m used in the Fourier representation of the single-valued part of the current potential.')
  call cdf_define(ncid, vn_xn_potential, xn_potential, dimname=mnmax_potential_dim)
  call cdf_setatt(ncid, vn_xm_potential, 'Values of toroidal mode number n used in the Fourier representation of the single-valued part of the current potential.')

  call cdf_define(ncid, vn_xm_plasma, xm_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_xm_plasma, 'Values of poloidal mode number m used in the Fourier representation of the plasma surface.')
  call cdf_define(ncid, vn_xn_plasma, xn_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_xm_plasma, 'Values of toroidal mode number n used in the Fourier representation of the plasma surface.')

  call cdf_define(ncid, vn_xm_coil, xm_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_xm_coil, 'Values of poloidal mode number m used in the Fourier representation of the coil winding surface.')
  call cdf_define(ncid, vn_xn_coil, xn_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_xm_coil, 'Values of toroidal mode number n used in the Fourier representation of the coil winding surface.')

  call cdf_define(ncid, vn_rmnc_plasma, rmnc_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_rmnc_plasma, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
       'R(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
  if (lasym) then
     call cdf_define(ncid, vn_rmns_plasma, rmns_plasma, dimname=mnmax_plasma_dim)
     call cdf_setatt(ncid, vn_rmns_plasma, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
          'R(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
     call cdf_define(ncid, vn_zmnc_plasma, zmnc_plasma, dimname=mnmax_plasma_dim)
     call cdf_setatt(ncid, vn_zmnc_plasma, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
          'Z(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')
  end if
  call cdf_define(ncid, vn_zmns_plasma, zmns_plasma, dimname=mnmax_plasma_dim)
  call cdf_setatt(ncid, vn_zmns_plasma, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
       'Z(theta,zeta) defining the plasma surface. The corresponding mode numbers (m,n) are stored in xm_plasma and xn_plasma.')

  call cdf_define(ncid, vn_rmnc_coil, rmnc_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_rmnc_coil, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
       'R(theta,zeta) defining the coil surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
  if (lasym) then
     call cdf_define(ncid, vn_rmns_coil, rmns_coil, dimname=mnmax_coil_dim)
     call cdf_setatt(ncid, vn_rmns_coil, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the cylindrical coordinate ' // &
          'R(theta,zeta) defining the coil surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
     call cdf_define(ncid, vn_zmnc_coil, zmnc_coil, dimname=mnmax_coil_dim)
     call cdf_setatt(ncid, vn_zmnc_coil, 'Amplitudes of the cosine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
          'Z(theta,zeta) defining the coil surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')
  end if
  call cdf_define(ncid, vn_zmns_coil, zmns_coil, dimname=mnmax_coil_dim)
  call cdf_setatt(ncid, vn_zmns_coil, 'Amplitudes of the sine(m*theta-n*zeta) terms in a Fourier expansion of the coordinate ' // &
       'Z(theta,zeta) defining the coil surface. The corresponding mode numbers (m,n) are stored in xm_coil and xn_coil.')

  call cdf_define(ncid, vn_h, h, dimname=nthetanzeta_plasma_dim)
  call cdf_define(ncid, vn_RHS_B, RHS_B, dimname=num_basis_functions_dim)
  call cdf_define(ncid, vn_RHS_regularization, RHS_regularization, dimname=num_basis_functions_dim)

  call cdf_define(ncid, vn_lambda, lambda(1:Nlambda), dimname=nlambda_dim)
  call cdf_setatt(ncid, vn_lambda, 'Values of the regularization parameter that were used, in SI units (Tesla^2 meter^2 / Ampere^2)')

  call cdf_define(ncid, vn_chi2_B, chi2_B(1:Nlambda), dimname=nlambda_dim)
  call cdf_setatt(ncid, vn_chi2_B, 'Values of chi^2_B (the area integral over the plasma surface of |B_normal|^2) that resulted for each value of lambda, in SI units (Tesla^2 meter^2)')

  call cdf_define(ncid, vn_chi2_K, chi2_K(1:Nlambda), dimname=nlambda_dim)
  call cdf_setatt(ncid, vn_chi2_K, 'Values of chi^2_K (the area integral over the coil winding surface of current density squared) that resulted for each value of lambda, in SI units (Ampere^2)')

  call cdf_define(ncid, vn_chi2_Laplace_Beltrami, chi2_Laplace_Beltrami(1:Nlambda), dimname=nlambda_dim)
  call cdf_setatt(ncid, vn_chi2_Laplace_Beltrami, '')

  call cdf_define(ncid, vn_max_Bnormal, max_Bnormal(1:Nlambda), dimname=nlambda_dim)
  call cdf_setatt(ncid, vn_max_Bnormal, 'Maximum (over the plasma surface) magnetic field normal to the target plasma shape that resulted for each value of lambda, in Tesla.')

  call cdf_define(ncid, vn_max_K, max_K(1:Nlambda), dimname=nlambda_dim) ! We only write elements 1:Nlambda in case of a lambda search.
  call cdf_setatt(ncid, vn_max_K, 'Maximum (over the coil surface) current density that resulted for each value of lambda, in Amperes/meter.')
  if (sensitivity_option > 1) then
    call cdf_define(ncid, vn_xn_sensitivity, xn_sensitivity, dimname=nomega_coil_dim)
    call cdf_define(ncid, vn_xm_sensitivity, xm_sensitivity, dimname=nomega_coil_dim)
    call cdf_define(ncid, vn_omega_coil, omega_coil, dimname=nomega_coil_dim)
    call cdf_define(ncid, vn_dvolume_coildomega, dvolume_coildomega, dimname=nomega_coil_dim)
    call cdf_define(ncid, vn_darea_coildomega, darea_coildomega)
    call cdf_define(ncid, vn_dcoil_plasma_dist_mindomega, dcoil_plasma_dist_mindomega, dimname=nomega_coil_dim)
  end if
  if (trim(target_option)==target_option_max_K_lse) then
    call cdf_define(ncid, vn_max_K_lse, max_K_lse, dimname=nlambda_dim)
  end if
  if (trim(target_option)==target_option_lp_norm_K) then
    call cdf_define(ncid, vn_lp_norm_K, lp_norm_K, dimname=nlambda_dim)
  end if

  ! Arrays with dimension 2

  call cdf_define(ncid, vn_norm_normal_plasma,  norm_normal_plasma,  dimname=ntheta_nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_norm_normal_plasma, '|N|, where N = (d r / d zeta) cross (d r / d theta) is a non-unit-length normal vector ' // &
       'and r is the posiiton vector, for the plasma surface. This quantity is the Jacobian appearing in area integrals: ' // &
       'int d^2a = int dtheta int dzeta |N|. Units = meters^2.')

  call cdf_define(ncid, vn_norm_normal_coil,  norm_normal_coil,  dimname=ntheta_nzeta_coil_dim)
  call cdf_setatt(ncid, vn_norm_normal_coil, '|N|, where N = (d r / d zeta) cross (d r / d theta) is a non-unit-length normal vector ' // &
       'and r is the posiiton vector, for the coil surface. This quantity is the Jacobian appearing in area integrals: ' // &
       'int d^2a = int dtheta int dzeta |N|. Units = meters^2.')

  call cdf_define(ncid, vn_Bnormal_from_plasma_current, Bnormal_from_plasma_current, dimname=ntheta_nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_Bnormal_from_plasma_current, 'Contribution to the magnetic field normal to the plasma surface from currents inside the plasma. ' // &
       'This is the quantity named B_{normal}^{plasma} in the 2017 Nuclear Fusion paper. Units = Tesla.')

  call cdf_define(ncid, vn_Bnormal_from_net_coil_currents, Bnormal_from_net_coil_currents, dimname=ntheta_nzeta_plasma_dim)
  call cdf_setatt(ncid, vn_Bnormal_from_plasma_current, 'Contribution to the magnetic field normal to the plasma surface from the secular (nonperiodic) ' // &
       'part of the current potential. This is the quantity named B_{normal}^{GI} in the 2017 Nuclear Fusion paper. Units = Tesla.')

  if (save_level<1) then
     call cdf_define(ncid, vn_inductance, inductance, dimname=nthetanzeta_plasma_nthetanzeta_coil_dim)
  end if
  if (save_level<2) then
     call cdf_define(ncid, vn_g, g, dimname=nthetanzeta_plasma_basis_dim)
  end if
  !call cdf_define(ncid, vn_matrix_B, matrix_B, dimname=basis_basis_dim)
  !call cdf_define(ncid, vn_matrix_K, matrix_K, dimname=basis_basis_dim)
  call cdf_define(ncid, vn_single_valued_current_potential_mn, single_valued_current_potential_mn(:,1:Nlambda), &
       dimname=basis_nlambda_dim)

  if (sensitivity_option > 1 .and. exit_code == 0) then
    call cdf_define(ncid, vn_dchi2domega, dchi2domega(:,1:Nlambda),dimname=nomega_coil_nlambda_dim)
  end if
  if (sensitivity_option > 2 .and. exit_code == 0) then
    call cdf_define(ncid, vn_dchi2Bdomega, dchi2Bdomega(:,1:Nlambda),dimname=nomega_coil_nlambda_dim)
    call cdf_define(ncid, vn_dchi2Kdomega, dchi2Kdomega(:,1:Nlambda),dimname=nomega_coil_nlambda_dim)
  end if

  ! Arrays with dimension 3

  call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
  call cdf_setatt(ncid, vn_r_plasma, 'Position vector describing the plasma boundary surface, in Cartesian components, with units of meters. ' // &
       'The dimension of size 3 corresponds to the (x,y,z) components.')

  call cdf_define(ncid, vn_r_coil,  r_coil,  dimname=xyz_ntheta_nzetal_coil_dim)
  call cdf_setatt(ncid, vn_r_coil,  'Position vector describing the plasma boundary surface, in Cartesian components, with units of meters. ' // &
       'The dimension of size 3 corresponds to the (x,y,z) components.')

  if (save_level < 3) then
     call cdf_define(ncid, vn_drdtheta_plasma,  drdtheta_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_drdtheta_coil,  drdtheta_coil,  dimname=xyz_ntheta_nzetal_coil_dim)
     
     call cdf_define(ncid, vn_drdzeta_plasma,  drdzeta_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_drdzeta_coil,  drdzeta_coil,  dimname=xyz_ntheta_nzetal_coil_dim)

     call cdf_define(ncid, vn_normal_plasma,  normal_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)
     call cdf_define(ncid, vn_normal_coil,  normal_coil,  dimname=xyz_ntheta_nzetal_coil_dim)

  end if

  call cdf_define(ncid, vn_single_valued_current_potential_thetazeta, single_valued_current_potential_thetazeta(:,:,1:Nlambda), &
       dimname=ntheta_nzeta_coil_nlambda_dim)
  call cdf_setatt(ncid, vn_single_valued_current_potential_thetazeta, 'Periodic (single-valued) term in the current potential on the coil winding surface, in units of Amperes, ' // &
       'for each value of the regularization parameter lambda considered.')

  call cdf_define(ncid, vn_current_potential, current_potential(:,:,1:Nlambda), dimname=ntheta_nzeta_coil_nlambda_dim)
  call cdf_setatt(ncid, vn_current_potential, 'Total (multiple-valued) current potential on the coil winding surface, in units of Amperes, ' // &
       'for each value of the regularization parameter lambda considered.')

  call cdf_define(ncid, vn_Bnormal_total, Bnormal_total(:,:,1:Nlambda), dimname=ntheta_nzeta_plasma_nlambda_dim)
  call cdf_setatt(ncid, vn_Bnormal_total, 'Residual magnetic field normal to the plasma surface, in units of Tesla, ' // &
       'for each value of the regularization parameter lambda considered.')

  call cdf_define(ncid, vn_K2, K2(:,:,1:Nlambda), dimname=ntheta_nzeta_coil_nlambda_dim)
  call cdf_setatt(ncid, vn_K2, 'Squared current density on the coil winding surface, in units of Amperes^2/meter^2, ' // &
       'for each value of the regularization parameter lambda considered.')

  call cdf_define(ncid, vn_Laplace_Beltrami2, Laplace_Beltrami2(:,:,1:Nlambda), dimname=ntheta_nzeta_coil_nlambda_dim)
  call cdf_setatt(ncid, vn_Laplace_Beltrami2, '')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_geometry_option_plasma, geometry_option_plasma)
  call cdf_write(ncid, vn_geometry_option_coil, geometry_option_coil)
  call cdf_write(ncid, vn_ntheta_plasma, ntheta_plasma)
  call cdf_write(ncid, vn_nzeta_plasma, nzeta_plasma)
  call cdf_write(ncid, vn_nzetal_plasma, nzetal_plasma)
  call cdf_write(ncid, vn_ntheta_coil, ntheta_coil)
  call cdf_write(ncid, vn_nzeta_coil, nzeta_coil)
  call cdf_write(ncid, vn_nzetal_coil, nzetal_coil)
  call cdf_write(ncid, vn_a_plasma, a_plasma)
  call cdf_write(ncid, vn_a_coil, a_coil)
  call cdf_write(ncid, vn_R0_plasma, R0_plasma)
  call cdf_write(ncid, vn_R0_coil, R0_coil)
  call cdf_write(ncid, vn_mpol_potential, mpol_potential)
  call cdf_write(ncid, vn_ntor_potential, ntor_potential)
  call cdf_write(ncid, vn_mnmax_potential, mnmax_potential)
  call cdf_write(ncid, vn_mnmax_plasma, mnmax_plasma)
  call cdf_write(ncid, vn_mnmax_coil, mnmax_coil)
  call cdf_write(ncid, vn_num_basis_functions, num_basis_functions)
  call cdf_write(ncid, vn_symmetry_option, symmetry_option)
  call cdf_write(ncid, vn_area_plasma, area_plasma)
  call cdf_write(ncid, vn_area_coil, area_coil)
  call cdf_write(ncid, vn_volume_plasma, volume_plasma)
  call cdf_write(ncid, vn_volume_coil, volume_coil)
  call cdf_write(ncid, vn_net_poloidal_current_Amperes, net_poloidal_current_Amperes)
  call cdf_write(ncid, vn_net_toroidal_current_Amperes, net_toroidal_current_Amperes)
  call cdf_write(ncid, vn_curpol, curpol)
  call cdf_write(ncid, vn_helicity_ratio, helicity_ratio)
  call cdf_write(ncid, vn_nlambda, nlambda)
  call cdf_write(ncid, vn_total_time, total_time)
  call cdf_write(ncid, vn_exit_code, exit_code)
  if (general_option==4 .or. general_option==5) call cdf_write(ncid, vn_chi2_B_target, chi2_B_target)
  call cdf_write(ncid, vn_sensitivity_option, sensitivity_option)
  if (sensitivity_option > 1) then
    call cdf_write(ncid, vn_mmax_sensitivity, mmax_sensitivity)
    call cdf_write(ncid, vn_nmax_sensitivity, nmax_sensitivity)
    call cdf_write(ncid, vn_mnmax_sensitivity, mnmax_sensitivity)
    call cdf_write(ncid, vn_nomega_coil, nomega_coil)
    call cdf_write(ncid, vn_sensitivity_symmetry_option, sensitivity_symmetry_option)
    call cdf_write(ncid, vn_fixed_norm_sensitivity_option, fixed_norm_sensitivity_option)
  endif
  if (sensitivity_option > 1) then
    call cdf_write(ncid, vn_coil_plasma_dist_min, coil_plasma_dist_min)
    call cdf_write(ncid, vn_coil_plasma_dist_max, coil_plasma_dist_max)
    call cdf_write(ncid, vn_coil_plasma_dist_lse_p, coil_plasma_dist_lse_p)
    call cdf_write(ncid, vn_coil_plasma_dist_min_lse, coil_plasma_dist_min_lse)
    call cdf_write(ncid, vn_coil_plasma_dist_max_lse, coil_plasma_dist_max_lse)
  end if

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_theta_plasma, theta_plasma)
  call cdf_write(ncid, vn_zeta_plasma, zeta_plasma)
  call cdf_write(ncid, vn_zetal_plasma, zetal_plasma)
  call cdf_write(ncid, vn_theta_coil, theta_coil)
  call cdf_write(ncid, vn_zeta_coil, zeta_coil)
  call cdf_write(ncid, vn_zetal_coil, zetal_coil)
  call cdf_write(ncid, vn_xm_potential, xm_potential)
  call cdf_write(ncid, vn_xn_potential, xn_potential)
  call cdf_write(ncid, vn_xm_plasma, xm_plasma)
  call cdf_write(ncid, vn_xn_plasma, xn_plasma)
  call cdf_write(ncid, vn_xm_coil, xm_coil)
  call cdf_write(ncid, vn_xn_coil, xn_coil)
  call cdf_write(ncid, vn_rmnc_plasma, rmnc_plasma)
  if (lasym) then
     call cdf_write(ncid, vn_rmns_plasma, rmns_plasma)
     call cdf_write(ncid, vn_zmnc_plasma, zmnc_plasma)
  end if
  call cdf_write(ncid, vn_zmns_plasma, zmns_plasma)
  call cdf_write(ncid, vn_rmnc_coil, rmnc_coil)
  if (lasym) then
     call cdf_write(ncid, vn_rmns_coil, rmns_coil)
     call cdf_write(ncid, vn_zmnc_coil, zmnc_coil)
  end if
  call cdf_write(ncid, vn_zmns_coil, zmns_coil)
  call cdf_write(ncid, vn_h, h)
  call cdf_write(ncid, vn_RHS_B, RHS_B)
  call cdf_write(ncid, vn_RHS_regularization, RHS_regularization)
  call cdf_write(ncid, vn_lambda, lambda(1:Nlambda))
  call cdf_write(ncid, vn_chi2_B, chi2_B(1:Nlambda))
  call cdf_write(ncid, vn_chi2_K, chi2_K(1:Nlambda))
  call cdf_write(ncid, vn_chi2_Laplace_Beltrami, chi2_Laplace_Beltrami(1:Nlambda))
  call cdf_write(ncid, vn_max_Bnormal, max_Bnormal(1:Nlambda))
  call cdf_write(ncid, vn_max_K, max_K(1:Nlambda))
  if (sensitivity_option > 1) then
    call cdf_write(ncid, vn_xn_sensitivity, xn_sensitivity)
    call cdf_write(ncid, vn_xm_sensitivity, xm_sensitivity)
    call cdf_write(ncid, vn_omega_coil, omega_coil)
    call cdf_write(ncid, vn_dvolume_coildomega, dvolume_coildomega)
    call cdf_write(ncid, vn_darea_coildomega, darea_coildomega)
    call cdf_write(ncid, vn_dcoil_plasma_dist_mindomega, dcoil_plasma_dist_mindomega)
  end if

  if (trim(target_option)==target_option_max_K_lse .and. exit_code==0) then
    call cdf_write(ncid, vn_max_K_lse, max_K_lse)
  end if
  if (trim(target_option)==target_option_lp_norm_K .and. exit_code==0) then
    call cdf_write(ncid, vn_lp_norm_K, lp_norm_K)
  end if

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_norm_normal_plasma,  norm_normal_plasma)
  call cdf_write(ncid, vn_norm_normal_coil,  norm_normal_coil)
  call cdf_write(ncid, vn_Bnormal_from_plasma_current, Bnormal_from_plasma_current)
  call cdf_write(ncid, vn_Bnormal_from_net_coil_currents, Bnormal_from_net_coil_currents)
  if (save_level<1) then
     call cdf_write(ncid, vn_inductance, inductance)
  end if
  if (save_level<2) then
     call cdf_write(ncid, vn_g, g)
  end if
  !call cdf_write(ncid, vn_matrix_B, matrix_B)
  !call cdf_write(ncid, vn_matrix_K, matrix_K)
  call cdf_write(ncid, vn_single_valued_current_potential_mn, single_valued_current_potential_mn(:,1:Nlambda))
  if (sensitivity_option > 1 .and. exit_code == 0) then
    call cdf_write(ncid, vn_dchi2domega, dchi2domega(:,1:Nlambda))
  end if
  if (sensitivity_option > 2 .and. exit_code == 0) then
    call cdf_write(ncid, vn_dchi2Kdomega, dchi2Kdomega(:,1:Nlambda))
    call cdf_write(ncid, vn_dchi2Bdomega, dchi2Bdomega(:,1:Nlambda))
    call cdf_write(ncid, vn_dRMSKdomega, dRMSKdomega(:,1:Nlambda))
  end if

  ! Arrays with dimension 3

  call cdf_write(ncid, vn_r_plasma, r_plasma)
  call cdf_write(ncid, vn_r_coil, r_coil)

  if (save_level < 3) then
     call cdf_write(ncid, vn_drdtheta_plasma, drdtheta_plasma)
     call cdf_write(ncid, vn_drdtheta_coil, drdtheta_coil)

     call cdf_write(ncid, vn_drdzeta_plasma, drdzeta_plasma)
     call cdf_write(ncid, vn_drdzeta_coil, drdzeta_coil)

     call cdf_write(ncid, vn_normal_plasma, normal_plasma)
     call cdf_write(ncid, vn_normal_coil, normal_coil)

  end if

  call cdf_write(ncid, vn_single_valued_current_potential_thetazeta, single_valued_current_potential_thetazeta(:,:,1:Nlambda))
  call cdf_write(ncid, vn_current_potential, current_potential(:,:,1:Nlambda))
  call cdf_write(ncid, vn_Bnormal_total, Bnormal_total(:,:,1:Nlambda))
  call cdf_write(ncid, vn_K2, K2(:,:,1:Nlambda))
  call cdf_write(ncid, vn_Laplace_Beltrami2, Laplace_Beltrami2(:,:,1:Nlambda))

  ! Finish up:
  call cdf_close(ncid)

end subroutine regcoil_write_output
