subroutine regcoil_diagnostics(prob, ilambda)

  use regcoil_variables, only: regcoil_t, target_option_lp_norm_K
  use stel_constants
  use stel_kinds

  implicit none


  type(regcoil_t), intent(inout) :: prob
  integer, intent(in) :: ilambda
  integer :: tic, toc, countrate
  real(dp) :: factor_theta, factor_zeta
  integer :: itheta, izeta

  associate ( &
       ntheta_plasma => prob%plasma%ntheta_plasma, &
       nzeta_plasma => prob%plasma%nzeta_plasma, &
       Bnormal_from_plasma_current => prob%plasma%Bnormal_from_plasma_current, &
       norm_normal_plasma => prob%plasma%norm_normal_plasma, &
       dtheta_plasma => prob%plasma%dtheta_plasma, &
       dzeta_plasma => prob%plasma%dzeta_plasma, &
       nfp => prob%plasma%nfp, &
       ntheta_coil => prob%coil%ntheta_coil, &
       nzeta_coil => prob%coil%nzeta_coil, &
       theta_coil => prob%coil%theta_coil, &
       zeta_coil => prob%coil%zeta_coil, &
       Bnormal_from_net_coil_currents => prob%coil%Bnormal_from_net_coil_currents, &
       norm_normal_coil => prob%coil%norm_normal_coil, &
       dtheta_coil => prob%coil%dtheta_coil, &
       dzeta_coil => prob%coil%dzeta_coil, &
       area_coil => prob%coil%area_coil, &
       verbose => prob%input%verbose, &
       net_poloidal_current_Amperes => prob%input%net_poloidal_current_Amperes, &
       net_toroidal_current_Amperes => prob%input%net_toroidal_current_Amperes, &
       target_option => prob%input%target_option, &
       target_option_p => prob%input%target_option_p, &
       single_valued_current_potential_mn => prob%output%single_valued_current_potential_mn, &
       single_valued_current_potential_thetazeta => prob%output%single_valued_current_potential_thetazeta, &
       current_potential => prob%output%current_potential, &
       Bnormal_total => prob%output%Bnormal_total, &
       K2 => prob%output%K2, &
       Laplace_Beltrami2 => prob%output%Laplace_Beltrami2, &
       chi2_B => prob%output%chi2_B, &
       chi2_K => prob%output%chi2_K, &
       max_Bnormal => prob%output%max_Bnormal, &
       max_K => prob%output%max_K, &
       chi2_Laplace_Beltrami => prob%output%chi2_Laplace_Beltrami, &
       lp_norm_K => prob%output%lp_norm_K, &
       max_K_lse => prob%output%max_K_lse, &
       g => prob%work%g, &
       f_x => prob%work%f_x, &
       f_y => prob%work%f_y, &
       f_z => prob%work%f_z, &
       f_Laplace_Beltrami => prob%work%f_Laplace_Beltrami, &
       d_x => prob%work%d_x, &
       d_y => prob%work%d_y, &
       d_z => prob%work%d_z, &
       d_Laplace_Beltrami => prob%work%d_Laplace_Beltrami, &
       basis_functions => prob%work%basis_functions, &
       solution => prob%work%solution, &
       this_current_potential => prob%work%this_current_potential, &
       KDifference_x => prob%work%KDifference_x, &
       KDifference_y => prob%work%KDifference_y, &
       KDifference_z => prob%work%KDifference_z, &
       KDifference_Laplace_Beltrami => prob%work%KDifference_Laplace_Beltrami, &
       this_K2_times_N => prob%work%this_K2_times_N, &
       this_Laplace_Beltrami2_times_N => prob%work%this_Laplace_Beltrami2_times_N &
       )

  call system_clock(tic,countrate)

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  single_valued_current_potential_mn(:,ilambda) = solution
  this_current_potential = reshape(matmul(basis_functions, solution), (/ ntheta_coil, nzeta_coil /)) ! Could use BLAS2 for this line for extra speed.
  single_valued_current_potential_thetazeta(:,:,ilambda) = this_current_potential
  do izeta = 1,nzeta_coil
     do itheta = 1,ntheta_coil
        this_current_potential(itheta,izeta) = this_current_potential(itheta,izeta) &
             + factor_zeta*zeta_coil(izeta) + factor_theta*theta_coil(itheta)
     end do
  end do
  current_potential(:,:,ilambda) = this_current_potential
  
  KDifference_x = d_x - matmul(f_x, solution)
  KDifference_y = d_y - matmul(f_y, solution)
  KDifference_z = d_z - matmul(f_z, solution)
  KDifference_Laplace_Beltrami = d_Laplace_Beltrami - matmul(f_Laplace_Beltrami, solution)
  this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
       / norm_normal_coil
  this_Laplace_Beltrami2_times_N = reshape(KDifference_Laplace_Beltrami*KDifference_Laplace_Beltrami, (/ ntheta_coil, nzeta_coil /)) &
       / norm_normal_coil
  chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)
  K2(:,:,ilambda) = this_K2_times_N / norm_normal_coil
  chi2_Laplace_Beltrami(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_Laplace_Beltrami2_times_N)
  Laplace_Beltrami2(:,:,ilambda) = this_Laplace_Beltrami2_times_N / norm_normal_coil
  
  Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
       + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents
  
  max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
  max_K(ilambda) = sqrt(maxval(K2(:,:,ilambda)))
  
  chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
       * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)

  if (trim(target_option)=='max_K_lse') then
    max_K_lse(ilambda) = (1/target_option_p)*log(sum(norm_normal_coil*dtheta_coil*dzeta_coil*nfp &
        /area_coil*exp(target_option_p &
        *(K2(:,:,ilambda)**(0.5) - max_K(ilambda))))) + max_K(ilambda)
  else if (trim(target_option)==target_option_lp_norm_K) then
    lp_norm_K(ilambda) = (dtheta_coil*dzeta_coil*nfp*sum(norm_normal_coil* &
        K2(:,:,ilambda)**(target_option_p/2.0))/area_coil)**(1.0/target_option_p)
  end if
  
  call system_clock(toc)
  if (verbose) print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
  if (verbose) print "(a,es10.3,a,es10.3,a,es10.3)","   chi2_B:",chi2_B(ilambda),",  chi2_K:",chi2_K(ilambda),",  chi2_Laplace_Beltrami:",chi2_Laplace_Beltrami(ilambda)
  if (verbose) print "(a,es10.3,a,es10.3,a,es10.3)","   max(B_n):",max_Bnormal(ilambda),",  max(K):",max_K(ilambda),",  rms K:",sqrt(chi2_K(ilambda)/area_coil)
  if (trim(target_option)=='max_K_lse') then
    if (verbose) print "(a,es10.3)","    max_K_lse:",max_K_lse(ilambda)
  end if
  if (trim(target_option)=='lp_norm_K') then
    if (verbose) print "(a,es10.3)","    lp_norm_K:", lp_norm_K(ilambda)
  end if


  end associate
end subroutine regcoil_diagnostics
