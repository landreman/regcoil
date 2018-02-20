module regcoil_init_coil_surface

contains
 
subroutine init_coil_surface(lscreen_optin)

  use regcoil_variables
  use init_surface_mod

  implicit none

  integer :: tic, toc, countrate
  integer :: which_surface

  ! variables to handle printing to the screen
  logical, optional :: lscreen_optin
  logical :: lscreen

  if(present(lscreen_optin)) then 
    lscreen = lscreen_optin
  else
    lscreen = .true.
  endif

  call system_clock(tic,countrate)
  !print "(a,l8)", "lscreen=",lscreen
  if (lscreen) print *,"Initializing coil surface."
  which_surface = 1
  call init_surface(ntheta_coil, nzeta_coil, nzetal_coil, theta_coil, zeta_coil, zetal_coil, &
       r_coil, drdtheta_coil, drdzeta_coil, &
       normal_coil, norm_normal_coil, area_coil, &
       geometry_option_coil, R0_coil, a_coil, separation, dtheta_coil, dzeta_coil, &
       nescin_filename, which_surface, lscreen)
  call system_clock(toc)
  if (lscreen) print *,"Done initializing coil surface. Took ",real(toc-tic)/countrate," sec."

!!$  call system_clock(tic)
!!$  print *,"Initializing outer surface."
!!$  which_surface = 2
!!$  call init_surface(nu_outer, nv_outer, nvl_outer, u_outer, v_outer, vl_outer, &
!!$       r_outer, drdu_outer, drdv_outer, &
!!$       d2rdu2_middle, d2rdudv_middle, d2rdv2_middle, & ! The 2nd derivatives of r are only set for the middle surface, so just re-use these arrays when initializing the middle surface.
!!$       normal_outer, norm_normal_outer, area_outer, &
!!$       geometry_option_outer, R0_outer, a_outer, separation_outer, du_outer, dv_outer, &
!!$       nescin_filename_outer, which_surface)
!!$  call system_clock(toc)
!!$  print *,"Done initializing outer surface. Took ",real(toc-tic)/countrate," sec."
  

end subroutine init_coil_surface

end module regcoil_init_coil_surface
