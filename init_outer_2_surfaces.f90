subroutine init_outer_2_surfaces

  use global_variables
  use init_surface_mod

  implicit none

  integer :: tic, toc, countrate
  integer :: which_surface

  call system_clock(tic,countrate)
  print *,"Initializing middle surface."
  which_surface = 1
  call init_surface(nu_middle, nv_middle, nvl_middle, u_middle, v_middle, vl_middle, &
       r_middle, drdu_middle, drdv_middle, &
       d2rdu2_middle, d2rdudv_middle, d2rdv2_middle, & 
       normal_middle, norm_normal_middle, area_middle, &
       geometry_option_middle, R0_middle, a_middle, separation_middle, du_middle, dv_middle, &
       nescin_filename_middle, which_surface)
  call system_clock(toc)
  print *,"Done initializing middle surface. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Initializing outer surface."
  which_surface = 2
  call init_surface(nu_outer, nv_outer, nvl_outer, u_outer, v_outer, vl_outer, &
       r_outer, drdu_outer, drdv_outer, &
       d2rdu2_middle, d2rdudv_middle, d2rdv2_middle, & ! The 2nd derivatives of r are only set for the middle surface, so just re-use these arrays when initializing the middle surface.
       normal_outer, norm_normal_outer, area_outer, &
       geometry_option_outer, R0_outer, a_outer, separation_outer, du_outer, dv_outer, &
       nescin_filename_outer, which_surface)
  call system_clock(toc)
  print *,"Done initializing outer surface. Took ",real(toc-tic)/countrate," sec."
  

end subroutine init_outer_2_surfaces
