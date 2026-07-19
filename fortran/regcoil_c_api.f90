! C-compatible API for the Python extension (Phase 7: stateless kernels).
!
! Every entry point below is a thin bind(C) wrapper around a pure Fortran
! kernel (regcoil_kernels_mod / regcoil_uniform_offset_surface_mod): it
! reshapes the caller's raw buffers via c_f_pointer, calls the kernel, and
! turns a nonzero `info` into the returned error code. No opaque handle, no
! module state -- see ADR-020. (The earlier Phase 4/5 opaque-handle API,
! regcoil_t-backed, is retired; the namelist-driven solve chain it wrapped is
! superseded per ADR-019.)

module regcoil_c_api
  use iso_c_binding
  use stel_kinds, only: dp
  use regcoil_kernels_mod, only: regcoil_build_inductance, regcoil_build_g_and_h
  use regcoil_uniform_offset_surface_mod, only: regcoil_uniform_offset_surface

  implicit none

contains

  function regcoil_c_build_inductance( &
       ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, &
       r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil, &
       net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil, &
       inductance, h) result(ierr) bind(C, name="regcoil_c_build_inductance")

    integer(c_int), value, intent(in) :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp
    type(c_ptr), value, intent(in) :: r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil
    real(c_double), value, intent(in) :: net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil
    type(c_ptr), value, intent(in) :: inductance, h
    integer(c_int) :: ierr

    real(c_double), pointer :: f_r_plasma(:,:,:), f_normal_plasma(:,:,:)
    real(c_double), pointer :: f_r_coil(:,:,:), f_normal_coil(:,:,:), f_drdtheta_coil(:,:,:), f_drdzeta_coil(:,:,:)
    real(c_double), pointer :: f_inductance(:,:), f_h(:)
    integer :: nt_p, nz_p, nt_c, nz_c, nfp_i, nzetal_c, info

    nt_p = int(ntheta_plasma); nz_p = int(nzeta_plasma)
    nt_c = int(ntheta_coil); nz_c = int(nzeta_coil); nfp_i = int(nfp)
    nzetal_c = nz_c * nfp_i

    call c_f_pointer(r_plasma, f_r_plasma, [3, nt_p, nz_p])
    call c_f_pointer(normal_plasma, f_normal_plasma, [3, nt_p, nz_p])
    call c_f_pointer(r_coil, f_r_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(normal_coil, f_normal_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(drdtheta_coil, f_drdtheta_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(drdzeta_coil, f_drdzeta_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(inductance, f_inductance, [nt_p*nz_p, nt_c*nz_c])
    call c_f_pointer(h, f_h, [nt_p*nz_p])

    call regcoil_build_inductance(nt_p, nz_p, nt_c, nz_c, nfp_i, &
         f_r_plasma, f_normal_plasma, f_r_coil, f_normal_coil, f_drdtheta_coil, f_drdzeta_coil, &
         real(net_poloidal_current, kind=dp), real(net_toroidal_current, kind=dp), &
         real(dtheta_coil, kind=dp), real(dzeta_coil, kind=dp), &
         f_inductance, f_h, info)
    ierr = int(info, kind=c_int)
  end function regcoil_c_build_inductance

  function regcoil_c_build_g_and_h( &
       ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf, &
       r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil, &
       basis_functions, net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil, &
       g, h) result(ierr) bind(C, name="regcoil_c_build_g_and_h")

    integer(c_int), value, intent(in) :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf
    type(c_ptr), value, intent(in) :: r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil
    type(c_ptr), value, intent(in) :: basis_functions
    real(c_double), value, intent(in) :: net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil
    type(c_ptr), value, intent(in) :: g, h
    integer(c_int) :: ierr

    real(c_double), pointer :: f_r_plasma(:,:,:), f_normal_plasma(:,:,:)
    real(c_double), pointer :: f_r_coil(:,:,:), f_normal_coil(:,:,:), f_drdtheta_coil(:,:,:), f_drdzeta_coil(:,:,:)
    real(c_double), pointer :: f_basis_functions(:,:), f_g(:,:), f_h(:)
    integer :: nt_p, nz_p, nt_c, nz_c, nfp_i, nbf_i, nzetal_c, info

    nt_p = int(ntheta_plasma); nz_p = int(nzeta_plasma)
    nt_c = int(ntheta_coil); nz_c = int(nzeta_coil); nfp_i = int(nfp); nbf_i = int(nbf)
    nzetal_c = nz_c * nfp_i

    call c_f_pointer(r_plasma, f_r_plasma, [3, nt_p, nz_p])
    call c_f_pointer(normal_plasma, f_normal_plasma, [3, nt_p, nz_p])
    call c_f_pointer(r_coil, f_r_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(normal_coil, f_normal_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(drdtheta_coil, f_drdtheta_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(drdzeta_coil, f_drdzeta_coil, [3, nt_c, nzetal_c])
    call c_f_pointer(basis_functions, f_basis_functions, [nt_c*nz_c, nbf_i])
    call c_f_pointer(g, f_g, [nt_p*nz_p, nbf_i])
    call c_f_pointer(h, f_h, [nt_p*nz_p])

    call regcoil_build_g_and_h(nt_p, nz_p, nt_c, nz_c, nfp_i, nbf_i, &
         f_r_plasma, f_normal_plasma, f_r_coil, f_normal_coil, f_drdtheta_coil, f_drdzeta_coil, &
         f_basis_functions, real(net_poloidal_current, kind=dp), real(net_toroidal_current, kind=dp), &
         real(dtheta_coil, kind=dp), real(dzeta_coil, kind=dp), &
         f_g, f_h, info)
    ierr = int(info, kind=c_int)
  end function regcoil_c_build_g_and_h

  function regcoil_c_uniform_offset_surface( &
       mnmax_in, xm_in, xn_in, rmnc_in, rmns_in, zmnc_in, zmns_in, lasym_flag, nfp, &
       separation, mpol_out, ntor_out, ntheta_transform, nzeta_transform, tol, &
       mnmax_out, xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out) result(ierr) &
       bind(C, name="regcoil_c_uniform_offset_surface")

    integer(c_int), value, intent(in) :: mnmax_in
    type(c_ptr), value, intent(in) :: xm_in, xn_in
    type(c_ptr), value, intent(in) :: rmnc_in, rmns_in, zmnc_in, zmns_in
    integer(c_int), value, intent(in) :: lasym_flag, nfp
    real(c_double), value, intent(in) :: separation
    integer(c_int), value, intent(in) :: mpol_out, ntor_out, ntheta_transform, nzeta_transform
    real(c_double), value, intent(in) :: tol
    integer(c_int), value, intent(in) :: mnmax_out
    type(c_ptr), value, intent(in) :: xm_out, xn_out
    type(c_ptr), value, intent(in) :: rmnc_out, rmns_out, zmnc_out, zmns_out
    integer(c_int) :: ierr

    integer(c_int), pointer :: c_xm_in(:), c_xn_in(:), c_xm_out(:), c_xn_out(:)
    real(c_double), pointer :: f_rmnc_in(:), f_rmns_in(:), f_zmnc_in(:), f_zmns_in(:)
    real(c_double), pointer :: f_rmnc_out(:), f_rmns_out(:), f_zmnc_out(:), f_zmns_out(:)
    integer, allocatable :: f_xm_in(:), f_xn_in(:), f_xm_out(:), f_xn_out(:)
    integer :: mnmax_in_i, mnmax_out_i, info

    mnmax_in_i = int(mnmax_in)
    mnmax_out_i = int(mnmax_out)

    call c_f_pointer(xm_in, c_xm_in, [mnmax_in_i])
    call c_f_pointer(xn_in, c_xn_in, [mnmax_in_i])
    call c_f_pointer(rmnc_in, f_rmnc_in, [mnmax_in_i])
    call c_f_pointer(rmns_in, f_rmns_in, [mnmax_in_i])
    call c_f_pointer(zmnc_in, f_zmnc_in, [mnmax_in_i])
    call c_f_pointer(zmns_in, f_zmns_in, [mnmax_in_i])
    call c_f_pointer(xm_out, c_xm_out, [mnmax_out_i])
    call c_f_pointer(xn_out, c_xn_out, [mnmax_out_i])
    call c_f_pointer(rmnc_out, f_rmnc_out, [mnmax_out_i])
    call c_f_pointer(rmns_out, f_rmns_out, [mnmax_out_i])
    call c_f_pointer(zmnc_out, f_zmnc_out, [mnmax_out_i])
    call c_f_pointer(zmns_out, f_zmns_out, [mnmax_out_i])

    allocate(f_xm_in(mnmax_in_i), f_xn_in(mnmax_in_i))
    f_xm_in = int(c_xm_in)
    f_xn_in = int(c_xn_in)
    allocate(f_xm_out(mnmax_out_i), f_xn_out(mnmax_out_i))

    call regcoil_uniform_offset_surface(mnmax_in_i, f_xm_in, f_xn_in, &
         f_rmnc_in, f_rmns_in, f_zmnc_in, f_zmns_in, (lasym_flag /= 0), int(nfp), &
         real(separation, kind=dp), int(mpol_out), int(ntor_out), &
         int(ntheta_transform), int(nzeta_transform), real(tol, kind=dp), &
         mnmax_out_i, f_xm_out, f_xn_out, &
         f_rmnc_out, f_rmns_out, f_zmnc_out, f_zmns_out, info)

    c_xm_out = int(f_xm_out, kind=c_int)
    c_xn_out = int(f_xn_out, kind=c_int)

    deallocate(f_xm_in, f_xn_in, f_xm_out, f_xn_out)

    ierr = int(info, kind=c_int)
  end function regcoil_c_uniform_offset_surface

end module regcoil_c_api
