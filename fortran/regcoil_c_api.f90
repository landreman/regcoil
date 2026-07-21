! C-compatible API for the Python extension (Phase 7: stateless kernels).
!
! Every entry point below is a thin bind(C) wrapper around a pure Fortran
! kernel (regcoil_kernels_mod): it
! reshapes the caller's raw buffers via c_f_pointer, calls the kernel, and
! turns a nonzero `info` into the returned error code.

module regcoil_c_api
  use iso_c_binding
  use regcoil_kinds_mod, only: dp
     use omp_lib, only: omp_get_max_threads
  use regcoil_kernels_mod, only: regcoil_build_inductance, regcoil_build_g_and_h

  implicit none

contains

     function regcoil_c_omp_max_threads() result(nthreads) bind(C, name="regcoil_c_omp_max_threads")
          integer(c_int) :: nthreads

          nthreads = int(omp_get_max_threads(), kind=c_int)
     end function regcoil_c_omp_max_threads

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


end module regcoil_c_api
