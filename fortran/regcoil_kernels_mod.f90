module regcoil_kernels_mod

  ! Stateless Fortran kernels for the Python extension.
  !
  ! Pure subroutines: intent(in)/intent(out) arrays with explicit extents,
  ! no module state, no `stop` (callers get an `info` code instead), and no
  ! use of `regcoil_variables`. Both routines are safe to call concurrently
  ! from multiple threads with different problem sizes since all working
  ! storage is local (stack/automatic), never module-level.

  use regcoil_kinds_mod, only: dp, pi, mu0

  implicit none

  private
  public :: regcoil_build_inductance, regcoil_build_g_and_h

contains

  ! Debug/regression entry point: returns the full (ntheta_plasma*nzeta_plasma)
  ! x (ntheta_coil*nzeta_coil) inductance matrix (already scaled by mu0/(4*pi)),
  ! plus h (the Bnormal contribution from the net coil currents, scaled by
  ! dtheta_coil*dzeta_coil*mu0/(8*pi**2), matching the legacy convention).
  subroutine regcoil_build_inductance( &
       ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, &
       r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil, &
       net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil, &
       inductance, h, info)

    integer,  intent(in)  :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp
    real(dp), intent(in)  :: r_plasma(3, ntheta_plasma, nzeta_plasma)
    real(dp), intent(in)  :: normal_plasma(3, ntheta_plasma, nzeta_plasma)
    real(dp), intent(in)  :: r_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: normal_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: drdtheta_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: drdzeta_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: net_poloidal_current, net_toroidal_current
    real(dp), intent(in)  :: dtheta_coil, dzeta_coil
    real(dp), intent(out) :: inductance(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil)
    real(dp), intent(out) :: h(ntheta_plasma*nzeta_plasma)
    integer,  intent(out) :: info

    real(dp), allocatable :: factor_for_h(:,:,:)
    integer :: izeta_plasma

    info = 0
    if (ntheta_plasma < 1 .or. nzeta_plasma < 1 .or. ntheta_coil < 1 .or. nzeta_coil < 1 .or. nfp < 1) then
       info = 1
       return
    end if

    allocate(factor_for_h(3, ntheta_coil, nzeta_coil*nfp))
    factor_for_h = net_poloidal_current * drdtheta_coil - net_toroidal_current * drdzeta_coil

    !$OMP PARALLEL DO SCHEDULE(static)
    do izeta_plasma = 1, nzeta_plasma
       call fill_plasma_row_block(izeta_plasma)
    end do
    !$OMP END PARALLEL DO

    deallocate(factor_for_h)

  contains

    subroutine fill_plasma_row_block(izeta_plasma)
      integer, intent(in) :: izeta_plasma
      integer :: itheta_plasma, itheta_coil, izeta_coil, izetal_coil, l_coil
      integer :: index_plasma, index_coil
      real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv
      real(dp) :: dx_norm2, dx_norm3, dy_norm1, dy_norm3, dz_norm1, dz_norm2, this_h

      do itheta_plasma = 1, ntheta_plasma
         index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
         x = r_plasma(1,itheta_plasma,izeta_plasma)
         y = r_plasma(2,itheta_plasma,izeta_plasma)
         z = r_plasma(3,itheta_plasma,izeta_plasma)
         h(index_plasma) = 0
         do izeta_coil = 1, nzeta_coil
            do itheta_coil = 1, ntheta_coil
               index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
               inductance(index_plasma, index_coil) = 0
               do l_coil = 0, nfp-1
                  izetal_coil = izeta_coil + l_coil*nzeta_coil
                  dx = x - r_coil(1,itheta_coil,izetal_coil)
                  dy = y - r_coil(2,itheta_coil,izetal_coil)
                  dz = z - r_coil(3,itheta_coil,izetal_coil)

                  dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                  dr32inv = dr2inv*sqrt(dr2inv)

                  dy_norm3 = dy * normal_plasma(3,itheta_plasma,izeta_plasma)
                  dz_norm1 = dz * normal_plasma(1,itheta_plasma,izeta_plasma)
                  dx_norm2 = dx * normal_plasma(2,itheta_plasma,izeta_plasma)
                  dy_norm1 = dy * normal_plasma(1,itheta_plasma,izeta_plasma)
                  dz_norm2 = dz * normal_plasma(2,itheta_plasma,izeta_plasma)
                  dx_norm3 = dx * normal_plasma(3,itheta_plasma,izeta_plasma)

                  this_h = (factor_for_h(1,itheta_coil,izetal_coil) * dy_norm3 + &
                       factor_for_h(2,itheta_coil,izetal_coil) * dz_norm1 + &
                       factor_for_h(3,itheta_coil,izetal_coil) * dx_norm2  &
                       - factor_for_h(3,itheta_coil,izetal_coil) * dy_norm1 &
                       - factor_for_h(1,itheta_coil,izetal_coil) * dz_norm2 &
                       - factor_for_h(2,itheta_coil,izetal_coil) * dx_norm3 ) * dr32inv

                  inductance(index_plasma, index_coil) = inductance(index_plasma, index_coil) + &
                       (normal_plasma(1,itheta_plasma,izeta_plasma)*normal_coil(1,itheta_coil,izetal_coil) &
                       +normal_plasma(2,itheta_plasma,izeta_plasma)*normal_coil(2,itheta_coil,izetal_coil) &
                       +normal_plasma(3,itheta_plasma,izeta_plasma)*normal_coil(3,itheta_coil,izetal_coil) &
                       - (3*dr2inv) * &
                       (normal_plasma(1,itheta_plasma,izeta_plasma)*dx &
                       + normal_plasma(2,itheta_plasma,izeta_plasma)*dy &
                       + normal_plasma(3,itheta_plasma,izeta_plasma)*dz) * &
                       (normal_coil(1,itheta_coil,izetal_coil)*dx &
                       +normal_coil(2,itheta_coil,izetal_coil)*dy &
                       +normal_coil(3,itheta_coil,izetal_coil)*dz)) * dr32inv

                  h(index_plasma) = h(index_plasma) + this_h
               end do
            end do
         end do
      end do

      h((izeta_plasma-1)*ntheta_plasma+1 : izeta_plasma*ntheta_plasma) = &
           h((izeta_plasma-1)*ntheta_plasma+1 : izeta_plasma*ntheta_plasma) * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi))
      inductance((izeta_plasma-1)*ntheta_plasma+1 : izeta_plasma*ntheta_plasma, :) = &
           inductance((izeta_plasma-1)*ntheta_plasma+1 : izeta_plasma*ntheta_plasma, :) * (mu0/(4*pi))

    end subroutine fill_plasma_row_block

  end subroutine regcoil_build_inductance

  ! Fused kernel: computes g = dtheta_coil*dzeta_coil * (inductance @ basis_functions)
  ! and h, without ever materializing the full inductance matrix. Blocked over
  ! plasma rows: one field-period-zeta's worth of plasma points (ntheta_plasma
  ! rows) is accumulated into a small per-thread buffer, then contracted
  ! against basis_functions with DGEMM (faster than the `matmul` intrinsic for
  ! this contraction), keeping peak memory at
  ! O(ntheta_plasma * ntheta_coil*nzeta_coil) rather than O(nplasma * ncoil).
  subroutine regcoil_build_g_and_h( &
       ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf, &
       r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil, &
       basis_functions, net_poloidal_current, net_toroidal_current, &
       dtheta_coil, dzeta_coil, g, h, info)

    integer,  intent(in)  :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf
    real(dp), intent(in)  :: r_plasma(3, ntheta_plasma, nzeta_plasma)
    real(dp), intent(in)  :: normal_plasma(3, ntheta_plasma, nzeta_plasma)
    real(dp), intent(in)  :: r_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: normal_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: drdtheta_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: drdzeta_coil(3, ntheta_coil, nzeta_coil*nfp)
    real(dp), intent(in)  :: basis_functions(ntheta_coil*nzeta_coil, nbf)
    real(dp), intent(in)  :: net_poloidal_current, net_toroidal_current
    real(dp), intent(in)  :: dtheta_coil, dzeta_coil
    real(dp), intent(out) :: g(ntheta_plasma*nzeta_plasma, nbf)
    real(dp), intent(out) :: h(ntheta_plasma*nzeta_plasma)
    integer,  intent(out) :: info

    real(dp), allocatable :: factor_for_h(:,:,:)
    integer :: izeta_plasma, ncoil

    info = 0
    if (ntheta_plasma < 1 .or. nzeta_plasma < 1 .or. ntheta_coil < 1 .or. nzeta_coil < 1 &
         .or. nfp < 1 .or. nbf < 1) then
       info = 1
       return
    end if

    ncoil = ntheta_coil*nzeta_coil

    allocate(factor_for_h(3, ntheta_coil, nzeta_coil*nfp))
    factor_for_h = net_poloidal_current * drdtheta_coil - net_toroidal_current * drdzeta_coil

    !$OMP PARALLEL DO SCHEDULE(static)
    do izeta_plasma = 1, nzeta_plasma
       call process_plasma_row_block(izeta_plasma)
    end do
    !$OMP END PARALLEL DO

    deallocate(factor_for_h)

  contains

    subroutine process_plasma_row_block(izeta_plasma)
      integer, intent(in) :: izeta_plasma
      integer :: itheta_plasma, itheta_coil, izeta_coil, izetal_coil, l_coil
      integer :: index_coil, row_start, row_end
      real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv
      real(dp) :: dx_norm2, dx_norm3, dy_norm1, dy_norm3, dz_norm1, dz_norm2, this_h
      real(dp) :: inductance_block(ntheta_plasma, ncoil)
      real(dp) :: h_block(ntheta_plasma)
      real(dp) :: g_block(ntheta_plasma, nbf)
      real(dp) :: blas_alpha, blas_beta

      inductance_block = 0
      h_block = 0

      do itheta_plasma = 1, ntheta_plasma
         x = r_plasma(1,itheta_plasma,izeta_plasma)
         y = r_plasma(2,itheta_plasma,izeta_plasma)
         z = r_plasma(3,itheta_plasma,izeta_plasma)
         do izeta_coil = 1, nzeta_coil
            do itheta_coil = 1, ntheta_coil
               index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
               do l_coil = 0, nfp-1
                  izetal_coil = izeta_coil + l_coil*nzeta_coil
                  dx = x - r_coil(1,itheta_coil,izetal_coil)
                  dy = y - r_coil(2,itheta_coil,izetal_coil)
                  dz = z - r_coil(3,itheta_coil,izetal_coil)

                  dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                  dr32inv = dr2inv*sqrt(dr2inv)

                  dy_norm3 = dy * normal_plasma(3,itheta_plasma,izeta_plasma)
                  dz_norm1 = dz * normal_plasma(1,itheta_plasma,izeta_plasma)
                  dx_norm2 = dx * normal_plasma(2,itheta_plasma,izeta_plasma)
                  dy_norm1 = dy * normal_plasma(1,itheta_plasma,izeta_plasma)
                  dz_norm2 = dz * normal_plasma(2,itheta_plasma,izeta_plasma)
                  dx_norm3 = dx * normal_plasma(3,itheta_plasma,izeta_plasma)

                  this_h = (factor_for_h(1,itheta_coil,izetal_coil) * dy_norm3 + &
                       factor_for_h(2,itheta_coil,izetal_coil) * dz_norm1 + &
                       factor_for_h(3,itheta_coil,izetal_coil) * dx_norm2  &
                       - factor_for_h(3,itheta_coil,izetal_coil) * dy_norm1 &
                       - factor_for_h(1,itheta_coil,izetal_coil) * dz_norm2 &
                       - factor_for_h(2,itheta_coil,izetal_coil) * dx_norm3 ) * dr32inv

                  inductance_block(itheta_plasma, index_coil) = inductance_block(itheta_plasma, index_coil) + &
                       (normal_plasma(1,itheta_plasma,izeta_plasma)*normal_coil(1,itheta_coil,izetal_coil) &
                       +normal_plasma(2,itheta_plasma,izeta_plasma)*normal_coil(2,itheta_coil,izetal_coil) &
                       +normal_plasma(3,itheta_plasma,izeta_plasma)*normal_coil(3,itheta_coil,izetal_coil) &
                       - (3*dr2inv) * &
                       (normal_plasma(1,itheta_plasma,izeta_plasma)*dx &
                       + normal_plasma(2,itheta_plasma,izeta_plasma)*dy &
                       + normal_plasma(3,itheta_plasma,izeta_plasma)*dz) * &
                       (normal_coil(1,itheta_coil,izetal_coil)*dx &
                       +normal_coil(2,itheta_coil,izetal_coil)*dy &
                       +normal_coil(3,itheta_coil,izetal_coil)*dz)) * dr32inv

                  h_block(itheta_plasma) = h_block(itheta_plasma) + this_h
               end do
            end do
         end do
      end do

      row_start = (izeta_plasma-1)*ntheta_plasma + 1
      row_end   = izeta_plasma*ntheta_plasma

      h(row_start:row_end) = h_block * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi))

      ! DGEMM on one ntheta_plasma-row chunk at a time: the full inductance
      ! matrix is never materialized, and DGEMM (threaded BLAS) contracts
      ! each chunk against basis_functions faster than the `matmul`
      ! intrinsic would.
      blas_alpha = dtheta_coil*dzeta_coil*mu0/(4*pi)
      blas_beta = 0
      g_block = 0
      call DGEMM('N', 'N', ntheta_plasma, nbf, ncoil, blas_alpha, &
           inductance_block, ntheta_plasma, basis_functions, ncoil, blas_beta, g_block, ntheta_plasma)

      g(row_start:row_end, :) = g_block

    end subroutine process_plasma_row_block

  end subroutine regcoil_build_g_and_h

end module regcoil_kernels_mod
