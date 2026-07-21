module regcoil_uniform_offset_surface_mod

  ! Stateless offset-surface kernel (Phase 7 / ADR-020): given a plasma
  ! surface's Fourier coefficients, returns the Fourier coefficients of the
  ! surface offset uniformly outward by `separation` along the plasma normal
  ! (legacy geometry_option_coil=2, without the constant-arclength theta
  ! iteration). No `regcoil_t`/module state crosses this boundary; every
  ! working array is local, so concurrent calls with different problem sizes
  ! do not interfere.

  use regcoil_kinds_mod, only: dp, pi, twopi

  implicit none

  private
  public :: regcoil_uniform_offset_surface

contains

  subroutine regcoil_uniform_offset_surface( &
       mnmax_in, xm_in, xn_in, rmnc_in, rmns_in, zmnc_in, zmns_in, lasym, nfp, &
       separation, mpol_out, ntor_out, ntheta_transform, nzeta_transform, tol, &
       mnmax_out, xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out, info)

    integer,  intent(in)  :: mnmax_in
    integer,  intent(in)  :: xm_in(mnmax_in), xn_in(mnmax_in)
    real(dp), intent(in)  :: rmnc_in(mnmax_in), rmns_in(mnmax_in), zmnc_in(mnmax_in), zmns_in(mnmax_in)
    logical,  intent(in)  :: lasym
    integer,  intent(in)  :: nfp
    real(dp), intent(in)  :: separation
    integer,  intent(in)  :: mpol_out, ntor_out, ntheta_transform, nzeta_transform
    real(dp), intent(in)  :: tol
    integer,  intent(in)  :: mnmax_out
    integer,  intent(out) :: xm_out(mnmax_out), xn_out(mnmax_out)
    real(dp), intent(out) :: rmnc_out(mnmax_out), rmns_out(mnmax_out), zmnc_out(mnmax_out), zmns_out(mnmax_out)
    integer,  intent(out) :: info

    integer :: itheta, izeta, jm, jn, j, index, mnmax_check, fail_flag
    real(dp) :: x_out, y_out, z_out, factor, factor2, angle, sinangle, cosangle
    integer :: iflag_local
    real(dp), allocatable :: theta_grid(:), zeta_grid(:), major_r(:,:), z_val(:,:)

    info = 0

    if (mnmax_in < 1 .or. nfp < 1 .or. ntheta_transform < 1 .or. nzeta_transform < 1 &
         .or. mpol_out < 0 .or. ntor_out < 0) then
       info = 1
       return
    end if

    mnmax_check = mpol_out*(2*ntor_out+1) + ntor_out + 1
    if (mnmax_check /= mnmax_out) then
       info = 2
       return
    end if

    allocate(theta_grid(ntheta_transform), zeta_grid(nzeta_transform))
    allocate(major_r(ntheta_transform, nzeta_transform), z_val(ntheta_transform, nzeta_transform))

    do itheta = 1, ntheta_transform
       theta_grid(itheta) = twopi*(itheta-1.0_dp)/ntheta_transform
    end do
    do izeta = 1, nzeta_transform
       zeta_grid(izeta) = (twopi/nfp)*(izeta-1.0_dp)/nzeta_transform
    end do

    fail_flag = 0

    !$OMP PARALLEL DO PRIVATE(itheta,x_out,y_out,z_out,iflag_local) SCHEDULE(static)
    do izeta = 1, nzeta_transform
       do itheta = 1, ntheta_transform
          call regcoil_offset_surface_point(theta_grid(itheta), zeta_grid(izeta), separation, &
               mnmax_in, xm_in, xn_in, rmnc_in, rmns_in, zmnc_in, zmns_in, lasym, tol, &
               x_out, y_out, z_out, iflag_local)
          if (iflag_local == 4) then
             !$OMP ATOMIC WRITE
             fail_flag = 1
          end if
          major_r(itheta, izeta) = sqrt(x_out*x_out + y_out*y_out)
          z_val(itheta, izeta) = z_out
       end do
    end do
    !$OMP END PARALLEL DO

    if (fail_flag /= 0) then
       info = 3
       deallocate(theta_grid, zeta_grid, major_r, z_val)
       return
    end if

    ! Generate the output Fourier mode numbers (same convention as
    ! regcoil_init_Fourier_modes_mod with include_00=.true.): xm >= 0, xn can
    ! be any sign, and xn ends up scaled by nfp (VMEC convention).
    index = 1
    xm_out(1) = 0
    xn_out(1) = 0
    do jn = 1, ntor_out
       index = index + 1
       xm_out(index) = 0
       xn_out(index) = jn
    end do
    do jm = 1, mpol_out
       do jn = -ntor_out, ntor_out
          index = index + 1
          xm_out(index) = jm
          xn_out(index) = jn
       end do
    end do
    xn_out = xn_out * nfp

    rmnc_out = 0
    rmns_out = 0
    zmnc_out = 0
    zmns_out = 0

    factor = 2.0_dp / (ntheta_transform * nzeta_transform)
    do izeta = 1, nzeta_transform
       do itheta = 1, ntheta_transform
          do j = 2, mnmax_out
             angle = xm_out(j) * theta_grid(itheta) - xn_out(j) * zeta_grid(izeta)
             sinangle = sin(angle)
             cosangle = cos(angle)
             factor2 = factor
             ! Halve the weight of the Nyquist mode so that
             ! inverse-transform(transform(.)) is the identity.
             if (mod(ntheta_transform,2) == 0 .and. xm_out(j) == (ntheta_transform/2)) factor2 = factor2 / 2
             if (mod(nzeta_transform,2) == 0 .and. abs(xn_out(j)) == nfp*(nzeta_transform/2)) factor2 = factor2 / 2
             rmnc_out(j) = rmnc_out(j) + major_r(itheta,izeta) * cosangle * factor2
             rmns_out(j) = rmns_out(j) + major_r(itheta,izeta) * sinangle * factor2
             zmnc_out(j) = zmnc_out(j) + z_val(itheta,izeta) * cosangle * factor2
             zmns_out(j) = zmns_out(j) + z_val(itheta,izeta) * sinangle * factor2
          end do
       end do
    end do
    rmnc_out(1) = sum(major_r) / (ntheta_transform * nzeta_transform)
    zmnc_out(1) = sum(z_val) / (ntheta_transform * nzeta_transform)
    rmns_out(1) = 0
    zmns_out(1) = 0

    if (.not. lasym) then
       rmns_out = 0
       zmnc_out = 0
    end if

    deallocate(theta_grid, zeta_grid, major_r, z_val)

  end subroutine regcoil_uniform_offset_surface

  ! Root-solves for the point on the offset surface at poloidal angle
  ! theta_in whose *toroidal Cartesian angle* atan2(y,x) equals
  ! zeta_target_in (the legacy "cosm" construction), then returns its
  ! Cartesian position. iflag_out follows regcoil_fzero's convention (2 =
  ! exact hit, 1/3 = converged, 4 = no sign change / failure, 5 = too many
  ! evaluations).
  subroutine regcoil_offset_surface_point(theta_in, zeta_target_in, separation_in, &
       mnmax_in, xm_in, xn_in, rmnc_in, rmns_in, zmnc_in, zmns_in, lasym, tol_in, &
       x_out, y_out, z_out, iflag_out)

    real(dp), intent(in)  :: theta_in, zeta_target_in, separation_in, tol_in
    integer,  intent(in)  :: mnmax_in
    integer,  intent(in)  :: xm_in(mnmax_in), xn_in(mnmax_in)
    real(dp), intent(in)  :: rmnc_in(mnmax_in), rmns_in(mnmax_in), zmnc_in(mnmax_in), zmns_in(mnmax_in)
    logical,  intent(in)  :: lasym
    real(dp), intent(out) :: x_out, y_out, z_out
    integer,  intent(out) :: iflag_out

    real(dp) :: zeta_min, zeta_max, zeta_solution

    zeta_min = zeta_target_in - 1.0_dp
    zeta_max = zeta_target_in + 1.0_dp
    call regcoil_fzero(residual, zeta_min, zeta_max, zeta_target_in, tol_in, tol_in, iflag_out)
    ! regcoil_fzero returns its answer in the (intent inout) first bound.
    zeta_solution = zeta_min

    call evaluate_point(theta_in, zeta_solution, x_out, y_out, z_out)

  contains

    function residual(zeta_test)
      real(dp) :: zeta_test, residual
      real(dp) :: x_t, y_t, z_t, zeta_new, zeta_err

      call evaluate_point(theta_in, zeta_test, x_t, y_t, z_t)
      zeta_new = atan2(y_t, x_t)
      zeta_err = zeta_new - zeta_target_in
      if (zeta_err < -pi) zeta_err = zeta_err + twopi
      if (zeta_err > pi) zeta_err = zeta_err - twopi
      residual = zeta_err
    end function residual

    ! Evaluates the plasma surface at (theta_local, zeta_local) and offsets
    ! by separation_in along the (unit) local normal -- a stateless
    ! equivalent of regcoil_expand_plasma_surface.
    subroutine evaluate_point(theta_local, zeta_local, x_l, y_l, z_l)
      real(dp), intent(in)  :: theta_local, zeta_local
      real(dp), intent(out) :: x_l, y_l, z_l

      integer :: imn
      real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
      real(dp) :: angle, sinangle, cosangle
      real(dp) :: dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
      real(dp) :: phi, sinphi, cosphi, dsinphidzeta, dcosphidzeta
      real(dp) :: normal_x, normal_y, normal_z, normal_abs

      x_l = 0; y_l = 0; z_l = 0
      dxdtheta = 0; dydtheta = 0; dzdtheta = 0
      dxdzeta = 0; dydzeta = 0; dzdzeta = 0

      phi = zeta_local
      sinphi = sin(phi)
      cosphi = cos(phi)
      dsinphidzeta = cosphi
      dcosphidzeta = -sinphi

      do imn = 1, mnmax_in
         angle = xm_in(imn)*theta_local - xn_in(imn)*zeta_local
         cosangle = cos(angle)
         sinangle = sin(angle)
         dsinangledtheta = cosangle*xm_in(imn)
         dcosangledtheta = -sinangle*xm_in(imn)
         dsinangledzeta = -cosangle*xn_in(imn)
         dcosangledzeta = sinangle*xn_in(imn)

         x_l = x_l + rmnc_in(imn) * cosangle * cosphi
         y_l = y_l + rmnc_in(imn) * cosangle * sinphi
         z_l = z_l + zmns_in(imn) * sinangle

         dxdtheta = dxdtheta + rmnc_in(imn) * dcosangledtheta * cosphi
         dydtheta = dydtheta + rmnc_in(imn) * dcosangledtheta * sinphi
         dzdtheta = dzdtheta + zmns_in(imn) * dsinangledtheta

         dxdzeta = dxdzeta + rmnc_in(imn) * (dcosangledzeta * cosphi + cosangle * dcosphidzeta)
         dydzeta = dydzeta + rmnc_in(imn) * (dcosangledzeta * sinphi + cosangle * dsinphidzeta)
         dzdzeta = dzdzeta + zmns_in(imn) * dsinangledzeta

         if (lasym) then
            x_l = x_l + rmns_in(imn) * sinangle * cosphi
            y_l = y_l + rmns_in(imn) * sinangle * sinphi
            z_l = z_l + zmnc_in(imn) * cosangle

            dxdtheta = dxdtheta + rmns_in(imn) * dsinangledtheta * cosphi
            dydtheta = dydtheta + rmns_in(imn) * dsinangledtheta * sinphi
            dzdtheta = dzdtheta + zmnc_in(imn) * dcosangledtheta

            dxdzeta = dxdzeta + rmns_in(imn) * (dsinangledzeta * cosphi + sinangle * dcosphidzeta)
            dydzeta = dydzeta + rmns_in(imn) * (dsinangledzeta * sinphi + sinangle * dsinphidzeta)
            dzdzeta = dzdzeta + zmnc_in(imn) * dcosangledzeta
         end if
      end do

      normal_x = dydzeta * dzdtheta - dydtheta * dzdzeta
      normal_y = dzdzeta * dxdtheta - dzdtheta * dxdzeta
      normal_z = dxdzeta * dydtheta - dxdtheta * dydzeta
      normal_abs = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
      normal_x = normal_x / normal_abs
      normal_y = normal_y / normal_abs
      normal_z = normal_z / normal_abs

      x_l = x_l + normal_x * separation_in
      y_l = y_l + normal_y * separation_in
      z_l = z_l + normal_z * separation_in

    end subroutine evaluate_point

  end subroutine regcoil_offset_surface_point

end module regcoil_uniform_offset_surface_mod
