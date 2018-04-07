module regcoil_nescin_utils_2

contains

subroutine initupdate_nescin_winding_surface(lscreen_optin)

  use stel_constants
  use stel_kinds
  use regcoil_variables, only: &
        ntheta_coil, nzeta_coil, nzetal_coil, theta_coil, zeta_coil, zetal_coil, &
        r_coil, drdtheta_coil, drdzeta_coil, &
        normal_coil, norm_normal_coil, area_coil, volume_coil, &
        R0_coil, a_coil, &
        area_coil, volume_coil, dtheta_coil, dzeta_coil, nfp, &
        rc_rmnc_stellopt, rc_rmns_stellopt, &
        rc_zmnc_stellopt, rc_zmns_stellopt, mpol_rcws, ntor_rcws
  implicit none

  integer :: iflag, which_surface, iunit  = 7
  integer  :: i, itheta, izeta, nummodes, xm_in, xn_in, istat, ii, jj, kk
  integer :: m, n, ntotal, k,  mr, nr, nonzeromodes

  real(dp), dimension(:,:), allocatable :: major_R_squared
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta, &
                dcosangledzeta, dsinangledzeta
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
  real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
  !real(dp) :: d2sinangle2dzeta2, d2cosangle2dzeta2
  !real(dp) :: d2sinangledtheta2,
  real(dp)  :: d2sinangledthetadzeta, d2sinangledzeta2
  real(dp)  :: d2cosangledthetadzeta, d2cosangledzeta2
!  real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_theta, delta_zeta, temp
  real(dp) :: rmnc, rmns, zmnc, zmns

  ! variables to handle printing to the screen
  logical, optional :: lscreen_optin
  logical :: lscreen
  
  logical :: compute_2nd_derivs = .false.

  if(present(lscreen_optin)) then 
    lscreen = lscreen_optin
  else
    lscreen = .true.
  endif

  ! zero out the 2-D arrays of r_coil, drdtheta_coil and drdzeta_coil
  r_coil=0
  drdtheta_coil=0
  drdzeta_coil=0

  if (compute_2nd_derivs) then
     !d2rdtheta2 = 0
     !d2rdthetadzeta = 0
     !d2rdzeta2 = 0
  end if

  nonzeromodes = 0
  do ii = -mpol_rcws,mpol_rcws
    do jj = -ntor_rcws,ntor_rcws
     m = ii
     n = jj
     rmnc = rc_rmnc_stellopt(ii,jj)
     rmns = rc_rmns_stellopt(ii,jj)
     zmnc = rc_zmnc_stellopt(ii,jj)
     zmns = rc_zmns_stellopt(ii,jj)
     if ( (rmnc .ne. 0) .or. (rmns .ne. 0) .or. (zmnc .ne. 0) .or. (zmns .ne. 0) ) then
       nonzeromodes = nonzeromodes + 1
       do izeta = 1,nzetal_coil
          angle2 = zetal_coil(izeta)
          sinangle2 = sin(angle2)
          cosangle2 = cos(angle2)
          dsinangle2dzeta = cosangle2
          dcosangle2dzeta = -sinangle2
          d2sinangle2dzeta2 = -sinangle2
          d2cosangle2dzeta2 = -cosangle2
          do itheta = 1,ntheta_coil
             angle = m*theta_coil(itheta) + n*nfp*zetal_coil(izeta)
             sinangle = sin(angle)
             cosangle = cos(angle)
             dsinangledtheta = cosangle*m
             dcosangledtheta = -sinangle*m
             dsinangledzeta = cosangle*n*nfp
             dcosangledzeta = -sinangle*n*nfp
             d2sinangledtheta2  = -m*m*sinangle
             d2sinangledthetadzeta = -m*n*nfp*sinangle
             d2sinangledzeta2  = -n*n*nfp*nfp*sinangle
             d2cosangledtheta2  = -m*m*cosangle
             d2cosangledthetadzeta = -m*n*nfp*cosangle
             d2cosangledzeta2  = -n*n*nfp*nfp*cosangle
             
             r_coil(1,itheta,izeta) = r_coil(1,itheta,izeta) + rmnc * cosangle * cosangle2 + rmns * sinangle * cosangle2
             r_coil(2,itheta,izeta) = r_coil(2,itheta,izeta) + rmnc * cosangle * sinangle2 + rmns * sinangle * sinangle2
             r_coil(3,itheta,izeta) = r_coil(3,itheta,izeta) + zmns * sinangle             + zmnc * cosangle
             
             drdtheta_coil(1,itheta,izeta) = drdtheta_coil(1,itheta,izeta) + rmnc * dcosangledtheta * cosangle2 + rmns * dsinangledtheta * cosangle2
             drdtheta_coil(2,itheta,izeta) = drdtheta_coil(2,itheta,izeta) + rmnc * dcosangledtheta * sinangle2 + rmns * dsinangledtheta * sinangle2
             drdtheta_coil(3,itheta,izeta) = drdtheta_coil(3,itheta,izeta) + zmns * dsinangledtheta + zmnc * dcosangledtheta
             
             drdzeta_coil(1,itheta,izeta) = drdzeta_coil(1,itheta,izeta) + rmnc * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta) &
                  + rmns * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
             drdzeta_coil(2,itheta,izeta) = drdzeta_coil(2,itheta,izeta) + rmnc * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta) &
                  + rmns * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
             drdzeta_coil(3,itheta,izeta) = drdzeta_coil(3,itheta,izeta) + zmns * dsinangledzeta + zmnc * dcosangledzeta
  
!             if (compute_2nd_derivs) then
!                d2rdtheta2(1,itheta,izeta) = d2rdtheta2(1,itheta,izeta) + rmnc * d2cosangledtheta2 * cosangle2 + rmns * d2sinangledtheta2 * cosangle2
!                d2rdtheta2(2,itheta,izeta) = d2rdtheta2(2,itheta,izeta) + rmnc * d2cosangledtheta2 * sinangle2 + rmns * d2sinangledtheta2 * sinangle2
!                d2rdtheta2(3,itheta,izeta) = d2rdtheta2(3,itheta,izeta) + zmns * d2sinangledtheta2 + zmnc * d2cosangledtheta2
!  
!                d2rdthetadzeta(1,itheta,izeta) = d2rdthetadzeta(1,itheta,izeta) + rmnc * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta) &
!                     + rmns * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
!                d2rdthetadzeta(2,itheta,izeta) = d2rdthetadzeta(2,itheta,izeta) + rmnc * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta) &
!                     + rmns * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
!                d2rdthetadzeta(3,itheta,izeta) = d2rdthetadzeta(3,itheta,izeta) + zmns * d2sinangledthetadzeta + zmnc * d2cosangledthetadzeta
!             
!                d2rdzeta2(1,itheta,izeta) = d2rdzeta2(1,itheta,izeta) + rmnc * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
!                     + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2) &
!                     + rmns * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
!                     + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
!                d2rdzeta2(2,itheta,izeta) = d2rdzeta2(2,itheta,izeta) + rmnc * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
!                     + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2) &
!                     + rmns * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
!                     + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
!                d2rdzeta2(3,itheta,izeta) = d2rdzeta2(3,itheta,izeta) + zmns * d2sinangledzeta2 + zmnc * d2cosangledzeta2
!             end if
          end do
       end do
     end if
   end do
 end do

  !print *,"<---REGCOIL found", nonzeromodes," non-zero Fourier modes"
  

      ! Evaluate cross product:
      normal_coil(1,:,:) = drdzeta_coil(2,:,:) * drdtheta_coil(3,:,:) - drdtheta_coil(2,:,:) * drdzeta_coil(3,:,:)
      normal_coil(2,:,:) = drdzeta_coil(3,:,:) * drdtheta_coil(1,:,:) - drdtheta_coil(3,:,:) * drdzeta_coil(1,:,:)
      normal_coil(3,:,:) = drdzeta_coil(1,:,:) * drdtheta_coil(2,:,:) - drdtheta_coil(1,:,:) * drdzeta_coil(2,:,:)


      norm_normal_coil = sqrt(normal_coil(1,:,1:nzeta_coil)**2 + normal_coil(2,:,1:nzeta_coil)**2 + normal_coil(3,:,1:nzeta_coil)**2)
      !print *,"<----norm_normal:", norm_normal_coil

      area_coil = nfp * dtheta_coil * dzeta_coil * sum(norm_normal_coil)

      ! Compute coil surface volume using \int (1/2) R^2 dZ dzeta.
      ! These quantities will be evaluated on the half theta grid, which is the natural grid for dZ,
      ! but we will need to interpolate R^2 from the full to half grid.
      allocate(major_R_squared(ntheta_coil,nzetal_coil))
      major_R_squared = r_coil(1,:,:)*r_coil(1,:,:) + r_coil(2,:,:)*r_coil(2,:,:)
      ! First handle the interior of the theta grid:
      volume_coil = sum((major_R_squared(1:ntheta_coil-1,:) + major_R_squared(2:ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r_coil(3,2:ntheta_coil,:)-r_coil(3,1:ntheta_coil-1,:))) ! dZ
      ! Add the contribution from the ends of the theta grid:
      volume_coil = volume_coil + sum((major_R_squared(1,:) + major_R_squared(ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r_coil(3,1,:)-r_coil(3,ntheta_coil,:))) ! dZ
      volume_coil = abs(volume_coil * dzeta_coil / 2) ! r includes all nfp periods already, so no factor of nfp needed.
      deallocate(major_R_squared)
      if (lscreen) print "(a,es10.3,a,es10.3,a)"," Coil surface area:",area_coil," m^2. Volume:",volume_coil," m^3."

      !print *,"<----Exiting regcoil initupdate_nescin_winding_surface."
 
end subroutine initupdate_nescin_winding_surface


end module regcoil_nescin_utils_2
