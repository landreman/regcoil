subroutine regcoil_evaluate_coil_surface()

  ! This subroutine takes the arrays rmnc_coil, zmns_coil, etc, and evaluates the position vector r_coil
  ! and its first and second derivatives with respect to theta and zeta.

  use regcoil_variables
  use stel_kinds
  
  implicit none
  
  integer :: imn, m, n, iflag, j
  real(dp) :: rmnc, rmns, zmnc, zmns
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta, dsinangledzeta, dcosangledzeta
  real(dp) :: d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2, d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta, d2sinangle2dzeta2, d2cosangle2dzeta2
  real(dp), dimension(:,:), allocatable :: major_R_squared
  integer :: itheta, izeta
  integer :: tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: sin_m_theta, cos_m_theta, sin_n_zeta, cos_n_zeta

  call system_clock(tic,countrate)

  allocate(sin_m_theta(mnmax_coil,ntheta_coil))
  allocate(cos_m_theta(mnmax_coil,ntheta_coil))
  allocate(sin_n_zeta( mnmax_coil,nzetal_coil))
  allocate(cos_n_zeta( mnmax_coil,nzetal_coil))
  do itheta = 1,ntheta_coil
     sin_m_theta(:,itheta) = sin(xm_coil * theta_coil(itheta))
     cos_m_theta(:,itheta) = cos(xm_coil * theta_coil(itheta))
  end do
  do izeta = 1,nzetal_coil
     sin_n_zeta(:,  izeta) = sin(xn_coil * zetal_coil(izeta))
     cos_n_zeta(:,  izeta) = cos(xn_coil * zetal_coil(izeta))
  end do

  r_coil = 0
  drdtheta_coil = 0
  drdzeta_coil = 0
  d2rdtheta2_coil = 0
  d2rdthetadzeta_coil = 0
  d2rdzeta2_coil = 0
  
  do izeta = 1,nzetal_coil
     angle2 = zetal_coil(izeta)
     sinangle2 = sin(angle2)
     cosangle2 = cos(angle2)
     dsinangle2dzeta = cosangle2
     dcosangle2dzeta = -sinangle2
     d2sinangle2dzeta2 = -sinangle2
     d2cosangle2dzeta2 = -cosangle2
     do imn = 1, mnmax_coil
        m = xm_coil(imn)
        n = xn_coil(imn)
        rmnc = rmnc_coil(imn)
        rmns = rmns_coil(imn)
        zmnc = zmnc_coil(imn)
        zmns = zmns_coil(imn)
        do itheta = 1,ntheta_coil
           !angle = m*theta_coil(itheta) - n*zetal_coil(izeta)
           !sinangle = sin(angle)
           !cosangle = cos(angle)
           ! Trig angle sum formulae for angle = m*theta_coil(itheta) - n*zetal_coil(izeta):
           sinangle = sin_m_theta(imn,itheta) * cos_n_zeta(imn,izeta) - cos_m_theta(imn,itheta) * sin_n_zeta(imn,izeta)
           cosangle = cos_m_theta(imn,itheta) * cos_n_zeta(imn,izeta) + sin_m_theta(imn,itheta) * sin_n_zeta(imn,izeta)
           !if (abs(sinangle - sin(angle)) > 1d-10) stop "Error sin"
           !if (abs(cosangle - cos(angle)) > 1d-10) stop "Error cos"
           dsinangledtheta = cosangle*m
           dcosangledtheta = -sinangle*m
           dsinangledzeta = -cosangle*n
           dcosangledzeta = sinangle*n
           
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

           ! 2nd derivatives are only constructed on 1 field period.
           if (izeta > nzeta_coil) cycle

           d2sinangledtheta2  = -m*m*sinangle
           d2sinangledthetadzeta = m*n*sinangle
           d2sinangledzeta2  = -n*n*sinangle
           d2cosangledtheta2  = -m*m*cosangle
           d2cosangledthetadzeta = m*n*cosangle
           d2cosangledzeta2  = -n*n*cosangle

           d2rdtheta2_coil(1,itheta,izeta) = d2rdtheta2_coil(1,itheta,izeta) + rmnc * d2cosangledtheta2 * cosangle2 + rmns * d2sinangledtheta2 * cosangle2
           d2rdtheta2_coil(2,itheta,izeta) = d2rdtheta2_coil(2,itheta,izeta) + rmnc * d2cosangledtheta2 * sinangle2 + rmns * d2sinangledtheta2 * sinangle2
           d2rdtheta2_coil(3,itheta,izeta) = d2rdtheta2_coil(3,itheta,izeta) + zmns * d2sinangledtheta2 + zmnc * d2cosangledtheta2
           
           d2rdthetadzeta_coil(1,itheta,izeta) = d2rdthetadzeta_coil(1,itheta,izeta) + rmnc * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta) &
                + rmns * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
           d2rdthetadzeta_coil(2,itheta,izeta) = d2rdthetadzeta_coil(2,itheta,izeta) + rmnc * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta) &
                + rmns * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
           d2rdthetadzeta_coil(3,itheta,izeta) = d2rdthetadzeta_coil(3,itheta,izeta) + zmns * d2sinangledthetadzeta + zmnc * d2cosangledthetadzeta
           
           d2rdzeta2_coil(1,itheta,izeta) = d2rdzeta2_coil(1,itheta,izeta) + rmnc * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
                + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2) &
                + rmns * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
                + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
           d2rdzeta2_coil(2,itheta,izeta) = d2rdzeta2_coil(2,itheta,izeta) + rmnc * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
                + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2) &
                + rmns * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
                + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
           d2rdzeta2_coil(3,itheta,izeta) = d2rdzeta2_coil(3,itheta,izeta) + zmns * d2sinangledzeta2 + zmnc * d2cosangledzeta2
        end do
     end do
  end do

  deallocate(sin_m_theta, cos_m_theta, sin_n_zeta, cos_n_zeta)

  call system_clock(toc)
  if (verbose) print *,"  Evaluating coil surface & derivatives:",real(toc-tic)/countrate," sec."

!!$  print *,"mnmax_coil:",mnmax_coil
!!$  print *,"xm_coil:",xm_coil
!!$  print *,"xn_coil:",xn_coil
!!$  print *,"rmnc_coil:",rmnc_coil
!!$  print *,"zmns_coil:",zmns_coil
!!$  print *,"r_coil(1,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(1,j,:)
!!$  end do
!!$  print *,"r_coil(2,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(2,j,:)
!!$  end do
!!$  print *,"r_coil(3,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(3,j,:)
!!$  end do
  
  ! Evaluate cross product:
  normal_coil(1,:,:) = drdzeta_coil(2,:,:) * drdtheta_coil(3,:,:) - drdtheta_coil(2,:,:) * drdzeta_coil(3,:,:)
  normal_coil(2,:,:) = drdzeta_coil(3,:,:) * drdtheta_coil(1,:,:) - drdtheta_coil(3,:,:) * drdzeta_coil(1,:,:)
  normal_coil(3,:,:) = drdzeta_coil(1,:,:) * drdtheta_coil(2,:,:) - drdtheta_coil(1,:,:) * drdzeta_coil(2,:,:)
  
  if (allocated(norm_normal_coil)) deallocate(norm_normal_coil)
  allocate(norm_normal_coil(ntheta_coil, nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 11'
  norm_normal_coil = sqrt(normal_coil(1,:,1:nzeta_coil)**2 + normal_coil(2,:,1:nzeta_coil)**2 &
       +  normal_coil(3,:,1:nzeta_coil)**2)
  
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
  if (verbose) print "(a,es10.3,a,es10.3,a)"," Coil surface area:",area_coil," m^2. Volume:",volume_coil," m^3."
  
end subroutine  regcoil_evaluate_coil_surface
