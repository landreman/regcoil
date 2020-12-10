subroutine regcoil_double_to_single()

    use regcoil_variables
    use regcoil_splines
    use read_wout_mod, only: raxis, zaxis

    use stel_kinds
    use stel_constants
    use omp_lib

    implicit none

    integer :: iflag
    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta, d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2, d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2
    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta, d2sinangle2dzeta2, d2cosangle2dzeta2
    real(dp) :: cos_norm, sin_norm, dtheta, dzeta
    integer :: i, j, imn, itheta, izeta
    integer :: tic, toc, countrate, tic1, toc1
    type (periodic_spline) :: theta_spline
    real(dp), allocatable, dimension(:) :: theta, zeta
    real(dp), allocatable, dimension(:,:) :: R, Z, B, R_0, Z_0, x, y
    real(dp), allocatable, dimension(:,:) :: l, vartheta_on_theta_grid, l_vartheta, arclength, l_arclength, omega, omega_on_arclength_grid, norm_normal
    real(dp), allocatable, dimension(:,:,:) :: drdtheta, drdzeta, normal, d2rdtheta2, d2rdthetadzeta, d2rdzeta2
    real(dp), allocatable, dimension(:) :: x1, y1, x2, y2
    real(dp), dimension(:,:), allocatable :: mean_curvature_temp, E_big, F_big, G_big, e_small, f_small, g_small, K_curvature, kappa1, kappa2
    real(dp), dimension(:,:), allocatable :: x_offset, y_offset, z_offset, R_offset
    real(dp) :: max_separation_arclength
    real(dp) :: dR, dZ, tot_length
    integer :: ntheta, nzeta
    real(dp) :: norm2l, norm2o, sum_amp2l, sum_amp2o
    real(dp) :: temp
    logical :: sorted

    if (geometry_option_plasma .ne. 8) then
        print *,"Invalid setting for geometry_option_plasma: ",geometry_option_plasma
        stop
    end if

    call system_clock(tic,countrate)
    if (verbose) print *,"Converting from double to single Fourier representation of plasma surface."

    ntheta = ntheta_plasma * theta_transform_refinement
    nzeta = nzeta_plasma * zeta_transform_refinement

    allocate(theta(ntheta))
    allocate(zeta(nzeta))

    do i = 1,ntheta
     theta(i) = twopi*(i-1.0_dp)/ntheta
    end do

    do i = 1,nzeta
     zeta(i) = twopi/nfp*(i-1.0_dp)/nzeta
    end do

    allocate(R(ntheta,nzeta))
    allocate(Z(ntheta,nzeta))
    allocate(B(ntheta,nzeta))
    allocate(norm_normal(ntheta,nzeta))
    allocate(x(ntheta,nzeta))
    allocate(y(ntheta,nzeta))

    allocate(drdtheta(3,ntheta,nzeta))
    allocate(drdzeta(3,ntheta,nzeta))
    allocate(normal(3,ntheta,nzeta))

    allocate(d2rdtheta2(3,ntheta,nzeta))
    allocate(d2rdthetadzeta(3,ntheta,nzeta))
    allocate(d2rdzeta2(3,ntheta,nzeta))

    allocate(R_0(ntheta,nzeta))
    allocate(Z_0(ntheta,nzeta))

    R=0
    Z=0
    B=0
    drdtheta=0
    drdzeta=0
    d2rdtheta2=0
    d2rdthetadzeta=0
    d2rdzeta2=0
    R_0=0
    Z_0=0
    do izeta = 1, nzeta
        angle2 = zeta(izeta)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dzeta = cosangle2
        dcosangle2dzeta = -sinangle2
        do itheta = 1, ntheta
            do imn = 1, mnmax_plasma
                angle = xm_plasma(imn) * theta(itheta) - xn_plasma(imn) * zeta(izeta)
                sinangle = sin(angle)
                cosangle = cos(angle)
                dsinangledtheta = cosangle*xm_plasma(imn)
                dcosangledtheta = -sinangle*xm_plasma(imn)
                dsinangledzeta = -cosangle*xn_plasma(imn)
                dcosangledzeta = sinangle*xn_plasma(imn)
                d2sinangledtheta2 = -sinangle*xm_plasma(imn)**2
                d2cosangledtheta2 = -cosangle*xm_plasma(imn)**2
                d2sinangledthetadzeta = sinangle*xm_plasma(imn)*xn_plasma(imn)
                d2cosangledthetadzeta = cosangle*xm_plasma(imn)*xn_plasma(imn)
                d2sinangledzeta2 = -sinangle*xn_plasma(imn)**2
                d2cosangledzeta2 = -cosangle*xn_plasma(imn)**2

                R(itheta,izeta) = R(itheta,izeta) + rmnc_plasma(imn)*cosangle
                Z(itheta,izeta) = Z(itheta,izeta) + zmns_plasma(imn)*sinangle
                B(itheta,izeta) = B(itheta,izeta) + bmnc_plasma(imn)*cosangle

                drdtheta(1,itheta,izeta) = drdtheta(1,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * cosangle2
                drdtheta(2,itheta,izeta) = drdtheta(2,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * sinangle2
                drdtheta(3,itheta,izeta) = drdtheta(3,itheta,izeta) + zmns_plasma(imn) * dsinangledtheta

                drdzeta(1,itheta,izeta) = drdzeta(1,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta)
                drdzeta(2,itheta,izeta) = drdzeta(2,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta)
                drdzeta(3,itheta,izeta) = drdzeta(3,itheta,izeta) + zmns_plasma(imn) * dsinangledzeta

                d2rdtheta2(1,itheta,izeta) = d2rdtheta2(1,itheta,izeta) + rmnc_plasma(imn) * d2cosangledtheta2 * cosangle2 + rmns_plasma(imn) * d2sinangledtheta2 * cosangle2
                d2rdtheta2(2,itheta,izeta) = d2rdtheta2(2,itheta,izeta) + rmnc_plasma(imn) * d2cosangledtheta2 * sinangle2 + rmns_plasma(imn) * d2sinangledtheta2 * sinangle2
                d2rdtheta2(3,itheta,izeta) = d2rdtheta2(3,itheta,izeta) + zmns_plasma(imn) * d2sinangledtheta2 + zmnc_plasma(imn) * d2cosangledtheta2
                
                d2rdthetadzeta(1,itheta,izeta) = d2rdthetadzeta(1,itheta,izeta) + rmnc_plasma(imn) * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta)
                d2rdthetadzeta(2,itheta,izeta) = d2rdthetadzeta(2,itheta,izeta) + rmnc_plasma(imn) * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta)
                d2rdthetadzeta(3,itheta,izeta) = d2rdthetadzeta(3,itheta,izeta) + zmns_plasma(imn) * d2sinangledthetadzeta
                
                d2rdzeta2(1,itheta,izeta) = d2rdzeta2(1,itheta,izeta) + rmnc_plasma(imn) * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
                     + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2)
                d2rdzeta2(2,itheta,izeta) = d2rdzeta2(2,itheta,izeta) + rmnc_plasma(imn) * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
                     + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2)
                d2rdzeta2(3,itheta,izeta) = d2rdzeta2(3,itheta,izeta) + zmns_plasma(imn) * d2sinangledzeta2
                if (lasym) then
                    R(itheta,izeta) = R(itheta,izeta) + rmns_plasma(imn)*sinangle
                    Z(itheta,izeta) = Z(itheta,izeta) + zmnc_plasma(imn)*cosangle
                    B(itheta,izeta) = B(itheta,izeta) + bmns_plasma(imn)*sinangle

                    drdtheta(1,itheta,izeta) = drdtheta(1,itheta,izeta) + rmns_plasma(imn) * dsinangledtheta * cosangle2
                    drdtheta(2,itheta,izeta) = drdtheta(2,itheta,izeta) + rmns_plasma(imn) * dsinangledtheta * sinangle2
                    drdtheta(3,itheta,izeta) = drdtheta(3,itheta,izeta) + zmnc_plasma(imn) * dcosangledtheta

                    drdzeta(1,itheta,izeta) = drdzeta(1,itheta,izeta) + rmns_plasma(imn) * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
                    drdzeta(2,itheta,izeta) = drdzeta(2,itheta,izeta) + rmns_plasma(imn) * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
                    drdzeta(3,itheta,izeta) = drdzeta(3,itheta,izeta) + zmnc_plasma(imn) * dcosangledzeta

                    d2rdtheta2(1,itheta,izeta) = d2rdtheta2(1,itheta,izeta) + rmns_plasma(imn) * d2sinangledtheta2 * cosangle2
                    d2rdtheta2(2,itheta,izeta) = d2rdtheta2(2,itheta,izeta) + rmns_plasma(imn) * d2sinangledtheta2 * sinangle2
                    d2rdtheta2(3,itheta,izeta) = d2rdtheta2(3,itheta,izeta)+ zmnc_plasma(imn) * d2cosangledtheta2
                    
                    d2rdthetadzeta(1,itheta,izeta) = d2rdthetadzeta(1,itheta,izeta) + rmns_plasma(imn) * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
                    d2rdthetadzeta(2,itheta,izeta) = d2rdthetadzeta(2,itheta,izeta) + rmns_plasma(imn) * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
                    d2rdthetadzeta(3,itheta,izeta) = d2rdthetadzeta(3,itheta,izeta) + zmnc_plasma(imn) * d2cosangledthetadzeta
                    
                    d2rdzeta2(1,itheta,izeta) = d2rdzeta2(1,itheta,izeta) + rmns_plasma(imn) * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
                         + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
                    d2rdzeta2(2,itheta,izeta) = d2rdzeta2(2,itheta,izeta) + rmns_plasma(imn) * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
                         + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
                    d2rdzeta2(3,itheta,izeta) = d2rdzeta2(3,itheta,izeta) + zmnc_plasma(imn) * d2cosangledzeta2
                end if
            end do
        end do
        x(:,izeta) = R(:,izeta) * cosangle2
        y(:,izeta) = R(:,izeta) * sinangle2
    end do

    ! Evaluate cross product
    normal(1,:,:) = drdzeta(2,:,:) * drdtheta(3,:,:) - drdtheta(2,:,:) * drdzeta(3,:,:)
    normal(2,:,:) = drdzeta(3,:,:) * drdtheta(1,:,:) - drdtheta(3,:,:) * drdzeta(1,:,:)
    normal(3,:,:) = drdzeta(1,:,:) * drdtheta(2,:,:) - drdtheta(1,:,:) * drdzeta(2,:,:)

    norm_normal = sqrt(normal(1,:,:)**2 &
        + normal(2,:,:)**2 &
        + normal(3,:,:)**2)

    dtheta = theta(2)-theta(1)
    dzeta = zeta(2)-zeta(1)

    B_0 = sqrt(sum(B**2 * norm_normal) / sum(norm_normal))
!print *,B_0
!stop

    if (allocated(mean_curvature_temp)) deallocate(mean_curvature_temp)
    allocate(mean_curvature_temp(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(E_big(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(F_big(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(G_big(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(e_small(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(f_small(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(g_small(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(K_curvature(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(kappa1(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(kappa2(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'

    E_big = drdtheta(1,:,:)**2 + drdtheta(2,:,:)**2 + drdtheta(3,:,:)**2
    F_big = drdtheta(1,:,:)*drdzeta(1,:,:) + drdtheta(2,:,:)*drdzeta(2,:,:) + drdtheta(3,:,:)*drdzeta(3,:,:)
    G_big = drdzeta(1,:,:)**2 + drdzeta(2,:,:)**2 + drdzeta(3,:,:)**2
    e_small = -(d2rdtheta2(1,:,:)*normal(1,:,:) + d2rdtheta2(2,:,:)*normal(2,:,:) + d2rdtheta2(3,:,:)*normal(3,:,:)) / norm_normal
    f_small = -(d2rdthetadzeta(1,:,:)*normal(1,:,:) + d2rdthetadzeta(2,:,:)*normal(2,:,:) + d2rdthetadzeta(3,:,:)*normal(3,:,:)) / norm_normal
    g_small = -(d2rdzeta2(1,:,:)*normal(1,:,:) + d2rdzeta2(2,:,:)*normal(2,:,:) + d2rdzeta2(3,:,:)*normal(3,:,:)) / norm_normal

    mean_curvature_temp = (e_small*G_big - 2*f_small*F_big + g_small*E_big)/(2*(E_big*G_big-F_big*F_big))
    K_curvature = (e_small*g_small - f_small*f_small)/(E_big*G_big-F_big*F_big)
    kappa1 = mean_curvature_temp + sqrt(mean_curvature_temp**2 - K_curvature)
    kappa2 = mean_curvature_temp - sqrt(mean_curvature_temp**2 - K_curvature)

    max_separation_arclength = 1/abs(minval(kappa2))

    if (arclength_separation >= max_separation_arclength) then
      arclength_separation = max_separation_arclength
      print "(a,es10.3,a)","Separation is too large. Setting to maximum separation = ",max_separation_arclength," m."
    end if

    ! Can deallocate derivatives of plasma surface and curvature information now

    allocate(x_offset(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(y_offset(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(z_offset(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'
    allocate(R_offset(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_double_to_single Allocation error!'

    x_offset = x + arclength_separation * normal(1,:,:) / norm_normal
    y_offset = y + arclength_separation * normal(2,:,:) / norm_normal
    z_offset = Z + arclength_separation * normal(3,:,:) / norm_normal
    R_offset = sqrt(x_offset**2 + y_offset**2);

    raxis_cc = raxis(:,1)
    zaxis_cs = zaxis(:,1)
    nmax_axis = size(raxis_cc)
    allocate(xn_axis(nmax_axis))
    do izeta = 1, nzeta
        do itheta = 1, ntheta
            do j=1,nmax_axis
                xn_axis(j) = nfp * (j-1)
                angle = -xn_axis(j)*zeta(izeta)

                R_0(itheta,izeta) = R_0(itheta,izeta) + raxis_cc(j)*cos(angle)
                Z_0(itheta,izeta) = Z_0(itheta,izeta) + zaxis_cs(j)*sin(angle)
            end do
        end do
    end do

    allocate(l(ntheta,nzeta))
    allocate(vartheta_on_theta_grid(ntheta,nzeta))
    l = sqrt((R-R_0)**2 + (Z-Z_0)**2)
    vartheta_on_theta_grid = atan2(Z-Z_0, R-R_0)
    do izeta = 1, nzeta
        do itheta = 1, ntheta
            if (vartheta_on_theta_grid(itheta,izeta)<0) then
                vartheta_on_theta_grid(itheta,izeta) = vartheta_on_theta_grid(itheta,izeta) + twopi
            end if
        end do
    end do

    if (use_arclength_angle) then
      allocate(arclength(ntheta,nzeta))
      do izeta = 1,nzeta
        do itheta = 2,ntheta
          dR = R_offset(itheta,izeta)-R_offset(itheta-1,izeta)
          dZ = z_offset(itheta,izeta)-z_offset(itheta-1,izeta)
          arclength(itheta,izeta) = arclength(itheta-1,izeta) + sqrt(dR**2 + dZ**2)
        end do
        ! Total length
        dR = R_offset(1,izeta)-R_offset(ntheta,izeta)
        dZ = z_offset(1,izeta)-z_offset(ntheta,izeta)
        tot_length = arclength(ntheta,izeta) + sqrt(dR**2 + dZ**2)
        ! Normalize angle
        arclength(:,izeta) = arclength(:,izeta) * twopi / tot_length
      end do

      ! Compute difference between arclength and parameterization angle
      allocate(omega(ntheta,nzeta))
      omega = arclength - vartheta_on_theta_grid
      do izeta = 1, nzeta
          do itheta = 1, ntheta
              if (omega(itheta,izeta)>=twopi/2) then
                  omega(itheta,izeta) = omega(itheta,izeta) - twopi
              else if (omega(itheta,izeta)<-twopi/2) then
                  omega(itheta,izeta) = omega(itheta,izeta) + twopi
              end if
          end do
      end do

      allocate(l_arclength(ntheta,nzeta))
      allocate(omega_on_arclength_grid(ntheta,nzeta))
      ! Perform splines
      do izeta = 1, nzeta
          call new_periodic_spline(ntheta, arclength(:,izeta), &
              l(:,izeta), twopi, theta_spline)
          do itheta = 1, ntheta
              l_arclength(itheta,izeta) = periodic_splint(theta(itheta), theta_spline)
          end do
          call delete_periodic_spline(theta_spline)
      end do
      
      do izeta = 1, nzeta
          call new_periodic_spline(ntheta, arclength(:,izeta), &
              omega(:,izeta), twopi, theta_spline)
          do itheta = 1, ntheta
              omega_on_arclength_grid(itheta,izeta) = periodic_splint(theta(itheta), theta_spline)
          end do
          call delete_periodic_spline(theta_spline)
      end do

      do izeta = 1, nzeta
          do itheta = 1, ntheta
              if (l_arclength(itheta,izeta) .ne. l_arclength(itheta,izeta)) then
                  print *, "NaN found in l_arclength!"
                  stop
              end if
              if (omega_on_arclength_grid(itheta,izeta) .ne. omega_on_arclength_grid(itheta,izeta)) then
                  print *, "NaN found in omega_on_arclength_grid!"
                  stop
              end if
          end do
      end do

      do izeta = 1, nzeta
          do itheta = 1, ntheta
              if (omega_on_arclength_grid(itheta,izeta)>=twopi/2) then
                  omega_on_arclength_grid(itheta,izeta) = omega_on_arclength_grid(itheta,izeta) - twopi
              else if (omega_on_arclength_grid(itheta,izeta)<-twopi/2) then
                  omega_on_arclength_grid(itheta,izeta) = omega_on_arclength_grid(itheta,izeta) + twopi
              end if
          end do
      end do

      norm2l = nfp * dtheta * dzeta * sum(l_arclength**2)
      norm2o = nfp * dtheta * dzeta * sum(omega_on_arclength_grid**2)

    else

      ! Sorting algorithm to ensure continuity when performing splines.
      allocate(x1(ntheta))
      allocate(y1(ntheta))
      allocate(x2(ntheta))
      allocate(y2(ntheta))
      do izeta = 1, nzeta
          x1 = vartheta_on_theta_grid(:,izeta)
          x2 = l(:,izeta)
          y1 = x1
          y2 = x2
          ! Simple swapping of slices of the vartheta array
          do itheta = 1, ntheta-1
              if (x1(itheta+1) < x1(itheta)) then
                  y1(1:ntheta-itheta) = x1(itheta+1:ntheta)
                  y2(1:ntheta-itheta) = x2(itheta+1:ntheta)
                  y1(ntheta-itheta+1:ntheta) = x1(1:itheta)
                  y2(ntheta-itheta+1:ntheta) = x2(1:itheta)
                  exit
              end if
          end do
          ! Bubble sort in case vartheta is still out of order
          do
              sorted = .true.
              do itheta = 1, ntheta-1
                  if (y1(itheta+1) < y1(itheta)) then
                      temp = y1(itheta)
                      y1(itheta) = y1(itheta+1)
                      y1(itheta+1) = temp
                      temp = y2(itheta)
                      y2(itheta) = y2(itheta+1)
                      y2(itheta+1) = temp
                      sorted = .false.
                  end if
              end do
              if (sorted) exit
          end do
          vartheta_on_theta_grid(:,izeta) = y1
          l(:,izeta) = y2
      end do
      deallocate(x1,x2,y1,y2)
      ! End of sorting algorithm

      allocate(l_vartheta(ntheta,nzeta))
      ! Perform splines
      do izeta = 1, nzeta
          call new_periodic_spline(ntheta, vartheta_on_theta_grid(:,izeta), &
              l(:,izeta), twopi, theta_spline)
          do itheta = 1, ntheta
              l_vartheta(itheta,izeta) = periodic_splint(theta(itheta), theta_spline)
          end do
          call delete_periodic_spline(theta_spline)
      end do

      do izeta = 1, nzeta
          do itheta = 1, ntheta
              if (l_vartheta(itheta,izeta) .ne. l_vartheta(itheta,izeta)) then
                  print *, "NaN found in l_vartheta!"
                  stop
              end if
          end do
      end do

      norm2l = nfp * dtheta * dzeta * sum(l_vartheta**2)
    
    end if

    ! Perform Fourier transforms
    mnmax_plasma = (m_max+1)*(2*n_max+1)
    if (allocated(xm_plasma)) deallocate(xm_plasma)
    allocate(xm_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error! regcoil_double_to_single 1'

    if (allocated(xn_plasma)) deallocate(xn_plasma)
    allocate(xn_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error! regcoil_double_to_single 2'

    allocate(lmnc(mnmax_plasma))
    lmnc = 0
    if (lasym) then
      allocate(lmns(mnmax_plasma))
      lmns = 0
    end if
    if (use_arclength_angle) then
      allocate(omns(mnmax_plasma))
      omns = 0
    end if

    cos_norm = twopi**2 / 2
    sin_norm = twopi**2 / 2
    do j=1,2*n_max+1
        do i=1,m_max+1
            imn = (m_max+1)*(j-1) + i
            xm_plasma(imn) = (i-1)
            xn_plasma(imn) = nfp * (j-n_max-1)
            if (xm_plasma(imn)==0 .and. xn_plasma(imn)==0) then
                do izeta = 1, nzeta
                    do itheta = 1, ntheta
                        if (use_arclength_angle) then
                          lmnc(imn) = lmnc(imn) + l_arclength(itheta,izeta) * nfp * (dtheta*dzeta) / twopi**2
                        else
                          lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * nfp * (dtheta*dzeta) / twopi**2
                        end if
                    end do
                end do
            end if
            if (xm_plasma(imn)==0 .and. xn_plasma(imn)<=0) then
                cycle
            end if
            do izeta = 1, nzeta
                do itheta = 1, ntheta
                    angle = xm_plasma(imn) * theta(itheta) - xn_plasma(imn) * zeta(izeta)
                    if (use_arclength_angle) then
                      ! Cosine terms
                      lmnc(imn) = lmnc(imn) + l_arclength(itheta,izeta) * cos(angle) * nfp * (dtheta*dzeta)
                      ! Sine terms
                      if (lasym) then
                          lmns(imn) = lmns(imn) + l_arclength(itheta,izeta) * sin(angle) * nfp * (dtheta*dzeta)
                      end if
                      ! Omega terms
                      omns(imn) = omns(imn) + omega_on_arclength_grid(itheta,izeta) * sin(angle) * nfp * (dtheta*dzeta)
                    else
                      ! Cosine terms
                      lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * cos(angle) * nfp * (dtheta*dzeta)
                      ! Sine terms
                      if (lasym) then
                          lmns(imn) = lmns(imn) + l_vartheta(itheta,izeta) * sin(angle) * nfp * (dtheta*dzeta)
                      end if
                    end if
                end do
            end do
            lmnc(imn) = lmnc(imn) / cos_norm
            if (lasym) then
                lmns(imn) = lmns(imn) / sin_norm
            end if
            if (use_arclength_angle) then
              omns(imn) = omns(imn) / sin_norm
            end if
        end do
    end do

    sum_amp2l = sum(lmnc**2 * cos_norm) + lmnc((m_max+1)*n_max + 1)**2 * cos_norm
    if (lasym) then
        sum_amp2l = sum_amp2l + sum(lmns**2 * sin_norm)
    end if
    if (use_arclength_angle) then
      sum_amp2o = sum(omns**2 * sin_norm)
    end if

    print *, "L2 norm^2 of l: ", norm2l
    print *, "Sum of squared amplitudes of l: ", sum_amp2l
    print *, "Difference: ", norm2l - sum_amp2l
    if (use_arclength_angle) then
      print *, "L2 norm^2 of omega: ", norm2o
      print *, "Sum of squared amplitudes of omega: ", sum_amp2o
      print *, "Difference: ", norm2o - sum_amp2o
    end if

    call system_clock(toc)
    if (verbose) print *,"Done converting to single Fourier coefficients. Took ",real(toc-tic)/countrate," sec."

end subroutine regcoil_double_to_single












