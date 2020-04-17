subroutine regcoil_double_to_single()

    use regcoil_variables
    use regcoil_splines
    use read_wout_mod, only: raxis, zaxis

    use stel_kinds
    use stel_constants
    use omp_lib

    implicit none

    integer :: iflag
    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
    real(dp) :: cos_norm, sin_norm, dtheta, dzeta
    integer :: i, j, imn, itheta, izeta
    integer :: tic, toc, countrate, tic1, toc1
    type (periodic_spline) :: theta_spline
    real(dp), allocatable, dimension(:) :: theta, zeta
    real(dp), allocatable, dimension(:,:) :: R, Z, B, R_0, Z_0
    real(dp), allocatable, dimension(:,:) :: l, vartheta_on_theta_grid, l_vartheta, norm_normal
    real(dp), allocatable, dimension(:,:,:) :: drdtheta, drdzeta, normal
    real(dp), allocatable, dimension(:) :: x1, y1, x2, y2
    integer :: ntheta, nzeta
    real(dp) :: norm2l, sum_amp2
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

    allocate(drdtheta(3,ntheta,nzeta))
    allocate(drdzeta(3,ntheta,nzeta))
    allocate(normal(3,ntheta,nzeta))

    allocate(R_0(ntheta,nzeta))
    allocate(Z_0(ntheta,nzeta))

    R=0
    Z=0
    B=0
    drdtheta=0
    drdzeta=0
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

                R(itheta,izeta) = R(itheta,izeta) + rmnc_plasma(imn)*cosangle
                Z(itheta,izeta) = Z(itheta,izeta) + zmns_plasma(imn)*sinangle
                B(itheta,izeta) = B(itheta,izeta) + bmnc_plasma(imn)*cosangle

                drdtheta(1,itheta,izeta) = drdtheta(1,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * cosangle2
                drdtheta(2,itheta,izeta) = drdtheta(2,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * sinangle2
                drdtheta(3,itheta,izeta) = drdtheta(3,itheta,izeta) + zmns_plasma(imn) * dsinangledtheta

                drdzeta(1,itheta,izeta) = drdzeta(1,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta)
                drdzeta(2,itheta,izeta) = drdzeta(2,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta)
                drdzeta(3,itheta,izeta) = drdzeta(3,itheta,izeta) + zmns_plasma(imn) * dsinangledzeta
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
                end if
            end do
        end do
    end do

    ! Evaluate cross product
    normal(1,:,:) = drdzeta(2,:,:) * drdtheta(3,:,:) - drdtheta(2,:,:) * drdzeta(3,:,:)
    normal(2,:,:) = drdzeta(3,:,:) * drdtheta(1,:,:) - drdtheta(3,:,:) * drdzeta(1,:,:)
    normal(3,:,:) = drdzeta(1,:,:) * drdtheta(2,:,:) - drdtheta(1,:,:) * drdzeta(2,:,:)

    norm_normal = sqrt(normal(1,:,1:nzeta)**2 &
        + normal(2,:,1:nzeta)**2 &
        + normal(3,:,1:nzeta)**2)

    dtheta = theta(2)-theta(1)
    dzeta = zeta(2)-zeta(1)

    B_0 = sqrt(sum(B**2 * norm_normal) / sum(norm_normal))
!print *,B_0
!stop

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
            if (l_vartheta(itheta,izeta)/=l_vartheta(itheta,izeta)) then
                print *, "NaN found in l_vartheta!"
                stop
            end if
        end do
    end do

    norm2l = nfp * dtheta * dzeta * sum(l_vartheta**2)

    ! Perform Fourier transforms
    mnmax_plasma = (m_max+1)*(2*n_max+1)
    if (allocated(xm_plasma)) deallocate(xm_plasma)
    allocate(xm_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error! regcoil_double_to_single 1'

    if (allocated(xn_plasma)) deallocate(xn_plasma)
    allocate(xn_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error! regcoil_double_to_single 2'

    allocate(lmnc(mnmax_plasma))
    allocate(lmns(mnmax_plasma))

    lmnc = 0
    lmns = 0
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
                        lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * nfp * (dtheta*dzeta) / twopi**2
                    end do
                end do
            end if
            if (xm_plasma(imn)==0 .and. xn_plasma(imn)<=0) then
                cycle
            end if
            do izeta = 1, nzeta
                do itheta = 1, ntheta
                    angle = xm_plasma(imn) * theta(itheta) - xn_plasma(imn) * zeta(izeta)

                    ! Cosine terms
                    lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * cos(angle) * nfp * (dtheta*dzeta)
                    ! Sine terms
                    if (lasym) then
                        lmns(imn) = lmns(imn) + l_vartheta(itheta,izeta) * sin(angle) * nfp * (dtheta*dzeta)
                    end if
                end do
            end do
            lmnc(imn) = lmnc(imn) / cos_norm
            if (lasym) then
                lmns(imn) = lmns(imn) / sin_norm
            end if
        end do
    end do

    sum_amp2 = sum(lmnc**2 * cos_norm) + lmnc((m_max+1)*n_max + 1)**2 * cos_norm
    if (lasym) then
        sum_amp2 = sum_amp2 + sum(lmns**2 * sin_norm)
    end if

    print *, "L2 norm^2 of l: ", norm2l
    print *, "Sum of squared amplitudes of l: ", sum_amp2
    print *, "Difference: ", norm2l - sum_amp2

    call system_clock(toc)
    if (verbose) print *,"Done converting to single Fourier coefficients. Took ",real(toc-tic)/countrate," sec."

end subroutine regcoil_double_to_single












