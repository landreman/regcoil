subroutine regcoil_double_to_single()

    !use regcoil_compute_offset_surface_mod
    use regcoil_variables
    !use regcoil_init_Fourier_modes_mod
    use regcoil_splines
    !use regcoil_init_plasma_mod
    use read_wout_mod, only: raxis, zaxis
    !use regcoil_write_single_Fourier

    use stel_kinds
    use stel_constants
    use omp_lib

    implicit none

    integer :: iflag
    !real(dp), dimension(:,:), allocatable :: major_R_coil, z_coil
    !real(dp) :: R0_to_use, relative_variation_in_dl
    !real(dp) :: angle, cos_norm, sin_norm, sinangle, cosangle, dsinangledtheta, dcosangledtheta
    !real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
    real(dp) :: cos_norm, sin_norm, dtheta, dzeta
    !real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
    integer :: i, j, imn, itheta, izeta!, iteration, max_iterations
    !real(dp) :: x_new, y_new, z_new, delta_theta, delta_zeta, temp, factor, factor2
    integer :: tic, toc, countrate, tic1, toc1
    !integer :: mpol_coil, ntor_coil
    !integer, allocatable, dimension(:) :: xm, xn
    !real(dp), allocatable, dimension(:) :: rmnc, rmns, zmnc, zmns
    !real(dp), allocatable, dimension(:) :: theta_plasma_corresponding_to_theta_coil, dl, R_slice, z_slice
    !real(dp), allocatable, dimension(:) :: constant_arclength_theta_on_old_theta_grid
    type (periodic_spline) :: theta_spline
    real(dp), allocatable, dimension(:) :: theta, zeta
    real(dp), allocatable, dimension(:,:) :: R, Z, B, R_0, Z_0
    real(dp), allocatable, dimension(:,:) :: l, vartheta_on_theta_grid, l_vartheta, norm_normal
    real(dp), allocatable, dimension(:,:,:) :: drdtheta, drdzeta, normal
    real(dp), allocatable, dimension(:) :: x1, y1, x2, y2
    real(dp) :: norm2l, sum_amp2
    real(dp) :: temp
    logical :: sorted

    if (geometry_option_plasma .ne. 8) then
        print *,"Invalid setting for geometry_option_plasma: ",geometry_option_plasma
        stop
    end if

    call system_clock(tic,countrate)
    if (verbose) print *,"Converting from double to single Fourier representation of plasma surface."

    ! In case STELLOPT or some other code calls subroutines of regcoil multiple times, make sure arrays are deallocated:
!    if(allocated(xm_coil)) deallocate(xm_coil)
!    if(allocated(xn_coil)) deallocate(xn_coil)
!    if(allocated(rmnc_coil)) deallocate(rmnc_coil)
!    if(allocated(rmns_coil)) deallocate(rmns_coil)
!    if(allocated(zmnc_coil)) deallocate(zmnc_coil)
!    if(allocated(zmns_coil)) deallocate(zmns_coil)

!    nzetal_plasma   = nzeta_plasma   * nfp

!    if (allocated(theta_plasma)) deallocate(theta_plasma)
!    allocate(theta_plasma(ntheta_plasma),stat=iflag)
!    if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 1'
!
!    if (allocated(zeta_plasma)) deallocate(zeta_plasma)
!    allocate(zeta_plasma(nzeta_plasma),stat=iflag)
!    if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 2'
!
!    if (allocated(zetal_plasma)) deallocate(zetal_plasma)
!    allocate(zetal_plasma(nzetal_plasma),stat=iflag)
!    if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 3'
!
    allocate(theta(ntheta_plasma))
    allocate(zeta(nzeta_plasma))

    do i = 1,ntheta_plasma
     theta(i) = twopi*(i-1.0_dp)/ntheta_plasma
    end do

    do i = 1,nzeta_plasma
     zeta(i) = twopi/nfp*(i-1.0_dp)/nzeta_plasma
    end do
!
!    do i = 1,nzetal_plasma
!     zetal_plasma(i) = twopi*(i-1.0_dp)/nzetal_plasma
!    end do

!    dtheta_plasma = theta_plasma(2)-theta_plasma(1)
!    dzeta_plasma  = zeta_plasma(2)-zeta_plasma(1)


    allocate(R(ntheta_plasma,nzeta_plasma))
    allocate(Z(ntheta_plasma,nzeta_plasma))
    allocate(B(ntheta_plasma,nzeta_plasma))
    allocate(norm_normal(ntheta_plasma,nzeta_plasma))

    allocate(drdtheta(3,ntheta_plasma,nzeta_plasma))
    allocate(drdzeta(3,ntheta_plasma,nzeta_plasma))
    allocate(normal(3,ntheta_plasma,nzeta_plasma))

    allocate(R_0(ntheta_plasma,nzeta_plasma))
    allocate(Z_0(ntheta_plasma,nzeta_plasma))
!    R=sqrt(r_plasma(1,:,1:nzeta_plasma)**2+r_plasma(2,:,1:nzeta_plasma)**2)
!    Z=r_plasma(3,:,1:nzeta_plasma)
!        print *,"R(1,1:10): ",r_plasma(1,1,1:10)
!        print *,"R(2,1:10): ",r_plasma(1,2,1:10)
!        print *,"R(3,1:10): ",r_plasma(1,3,1:10)
!        stop
    R=0
    Z=0
    B=0
    drdtheta=0
    drdzeta=0
    R_0=0
    Z_0=0
    do izeta = 1, nzeta_plasma
        angle2 = zeta(izeta)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dzeta = cosangle2
        dcosangle2dzeta = -sinangle2
        do itheta = 1, ntheta_plasma
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

    norm_normal = sqrt(normal(1,:,1:nzeta_plasma)**2 &
        + normal(2,:,1:nzeta_plasma)**2 &
        + normal(3,:,1:nzeta_plasma)**2)

    dtheta = theta(2)-theta(1)
    dzeta = zeta(2)-zeta(1)

    B_0 = sqrt(nfp * dtheta * dzeta * sum(B**2 * norm_normal))

    raxis_cc = raxis(:,1)
    zaxis_cs = zaxis(:,1)
    nmax_axis = size(raxis_cc)
    allocate(xn_axis(nmax_axis))
    do izeta = 1, nzeta_plasma
        do itheta = 1, ntheta_plasma
            do j=1,nmax_axis
                xn_axis(j) = nfp * (j-1)
                angle = -xn_axis(j)*zeta(izeta)

                R_0(itheta,izeta) = R_0(itheta,izeta) + raxis_cc(j)*cos(angle)
                Z_0(itheta,izeta) = Z_0(itheta,izeta) + zaxis_cs(j)*sin(angle)
            end do
        end do
    end do
!        print *,"R_0: ",raxis_cc
!        print *,"Z_0: ",zaxis_cs
!        stop

    allocate(l(ntheta_plasma,nzeta_plasma))
    allocate(vartheta_on_theta_grid(ntheta_plasma,nzeta_plasma))
    l = sqrt((R-R_0)**2 + (Z-Z_0)**2)
    norm2l = nfp * dtheta * dzeta * sum(l**2)
!        print *,"l(1,1:10): ",l(1,1:10)
!        print *,"l(2,1:10): ",l(2,1:10)
!        print *,"l(3,1:10): ",l(3,1:10)
!        stop
    vartheta_on_theta_grid = atan2(Z-Z_0, R-R_0)
        !vartheta_on_theta_grid = vartheta_on_theta_grid + twopi * (vartheta_on_theta_grid<0)
    do izeta = 1, nzeta_plasma
        do itheta = 1, ntheta_plasma
            !l(itheta,izetal)=sqrt((R(itheta,izetal)-R_0(itheta,izetal))^2+(Z(itheta,izetal)-Z_0(itheta,izetal))^2)
            !vartheta_on_theta_grid(itheta,izetal)=atan2((Z(itheta,izetal)-Z_0(itheta,izetal)),(R(itheta,izetal)-R_0(itheta,izetal)))
            if (vartheta_on_theta_grid(itheta,izeta)<0) then
                vartheta_on_theta_grid(itheta,izeta) = vartheta_on_theta_grid(itheta,izeta) + twopi
            end if
        end do
    end do
!        print *,"vartheta_on_theta_grid(1,1:10): ",vartheta_on_theta_grid(1,1:10)
!        print *,"vartheta_on_theta_grid(2,1:10): ",vartheta_on_theta_grid(2,1:10)
!        print *,"vartheta_on_theta_grid(3,1:10): ",vartheta_on_theta_grid(3,1:10)
!        stop

        ! Sorting algorithm to ensure continuity when performing splines.
    allocate(x1(ntheta_plasma))
    allocate(y1(ntheta_plasma))
    allocate(x2(ntheta_plasma))
    allocate(y2(ntheta_plasma))
    do izeta = 1, nzeta_plasma
        x1 = vartheta_on_theta_grid(:,izeta)
        x2 = l(:,izeta)
        y1 = x1
        y2 = x2
        ! Simple swapping of slices of the vartheta array
        do itheta = 1, ntheta_plasma-1
            if (x1(itheta+1) < x1(itheta)) then
                y1(1:ntheta_plasma-itheta) = x1(itheta+1:ntheta_plasma)
                y2(1:ntheta_plasma-itheta) = x2(itheta+1:ntheta_plasma)
                y1(ntheta_plasma-itheta+1:ntheta_plasma) = x1(1:itheta)
                y2(ntheta_plasma-itheta+1:ntheta_plasma) = x2(1:itheta)
                exit
            end if
        end do
        ! Bubble sort in case vartheta is still out of order
        do
            sorted = .true.
            do itheta = 1, ntheta_plasma-1
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

!        print *,"vartheta_on_theta_grid(1:20,1): ",vartheta_on_theta_grid(1:20,1)
!        print *,"vartheta_on_theta_grid(1:20,2): ",vartheta_on_theta_grid(1:20,2)
!        print *,"vartheta_on_theta_grid(1:20,3): ",vartheta_on_theta_grid(1:20,3)
!        stop

    allocate(l_vartheta(ntheta_plasma,nzeta_plasma))
    ! Perform splines
    do izeta = 1, nzeta_plasma
        call new_periodic_spline(ntheta_plasma, vartheta_on_theta_grid(:,izeta), &
            l(:,izeta), twopi, theta_spline)
        do itheta = 1, ntheta_plasma
            l_vartheta(itheta,izeta) = periodic_splint(theta(itheta), theta_spline)
        end do
        call delete_periodic_spline(theta_spline)
    end do

    do izeta = 1, nzeta_plasma
        do itheta = 1, ntheta_plasma
            if (l_vartheta(itheta,izeta)/=l_vartheta(itheta,izeta)) then
                print *, "NaN found in l_vartheta!"
                stop
            end if
        end do
    end do

!        print *,"l_vartheta(10:30,1): ",l_vartheta(10:30,1)
!        print *,"l_vartheta(10:30,2): ",l_vartheta(10:30,2)
!        print *,"l_vartheta(10:30,3): ",l_vartheta(:,3)
!        stop

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
!        allocate(cos_norm(m_max+1,2*n_max+1))
!        allocate(sin_norm(m_max+1,2*n_max+1))
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
                do izeta = 1, nzeta_plasma
                    do itheta = 1, ntheta_plasma
                        lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * nfp * (dtheta*dzeta) / twopi**2
                    end do
                end do
            end if
            if (xm_plasma(imn)==0 .and. xn_plasma(imn)<=0) then
                cycle
            end if
            do izeta = 1, nzeta_plasma
                do itheta = 1, ntheta_plasma
                    angle = xm_plasma(imn) * theta(itheta) - xn_plasma(imn) * zeta(izeta)

                    ! Cosine terms
                    lmnc(imn) = lmnc(imn) + l_vartheta(itheta,izeta) * cos(angle) * nfp * (dtheta*dzeta)
                    !cos_norm = cos_norm + (cos(angle))**2 * (dtheta_plasma*dzeta_plasma)
                    ! Sine terms
                    lmns(imn) = lmns(imn) + l_vartheta(itheta,izeta) * sin(angle) * nfp * (dtheta*dzeta)
                    !sin_norm = sin_norm + (sin(angle))**2 * (dtheta_plasma*dzeta_plasma)
                end do
            end do
            lmnc(imn) = lmnc(imn) / cos_norm
            lmns(imn) = lmns(imn) / sin_norm
        end do
    end do

    sum_amp2 = sum(lmnc**2 * cos_norm + lmns**2 * sin_norm) + lmnc((m_max+1)*n_max + 1)**2 * cos_norm

    print *, "L2 norm of l: ", norm2l
    print *, "Sum of squared amplitudes of l: ", sum_amp2
    ! Handle the 0th cosine mode
!    do izeta = 1, nzeta_plasma
!        do itheta = 1, ntheta_plasma
!            lmnc(1,n_max+1) = lmnc(1,n_max+1) + l_vartheta(itheta,izeta) * nfp * (dtheta_plasma*dzeta_plasma) / twopi**2
!        end do
!    end do

    deallocate(R, Z, R_0, Z_0, l, vartheta_on_theta_grid, l_vartheta)

!    print *,"lmnc(1,n_max+1:2*n_max+1): ", lmnc(1,n_max+1:2*n_max+1)
!    print *,"lmns(1,n_max+1:2*n_max+1): ", lmns(1,n_max+1:2*n_max+1)
!    stop

    call regcoil_write_single_Fourier()

    call system_clock(toc)
    if (verbose) print *,"Done converting to single Fourier coefficients. Took ",real(toc-tic)/countrate," sec."

end subroutine regcoil_double_to_single












