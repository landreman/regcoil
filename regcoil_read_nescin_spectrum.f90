module regcoil_nescin_utils_1

contains

  subroutine read_nescin_spectrum(nescin_filename, lscreen_optin)
  ! Modification of the original read_nescin/init_surface_mod subroutines.
  ! The spectrum is read in and returned via the passed-in variables, which are
  ! allocated during this call.

    use safe_open_mod
    use stel_constants
    use stel_kinds
    use regcoil_variables, only:  &
         ntheta_coil, nzeta_coil, nzetal_coil, theta_coil, zeta_coil, zetal_coil, &
         r_coil, drdtheta_coil, drdzeta_coil, &
         normal_coil, norm_normal_coil, area_coil, volume_coil, &
         R0_coil, a_coil, &
         dtheta_coil, dzeta_coil, nfp, &
         rc_rmnc_stellopt, rc_rmns_stellopt, &
         rc_zmnc_stellopt, rc_zmns_stellopt
    implicit none

    character(*) :: nescin_filename
    integer :: itheta, izeta, nummodes, xm_in, xn_in, istat, ii, jj, kk
    integer :: iflag, which_surface, iunit  = 7
    !real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta
    !real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
    !real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
    !real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_theta, delta_zeta, temp
    real(dp) :: rmnc_in, zmns_in, rmns_in, zmnc_in

    character(300) :: myline
    character(*), parameter :: matchString = "------ Current Surface"

    ! variables to handle printing to the screen
    logical, optional :: lscreen_optin
    logical :: lscreen
    
    if(present(lscreen_optin)) then 
      lscreen = lscreen_optin
    else
      lscreen = .true.
    endif

    call safe_open(iunit, istat, trim(nescin_filename), 'old', 'formatted')
    if (istat .ne. 0) then
       stop 'Error opening nescin file'
    endif

    ! Skip down to the line in the file that begins with matchString
    do
      read (iunit,"(a)") myline
      if (myline(:len(matchString)) == matchString) then
        exit
      end if
    end do
    read (iunit, *)

    read (iunit, *) nummodes

    ! STELLOPT may not have called regcoil_init_plasma yet.
    nzetal_coil = nzeta_coil * nfp

    read (iunit, *)
    read (iunit, *)
    do ii = 1, nummodes
       ! (mu + nv) in NESCOIL is (xm * theta + xn * zeta) in REGCOIL
       ! Also note the order: m,n, RmnC, ZmnS, RmnS, ZmnC
       read (iunit, *) xm_in, xn_in, rmnc_in, zmns_in, rmns_in, zmnc_in
       rc_rmnc_stellopt(xm_in, xn_in) = rmnc_in
       rc_zmns_stellopt(xm_in, xn_in) = zmns_in
       rc_rmns_stellopt(xm_in, xn_in) = rmns_in
       rc_zmnc_stellopt(xm_in, xn_in) = zmnc_in
    end do

    close(iunit)

        if (allocated(theta_coil)) then
          if (lscreen) print *,"<----Deallocating theta_coil"
          deallocate(theta_coil)
        end if
        allocate(theta_coil(ntheta_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin_spectrum 1'

        if (allocated(zeta_coil)) then
          if (lscreen) print *,"<----Deallocating zeta_coil"
          deallocate(zeta_coil)
        end if
        allocate(zeta_coil(nzeta_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin_spectrum 2'

        if (allocated(zetal_coil)) then
          if (lscreen) print *,"<----Deallocating zetal_coil"
          deallocate(zetal_coil)
        end if
        allocate(zetal_coil(nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin_spectrum 3'

        if (allocated(r_coil)) then
          if (lscreen) print *,"<----Deallocating r_coil"
          deallocate(r_coil)
        end if
        allocate(r_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin_spectrum 0'

        if (allocated(normal_coil)) then
          if (lscreen) print *,"<----Deallocating normal_coil"
          deallocate(normal_coil)
        end if
        allocate(normal_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin_spectrum 7'

        if (allocated(drdtheta_coil)) deallocate(drdtheta_coil)
        allocate(drdtheta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin 51'

        if (allocated(drdzeta_coil)) deallocate(drdzeta_coil)
        allocate(drdzeta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin 61'

        if (allocated(norm_normal_coil)) deallocate(norm_normal_coil)
        allocate(norm_normal_coil(ntheta_coil,nzeta_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error! regcoil_read_nescin 71'

        norm_normal_coil = 0
        r_coil = 0
        normal_coil = 0
        drdtheta_coil = 0
        drdzeta_coil = 0

        do ii = 1,ntheta_coil
           theta_coil(ii) = twopi*(ii-1.0_dp)/ntheta_coil
        end do
        do ii = 1,nzeta_coil
           zeta_coil(ii) = twopi/nfp*(ii-1.0_dp)/nzeta_coil
        end do

        do ii = 1,nzetal_coil
           zetal_coil(ii) = twopi*(ii-1.0_dp)/nzetal_coil
        end do

        dtheta_coil = theta_coil(2)-theta_coil(1)
        dzeta_coil  = zeta_coil(2)-zeta_coil(1)
  
  end subroutine read_nescin_spectrum
  
  
end module regcoil_nescin_utils_1
