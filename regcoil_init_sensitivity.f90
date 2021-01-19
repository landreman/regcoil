subroutine regcoil_init_sensitivity()

  use regcoil_variables
  use stel_constants
  use regcoil_init_Fourier_modes_mod
  use omp_lib

  implicit none

  integer :: iomega, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry, &
    tic, toc, countrate,indexl_coil
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta, &
    d2xdtheta2, d2xdthetadzeta, d2xdzeta2, d2ydtheta2, d2ydthetadzeta, d2ydzeta2, d2zdtheta2, d2zdthetadzeta, d2zdzeta2, &
    major_R_squared_half_grid, dzdtheta_coil_half_grid, sum_exp_min, sum_exp_max
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta, d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2, d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta, d2sinangle2dzeta2, d2cosangle2dzeta2
  real(dp) :: angle3, sinangle3, cosangle3, dsinangle3dtheta, dsinangle3dzeta, dcosangle3dtheta, dcosangle3dzeta, d2sinangle3dtheta2, d2sinangle3dthetadzeta, d2sinangle3dzeta2, d2cosangle3dtheta2, d2cosangle3dthetadzeta, d2cosangle3dzeta2
  integer :: izeta_plasma, itheta_plasma, izeta_coil, itheta_coil, index_coil, &
    l_coil, l_plasma, izetal_coil, izetal_plasma, index_plasma, index_coil_first, index_coil_last
  real(dp), dimension(:), allocatable :: dR_squared_domega_half_grid, dsum_exp_mindomega
  real(dp), dimension(:,:,:), allocatable :: dmajor_R_squareddomega
  real(dp), dimension(:,:,:,:), allocatable :: dist,ones,normal_coil_fourd,normal_plasma_fourd


  if (sensitivity_option == 6) then !----- Plasma derivatives -----

      ! Initialize Fourier arrays - need to include m = 0, n = 0 mode
      call regcoil_init_Fourier_modes_sensitivity(mmin_sensitivity, nmin_sensitivity,mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, nomega_plasma, &
        xm_sensitivity, xn_sensitivity, omega_plasma,sensitivity_symmetry_option,sensitivity_option)

      allocate(darea_plasmadomega(nomega_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdomega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(d2rdthetadomega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(d2rdzetadomega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      allocate(dnorm_normaldomega(nomega_plasma,ntheta_plasma,nzeta_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      
      allocate(dnormxdomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnormydomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnormzdomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      
      if (geometry_option_coil == 5) then
        allocate(d3rdtheta2domega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d3rdzeta2domega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d3rdthetadzetadomega(3,ntheta_plasma*nzeta_plasma,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'

        allocate(d2normxdthetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2normydthetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2normzdthetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2normxdzetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2normydzetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2normzdzetadomega(nomega_plasma,ntheta_plasma*nzeta_plasma,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'

        allocate(d2norm_normaldthetadomega(nomega_plasma,ntheta_plasma,nzeta_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2norm_normaldzetadomega(nomega_plasma,ntheta_plasma,nzeta_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'

        allocate(darea_coildomega(nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(drdomega_coil(3,ntheta_coil*nzeta_coil,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2rdthetadomega_coil(3,ntheta_coil*nzeta_coil,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(d2rdzetadomega_coil(3,ntheta_coil*nzeta_coil,nfp,nomega_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'

        allocate(dnorm_normaldomega_coil(nomega_plasma,ntheta_coil,nzeta_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        
        allocate(dnormxdomega_coil(nomega_plasma,ntheta_coil*nzeta_coil,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(dnormydomega_coil(nomega_plasma,ntheta_coil*nzeta_coil,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(dnormzdomega_coil(nomega_plasma,ntheta_coil*nzeta_coil,nfp),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'

        allocate(dddomega(3, nomega_plasma,ntheta_coil*nzetal_coil),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
      end if

  else !----- Coil derivatives -----

      ! Initialize Fourier arrays - need to include m = 0, n = 0 mode
      call regcoil_init_Fourier_modes_sensitivity(mmin_sensitivity, nmin_sensitivity,mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, nomega_coil, &
        xm_sensitivity, xn_sensitivity, omega_coil,sensitivity_symmetry_option,sensitivity_option)

      allocate(darea_coildomega(nomega_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdomega(3,ntheta_coil*nzeta_coil,nfp,nomega_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnorm_normaldomega(nomega_coil,ntheta_coil,nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dddomega(3, nomega_coil,ntheta_coil*nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnormxdomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnormydomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dnormzdomega(nomega_coil,ntheta_coil*nzeta_coil,nfp),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadxdtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadxdzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadydtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadydzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadzdtheta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(domegadzdzeta(nomega_coil,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

  end if

  call system_clock(tic,countrate)

  if (sensitivity_option == 6) then !----- Plasma derivatives -----

      do izeta_plasma = 1,nzeta_plasma
        do itheta_plasma = 1,ntheta_plasma
          index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
          do l_plasma = 0, (nfp-1)
            izetal_plasma = izeta_plasma + l_plasma*nzeta_plasma

            angle2 = zetal_plasma(izetal_plasma)
            sinangle2 = sin(angle2)
            cosangle2 = cos(angle2)
            dsinangle2dzeta = cosangle2
            dcosangle2dzeta = -sinangle2
            d2sinangle2dzeta2 = -sinangle2
            d2cosangle2dzeta2 = -cosangle2
            if (use_arclength_angle) then
              angle3 = theta_plasma(itheta_plasma) - omega_arclength(itheta_plasma,izetal_plasma)
              sinangle3 = sin(angle3)
              cosangle3 = cos(angle3)
              dsinangle3dtheta = cosangle3 * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma))
              dcosangle3dtheta = -sinangle3 * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma))
              dsinangle3dzeta = cosangle3 * (- domegadzeta_arclength(itheta_plasma,izetal_plasma))
              dcosangle3dzeta = -sinangle3 * (- domegadzeta_arclength(itheta_plasma,izetal_plasma))
              if (geometry_option_coil == 5) then
                d2sinangle3dtheta2 = dcosangle3dtheta * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma)) + cosangle3 * (- d2omegadtheta2_arclength(itheta_plasma,izetal_plasma))
                d2cosangle3dtheta2 = -dsinangle3dtheta * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma)) - sinangle3 * (- d2omegadtheta2_arclength(itheta_plasma,izetal_plasma))
                d2sinangle3dthetadzeta = dcosangle3dzeta * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma)) + cosangle3 * (- d2omegadthetadzeta_arclength(itheta_plasma,izetal_plasma))
                d2cosangle3dthetadzeta = -dsinangle3dzeta * (1 - domegadtheta_arclength(itheta_plasma,izetal_plasma)) - sinangle3 * (- d2omegadthetadzeta_arclength(itheta_plasma,izetal_plasma))
                d2sinangle3dzeta2 = dcosangle3dzeta * (- domegadzeta_arclength(itheta_plasma,izetal_plasma)) + cosangle3 * (- d2omegadzeta2_arclength(itheta_plasma,izetal_plasma))
                d2cosangle3dzeta2 = -dsinangle3dzeta * (- domegadzeta_arclength(itheta_plasma,izetal_plasma)) - sinangle3 * (- d2omegadzeta2_arclength(itheta_plasma,izetal_plasma))
              end if
            else
              angle3 = theta_plasma(itheta_plasma)
              sinangle3 = sin(angle3)
              cosangle3 = cos(angle3)
              dsinangle3dtheta = cosangle3
              dcosangle3dtheta = -sinangle3
              dsinangle3dzeta = 0
              dcosangle3dzeta = 0
              d2sinangle3dtheta2 = dcosangle3dtheta
              d2cosangle3dtheta2 = -dsinangle3dtheta
              d2sinangle3dthetadzeta = 0
              d2cosangle3dthetadzeta = 0
              d2sinangle3dzeta2 = 0
              d2cosangle3dzeta2 = 0
            end if

            do iomega = 1,nomega_plasma
              ! Opposite convention for angle as used by nescin input file
              angle = xm_sensitivity(iomega)*theta_plasma(itheta_plasma) - nfp*xn_sensitivity(iomega)*zetal_plasma(izetal_plasma)
              sinangle = sin(angle)
              cosangle = cos(angle)
              dsinangledtheta = xm_sensitivity(iomega)*cosangle
              dsinangledzeta  = -nfp*xn_sensitivity(iomega)*cosangle
              dcosangledtheta = -xm_sensitivity(iomega)*sinangle
              dcosangledzeta  = nfp*xn_sensitivity(iomega)*sinangle
              d2sinangledtheta2 = -sinangle*xm_sensitivity(iomega)**2
              d2cosangledtheta2 = -cosangle*xm_sensitivity(iomega)**2
              d2sinangledthetadzeta = sinangle*xm_sensitivity(iomega)*nfp*xn_sensitivity(iomega)
              d2cosangledthetadzeta = cosangle*xm_sensitivity(iomega)*nfp*xn_sensitivity(iomega)
              d2sinangledzeta2 = -sinangle*(nfp*xn_sensitivity(iomega))**2
              d2cosangledzeta2 = -cosangle*(nfp*xn_sensitivity(iomega))**2

              ! omega = lmnc / rmnc
              
              drdomega(1,index_plasma,l_plasma+1,iomega) = cosangle * cosangle3 * cosangle2
              drdomega(2,index_plasma,l_plasma+1,iomega) = cosangle * cosangle3 * sinangle2
              drdomega(3,index_plasma,l_plasma+1,iomega) = cosangle * sinangle3

              d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) = (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * cosangle2
              d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) = (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * sinangle2
              d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) = (dcosangledtheta * sinangle3 + cosangle * dsinangle3dtheta)
              
              d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) = (dcosangledzeta * cosangle3 * cosangle2 + cosangle * dcosangle3dzeta * cosangle2 + cosangle * cosangle3 * dcosangle2dzeta)
              d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) = (dcosangledzeta * cosangle3 * sinangle2 + cosangle * dcosangle3dzeta * sinangle2 + cosangle * cosangle3 * dsinangle2dzeta)
              d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) = (dcosangledzeta * sinangle3 + cosangle * dsinangle3dzeta)

              if (geometry_option_coil == 5) then
                d3rdtheta2domega(1,index_plasma,l_plasma+1,iomega) = (d2cosangledtheta2 * cosangle3 + 2*dcosangledtheta * dcosangle3dtheta + cosangle * d2cosangle3dtheta2) * cosangle2
                d3rdtheta2domega(2,index_plasma,l_plasma+1,iomega) = (d2cosangledtheta2 * cosangle3 + 2*dcosangledtheta * dcosangle3dtheta + cosangle * d2cosangle3dtheta2) * sinangle2
                d3rdtheta2domega(3,index_plasma,l_plasma+1,iomega) = (d2cosangledtheta2 * sinangle3 + 2*dcosangledtheta * dsinangle3dtheta + cosangle * d2sinangle3dtheta2)

                d3rdthetadzetadomega(1,index_plasma,l_plasma+1,iomega) = ( (d2cosangledthetadzeta * cosangle3 + dcosangledtheta * dcosangle3dzeta + dcosangledzeta * dcosangle3dtheta + cosangle * d2cosangle3dthetadzeta) * cosangle2 &
                 + (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * dcosangle2dzeta )
                d3rdthetadzetadomega(2,index_plasma,l_plasma+1,iomega) = ( (d2cosangledthetadzeta * cosangle3 + dcosangledtheta * dcosangle3dzeta + dcosangledzeta * dcosangle3dtheta + cosangle * d2cosangle3dthetadzeta) * sinangle2 &
                 + (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * dsinangle2dzeta )
                d3rdthetadzetadomega(3,index_plasma,l_plasma+1,iomega) = (d2cosangledthetadzeta * sinangle3 + dcosangledtheta * dsinangle3dzeta + dcosangledzeta * dsinangle3dtheta + cosangle * d2sinangle3dthetadzeta)

                d3rdzeta2domega(1,index_plasma,l_plasma+1,iomega) = ( (d2cosangledzeta2 * cosangle3 + 2*dcosangledzeta * dcosangle3dzeta + cosangle * d2cosangle3dzeta2) * cosangle2 + 2*(dcosangledzeta * cosangle3 + cosangle * dcosangle3dzeta)*dcosangle2dzeta + cosangle * cosangle3 * d2cosangle2dzeta2 )
                d3rdzeta2domega(2,index_plasma,l_plasma+1,iomega) = ( (d2cosangledzeta2 * cosangle3 + 2*dcosangledzeta * dcosangle3dzeta + cosangle * d2cosangle3dzeta2) * sinangle2 + 2*(dcosangledzeta * cosangle3 + cosangle * dcosangle3dzeta)*dsinangle2dzeta + cosangle * cosangle3 * d2sinangle2dzeta2 )
                d3rdzeta2domega(3,index_plasma,l_plasma+1,iomega) = (d2cosangledzeta2 * sinangle3 + 2*dcosangledzeta * dsinangle3dzeta + cosangle * d2sinangle3dzeta2)
              end if

              dxdtheta = drdtheta_plasma(1,itheta_plasma,izetal_plasma)
              dydtheta = drdtheta_plasma(2,itheta_plasma,izetal_plasma)
              dzdtheta = drdtheta_plasma(3,itheta_plasma,izetal_plasma)
              dxdzeta = drdzeta_plasma(1,itheta_plasma,izetal_plasma)
              dydzeta = drdzeta_plasma(2,itheta_plasma,izetal_plasma)
              dzdzeta = drdzeta_plasma(3,itheta_plasma,izetal_plasma)
              if (geometry_option_coil == 5) then
                d2xdtheta2 = d2rdtheta2_plasma(1,itheta_plasma,izetal_plasma)
                d2ydtheta2 = d2rdtheta2_plasma(2,itheta_plasma,izetal_plasma)
                d2zdtheta2 = d2rdtheta2_plasma(3,itheta_plasma,izetal_plasma)
                d2xdthetadzeta = d2rdthetadzeta_plasma(1,itheta_plasma,izetal_plasma)
                d2ydthetadzeta = d2rdthetadzeta_plasma(2,itheta_plasma,izetal_plasma)
                d2zdthetadzeta = d2rdthetadzeta_plasma(3,itheta_plasma,izetal_plasma)
                d2xdzeta2 = d2rdzeta2_plasma(1,itheta_plasma,izetal_plasma)
                d2ydzeta2 = d2rdzeta2_plasma(2,itheta_plasma,izetal_plasma)
                d2zdzeta2 = d2rdzeta2_plasma(3,itheta_plasma,izetal_plasma)
              end if
              dnormxdomega(iomega,index_plasma,l_plasma+1) = &
                  d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) * dzdtheta + dydzeta * d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) &
                - d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) * dzdzeta - dydtheta * d2rdzetadomega(3,index_plasma,l_plasma+1,iomega)
              dnormydomega(iomega,index_plasma,l_plasma+1) = &
                  d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) * dxdtheta + dzdzeta * d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) &
                - d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) * dxdzeta - dzdtheta * d2rdzetadomega(1,index_plasma,l_plasma+1,iomega)
              dnormzdomega(iomega,index_plasma,l_plasma+1) = &
                  d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) * dydtheta + dxdzeta * d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) &
                - d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) * dydzeta - dxdtheta * d2rdzetadomega(2,index_plasma,l_plasma+1,iomega)

              if (geometry_option_coil == 5) then
                d2normxdthetadomega(iomega,index_plasma,l_plasma+1) = &
                    d3rdthetadzetadomega(2,index_plasma,l_plasma+1,iomega) * dzdtheta + d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) * d2zdtheta2 &
                  - d3rdtheta2domega(2,index_plasma,l_plasma+1,iomega) * dzdzeta - d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) * d2zdthetadzeta &
                  + d2ydthetadzeta * d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) + dydzeta * d3rdtheta2domega(3,index_plasma,l_plasma+1,iomega) &
                  - d2ydtheta2 * d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) - dydtheta * d3rdthetadzetadomega(3,index_plasma,l_plasma+1,iomega)
                d2normydthetadomega(iomega,index_plasma,l_plasma+1) = &
                   d3rdthetadzetadomega(3,index_plasma,l_plasma+1,iomega) * dxdtheta + d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) * d2xdtheta2 &
                 - d3rdtheta2domega(3,index_plasma,l_plasma+1,iomega) * dxdzeta - d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) * d2xdthetadzeta &
                 + d2zdthetadzeta * d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) + dzdzeta * d3rdtheta2domega(1,index_plasma,l_plasma+1,iomega) &
                 - d2zdtheta2 * d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) - dzdtheta * d3rdthetadzetadomega(1,index_plasma,l_plasma+1,iomega)
                d2normzdthetadomega(iomega,index_plasma,l_plasma+1) = &
                   d3rdthetadzetadomega(1,index_plasma,l_plasma+1,iomega) * dydtheta + d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) * d2ydtheta2 &
                 - d3rdtheta2domega(1,index_plasma,l_plasma+1,iomega) * dydzeta - d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) * d2ydthetadzeta &
                 + d2xdthetadzeta * d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) + dxdzeta * d3rdtheta2domega(2,index_plasma,l_plasma+1,iomega) &
                 - d2xdtheta2 * d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) - dxdtheta * d3rdthetadzetadomega(2,index_plasma,l_plasma+1,iomega)
                
                d2normxdzetadomega(iomega,index_plasma,l_plasma+1) = &
                    d3rdzeta2domega(2,index_plasma,l_plasma+1,iomega) * dzdtheta + d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) * d2zdthetadzeta &
                  - d3rdthetadzetadomega(2,index_plasma,l_plasma+1,iomega) * dzdzeta - d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) * d2zdzeta2 &
                  + d2ydzeta2 * d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) + dydzeta * d3rdthetadzetadomega(3,index_plasma,l_plasma+1,iomega) &
                  - d2ydthetadzeta * d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) - dydtheta * d3rdzeta2domega(3,index_plasma,l_plasma+1,iomega)
                d2normydzetadomega(iomega,index_plasma,l_plasma+1) = &
                    d3rdzeta2domega(3,index_plasma,l_plasma+1,iomega) * dxdtheta + d2rdzetadomega(3,index_plasma,l_plasma+1,iomega) * d2xdthetadzeta &
                  - d3rdthetadzetadomega(3,index_plasma,l_plasma+1,iomega) * dxdzeta - d2rdthetadomega(3,index_plasma,l_plasma+1,iomega) * d2xdzeta2 &
                  + d2zdzeta2 * d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) + dzdzeta * d3rdthetadzetadomega(1,index_plasma,l_plasma+1,iomega) &
                  - d2zdthetadzeta * d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) - dzdtheta * d3rdzeta2domega(1,index_plasma,l_plasma+1,iomega)
                d2normzdzetadomega(iomega,index_plasma,l_plasma+1) = &
                    d3rdzeta2domega(1,index_plasma,l_plasma+1,iomega) * dydtheta + d2rdzetadomega(1,index_plasma,l_plasma+1,iomega) * d2ydthetadzeta &
                  - d3rdthetadzetadomega(1,index_plasma,l_plasma+1,iomega) * dydzeta - d2rdthetadomega(1,index_plasma,l_plasma+1,iomega) * d2ydzeta2 &
                  + d2xdzeta2 * d2rdthetadomega(2,index_plasma,l_plasma+1,iomega) + dxdzeta * d3rdthetadzetadomega(2,index_plasma,l_plasma+1,iomega) &
                  - d2xdthetadzeta * d2rdzetadomega(2,index_plasma,l_plasma+1,iomega) - dxdtheta * d3rdzeta2domega(2,index_plasma,l_plasma+1,iomega)
              end if
            enddo
          enddo
          dnorm_normaldomega(:,itheta_plasma,izeta_plasma) = &
             (normal_plasma(1, itheta_plasma, izeta_plasma)*dnormxdomega(:,index_plasma,1) &
            + normal_plasma(2, itheta_plasma, izeta_plasma)*dnormydomega(:,index_plasma,1) &
            + normal_plasma(3, itheta_plasma, izeta_plasma)*dnormzdomega(:,index_plasma,1)) &
              /norm_normal_plasma(itheta_plasma,izeta_plasma)
          
          if (geometry_option_coil == 5) then
            d2norm_normaldthetadomega(:,itheta_plasma,izeta_plasma) = ( &
                d2normxdthetadomega(:,index_plasma,1) * normal_plasma(1,itheta_plasma,izeta_plasma) &
              + d2normydthetadomega(:,index_plasma,1) * normal_plasma(2,itheta_plasma,izeta_plasma) &
              + d2normzdthetadomega(:,index_plasma,1) * normal_plasma(3,itheta_plasma,izeta_plasma) &
              + dnormaldtheta_plasma(1,itheta_plasma,izeta_plasma) * dnormxdomega(:,index_plasma,1) &
              + dnormaldtheta_plasma(2,itheta_plasma,izeta_plasma) * dnormydomega(:,index_plasma,1) &
              + dnormaldtheta_plasma(3,itheta_plasma,izeta_plasma) * dnormzdomega(:,index_plasma,1) &
              - dnorm_normaldtheta_plasma(itheta_plasma,izeta_plasma) * dnorm_normaldomega(:,itheta_plasma,izeta_plasma)  )/norm_normal_plasma(itheta_plasma,izeta_plasma)
            d2norm_normaldzetadomega(:,itheta_plasma,izeta_plasma) = ( &
                d2normxdzetadomega(:,index_plasma,1) * normal_plasma(1,itheta_plasma,izeta_plasma) &
              + d2normydzetadomega(:,index_plasma,1) * normal_plasma(2,itheta_plasma,izeta_plasma) &
              + d2normzdzetadomega(:,index_plasma,1) * normal_plasma(3,itheta_plasma,izeta_plasma) &
              + dnormaldzeta_plasma(1,itheta_plasma,izeta_plasma) * dnormxdomega(:,index_plasma,1) &
              + dnormaldzeta_plasma(2,itheta_plasma,izeta_plasma) * dnormydomega(:,index_plasma,1) &
              + dnormaldzeta_plasma(3,itheta_plasma,izeta_plasma) * dnormzdomega(:,index_plasma,1) &
              - dnorm_normaldzeta_plasma(itheta_plasma,izeta_plasma) * dnorm_normaldomega(:,itheta_plasma,izeta_plasma)  )/norm_normal_plasma(itheta_plasma,izeta_plasma)
          end if
        enddo
      enddo
      do iomega = 1, nomega_plasma
        darea_plasmadomega(iomega) = dtheta_plasma*dzeta_plasma*nfp*sum(dnorm_normaldomega(iomega,:,:))
      end do
      call system_clock(toc)
      if (verbose) then
        print *,"Loop over plasma geom in regcoil_init_sensitivity:",real(toc-tic)/countrate," sec."
      end if

  else !----- Coil derivatives -----

      do izeta_coil = 1,nzeta_coil
        do itheta_coil = 1,ntheta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
          do l_coil = 0, (nfp-1)
            izetal_coil = izeta_coil + l_coil*nzeta_coil

            angle2 = zetal_coil(izetal_coil)
            sinangle2 = sin(angle2)
            cosangle2 = cos(angle2)

            do iomega = 1,nomega_coil
              ! Same convention for angle as used by nescin input file
              angle = xm_sensitivity(iomega)*theta_coil(itheta_coil) + nfp*xn_sensitivity(iomega)*zetal_coil(izetal_coil)
              sinangle = sin(angle)
              cosangle = cos(angle)

              ! omega = rmnc
              if (omega_coil(iomega) == 1) then
                domegadxdtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle*cosangle2
                domegadxdzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle*cosangle2 - cosangle*sinangle2
                domegadydtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle*sinangle2
                domegadydzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle*sinangle2 + cosangle*cosangle2
                domegadzdtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadzdzeta(iomega,itheta_coil,izetal_coil) = 0
                drdomega(1,index_coil,l_coil+1,iomega) = cosangle*cosangle2
                drdomega(2,index_coil,l_coil+1,iomega) = cosangle*sinangle2
                drdomega(3,index_coil,l_coil+1,iomega) = 0
              endif
              ! omega = zmns
              if (omega_coil(iomega) == 2) then
                domegadxdtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadxdzeta(iomega,itheta_coil,izetal_coil) = 0
                domegadydtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadydzeta(iomega,itheta_coil,izetal_coil) = 0
                domegadzdtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle
                domegadzdzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle
                drdomega(1,index_coil,l_coil+1,iomega) = 0
                drdomega(2,index_coil,l_coil+1,iomega) = 0
                drdomega(3,index_coil,l_coil+1,iomega) = sinangle
              endif
              ! omega = rmns
              if (omega_coil(iomega) == 3) then
                domegadxdtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle*cosangle2
                domegadxdzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle*cosangle2 - sinangle*sinangle2
                domegadydtheta(iomega,itheta_coil,izetal_coil) = xm_sensitivity(iomega)*cosangle*sinangle2
                domegadydzeta(iomega,itheta_coil,izetal_coil) = nfp*xn_sensitivity(iomega)*cosangle*sinangle2 + sinangle*cosangle2
                domegadzdtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadzdzeta(iomega,itheta_coil,izetal_coil) = 0
                drdomega(1,index_coil,l_coil+1,iomega) = sinangle*cosangle2
                drdomega(2,index_coil,l_coil+1,iomega) = sinangle*sinangle2
                drdomega(3,index_coil,l_coil+1,iomega) = 0
              endif
              ! omega = zmnc
              if (omega_coil(iomega) == 4) then
                domegadxdtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadxdzeta(iomega,itheta_coil,izetal_coil) = 0
                domegadydtheta(iomega,itheta_coil,izetal_coil) = 0
                domegadydzeta(iomega,itheta_coil,izetal_coil) = 0
                domegadzdtheta(iomega,itheta_coil,izetal_coil) = -xm_sensitivity(iomega)*sinangle
                domegadzdzeta(iomega,itheta_coil,izetal_coil) = -nfp*xn_sensitivity(iomega)*sinangle
                drdomega(1,index_coil,l_coil+1,iomega) = 0
                drdomega(2,index_coil,l_coil+1,iomega) = 0
                drdomega(3,index_coil,l_coil+1,iomega) = cosangle
              endif
              dxdtheta = drdtheta_coil(1,itheta_coil,izetal_coil)
              dydtheta = drdtheta_coil(2,itheta_coil,izetal_coil)
              dzdtheta = drdtheta_coil(3,itheta_coil,izetal_coil)
              dxdzeta = drdzeta_coil(1,itheta_coil,izetal_coil)
              dydzeta = drdzeta_coil(2,itheta_coil,izetal_coil)
              dzdzeta = drdzeta_coil(3,itheta_coil,izetal_coil)
              dnormxdomega(iomega,index_coil,l_coil+1) = &
                  domegadydzeta(iomega,itheta_coil,izetal_coil)*dzdtheta &
                + domegadzdtheta(iomega,itheta_coil,izetal_coil)*dydzeta &
                - domegadydtheta(iomega,itheta_coil,izetal_coil)*dzdzeta &
                - domegadzdzeta(iomega,itheta_coil,izetal_coil)*dydtheta
              dnormydomega(iomega,index_coil,l_coil+1) = &
                  domegadzdzeta(iomega,itheta_coil,izetal_coil)*dxdtheta &
                + domegadxdtheta(iomega,itheta_coil,izetal_coil)*dzdzeta &
                - domegadzdtheta(iomega,itheta_coil,izetal_coil)*dxdzeta &
                - domegadxdzeta(iomega,itheta_coil,izetal_coil)*dzdtheta
              dnormzdomega(iomega,index_coil,l_coil+1) = &
                  domegadxdzeta(iomega,itheta_coil,izetal_coil)*dydtheta &
                + domegadydtheta(iomega,itheta_coil,izetal_coil)*dxdzeta &
                - domegadxdtheta(iomega,itheta_coil,izetal_coil)*dydzeta &
                - domegadydzeta(iomega,itheta_coil,izetal_coil)*dxdtheta
            enddo
          enddo
        enddo
      enddo
      call system_clock(toc)
      if (verbose) then
        print *,"Loop over coil geom in regcoil_init_sensitivity:",real(toc-tic)/countrate," sec."
      end if

  end if

  if (sensitivity_option == 6) then !----- Plasma derivatives -----

      if (geometry_option_coil == 5) then
        do izeta_coil = 1,nzeta_coil
          do itheta_coil = 1,ntheta_coil
            index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
            do l_coil = 0, (nfp-1)
              izetal_coil = izeta_coil + l_coil*nzeta_coil
              indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil

              do iomega = 1,nomega_plasma
                drdomega_coil(1,index_coil,l_coil+1,iomega) = drdomega(1,index_coil,l_coil+1,iomega) + separation * ( &
                    dnormxdomega(iomega,index_coil,l_coil+1) - dnorm_normaldomega(iomega,itheta_coil,izeta_coil) &
                  * normal_plasma(1,itheta_coil,izetal_coil) / norm_normal_plasma(itheta_coil,izeta_coil) ) &
                  / norm_normal_plasma(itheta_coil,izeta_coil)
                drdomega_coil(2,index_coil,l_coil+1,iomega) = drdomega(2,index_coil,l_coil+1,iomega) + separation * ( &
                    dnormydomega(iomega,index_coil,l_coil+1) - dnorm_normaldomega(iomega,itheta_coil,izeta_coil) &
                  * normal_plasma(2,itheta_coil,izetal_coil) / norm_normal_plasma(itheta_coil,izeta_coil) ) &
                  / norm_normal_plasma(itheta_coil,izeta_coil)
                drdomega_coil(3,index_coil,l_coil+1,iomega) = drdomega(3,index_coil,l_coil+1,iomega) + separation * ( &
                    dnormzdomega(iomega,index_coil,l_coil+1) - dnorm_normaldomega(iomega,itheta_coil,izeta_coil) &
                  * normal_plasma(3,itheta_coil,izetal_coil) / norm_normal_plasma(itheta_coil,izeta_coil) ) &
                  / norm_normal_plasma(itheta_coil,izeta_coil)

                d2rdthetadomega_coil(1,index_coil,l_coil+1,iomega) = d2rdthetadomega(1,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normxdthetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldthetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(1,itheta_coil,izetal_coil) &
                  - (dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnormxdomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldtheta_plasma(1,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(1,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)
                d2rdthetadomega_coil(2,index_coil,l_coil+1,iomega) = d2rdthetadomega(2,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normydthetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldthetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(2,itheta_coil,izetal_coil) &
                  - (dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnormydomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldtheta_plasma(2,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(2,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)
                d2rdthetadomega_coil(3,index_coil,l_coil+1,iomega) = d2rdthetadomega(3,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normzdthetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldthetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(3,itheta_coil,izetal_coil) &
                  - (dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnormzdomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldtheta_plasma(3,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldtheta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(3,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)
                
                d2rdzetadomega_coil(1,index_coil,l_coil+1,iomega) = d2rdzetadomega(1,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normxdzetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldzetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(1,itheta_coil,izetal_coil) &
                  - (dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnormxdomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldzeta_plasma(1,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(1,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)
                d2rdzetadomega_coil(2,index_coil,l_coil+1,iomega) = d2rdzetadomega(2,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normydzetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldzetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(2,itheta_coil,izetal_coil) &
                  - (dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnormydomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldzeta_plasma(2,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(2,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)
                d2rdzetadomega_coil(3,index_coil,l_coil+1,iomega) = d2rdzetadomega(3,index_coil,l_coil+1,iomega) + separation * ( &
                    norm_normal_plasma(itheta_coil,izeta_coil) * d2normzdzetadomega(iomega,index_coil,l_coil+1) &
                  - d2norm_normaldzetadomega(iomega,itheta_coil,izeta_coil) * normal_plasma(3,itheta_coil,izetal_coil) &
                  - (dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnormzdomega(iomega,index_coil,l_coil+1) + dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * dnormaldzeta_plasma(3,itheta_coil,izetal_coil)) &
                  + 2*dnorm_normaldzeta_plasma(itheta_coil,izeta_coil) * dnorm_normaldomega(iomega,itheta_coil,izeta_coil) * normal_plasma(3,itheta_coil,izetal_coil)/norm_normal_plasma(itheta_coil,izeta_coil)) &
                  /(norm_normal_plasma(itheta_coil,izeta_coil)**2)

                dxdtheta = drdtheta_coil(1,itheta_coil,izetal_coil)
                dydtheta = drdtheta_coil(2,itheta_coil,izetal_coil)
                dzdtheta = drdtheta_coil(3,itheta_coil,izetal_coil)
                dxdzeta = drdzeta_coil(1,itheta_coil,izetal_coil)
                dydzeta = drdzeta_coil(2,itheta_coil,izetal_coil)
                dzdzeta = drdzeta_coil(3,itheta_coil,izetal_coil)

                dnormxdomega_coil(iomega,index_coil,l_coil+1) = &
                    d2rdzetadomega_coil(2,index_coil,l_coil+1,iomega) * dzdtheta + dydzeta * d2rdthetadomega_coil(3,index_coil,l_coil+1,iomega) &
                  - d2rdthetadomega_coil(2,index_coil,l_coil+1,iomega) * dzdzeta - dydtheta * d2rdzetadomega_coil(3,index_coil,l_coil+1,iomega)
                dnormydomega_coil(iomega,index_coil,l_coil+1) = &
                    d2rdzetadomega_coil(3,index_coil,l_coil+1,iomega) * dxdtheta + dzdzeta * d2rdthetadomega_coil(1,index_coil,l_coil+1,iomega) &
                  - d2rdthetadomega_coil(3,index_coil,l_coil+1,iomega) * dxdzeta - dzdtheta * d2rdzetadomega_coil(1,index_coil,l_coil+1,iomega)
                dnormzdomega_coil(iomega,index_coil,l_coil+1) = &
                    d2rdzetadomega_coil(1,index_coil,l_coil+1,iomega) * dydtheta + dxdzeta * d2rdthetadomega_coil(2,index_coil,l_coil+1,iomega) &
                  - d2rdthetadomega_coil(1,index_coil,l_coil+1,iomega) * dydzeta - dxdtheta * d2rdzetadomega_coil(2,index_coil,l_coil+1,iomega)
              end do
              dddomega(1, :, indexl_coil) = (net_poloidal_current_Amperes*d2rdthetadomega_coil(1,index_coil,l_coil+1,:) &
                - net_toroidal_current_Amperes*d2rdzetadomega_coil(1,index_coil,l_coil+1,:))/twopi
              dddomega(2, :, indexl_coil) = (net_poloidal_current_Amperes*d2rdthetadomega_coil(2,index_coil,l_coil+1,:) &
                - net_toroidal_current_Amperes*d2rdzetadomega_coil(2,index_coil,l_coil+1,:))/twopi
              dddomega(3, :, indexl_coil) = (net_poloidal_current_Amperes*d2rdthetadomega_coil(3,index_coil,l_coil+1,:) &
                - net_toroidal_current_Amperes*d2rdzetadomega_coil(3,index_coil,l_coil+1,:))/twopi
            end do
            dnorm_normaldomega_coil(:,itheta_coil,izeta_coil) = &
               (normal_coil(1, itheta_coil, izeta_coil)*dnormxdomega_coil(:,index_coil,1) &
              + normal_coil(2, itheta_coil, izeta_coil)*dnormydomega_coil(:,index_coil,1) &
              + normal_coil(3, itheta_coil, izeta_coil)*dnormzdomega_coil(:,index_coil,1)) &
                /norm_normal_coil(itheta_coil,izeta_coil)
          end do
        end do

        do iomega = 1, nomega_plasma
          darea_coildomega(iomega) = dtheta_coil*dzeta_coil*nfp*sum(dnorm_normaldomega_coil(iomega,:,:))
        end do

      end if

      if (verbose) then
        print *,"Coil surface and d sensitivity loop in regcoil_init_sensitivity: ",real(toc-tic)/countrate," sec."
        print *,"Init sensitivity complete."
      end if

  else !----- Coil derivatives -----

      call system_clock(tic,countrate)

      allocate(dR_squared_domega_half_grid(nomega_coil))
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dvolume_coildomega(nomega_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      dvolume_coildomega = 0

      do itheta_coil=1,ntheta_coil-1
        do izeta_coil=1,nzeta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
          do l_coil = 0, (nfp-1)
            izetal_coil = izeta_coil + l_coil*nzeta_coil

            major_R_squared_half_grid = &
               (r_coil(1,itheta_coil,izetal_coil)*r_coil(1,itheta_coil,izetal_coil) &
              + r_coil(2,itheta_coil,izetal_coil)*r_coil(2,itheta_coil,izetal_coil) &
              + r_coil(1,itheta_coil+1,izetal_coil)*r_coil(1,itheta_coil+1,izetal_coil) &
              + r_coil(2,itheta_coil+1,izetal_coil)*r_coil(2,itheta_coil+1,izetal_coil)) * (0.5d+0)

            dR_squared_domega_half_grid = &
              drdomega(1,index_coil,l_coil+1,:)*r_coil(1,itheta_coil,izetal_coil) &
              + drdomega(2,index_coil,l_coil+1,:)*r_coil(2,itheta_coil,izetal_coil) &
              + drdomega(1,index_coil+1,l_coil+1,:)*r_coil(1,itheta_coil+1,izetal_coil) &
              + drdomega(2,index_coil+1,l_coil+1,:)*r_coil(2,itheta_coil+1,izetal_coil)

            dvolume_coildomega = dvolume_coildomega &
              + dR_squared_domega_half_grid &
              * (r_coil(3,itheta_coil+1,izetal_coil) - r_coil(3,itheta_coil,izetal_coil)) &
              + major_R_squared_half_grid &
              * (drdomega(3,index_coil+1,l_coil+1,:) - drdomega(3,index_coil,l_coil+1,:))
          end do
        end do
      end do
      do izeta_coil=1,nzeta_coil
        index_coil_first = (izeta_coil-1)*ntheta_coil + 1
        index_coil_last = (izeta_coil-1)*ntheta_coil + ntheta_coil
        do l_coil = 0, (nfp-1)
          izetal_coil = izeta_coil + l_coil*nzeta_coil

          major_R_squared_half_grid = &
             (r_coil(1,ntheta_coil,izetal_coil)*r_coil(1,ntheta_coil,izetal_coil) &
            + r_coil(2,ntheta_coil,izetal_coil)*r_coil(2,ntheta_coil,izetal_coil) &
            + r_coil(1,1,izetal_coil)*r_coil(1,1,izetal_coil) &
            + r_coil(2,1,izetal_coil)*r_coil(2,1,izetal_coil)) * (0.5d+0)

          dR_squared_domega_half_grid = &
            drdomega(1,index_coil_first,l_coil+1,:)*r_coil(1,1,izetal_coil) &
            + drdomega(2,index_coil_first,l_coil+1,:)*r_coil(2,1,izetal_coil) &
            + drdomega(1,index_coil_last,l_coil+1,:)*r_coil(1,ntheta_coil,izetal_coil) &
            + drdomega(2,index_coil_last,l_coil+1,:)*r_coil(2,ntheta_coil,izetal_coil)

          dvolume_coildomega = dvolume_coildomega &
            + dR_squared_domega_half_grid &
            * (r_coil(3,1,izetal_coil) - r_coil(3,ntheta_coil,izetal_coil)) &
            + major_R_squared_half_grid &
            * (drdomega(3,index_coil_first,l_coil+1,:) - drdomega(3,index_coil_last,l_coil+1,:))
        end do
      end do
      dvolume_coildomega = dvolume_coildomega*(0.5d+0)*dzeta_coil
      call system_clock(toc)
      if (verbose) then
        print *,"Coil volume computation in regcoil_init_sensitivity :",real(toc-tic)/countrate," sec."
      end if
      deallocate(dR_squared_domega_half_grid)

      call system_clock(tic)

      do izeta_coil = 1, nzeta_coil
        do itheta_coil = 1, ntheta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
            do l_coil = 0, (nfp-1)
              izetal_coil = izeta_coil + l_coil*nzeta_coil
              indexl_coil = (izetal_coil-1)*ntheta_coil + itheta_coil
              dnorm_normaldomega(:,itheta_coil,izeta_coil) = &
                 (normal_coil(1, itheta_coil, izeta_coil)*dnormxdomega(:,index_coil,1) &
                + normal_coil(2, itheta_coil, izeta_coil)*dnormydomega(:,index_coil,1) &
                + normal_coil(3, itheta_coil, izeta_coil)*dnormzdomega(:,index_coil,1)) &
                  /norm_normal_coil(itheta_coil,izeta_coil)
              dddomega(1, :, indexl_coil) = (net_poloidal_current_Amperes*domegadxdtheta(:,itheta_coil,izetal_coil) &
                - net_toroidal_current_Amperes*domegadxdzeta(:,itheta_coil,izetal_coil))/twopi
              dddomega(2, :, indexl_coil) = (net_poloidal_current_Amperes*domegadydtheta(:,itheta_coil,izetal_coil) &
                - net_toroidal_current_Amperes*domegadydzeta(:,itheta_coil,izetal_coil))/twopi
              dddomega(3, :, indexl_coil) = (net_poloidal_current_Amperes*domegadzdtheta(:,itheta_coil,izetal_coil) &
                - net_toroidal_current_Amperes*domegadzdzeta(:,itheta_coil,izetal_coil))/twopi
          enddo
        enddo
      enddo
      call system_clock(toc)
      if (verbose) then
        print *,"d sensitivity loop in regcoil_init_sensitivity: ",real(toc-tic)/countrate," sec."
      end if

      do iomega = 1, nomega_coil
        darea_coildomega(iomega) = dtheta_coil*dzeta_coil*nfp*sum(dnorm_normaldomega(iomega,:,:))
      end do

      call system_clock(tic,countrate)
      if (verbose) then
        print *,"Beginning calculations of d_min in regcoil_init_sensitivity."
      end if

      allocate(dist(ntheta_coil,nzeta_coil,ntheta_plasma,nzeta_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(normal_coil_fourd(ntheta_coil,nzeta_coil,ntheta_plasma,nzeta_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(normal_plasma_fourd(ntheta_coil,nzeta_coil,ntheta_plasma,nzeta_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(ones(ntheta_coil,nzeta_coil,ntheta_plasma,nzeta_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dcoil_plasma_dist_mindomega(nomega_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dsum_exp_mindomega(nomega_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      ones = 1

      do itheta_coil = 1,ntheta_coil
        do izeta_coil = 1,nzeta_coil
          do itheta_plasma = 1,ntheta_plasma
            do izeta_plasma = 1,nzeta_plasma
              dist(itheta_coil,izeta_coil,itheta_plasma,izeta_plasma) = &
             sqrt((r_coil(1,itheta_coil,izeta_coil)-r_plasma(1,itheta_plasma,izeta_plasma))**2 &
                + (r_coil(2,itheta_coil,izeta_coil)-r_plasma(2,itheta_plasma,izeta_plasma))**2 &
                + (r_coil(3,itheta_coil,izeta_coil)-r_plasma(3,itheta_plasma,izeta_plasma))**2)
            end do
          end do
          normal_coil_fourd(itheta_coil,izeta_coil,:,:) = norm_normal_coil(itheta_coil,izeta_coil)
        end do
      end do

      do itheta_plasma=1,ntheta_plasma
        do izeta_plasma=1,nzeta_plasma
          normal_plasma_fourd(:,:,itheta_plasma,izeta_plasma) = norm_normal_plasma(itheta_plasma,izeta_plasma)
        end do
      end do

      coil_plasma_dist_min = minval(dist)

      sum_exp_min = sum(normal_plasma_fourd*normal_coil_fourd*exp(-coil_plasma_dist_lse_p* &
        (dist-coil_plasma_dist_min*ones)))
      coil_plasma_dist_min_lse = -log(sum_exp_min/(sum(normal_plasma_fourd*normal_coil_fourd)))/coil_plasma_dist_lse_p &
        + coil_plasma_dist_min

      if (verbose) then
        print *,"coil_plasma_dist_min: ", coil_plasma_dist_min
        print *,"coil_plasma_dist_min_lse: ", coil_plasma_dist_min_lse
      end if

      dcoil_plasma_dist_mindomega = 0
      dsum_exp_mindomega = 0
      do itheta_coil=1,ntheta_coil
        do izeta_coil=1,nzeta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
          do itheta_plasma=1,ntheta_plasma
            do izeta_plasma=1,nzeta_plasma
              dsum_exp_mindomega = dsum_exp_mindomega + &
                norm_normal_plasma(itheta_plasma,izeta_plasma)*norm_normal_coil(itheta_coil,izeta_coil)* &
                (dnorm_normaldomega(:,itheta_coil,izeta_coil)/norm_normal_coil(itheta_coil,izeta_coil) &
              - coil_plasma_dist_lse_p*(drdomega(1,index_coil,1,:)*(r_coil(1,itheta_coil,izeta_coil)-r_plasma(1,itheta_plasma,izeta_plasma)) &
              + drdomega(2,index_coil,1,:)*(r_coil(2,itheta_coil,izeta_coil)-r_plasma(2,itheta_plasma,izeta_plasma)) &
              + drdomega(3,index_coil,1,:)*(r_coil(3,itheta_coil,izeta_coil)-r_plasma(3,itheta_plasma,izeta_plasma)))/dist(itheta_coil,izeta_coil,itheta_plasma,izeta_plasma))* &
                exp(-coil_plasma_dist_lse_p*(dist(itheta_coil,izeta_coil,itheta_plasma,izeta_plasma)-coil_plasma_dist_min))
            end do
          end do
        end do
      end do

      dcoil_plasma_dist_mindomega = -dsum_exp_mindomega/(sum_exp_min*coil_plasma_dist_lse_p) &
        + darea_coildomega/(sum(norm_normal_plasma)*area_coil*coil_plasma_dist_lse_p)

      call system_clock(toc)
      if (verbose) then
        print *,"Coil-plasma distance calculation in regcoil_init_sensitivity:",real(toc-tic)/countrate," sec."
        print *,"Init sensitivity complete."
      end if

  end if


end subroutine regcoil_init_sensitivity





