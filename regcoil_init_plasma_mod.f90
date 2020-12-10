module regcoil_init_plasma_mod


  use stel_kinds
  use omp_lib

  implicit none

  private

  public :: regcoil_init_plasma

  real(dp) :: theta_rootSolve_target, zeta

contains

  subroutine regcoil_init_plasma()

    use regcoil_variables
    use regcoil_read_efit_mod
    use read_wout_mod, only: nfp_vmec => nfp, xm_vmec => xm, xn_vmec => xn, &
         rmnc_vmec => rmnc, zmns_vmec => zmns, rmns_vmec => rmns, zmnc_vmec => zmnc, &
         bmnc_vmec => bmnc, bmns_vmec => bmns, lasym_vmec => lasym, mnmax_vmec => mnmax, ns, Rmajor, read_wout_file, &
         mpol_vmec => mpol, ntor_vmec => ntor, bvco, bsubvmnc
    use safe_open_mod
    use stel_constants
!    use regcoil_double_to_single
!    use regcoil_read_single_Fourier

    implicit none

    integer :: i, itheta, izeta, imn, tic, toc, countrate, iflag, ierr, iopen, tic1, toc1, iunit
    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta, d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2, d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2
    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta, d2sinangle2dzeta2, d2cosangle2dzeta2
    real(dp) :: angle3, sinangle3, cosangle3, dsinangle3dtheta, dsinangle3dzeta, dcosangle3dtheta, dcosangle3dzeta, d2sinangle3dtheta2, d2sinangle3dthetadzeta, d2sinangle3dzeta2, d2cosangle3dtheta2, d2cosangle3dthetadzeta, d2cosangle3dzeta2
    real(dp) :: weight1, weight2, weight3, weight4, theta, r_temp, z_temp, dnorm
    integer :: ntheta_coordTransform, nzeta_coordTransform
    real(dp), dimension(:,:), allocatable :: r_coordTransform, z_coordTransform, major_R_squared
    real(dp), dimension(:), allocatable :: rmnc_vmecLast, zmns_vmecLast
    real(dp) :: rootSolve_abserr, rootSolve_relerr, theta_rootSolve_min, theta_rootSolve_max, theta_rootSolve_soln
    integer :: fzeroFlag, mpol, ntor, jm, jn, index
    real(dp), dimension(:,:), allocatable :: E_big, F_big, G_big, e_small, f_small, g_small, K_curvature, kappa1, kappa2
    real(dp), dimension(:,:,:), allocatable :: f_test
    
    call system_clock(tic, countrate)
    if (verbose) print *,"Initializing plasma surface."
    
    select case (geometry_option_plasma)
    case (0,1)
       ! Plain circular torus
       if (verbose) print *,"  Building a plain circular torus."
       
       nfp = nfp_imposed
       mnmax_plasma = 2
       lasym = .false.
       
       call regcoil_allocate_plasma_surface_arrays()

       xm_plasma = (/0,1/)
       xn_plasma = (/0,0/)
       rmnc_plasma = (/ R0_plasma, a_plasma /)
       zmns_plasma = (/ 0.0d+0, a_plasma /)
       
    case(6)
       ! Read in an ASCII table
       call safe_open(iunit, ierr, trim(shape_filename_plasma), 'old', 'formatted')
       if (ierr .ne. 0) then
          stop 'Error opening nescin file'
       endif
       
       ! Skip first line
       read (iunit, *)
       read (iunit, *) mnmax_plasma
       
       call regcoil_allocate_plasma_surface_arrays()

       ! Skip a line
       read (iunit, *)
       do i = 1, mnmax_plasma
          read (iunit, *) xm_plasma(i), xn_plasma(i), rmnc_plasma(i), zmns_plasma(i), rmns_plasma(i), zmnc_plasma(i)
       end do
       
       close(iunit)
       
       nfp = nfp_imposed
       lasym = .true.
       
    case(7)
       ! Read in FOCUS format plasma boundary, for more information, please check
       ! https://princetonuniversity.github.io/FOCUS/rdsurf.pdf
       call safe_open(iunit, ierr, trim(shape_filename_plasma), 'old', 'formatted')
       if (ierr .ne. 0) then
          stop 'Error opening FOCUS file'
       endif
       if (verbose) print *,"Reading FOCUS format data from file:", trim(shape_filename_plasma)
       
       ! Skip first line
       read (iunit, *)
       read (iunit, *) mnmax_plasma, nfp, nbf
       
       if (allocated(xm_plasma)) deallocate(xm_plasma)
       allocate(xm_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error! 1"

       if (allocated(xn_plasma)) deallocate(xn_plasma)
       allocate(xn_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  2"

       if (allocated(rmnc_plasma)) deallocate(rmnc_plasma)
       allocate(rmnc_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  3"

       if (allocated(zmns_plasma)) deallocate(zmns_plasma)
       allocate(zmns_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  4"

       if (allocated(rmns_plasma)) deallocate(rmns_plasma)
       allocate(rmns_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  5"

       if (allocated(zmnc_plasma)) deallocate(zmnc_plasma)
       allocate(zmnc_plasma(mnmax_plasma),stat=iflag)
       if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  6"
       
       ! Skip two lines
       read (iunit, *)
       read (iunit, *)
       do i = 1, mnmax_plasma
          read (iunit, *) xn_plasma(i), xm_plasma(i), rmnc_plasma(i), rmns_plasma(i), zmnc_plasma(i), zmns_plasma(i)
       end do      

       xn_plasma = xn_plasma * nfp ! include nfp

       ! read bnorm coefficients if available
       if ( nbf > 0 ) then
          if (allocated(bfm)) deallocate(bfs)
          allocate(bfm(nbf),stat=iflag)
          if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  7"          

          if (allocated(bfn)) deallocate(bfc)
          allocate(bfn(nbf),stat=iflag)
          if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  8"  

          if (allocated(bfs)) deallocate(bfs)
          allocate(bfs(nbf),stat=iflag)
          if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  9"          

          if (allocated(bfc)) deallocate(bfc)
          allocate(bfc(nbf),stat=iflag)
          if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error! 10"  

          ! Skip two lines
          read (iunit, *)
          read (iunit, *) !empty line
          do i = 1, nbf
             read (iunit, *) bfn(i), bfm(i), bfc(i), bfs(i)
          enddo

          bfn = bfn * nfp ! include nfp

          if (verbose) print *,"Number of modes for Bnormal read from FOCUS file:", nbf
       endif

       close(iunit)
       
       lasym = .true.
       
    case (2,3,8)
       ! VMEC, "original" theta coordinate which is not a straight-field-line coordinate
       ! 8 switches to a Single Fourier series representation on a theta coordinate with respect to axis
       call read_wout_file(wout_filename, ierr, iopen)
       if (iopen .ne. 0) stop 'error opening wout file'
       if (ierr .ne. 0) stop 'error reading wout file'
       if (verbose) print *,"  Successfully read VMEC data from ",trim(wout_filename)
       
       if (geometry_option_plasma == 2 .or. geometry_option_plasma == 8) then
          ! Only use the outermost point in the full radial mesh:
          weight1 = 0
          weight2 = 1
          weight3 = 1.5
          weight4 = -0.5
          if (verbose) print *,"  Using outermost grid point in VMEC's FULL radial grid."
       else
          ! Average the two outermost points in the full radial mesh 
          ! to get a value on the outermost point of the half radial mesh:
          weight1 = 0.5_dp
          weight2 = 0.5_dp
          weight3 = 1.0
          weight4 = 0
          if (verbose) print *,"  Using outermost grid point in VMEC's HALF radial grid."
       end if
       
       nfp = nfp_vmec
       mnmax_plasma = mnmax_vmec
       lasym = lasym_vmec
       R0_plasma = Rmajor
      
       call regcoil_allocate_plasma_surface_arrays()
       
       xm_plasma = xm_vmec
       xn_plasma = xn_vmec
       if (verbose) print *,"size of rmnc_vmec:",size(rmnc_vmec,1),size(rmnc_vmec,2)
       rmnc_plasma = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
       zmns_plasma = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2
       if (sensitivity_option == 6) then
         bmnc_plasma = bmnc_vmec(:,ns-1) * weight4 + bmnc_vmec(:,ns) * weight3
       end if
       if (lasym) then
          rmns_plasma = rmns_vmec(:,ns-1) * weight1 + rmns_vmec(:,ns) * weight2
          zmnc_plasma = zmnc_vmec(:,ns-1) * weight1 + zmnc_vmec(:,ns) * weight2
          if (sensitivity_option == 6) then
            bmns_plasma = bmns_vmec(:,ns-1) * weight4 + bmns_vmec(:,ns) * weight3
          end if
       end if

       if (geometry_option_plasma .eq. 8) then
          call regcoil_double_to_single()
       end if

    case (9)
       ! Single Fourier file, theta coordinate with respect to axis
       call regcoil_read_single_Fourier()

    case (4)
       ! VMEC, straight-field-line poloidal coordinate
       call read_wout_file(wout_filename, ierr, iopen)
       if (iopen .ne. 0) stop 'error opening wout file'
       if (ierr .ne. 0) stop 'error reading wout file'
       if (verbose) print *,"  Successfully read VMEC data from ",trim(wout_filename)
       
       
       nfp = nfp_vmec
       lasym = lasym_vmec
       if (lasym) then
          stop "Error! geometry_option_plasma=4 is not yet implemented for lasym=true"
       end if
       R0_plasma = Rmajor
       
       ! Average R and Z from the outermost 2 grid points in vmec's full mesh
       ! to get R and Z on the outermost point of vmec's half mesh:
       weight1 = 0.5_dp
       weight2 = 0.5_dp
       weight3 = 1.0
       weight4 = 0

       if (allocated(rmnc_vmecLast)) deallocate(rmnc_vmecLast)
       allocate(rmnc_vmecLast(mnmax_vmec),stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

       if (allocated(zmns_vmecLast)) deallocate(zmns_vmecLast)
       allocate(zmns_vmecLast(mnmax_vmec),stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

       rmnc_vmecLast = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
       zmns_vmecLast = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2
       if (sensitivity_option==6) then
         bmnc_plasma = bmnc_vmec(:,ns-1) * weight4 + bmnc_vmec(:,ns) * weight3
         if (lasym) then
           bmns_plasma = bmns_vmec(:,ns-1) * weight4 + bmnc_vmec(:,ns) * weight3
         end if
       end if

       ! Since the "original" vmec poloidal angle is chosen to have a very condensed
       ! Fourier spectrum, we probably need more Fourier modes to represent the surface using the
       ! straight-field-line coordinate.
       mpol = mpol_vmec*mpol_transform_refinement
       ntor = ntor_vmec*ntor_transform_refinement
       
       ! Beginning of coordinate transformation.
       ! Set up high-resolution grid in the "new" theta coordinate:
       ntheta_coordTransform = mpol * 2 
       nzeta_coordTransform = ntor * 2

       if (allocated(r_coordTransform)) deallocate(r_coordTransform)
       allocate(r_coordTransform(ntheta_coordTransform, nzeta_coordTransform), stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

       if (allocated(z_coordTransform)) deallocate(z_coordTransform)
       allocate(z_coordTransform(ntheta_coordTransform, nzeta_coordTransform), stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
       r_coordTransform = 0
       z_coordTransform = 0
       
       call system_clock(tic1)
       rootSolve_abserr = 1.0e-10_dp
       rootSolve_relerr = 1.0e-10_dp
       !open(unit=5,file="testStraightFieldLines",status='new',form='formatted')
       !write (5,*) ntheta_coordTransform, nzeta_coordTransform
       do izeta = 1,nzeta_coordTransform
          zeta = (izeta-1.0_dp)/nzeta_coordTransform
          do itheta = 1,ntheta_coordTransform
             ! For each value of the new coordinates, solve for the old theta:
             theta_rootSolve_target = (itheta-1.0_dp)/ntheta_coordTransform
             theta_rootSolve_min = theta_rootSolve_target - 0.3
             theta_rootSolve_max = theta_rootSolve_target + 0.3
             
             call regcoil_fzero(fzero_residual, theta_rootSolve_min, theta_rootSolve_max, theta_rootSolve_target, &
                  rootSolve_relerr, rootSolve_abserr, fzeroFlag)
             ! Note: fzero returns its answer in theta_rootSolve_min
             theta_rootSolve_soln = theta_rootSolve_min
             if (fzeroFlag == 4) then
                stop "ERROR: fzero returned error 4: no sign change in residual"
             else if (fzeroFlag > 2) then
                print *,"WARNING in irp: fzero returned an error code:",fzeroFlag
             end if
             ! Now that we have the old theta, evaluate r and z:
             r_temp = 0
             z_temp = 0
             do imn = 1, mnmax_vmec
                !angle = twopi*(xm_vmec(imn)*theta_rootSolve_soln - xn_vmec(imn)*zeta/nfp)
                angle = xm_vmec(imn)*theta_rootSolve_soln - xn_vmec(imn)*zeta
                r_temp = r_temp + rmnc_vmecLast(imn)*cos(angle)
                z_temp = z_temp + zmns_vmecLast(imn)*sin(angle)
             end do
             r_coordTransform(itheta,izeta) = r_temp
             z_coordTransform(itheta,izeta) = z_temp
             !write(5,*) theta_rootSolve_soln
          end do
       end do
       !close(unit=5)
       call system_clock(toc1)
       if (verbose) print *,"  Time for root solving:",real(toc1-tic1)/countrate
       
       ! Now that we have R and Z on a grid in the new coordinates, Fourier transform the results.
       
       ! The next bit of code is much like initFourierModesMod, but with 2 differences: 
       ! 1. We need to keep the m=n=0 mode.
       ! 2. We follow VMEC convention that n includes the factor of nfp.
       
       ! xm is nonnegative.
       ! xn can be negative, zero, or positive.
       ! When xm is 0, xn must be positive.
       mnmax_plasma = mpol*(ntor*2+1) + ntor + 1
       call regcoil_allocate_plasma_surface_arrays()
       
       ! Handle the xm=0 modes:
       xm_plasma=0
       do jn=0,ntor
          xn_plasma(jn+1)=jn*nfp
       end do
       
       ! Handle the xm>0 modes:
       index = ntor + 1
       do jm = 1,mpol
          do jn = -ntor, ntor
             index = index + 1
             xn_plasma(index) = jn*nfp
             xm_plasma(index) = jm
          end do
       end do
       ! Initialization of xm and xn is now complete.
       
       call system_clock(tic1)
       do imn = 1, mnmax_plasma
          dnorm = (1.0_dp)/(ntheta_coordTransform*nzeta_coordTransform)
          if (xm_plasma(imn).ne.0 .or. xn_plasma(imn).ne.0) dnorm = 2*dnorm
          r_temp = 0
          z_temp = 0
          do izeta = 1, nzeta_coordTransform
             zeta = (izeta-1.0_dp)/nzeta_coordTransform
             do itheta = 1, ntheta_coordTransform
                theta = (itheta-1.0_dp)/ntheta_coordTransform
                angle = xm_plasma(imn)*theta-xn_plasma(imn)*zeta
                cosangle = cos(angle)
                sinangle = sin(angle)
                r_temp = r_temp + r_coordTransform(itheta,izeta) * cosangle
                z_temp = z_temp + z_coordTransform(itheta,izeta) * sinangle
             end do
          end do
          rmnc_plasma(imn) = r_temp*dnorm
          zmns_plasma(imn) = z_temp*dnorm
       end do
       call system_clock(toc1)
       if (verbose) print *,"  Time for Fourier transform:",real(toc1-tic1)/countrate
       
    case (5)
       ! EFIT
       
       lasym = .true.
       nfp = nfp_imposed
       mnmax_plasma = efit_num_modes

       call regcoil_allocate_plasma_surface_arrays()

       call regcoil_read_efit(efit_filename, efit_psiN, efit_num_modes, rmnc_plasma, zmns_plasma, rmns_plasma, zmnc_plasma)
       
       ! Set major radius equal to the zero-frequency component of R(theta)
       R0_plasma = rmnc_plasma(1)
       
       xn_plasma = 0
       do i=1,efit_num_modes
          xm_plasma(i) = i-1
       end do
       
    case default
       print *,"Error! Invalid setting for geometry_option_plasma:",geometry_option_plasma
       stop
    end select

    ! ---------------------------
    ! End of the parts of code specific to each geometry_option_plasma.
    ! Now comes the code that applies to all values of geometry_option_plasma.

    nzetal_plasma = nzeta_plasma * nfp
    
    if (allocated(theta_plasma)) deallocate(theta_plasma)
    allocate(theta_plasma(ntheta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (allocated(zeta_plasma)) deallocate(zeta_plasma)
    allocate(zeta_plasma(nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (allocated(zetal_plasma)) deallocate(zetal_plasma)
    allocate(zetal_plasma(nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
    
    do i=1,ntheta_plasma
       theta_plasma(i) = twopi*(i-1.0_dp)/ntheta_plasma
    end do
    
    do i=1,nzeta_plasma
       zeta_plasma(i) = twopi/nfp*(i-1.0_dp)/nzeta_plasma
    end do
    
    do i=1,nzetal_plasma
       zetal_plasma(i) = twopi*(i-1.0_dp)/nzetal_plasma
    end do
    
    ! First coordinate is the Cartesian component x, y, or z

    if (allocated(r_plasma)) deallocate(r_plasma)
    allocate(r_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (allocated(drdtheta_plasma)) deallocate(drdtheta_plasma)
    allocate(drdtheta_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (allocated(drdzeta_plasma)) deallocate(drdzeta_plasma)
    allocate(drdzeta_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (allocated(normal_plasma)) deallocate(normal_plasma)
    allocate(normal_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

    if (use_arclength_angle) then
      if (allocated(omega_arclength)) deallocate(omega_arclength)
      allocate(omega_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(domegadtheta_arclength)) deallocate(domegadtheta_arclength)
      allocate(domegadtheta_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(domegadzeta_arclength)) deallocate(domegadzeta_arclength)
      allocate(domegadzeta_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (geometry_option_coil == 5) then
        if (allocated(d2omegadtheta2_arclength)) deallocate(d2omegadtheta2_arclength)
        allocate(d2omegadtheta2_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

        if (allocated(d2omegadthetadzeta_arclength)) deallocate(d2omegadthetadzeta_arclength)
        allocate(d2omegadthetadzeta_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

        if (allocated(d2omegadzeta2_arclength)) deallocate(d2omegadzeta2_arclength)
        allocate(d2omegadzeta2_arclength(ntheta_plasma,nzetal_plasma),stat=iflag)
        if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      end if
    end if

    if (geometry_option_coil == 5) then
      if (allocated(d2rdtheta2_plasma)) deallocate(d2rdtheta2_plasma)
      allocate(d2rdtheta2_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(d2rdthetadzeta_plasma)) deallocate(d2rdthetadzeta_plasma)
      allocate(d2rdthetadzeta_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(d2rdzeta2_plasma)) deallocate(d2rdzeta2_plasma)
      allocate(d2rdzeta2_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(dnormaldtheta_plasma)) deallocate(dnormaldtheta_plasma)
      allocate(dnormaldtheta_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      if (allocated(dnormaldzeta_plasma)) deallocate(dnormaldzeta_plasma)
      allocate(dnormaldzeta_plasma(3,ntheta_plasma,nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
    end if

    r_plasma=0
    drdtheta_plasma=0
    drdzeta_plasma=0
    if (geometry_option_coil == 5) then
      d2rdtheta2_plasma=0
      d2rdthetadzeta_plasma=0
      d2rdzeta2_plasma=0
    end if
    if (use_arclength_angle) then
      omega_arclength=0
      domegadtheta_arclength=0
      domegadzeta_arclength=0
      if (geometry_option_coil == 5) then
        d2omegadtheta2_arclength=0
        d2omegadthetadzeta_arclength=0
        d2omegadzeta2_arclength=0
      end if
    end if

    if (geometry_option_plasma==8 .or. geometry_option_plasma==9) then
        !$OMP PARALLEL
        !$OMP MASTER
        if (verbose) then
          print *,"  Number of OpenMP threads:",omp_get_num_threads()
        end if
        !$OMP END MASTER
        !$OMP DO PRIVATE(angle,angle2,sinangle,cosangle,sinangle2,cosangle2,dsinangle2dzeta,dcosangle2dzeta,d2sinangle2dzeta2,d2cosangle2dzeta2,dsinangledzeta,dcosangledzeta,d2sinangledzeta2,d2cosangledzeta2)
        do izeta = 1, nzetal_plasma
            angle2 = zetal_plasma(izeta)
            sinangle2 = sin(angle2)
            cosangle2 = cos(angle2)
            dsinangle2dzeta = cosangle2
            dcosangle2dzeta = -sinangle2
            d2sinangle2dzeta2 = -sinangle2
            d2cosangle2dzeta2 = -cosangle2
            do itheta = 1, ntheta_plasma
                do i=1,nmax_axis
                    angle = -xn_axis(i)*zetal_plasma(izeta)
                    sinangle = sin(angle)
                    cosangle = cos(angle)
                    dsinangledzeta = -cosangle*xn_axis(i)
                    dcosangledzeta = sinangle*xn_axis(i)
                    d2sinangledzeta2 = -sinangle*xn_axis(i)**2
                    d2cosangledzeta2 = -cosangle*xn_axis(i)**2

                    r_plasma(1,itheta,izeta) = r_plasma(1,itheta,izeta) + raxis_cc(i) * cosangle * cosangle2
                    r_plasma(2,itheta,izeta) = r_plasma(2,itheta,izeta) + raxis_cc(i) * cosangle * sinangle2
                    r_plasma(3,itheta,izeta) = r_plasma(3,itheta,izeta) + zaxis_cs(i) * sinangle

                    drdzeta_plasma(1,itheta,izeta) = drdzeta_plasma(1,itheta,izeta) + raxis_cc(i) * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta)
                    drdzeta_plasma(2,itheta,izeta) = drdzeta_plasma(2,itheta,izeta) + raxis_cc(i) * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta)
                    drdzeta_plasma(3,itheta,izeta) = drdzeta_plasma(3,itheta,izeta) + zaxis_cs(i) * dsinangledzeta
                    
                    if (geometry_option_coil == 5) then
                      d2rdzeta2_plasma(1,itheta,izeta) = d2rdzeta2_plasma(1,itheta,izeta) &
                        + raxis_cc(i) * (d2cosangledzeta2 * cosangle2 + 2*dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2)
                      d2rdzeta2_plasma(2,itheta,izeta) = d2rdzeta2_plasma(2,itheta,izeta) &
                        + raxis_cc(i) * (d2cosangledzeta2 * sinangle2 + 2*dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2)
                      d2rdzeta2_plasma(3,itheta,izeta) = d2rdzeta2_plasma(3,itheta,izeta) + zaxis_cs(i) * d2sinangledzeta2
                    endif
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    allocate(f_test(mnmax_plasma,ntheta_plasma,nzetal_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
    f_test = 0

    if (use_arclength_angle) then
        !$OMP PARALLEL
        !$OMP MASTER
        if (verbose) then
          print *,"  Number of OpenMP threads:",omp_get_num_threads()
        end if
        !$OMP END MASTER
        !$OMP DO PRIVATE(angle,sinangle,cosangle,dsinangledtheta,dsinangledzeta,d2sinangledtheta2,d2sinangledthetadzeta,d2sinangledzeta2)
        do izeta = 1,nzetal_plasma
          do itheta = 1,ntheta_plasma
            do imn = 1,mnmax_plasma
              angle = xm_plasma(imn)*theta_plasma(itheta) - xn_plasma(imn)*zetal_plasma(izeta)
              sinangle = sin(angle)
              cosangle = cos(angle)
              dsinangledtheta = cosangle*xm_plasma(imn)
              dsinangledzeta = -cosangle*xn_plasma(imn)
              d2sinangledtheta2 = -sinangle*xm_plasma(imn)**2
              d2sinangledthetadzeta = sinangle*xm_plasma(imn)*xn_plasma(imn)
              d2sinangledzeta2 = -sinangle*xn_plasma(imn)**2

!              f_test(imn,itheta,izeta) = dsinangledtheta

              omega_arclength(itheta,izeta) = omega_arclength(itheta,izeta) + omns(imn) * sinangle

              domegadtheta_arclength(itheta,izeta) = domegadtheta_arclength(itheta,izeta) + omns(imn) * dsinangledtheta
              domegadzeta_arclength(itheta,izeta) = domegadzeta_arclength(itheta,izeta) + omns(imn) * dsinangledzeta

              if (geometry_option_coil == 5) then
                d2omegadtheta2_arclength(itheta,izeta) = d2omegadtheta2_arclength(itheta,izeta) + omns(imn) * d2sinangledtheta2
                d2omegadthetadzeta_arclength(itheta,izeta) = d2omegadthetadzeta_arclength(itheta,izeta) + omns(imn) * d2sinangledthetadzeta
                d2omegadzeta2_arclength(itheta,izeta) = d2omegadzeta2_arclength(itheta,izeta) + omns(imn) * d2sinangledzeta2
              end if
            end do
          end do
        end do

        !$OMP END DO
        !$OMP END PARALLEL
    end if

    !$OMP PARALLEL
    !$OMP MASTER
    if (verbose) then
      print *,"  Number of OpenMP threads:",omp_get_num_threads()
    end if
    !$OMP END MASTER
    !$OMP DO PRIVATE(angle,angle2,angle3,sinangle,cosangle,sinangle2,cosangle2,sinangle3,cosangle3,dsinangledtheta,dcosangledtheta,dsinangledzeta,dcosangledzeta,dsinangle2dzeta,dcosangle2dzeta,d2sinangle2dzeta2,d2cosangle2dzeta2,d2sinangledtheta2,d2cosangledtheta2,d2sinangledthetadzeta,d2cosangledthetadzeta,d2sinangledzeta2,d2cosangledzeta2,dsinangle3dtheta,dcosangle3dtheta,dsinangle3dzeta,dcosangle3dzeta,d2sinangle3dtheta2,d2cosangle3dtheta2,d2sinangle3dthetadzeta,d2cosangle3dthetadzeta,d2sinangle3dzeta2,d2cosangle3dzeta2)
    do izeta = 1,nzetal_plasma
       angle2 = zetal_plasma(izeta)
       sinangle2 = sin(angle2)
       cosangle2 = cos(angle2)
       dsinangle2dzeta = cosangle2
       dcosangle2dzeta = -sinangle2
       d2sinangle2dzeta2 = -sinangle2
       d2cosangle2dzeta2 = -cosangle2
       do itheta = 1,ntheta_plasma
          if (use_arclength_angle) then
            angle3 = theta_plasma(itheta) - omega_arclength(itheta,izeta)
            sinangle3 = sin(angle3)
            cosangle3 = cos(angle3)
            dsinangle3dtheta = cosangle3 * (1 - domegadtheta_arclength(itheta,izeta))!
            dcosangle3dtheta = -sinangle3 * (1 - domegadtheta_arclength(itheta,izeta))!
            dsinangle3dzeta = cosangle3 * (- domegadzeta_arclength(itheta,izeta))!
            dcosangle3dzeta = -sinangle3 * (- domegadzeta_arclength(itheta,izeta))!
            if (geometry_option_coil == 5) then
              d2sinangle3dtheta2 = dcosangle3dtheta * (1 - domegadtheta_arclength(itheta,izeta)) + cosangle3 * (- d2omegadtheta2_arclength(itheta,izeta))
              d2cosangle3dtheta2 = -dsinangle3dtheta * (1 - domegadtheta_arclength(itheta,izeta)) - sinangle3 * (- d2omegadtheta2_arclength(itheta,izeta))
              d2sinangle3dthetadzeta = dcosangle3dzeta * (1 - domegadtheta_arclength(itheta,izeta)) + cosangle3 * (- d2omegadthetadzeta_arclength(itheta,izeta))
              d2cosangle3dthetadzeta = -dsinangle3dzeta * (1 - domegadtheta_arclength(itheta,izeta)) - sinangle3 * (- d2omegadthetadzeta_arclength(itheta,izeta))
              d2sinangle3dzeta2 = dcosangle3dzeta * (- domegadzeta_arclength(itheta,izeta)) + cosangle3 * (- d2omegadzeta2_arclength(itheta,izeta))
              d2cosangle3dzeta2 = -dsinangle3dzeta * (- domegadzeta_arclength(itheta,izeta)) - sinangle3 * (- d2omegadzeta2_arclength(itheta,izeta))
            end if
          else
            angle3 = theta_plasma(itheta)
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
          do imn = 1,mnmax_plasma
             angle = xm_plasma(imn)*theta_plasma(itheta) - xn_plasma(imn)*zetal_plasma(izeta)
             sinangle = sin(angle)
             cosangle = cos(angle)
             dsinangledtheta = cosangle*xm_plasma(imn)
             dcosangledtheta = -sinangle*xm_plasma(imn)
             dsinangledzeta = -cosangle*xn_plasma(imn)
             dcosangledzeta = sinangle*xn_plasma(imn)
             d2sinangledtheta2 = -sinangle*xm_plasma(imn)**2
             d2cosangledtheta2 = -cosangle*xm_plasma(imn)**2!!
             d2sinangledthetadzeta = sinangle*xm_plasma(imn)*xn_plasma(imn)
             d2cosangledthetadzeta = cosangle*xm_plasma(imn)*xn_plasma(imn)
             d2sinangledzeta2 = -sinangle*xn_plasma(imn)**2
             d2cosangledzeta2 = -cosangle*xn_plasma(imn)**2

             select case (geometry_option_plasma)
             case (8,9)
!                  if ((xm_plasma(imn) > 10) .or. (abs(xn_plasma(imn)) > 10*nfp)) cycle
!                  f_test(imn,itheta,izeta) = (dcosangle3dtheta)
                 r_plasma(1,itheta,izeta) = r_plasma(1,itheta,izeta) + lmnc(imn) * cosangle * cosangle3 * cosangle2
                 r_plasma(2,itheta,izeta) = r_plasma(2,itheta,izeta) + lmnc(imn) * cosangle * cosangle3 * sinangle2
                 r_plasma(3,itheta,izeta) = r_plasma(3,itheta,izeta) + lmnc(imn) * cosangle * sinangle3

                 drdtheta_plasma(1,itheta,izeta) = drdtheta_plasma(1,itheta,izeta) + lmnc(imn) * (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * cosangle2
                 drdtheta_plasma(2,itheta,izeta) = drdtheta_plasma(2,itheta,izeta) + lmnc(imn) * (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * sinangle2
                 drdtheta_plasma(3,itheta,izeta) = drdtheta_plasma(3,itheta,izeta) + lmnc(imn) * (dcosangledtheta * sinangle3 + cosangle * dsinangle3dtheta)!
                 
                 drdzeta_plasma(1,itheta,izeta) = drdzeta_plasma(1,itheta,izeta) + lmnc(imn) * (dcosangledzeta * cosangle3 * cosangle2 + cosangle * dcosangle3dzeta * cosangle2 + cosangle * cosangle3 * dcosangle2dzeta)
                 drdzeta_plasma(2,itheta,izeta) = drdzeta_plasma(2,itheta,izeta) + lmnc(imn) * (dcosangledzeta * cosangle3 * sinangle2 + cosangle * dcosangle3dzeta * sinangle2 + cosangle * cosangle3 * dsinangle2dzeta)
                 drdzeta_plasma(3,itheta,izeta) = drdzeta_plasma(3,itheta,izeta) + lmnc(imn) * (dcosangledzeta * sinangle3 + cosangle * dsinangle3dzeta)!

                 if (geometry_option_coil == 5) then
                    d2rdtheta2_plasma(1,itheta,izeta) = d2rdtheta2_plasma(1,itheta,izeta) + lmnc(imn) * (d2cosangledtheta2 * cosangle3 + 2*dcosangledtheta * dcosangle3dtheta + cosangle * d2cosangle3dtheta2) * cosangle2!!
                    d2rdtheta2_plasma(2,itheta,izeta) = d2rdtheta2_plasma(2,itheta,izeta) + lmnc(imn) * (d2cosangledtheta2 * cosangle3 + 2*dcosangledtheta * dcosangle3dtheta + cosangle * d2cosangle3dtheta2) * sinangle2
                    d2rdtheta2_plasma(3,itheta,izeta) = d2rdtheta2_plasma(3,itheta,izeta) + lmnc(imn) * (d2cosangledtheta2 * sinangle3 + 2*dcosangledtheta * dsinangle3dtheta + cosangle * d2sinangle3dtheta2)

                    d2rdthetadzeta_plasma(1,itheta,izeta) = d2rdthetadzeta_plasma(1,itheta,izeta) + lmnc(imn) * ( (d2cosangledthetadzeta * cosangle3 + dcosangledtheta * dcosangle3dzeta + dcosangledzeta * dcosangle3dtheta + cosangle * d2cosangle3dthetadzeta) * cosangle2 &
                      + (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * dcosangle2dzeta )!
                    d2rdthetadzeta_plasma(2,itheta,izeta) = d2rdthetadzeta_plasma(2,itheta,izeta) + lmnc(imn) * ( (d2cosangledthetadzeta * cosangle3 + dcosangledtheta * dcosangle3dzeta + dcosangledzeta * dcosangle3dtheta + cosangle * d2cosangle3dthetadzeta) * sinangle2 &
                      + (dcosangledtheta * cosangle3 + cosangle * dcosangle3dtheta) * dsinangle2dzeta )
                    d2rdthetadzeta_plasma(3,itheta,izeta) = d2rdthetadzeta_plasma(3,itheta,izeta) + lmnc(imn) * (d2cosangledthetadzeta * sinangle3 + dcosangledtheta * dsinangle3dzeta + dcosangledzeta * dsinangle3dtheta + cosangle * d2sinangle3dthetadzeta)

                    d2rdzeta2_plasma(1,itheta,izeta) = d2rdzeta2_plasma(1,itheta,izeta) + lmnc(imn) * ( (d2cosangledzeta2 * cosangle3 + 2*dcosangledzeta * dcosangle3dzeta + cosangle * d2cosangle3dzeta2) * cosangle2 + 2*(dcosangledzeta * cosangle3 + cosangle * dcosangle3dzeta)*dcosangle2dzeta + cosangle * cosangle3 * d2cosangle2dzeta2 )
                    d2rdzeta2_plasma(2,itheta,izeta) = d2rdzeta2_plasma(2,itheta,izeta) + lmnc(imn) * ( (d2cosangledzeta2 * cosangle3 + 2*dcosangledzeta * dcosangle3dzeta + cosangle * d2cosangle3dzeta2) * sinangle2 + 2*(dcosangledzeta * cosangle3 + cosangle * dcosangle3dzeta)*dsinangle2dzeta + cosangle * cosangle3 * d2sinangle2dzeta2 )
                    d2rdzeta2_plasma(3,itheta,izeta) = d2rdzeta2_plasma(3,itheta,izeta) + lmnc(imn) * (d2cosangledzeta2 * sinangle3 + 2*dcosangledzeta * dsinangle3dzeta + cosangle * d2sinangle3dzeta2)!!
                    f_test(imn,itheta,izeta) = d2sinangle3dzeta2
                 end if

                 if (lasym) then
                    r_plasma(1,itheta,izeta) = r_plasma(1,itheta,izeta) + lmns(imn) * sinangle * cosangle3 * cosangle2
                    r_plasma(2,itheta,izeta) = r_plasma(2,itheta,izeta) + lmns(imn) * sinangle * cosangle3 * sinangle2
                    r_plasma(3,itheta,izeta) = r_plasma(3,itheta,izeta) + lmns(imn) * sinangle * sinangle3

                    drdtheta_plasma(1,itheta,izeta) = drdtheta_plasma(1,itheta,izeta) + lmns(imn) * (dsinangledtheta * cosangle3 + sinangle * dcosangle3dtheta) * cosangle2
                    drdtheta_plasma(2,itheta,izeta) = drdtheta_plasma(2,itheta,izeta) + lmns(imn) * (dsinangledtheta * cosangle3 + sinangle * dcosangle3dtheta) * sinangle2
                    drdtheta_plasma(3,itheta,izeta) = drdtheta_plasma(3,itheta,izeta) + lmns(imn) * (dsinangledtheta * sinangle3 + sinangle * dsinangle3dtheta)

                    drdzeta_plasma(1,itheta,izeta) = drdzeta_plasma(1,itheta,izeta) + lmns(imn) * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta) * cosangle3
                    drdzeta_plasma(2,itheta,izeta) = drdzeta_plasma(2,itheta,izeta) + lmns(imn) * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta) * cosangle3
                    drdzeta_plasma(3,itheta,izeta) = drdzeta_plasma(3,itheta,izeta) + lmns(imn) * dsinangledzeta * sinangle3
                 end if

             case default
                 r_plasma(1,itheta,izeta) = r_plasma(1,itheta,izeta) + rmnc_plasma(imn) * cosangle * cosangle2
                 r_plasma(2,itheta,izeta) = r_plasma(2,itheta,izeta) + rmnc_plasma(imn) * cosangle * sinangle2
                 r_plasma(3,itheta,izeta) = r_plasma(3,itheta,izeta) + zmns_plasma(imn) * sinangle
                 
                 drdtheta_plasma(1,itheta,izeta) = drdtheta_plasma(1,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * cosangle2
                 drdtheta_plasma(2,itheta,izeta) = drdtheta_plasma(2,itheta,izeta) + rmnc_plasma(imn) * dcosangledtheta * sinangle2
                 drdtheta_plasma(3,itheta,izeta) = drdtheta_plasma(3,itheta,izeta) + zmns_plasma(imn) * dsinangledtheta
                 
                 drdzeta_plasma(1,itheta,izeta) = drdzeta_plasma(1,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta)
                 drdzeta_plasma(2,itheta,izeta) = drdzeta_plasma(2,itheta,izeta) + rmnc_plasma(imn) * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta)
                 drdzeta_plasma(3,itheta,izeta) = drdzeta_plasma(3,itheta,izeta) + zmns_plasma(imn) * dsinangledzeta

                 if (lasym) then
                    r_plasma(1,itheta,izeta) = r_plasma(1,itheta,izeta) + rmns_plasma(imn) * sinangle * cosangle2
                    r_plasma(2,itheta,izeta) = r_plasma(2,itheta,izeta) + rmns_plasma(imn) * sinangle * sinangle2
                    r_plasma(3,itheta,izeta) = r_plasma(3,itheta,izeta) + zmnc_plasma(imn) * cosangle
                    
                    drdtheta_plasma(1,itheta,izeta) = drdtheta_plasma(1,itheta,izeta) + rmns_plasma(imn) * dsinangledtheta * cosangle2
                    drdtheta_plasma(2,itheta,izeta) = drdtheta_plasma(2,itheta,izeta) + rmns_plasma(imn) * dsinangledtheta * sinangle2
                    drdtheta_plasma(3,itheta,izeta) = drdtheta_plasma(3,itheta,izeta) + zmnc_plasma(imn) * dcosangledtheta
                    
                    drdzeta_plasma(1,itheta,izeta) = drdzeta_plasma(1,itheta,izeta) + rmns_plasma(imn) * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
                    drdzeta_plasma(2,itheta,izeta) = drdzeta_plasma(2,itheta,izeta) + rmns_plasma(imn) * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
                    drdzeta_plasma(3,itheta,izeta) = drdzeta_plasma(3,itheta,izeta) + zmnc_plasma(imn) * dcosangledzeta
                 end if
             end select
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! Evaluate cross product
    normal_plasma(1,:,:) = drdzeta_plasma(2,:,:) * drdtheta_plasma(3,:,:) - drdtheta_plasma(2,:,:) * drdzeta_plasma(3,:,:)
    normal_plasma(2,:,:) = drdzeta_plasma(3,:,:) * drdtheta_plasma(1,:,:) - drdtheta_plasma(3,:,:) * drdzeta_plasma(1,:,:)
    normal_plasma(3,:,:) = drdzeta_plasma(1,:,:) * drdtheta_plasma(2,:,:) - drdtheta_plasma(1,:,:) * drdzeta_plasma(2,:,:)

    if (geometry_option_coil == 5) then
      dnormaldtheta_plasma(1,:,:) = d2rdthetadzeta_plasma(2,:,:) * drdtheta_plasma(3,:,:) + drdzeta_plasma(2,:,:) * d2rdtheta2_plasma(3,:,:) &
        - d2rdtheta2_plasma(2,:,:) * drdzeta_plasma(3,:,:) - drdtheta_plasma(2,:,:) * d2rdthetadzeta_plasma(3,:,:)
      dnormaldtheta_plasma(2,:,:) = d2rdthetadzeta_plasma(3,:,:) * drdtheta_plasma(1,:,:) + drdzeta_plasma(3,:,:) * d2rdtheta2_plasma(1,:,:) &
        - d2rdtheta2_plasma(3,:,:) * drdzeta_plasma(1,:,:) - drdtheta_plasma(3,:,:) * d2rdthetadzeta_plasma(1,:,:)
      dnormaldtheta_plasma(3,:,:) = d2rdthetadzeta_plasma(1,:,:) * drdtheta_plasma(2,:,:) + drdzeta_plasma(1,:,:) * d2rdtheta2_plasma(2,:,:) &
        - d2rdtheta2_plasma(1,:,:) * drdzeta_plasma(2,:,:) - drdtheta_plasma(1,:,:) * d2rdthetadzeta_plasma(2,:,:)

      dnormaldzeta_plasma(1,:,:) = d2rdzeta2_plasma(2,:,:) * drdtheta_plasma(3,:,:) + drdzeta_plasma(2,:,:) * d2rdthetadzeta_plasma(3,:,:) &
        - d2rdthetadzeta_plasma(2,:,:) * drdzeta_plasma(3,:,:) - drdtheta_plasma(2,:,:) * d2rdzeta2_plasma(3,:,:)
      dnormaldzeta_plasma(2,:,:) = d2rdzeta2_plasma(3,:,:) * drdtheta_plasma(1,:,:) + drdzeta_plasma(3,:,:) * d2rdthetadzeta_plasma(1,:,:) &
        - d2rdthetadzeta_plasma(3,:,:) * drdzeta_plasma(1,:,:) - drdtheta_plasma(3,:,:) * d2rdzeta2_plasma(1,:,:)
      dnormaldzeta_plasma(3,:,:) = d2rdzeta2_plasma(1,:,:) * drdtheta_plasma(2,:,:) + drdzeta_plasma(1,:,:) * d2rdthetadzeta_plasma(2,:,:) &
        - d2rdthetadzeta_plasma(1,:,:) * drdzeta_plasma(2,:,:) - drdtheta_plasma(1,:,:) * d2rdzeta2_plasma(2,:,:)
    end if
   
!do imn = 1,mnmax_plasma
!  if (maxval(abs(sum(f_test(imn,:,:),2))) > 1e-10) then
!    print *, xm_plasma(imn), xn_plasma(imn), maxval(abs(sum(f_test(imn,:,:),2)))
!  end if
!end do
!do imn = 1,mnmax_plasma
!  print *, xm_plasma(imn), xn_plasma(imn), maxval(abs(f_test(imn,:,:)))
!end do
!print *, ""
!stop
!print *, maxval(sum(domegadtheta_arclength,1))
!print *, maxval(sum(domegadzeta_arclength,2))
!print *, maxval(sum(d2omegadtheta2_arclength,1))
!print *, maxval(sum(d2omegadzeta2_arclength,2))
!print *, maxval(sum(d2omegadthetadzeta_arclength,1))
!print *, maxval(sum(d2omegadthetadzeta_arclength,2))
!print *, ""
!print *, maxval(abs(sum(drdtheta_plasma(1,:,:),1)))!
!print *, maxval(abs(sum(drdtheta_plasma(2,:,:),1)))!
!print *, maxval(abs(sum(drdtheta_plasma(3,:,:),1)))!
!print *, ""
!print *, maxval(abs(sum(drdzeta_plasma(1,:,:),2)))
!print *, maxval(abs(sum(drdzeta_plasma(2,:,:),2)))
!print *, maxval(abs(sum(drdzeta_plasma(3,:,:),2)))!
!print *, ""
!print *, maxval(abs(sum(d2rdtheta2_plasma(1,:,:),1)))!!
!print *, maxval(abs(sum(d2rdtheta2_plasma(2,:,:),1)))
!print *, maxval(abs(sum(d2rdtheta2_plasma(3,:,:),1)))
!print *, ""
!print *, maxval(abs(sum(d2rdzeta2_plasma(1,:,:),2)))
!print *, maxval(abs(sum(d2rdzeta2_plasma(2,:,:),2)))
!print *, maxval(abs(sum(d2rdzeta2_plasma(3,:,:),2)))!!!!!!
!print *, ""
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(1,:,:),1)))!
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(2,:,:),1)))
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(3,:,:),1)))
!print *,""
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(1,:,:),2)))
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(2,:,:),2)))
!print *, maxval(abs(sum(d2rdthetadzeta_plasma(3,:,:),2)))
!print *,""
!print *, maxval(abs(sum(dnormaldtheta_plasma(1,:,:),1)))
!print *, maxval(abs(sum(dnormaldtheta_plasma(2,:,:),1)))!!
!print *, maxval(abs(sum(dnormaldtheta_plasma(3,:,:),1)))!!
!print *, maxval(abs(sum(dnormaldzeta_plasma(1,:,:),2)))
!print *, maxval(abs(sum(dnormaldzeta_plasma(2,:,:),2)))
!print *, maxval(abs(sum(dnormaldzeta_plasma(3,:,:),2)))!
!print *,""
!stop
    if (allocated(norm_normal_plasma)) deallocate(norm_normal_plasma)
    allocate(norm_normal_plasma(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
    norm_normal_plasma = sqrt(normal_plasma(1,:,1:nzeta_plasma)**2 &
         + normal_plasma(2,:,1:nzeta_plasma)**2 &
         + normal_plasma(3,:,1:nzeta_plasma)**2)

    if (geometry_option_coil == 5) then
      if (allocated(norm_normal_plasma_full)) deallocate(norm_normal_plasma_full)
      allocate(norm_normal_plasma_full(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      norm_normal_plasma_full = sqrt(normal_plasma(1,:,:)**2 &
        + normal_plasma(2,:,:)**2 &
        + normal_plasma(3,:,:)**2)
    end if

    if (geometry_option_coil == 5) then
      if (allocated(dnorm_normaldtheta_plasma)) deallocate(dnorm_normaldtheta_plasma)
      allocate(dnorm_normaldtheta_plasma(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      if (allocated(dnorm_normaldzeta_plasma)) deallocate(dnorm_normaldzeta_plasma)
      allocate(dnorm_normaldzeta_plasma(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      dnorm_normaldtheta_plasma = (dnormaldtheta_plasma(1,:,:) * normal_plasma(1,:,:) &
        + dnormaldtheta_plasma(2,:,:) * normal_plasma(2,:,:) &
        + dnormaldtheta_plasma(3,:,:) * normal_plasma(3,:,:) ) / norm_normal_plasma_full
      dnorm_normaldzeta_plasma = (dnormaldzeta_plasma(1,:,:) * normal_plasma(1,:,:) &
        + dnormaldzeta_plasma(2,:,:) * normal_plasma(2,:,:) &
        + dnormaldzeta_plasma(3,:,:) * normal_plasma(3,:,:) ) / norm_normal_plasma_full
    end if

!print *, maxval(abs(sum(dnorm_normaldtheta_plasma,1)))!
!print *, maxval(abs(sum(dnorm_normaldzeta_plasma,2)))!
!print *,""
!print *, sum(abs(d2omegadtheta2_arclength))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2omegadthetadzeta_arclength))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2omegadzeta2_arclength))/(ntheta_plasma*nzetal_plasma)
!print *,""
!print *, sum(abs(dnormaldtheta_plasma(1,:,:) * normal_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldtheta_plasma(2,:,:) * normal_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldtheta_plasma(3,:,:) * normal_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldzeta_plasma(1,:,:) * normal_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldzeta_plasma(2,:,:) * normal_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldzeta_plasma(3,:,:) * normal_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *,""
!stop
!print *, sum(abs(drdtheta_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(drdtheta_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(drdtheta_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(drdzeta_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(drdzeta_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(drdzeta_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *,""
!stop
!print *, sum(abs(d2rdtheta2_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdtheta2_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdtheta2_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdthetadzeta_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdthetadzeta_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdthetadzeta_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdzeta2_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdzeta2_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(d2rdzeta2_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)!
!print *,""
!stop
!print *, sum(abs(dnormaldtheta_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldtheta_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldtheta_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *, sum(abs(dnormaldzeta_plasma(1,:,:)))/(ntheta_plasma*nzetal_plasma)!
!print *, sum(abs(dnormaldzeta_plasma(2,:,:)))/(ntheta_plasma*nzetal_plasma)!
!print *, sum(abs(dnormaldzeta_plasma(3,:,:)))/(ntheta_plasma*nzetal_plasma)
!print *,""
!stop
    
    dtheta_plasma = theta_plasma(2)-theta_plasma(1)
    dzeta_plasma = zeta_plasma(2)-zeta_plasma(1)
    
    area_plasma = nfp * dtheta_plasma * dzeta_plasma * sum(norm_normal_plasma)
!print *, area_plasma
!print *, dtheta_plasma * dzeta_plasma * sum(norm_normal_plasma_full)

    ! Compute plasma volume using \int (1/2) R^2 dZ dzeta.
    ! These quantities will be evaluated on the half theta grid, which is the natural grid for dZ,
    ! but we will need to interpolate R^2 from the full to half grid.
    allocate(major_R_squared(ntheta_plasma,nzetal_plasma))
    major_R_squared = r_plasma(1,:,:)*r_plasma(1,:,:) + r_plasma(2,:,:)*r_plasma(2,:,:)
    ! First handle the interior of the theta grid:
    volume_plasma = sum((major_R_squared(1:ntheta_plasma-1,:) + major_R_squared(2:ntheta_plasma,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
         * (r_plasma(3,2:ntheta_plasma,:)-r_plasma(3,1:ntheta_plasma-1,:))) ! dZ
    ! Add the contribution from the ends of the theta grid:
    volume_plasma = volume_plasma + sum((major_R_squared(1,:) + major_R_squared(ntheta_plasma,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
         * (r_plasma(3,1,:)-r_plasma(3,ntheta_plasma,:))) ! dZ
    volume_plasma = abs(volume_plasma * dzeta_plasma / 2) ! r_plasma includes all nfp periods already, so no factor of nfp needed.
    deallocate(major_R_squared)
    if (verbose) print "(a,es10.3,a,es10.3,a)"," Plasma surface area:",area_plasma," m^2. Volume:",volume_plasma," m^3."
    
    if (geometry_option_coil == 5) then
      if (allocated(mean_curvature)) deallocate(mean_curvature)
      allocate(mean_curvature(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      
      allocate(E_big(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(F_big(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(G_big(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(e_small(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(f_small(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(g_small(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(K_curvature(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(kappa1(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
      allocate(kappa2(ntheta_plasma, nzetal_plasma),stat=iflag)
      if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

      E_big = drdtheta_plasma(1,:,:)**2 + drdtheta_plasma(2,:,:)**2 + drdtheta_plasma(3,:,:)**2
      F_big = drdtheta_plasma(1,:,:)*drdzeta_plasma(1,:,:) + drdtheta_plasma(2,:,:)*drdzeta_plasma(2,:,:) + drdtheta_plasma(3,:,:)*drdzeta_plasma(3,:,:)
      G_big = drdzeta_plasma(1,:,:)**2 + drdzeta_plasma(2,:,:)**2 + drdzeta_plasma(3,:,:)**2
      e_small = -(d2rdtheta2_plasma(1,:,:)*normal_plasma(1,:,:) + d2rdtheta2_plasma(2,:,:)*normal_plasma(2,:,:) + d2rdtheta2_plasma(3,:,:)*normal_plasma(3,:,:)) / norm_normal_plasma_full
      f_small = -(d2rdthetadzeta_plasma(1,:,:)*normal_plasma(1,:,:) + d2rdthetadzeta_plasma(2,:,:)*normal_plasma(2,:,:) + d2rdthetadzeta_plasma(3,:,:)*normal_plasma(3,:,:)) / norm_normal_plasma_full
      g_small = -(d2rdzeta2_plasma(1,:,:)*normal_plasma(1,:,:) + d2rdzeta2_plasma(2,:,:)*normal_plasma(2,:,:) + d2rdzeta2_plasma(3,:,:)*normal_plasma(3,:,:)) / norm_normal_plasma_full

      mean_curvature = (e_small*G_big - 2*f_small*F_big + g_small*E_big)/(2*(E_big*G_big-F_big*F_big))
      K_curvature = (e_small*g_small - f_small*f_small)/(E_big*G_big-F_big*F_big)
      kappa1 = mean_curvature + sqrt(mean_curvature**2 - K_curvature)
      kappa2 = mean_curvature - sqrt(mean_curvature**2 - K_curvature)

      max_separation = 1/abs(minval(kappa2(:,1:nzeta_plasma)))
      
      mean_curvature = -mean_curvature
    end if

    select case (geometry_option_plasma)
    case (2,3,4,8)
       ! A VMEC wout file is available
       ! VMEC stores the toroidal Boozer component B_zeta as "bvco", using the HALF mesh
       net_poloidal_current_Amperes = 2*pi/mu0*(1.5_dp*bvco(ns)-0.5_dp*bvco(ns-1))
       ! curpol is a number which multiplies the data in the bnorm file.
       curpol = (2*pi/nfp)*(1.5_dp*bsubvmnc(1,ns) - 0.5_dp*bsubvmnc(1,ns-1))

       if (verbose) print *,"Overriding net_poloidal_current_Amperes with value from the VMEC wout file."
       if (verbose) print *,"G = ", net_poloidal_current_Amperes, " ; curpol = ", curpol
    case (9)
       if (verbose) print *,"Overriding net_poloidal_current_Amperes with value from the single Fourier file."
       if (verbose) print *,"G = ", net_poloidal_current_Amperes, " ; curpol = ", curpol
    case default
       if (abs(net_poloidal_current_Amperes-1)<1e-12) then
          if (verbose) print *,"No VMEC file is available, and the default value of net_poloidal_current_Amperes (=1) will be used."
       else
          if (verbose) print *,"No VMEC file is available, so net_poloidal_current_Amperes will be taken from the bdistrib input file."
       end if
    end select

    !Save the necessary quantities to describe the single Fourier representation
    if (geometry_option_plasma .eq. 8) then
        call regcoil_write_single_Fourier()
    end if

    call system_clock(toc)
    if (verbose) print *,"Done initializing plasma surface. Took ",real(toc-tic)/countrate," sec."

  end subroutine regcoil_init_plasma

  ! --------------------------------------------------------------------------

  function fzero_residual(theta_old)

    use regcoil_variables, only: nfp
    use read_wout_mod, only: xm_vmec => xm, xn_vmec => xn, mnmax_vmec => mnmax, lmns, ns
    use stel_constants
    
    implicit none
    
    real(dp) :: theta_old, fzero_residual
    integer :: imn

    ! residual = twopi*(u_new - u_new_target) = (twopi*u_old + lambda) - u_new_target*twopi
    fzero_residual = theta_old - theta_rootSolve_target

    do imn = 1, mnmax_vmec
       fzero_residual = fzero_residual + lmns(imn,ns)*sin(xm_vmec(imn)*theta_old - xn_vmec(imn)*zeta)
    end do

  end function fzero_residual

  ! --------------------------------------------------------------------------

  subroutine regcoil_allocate_plasma_surface_arrays()

    use regcoil_variables

    implicit none

    integer :: iflag

    if (allocated(xm_plasma)) deallocate(xm_plasma)
    allocate(xm_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  1"

    if (allocated(xn_plasma)) deallocate(xn_plasma)
    allocate(xn_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  2"

    if (allocated(rmnc_plasma)) deallocate(rmnc_plasma)
    allocate(rmnc_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  3"

    if (allocated(rmns_plasma)) deallocate(rmns_plasma)
    allocate(rmns_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  4"

    if (allocated(zmnc_plasma)) deallocate(zmnc_plasma)
    allocate(zmnc_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  5"

    if (allocated(zmns_plasma)) deallocate(zmns_plasma)
    allocate(zmns_plasma(mnmax_plasma),stat=iflag)
    if (iflag .ne. 0) stop "regcoil_init_plasma Allocation error!  6"
    
    xm_plasma = 0
    xn_plasma = 0
    rmnc_plasma = 0
    rmns_plasma = 0
    zmnc_plasma = 0
    zmns_plasma = 0

  end subroutine regcoil_allocate_plasma_surface_arrays

end module regcoil_init_plasma_mod
