module regcoil_init_plasma_mod


  use stel_kinds

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
         lasym_vmec => lasym, mnmax_vmec => mnmax, ns, Rmajor, read_wout_file, &
         mpol_vmec => mpol, ntor_vmec => ntor, bvco, bsubvmnc
    use safe_open_mod
    use stel_constants
    
    implicit none

    integer :: i, itheta, izeta, imn, tic, toc, countrate, iflag, ierr, iopen, tic1, toc1, iunit
    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
    real(dp) :: weight1, weight2, theta, r_temp, z_temp, dnorm
    integer :: ntheta_coordTransform, nzeta_coordTransform
    real(dp), dimension(:,:), allocatable :: r_coordTransform, z_coordTransform, major_R_squared
    real(dp), dimension(:), allocatable :: rmnc_vmecLast, zmns_vmecLast
    real(dp) :: rootSolve_abserr, rootSolve_relerr, theta_rootSolve_min, theta_rootSolve_max, theta_rootSolve_soln
    integer :: fzeroFlag, mpol, ntor, jm, jn, index
    
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
       
    case (2,3)
       ! VMEC, "original" theta coordinate which is not a straight-field-line coordinate
       call read_wout_file(wout_filename, ierr, iopen)
       if (iopen .ne. 0) stop 'error opening wout file'
       if (ierr .ne. 0) stop 'error reading wout file'
       if (verbose) print *,"  Successfully read VMEC data from ",trim(wout_filename)
       
       if (geometry_option_plasma == 2) then
          ! Only use the outermost point in the full radial mesh:
          weight1 = 0
          weight2 = 1
          if (verbose) print *,"  Using outermost grid point in VMEC's FULL radial grid."
       else
          ! Average the two outermost points in the full radial mesh 
          ! to get a value on the outermost point of the half radial mesh:
          weight1 = 0.5_dp
          weight2 = 0.5_dp
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
       if (lasym) then
          rmns_plasma = rmns_vmec(:,ns-1) * weight1 + rmns_vmec(:,ns) * weight2
          zmnc_plasma = zmnc_vmec(:,ns-1) * weight1 + zmnc_vmec(:,ns) * weight2
       end if
       
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

       if (allocated(rmnc_vmecLast)) deallocate(rmnc_vmecLast)
       allocate(rmnc_vmecLast(mnmax_vmec),stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

       if (allocated(zmns_vmecLast)) deallocate(zmns_vmecLast)
       allocate(zmns_vmecLast(mnmax_vmec),stat=iflag)
       if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'

       rmnc_vmecLast = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
       zmns_vmecLast = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2
       
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
    
    r_plasma=0
    drdtheta_plasma=0
    drdzeta_plasma=0
    
    do izeta = 1,nzetal_plasma
       angle2 = zetal_plasma(izeta)
       sinangle2 = sin(angle2)
       cosangle2 = cos(angle2)
       dsinangle2dzeta = cosangle2
       dcosangle2dzeta = -sinangle2
       do itheta = 1,ntheta_plasma
          do imn = 1,mnmax_plasma
             angle = xm_plasma(imn)*theta_plasma(itheta) - xn_plasma(imn)*zetal_plasma(izeta)
             sinangle = sin(angle)
             cosangle = cos(angle)
             dsinangledtheta = cosangle*xm_plasma(imn)
             dcosangledtheta = -sinangle*xm_plasma(imn)
             dsinangledzeta = -cosangle*xn_plasma(imn)
             dcosangledzeta = sinangle*xn_plasma(imn)
             
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
          end do
       end do
    end do
    
    ! Evaluate cross product
    normal_plasma(1,:,:) = drdzeta_plasma(2,:,:) * drdtheta_plasma(3,:,:) - drdtheta_plasma(2,:,:) * drdzeta_plasma(3,:,:)
    normal_plasma(2,:,:) = drdzeta_plasma(3,:,:) * drdtheta_plasma(1,:,:) - drdtheta_plasma(3,:,:) * drdzeta_plasma(1,:,:)
    normal_plasma(3,:,:) = drdzeta_plasma(1,:,:) * drdtheta_plasma(2,:,:) - drdtheta_plasma(1,:,:) * drdzeta_plasma(2,:,:)
    
    if (allocated(norm_normal_plasma)) deallocate(norm_normal_plasma)
    allocate(norm_normal_plasma(ntheta_plasma, nzeta_plasma),stat=iflag)
    if (iflag .ne. 0) stop 'regcoil_init_plasma Allocation error!'
    norm_normal_plasma = sqrt(normal_plasma(1,:,1:nzeta_plasma)**2 &
         + normal_plasma(2,:,1:nzeta_plasma)**2 &
         + normal_plasma(3,:,1:nzeta_plasma)**2)
    
    dtheta_plasma = theta_plasma(2)-theta_plasma(1)
    dzeta_plasma = zeta_plasma(2)-zeta_plasma(1)
    
    area_plasma = nfp * dtheta_plasma * dzeta_plasma * sum(norm_normal_plasma)

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
    
    select case (geometry_option_plasma)
    case (2,3,4)
       ! A VMEC wout file is available
       ! VMEC stores the toroidal Boozer component B_zeta as "bvco", using the HALF mesh
       net_poloidal_current_Amperes = 2*pi/mu0*(1.5_dp*bvco(ns)-0.5_dp*bvco(ns-1)) - TF_coil_pol_current
       ! subtract the current that is already taken care of by the TF coils?
       if (helicity_ratio .ne. 0) then
          net_toroidal_current_Amperes = net_poloidal_current_Amperes / helicity_ratio / nfp
       else
          net_toroidal_current_Amperes = 0
       end if 
       ! curpol is a number which multiplies the data in the bnorm file.
       curpol = 1!I am supplying it in T already so should be ok... (2*pi/nfp)*(1.5_dp*bsubvmnc(1,ns) - 0.5_dp*bsubvmnc(1,ns-1))

       if (verbose) print *,"Overriding net_poloidal_current_Amperes with value from the VMEC wout file."
       if (verbose) print *,"G = ", net_poloidal_current_Amperes, " ; curpol = ", curpol
       if (verbose) print *,"I = ", net_toroidal_current_Amperes, " ; helicity_ratio = ", helicity_ratio

    case default
       if (abs(net_poloidal_current_Amperes-1)<1e-12) then
          if (verbose) print *,"No VMEC file is available, and the default value of net_poloidal_current_Amperes (=1) will be used."
       else
          if (verbose) print *,"No VMEC file is available, so net_poloidal_current_Amperes will be taken from the bdistrib input file."
       end if
    end select

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
