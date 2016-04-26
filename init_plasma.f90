subroutine init_plasma

  use global_variables
  use read_wout_mod, only: nfp_vmec => nfp, xm_vmec => xm, xn_vmec => xn, &
       rmnc_vmec => rmnc, zmns_vmec => zmns, rmns_vmec => rmns, zmnc_vmec => zmnc, &
       lasym_vmec => lasym, mnmax_vmec => mnmax, ns, Rmajor, read_wout_file, lmns, &
       mpol_vmec => mpol, ntor_vmec => ntor, bvco, bsubvmnc
  use safe_open_mod
  use stel_constants

  implicit none

  integer :: i, iu, iv, imn, tic, toc, countrate, iflag, ierr, iopen, tic1, toc1, iunit
  real(dp) :: angle, sinangle, cosangle, dsinangledu, dsinangledv, dcosangledu, dcosangledv
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv
  real(dp) :: weight1, weight2, u, v, r_temp, z_temp, dnorm
  integer :: nu_coordTransform, nv_coordTransform
  real(dp), dimension(:,:), allocatable :: r_coordTransform, z_coordTransform
  real(dp), dimension(:), allocatable :: rmnc_vmecLast, zmns_vmecLast
  real(dp) :: rootSolve_abserr, rootSolve_relerr, u_rootSolve_min, u_rootSolve_max, u_rootSolve_target, u_rootSolve_soln
  integer :: fzeroFlag, mpol, ntor, jm, jn, index

  call system_clock(tic, countrate)
  print *,"Initializing plasma surface."

  select case (geometry_option_plasma)
  case (0,1)
     ! Plain circular torus
     print *,"  Building a plain circular torus."

     nfp = nfp_imposed
     mnmax = 2
     lasym = .false.

     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmnc(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     xm = (/0,1/)
     xn = (/0,0/)
     rmnc(1) = R0_plasma
     rmnc(2) = a_plasma
     zmns(1) = 0
     zmns(2) = a_plasma
     
  case(6)
     ! Read in an ASCII table
     call safe_open(iunit, ierr, trim(shape_filename_plasma), 'old', 'formatted')
     if (ierr .ne. 0) then
        stop 'Error opening nescin file'
     endif

     ! Skip first line
     read (iunit, *)
     read (iunit, *) mnmax

     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! 1"
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error!  2"
     allocate(rmnc(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error!  3"
     allocate(zmns(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error!  4"
     allocate(rmns(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error!  5"
     allocate(zmnc(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error!  6"

     ! Skip a line
     read (iunit, *)
     do i = 1, mnmax
        read (iunit, *) xm(i), xn(i), rmnc(i), zmns(i), rmns(i), zmnc(i)
     end do

     close(iunit)

     nfp = nfp_imposed
     lasym = .true.

  case (2,3)
     ! VMEC, "original" theta coordinate which is not a straight-field-line coordinate
     call read_wout_file(wout_filename, ierr, iopen)
     if (iopen .ne. 0) stop 'error opening wout file'
     if (ierr .ne. 0) stop 'error reading wout file'
     print *,"  Successfully read VMEC data from ",trim(wout_filename)

     if (geometry_option_plasma == 2) then
        ! Only use the outermost point in the full radial mesh:
        weight1 = 0
        weight2 = 1
        print *,"  Using outermost grid point in VMEC's FULL radial grid."
     else
        ! Average the two outermost points in the full radial mesh 
        ! to get a value on the outermost point of the half radial mesh:
        weight1 = 0.5_dp
        weight2 = 0.5_dp
        print *,"  Using outermost grid point in VMEC's HALF radial grid."
     end if

     nfp = nfp_vmec
     mnmax = mnmax_vmec
     lasym = lasym_vmec
     R0_plasma = Rmajor

     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmnc(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     
     xm = xm_vmec
     xn = xn_vmec
     print *,"size of rmnc_vmec:",size(rmnc_vmec,1),size(rmnc_vmec,2)
     rmnc = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
     zmns = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2
     if (lasym) then
        allocate(rmns(mnmax),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        allocate(zmnc(mnmax),stat=iflag)
        if (iflag .ne. 0) stop 'Allocation error!'
        rmns = rmns_vmec(:,ns-1) * weight1 + rmns_vmec(:,ns) * weight2
        zmnc = zmnc_vmec(:,ns-1) * weight1 + zmnc_vmec(:,ns) * weight2
     end if

  case (4)
     ! VMEC, straight-field-line poloidal coordinate
     call read_wout_file(wout_filename, ierr, iopen)
     if (iopen .ne. 0) stop 'error opening wout file'
     if (ierr .ne. 0) stop 'error reading wout file'
     print *,"  Successfully read VMEC data from ",trim(wout_filename)


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
     allocate(rmnc_vmecLast(mnmax_vmec),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns_vmecLast(mnmax_vmec),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     rmnc_vmecLast = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
     zmns_vmecLast = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2

     ! Since the "original" vmec poloidal angle is chosen to have a very condensed
     ! Fourier spectrum, we probably need more Fourier modes to represent the surface using the
     ! straight-field-line coordinate.
     mpol = mpol_vmec*mpol_transform_refinement
     ntor = ntor_vmec*ntor_transform_refinement

     ! Beginning of coordinate transformation.
     ! Set up high-resolution grid in the "new" theta coordinate:
     nu_coordTransform = mpol * 2 
     nv_coordTransform = ntor * 2
     allocate(r_coordTransform(nu_coordTransform, nv_coordTransform), stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(z_coordTransform(nu_coordTransform, nv_coordTransform), stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     r_coordTransform = 0
     z_coordTransform = 0

     call system_clock(tic1)
     rootSolve_abserr = 1.0e-10_dp
     rootSolve_relerr = 1.0e-10_dp
     !open(unit=5,file="testStraightFieldLines",status='new',form='formatted')
     !write (5,*) nu_coordTransform, nv_coordTransform
     do iv = 1,nv_coordTransform
        v = (iv-1.0_dp)/nv_coordTransform
        do iu = 1,nu_coordTransform
           ! For each value of the new coordinates, solve for the old theta:
           u_rootSolve_target = (iu-1.0_dp)/nu_coordTransform
           u_rootSolve_min = u_rootSolve_target - 0.3
           u_rootSolve_max = u_rootSolve_target + 0.3

           call fzero(fzero_residual, u_rootSolve_min, u_rootSolve_max, u_rootSolve_target, &
                rootSolve_relerr, rootSolve_abserr, fzeroFlag)
           ! Note: fzero returns its answer in u_rootSolve_min
           u_rootSolve_soln = u_rootSolve_min
           if (fzeroFlag == 4) then
              stop "ERROR: fzero returned error 4: no sign change in residual"
           else if (fzeroFlag > 2) then
              print *,"WARNING: fzero returned an error code:",fzeroFlag
           end if
           ! Now that we have the old theta, evaluate r and z:
           r_temp = 0
           z_temp = 0
           do imn = 1, mnmax_vmec
              angle = twopi*(xm_vmec(imn)*u_rootSolve_soln - xn_vmec(imn)*v/nfp)
              r_temp = r_temp + rmnc_vmecLast(imn)*cos(angle)
              z_temp = z_temp + zmns_vmecLast(imn)*sin(angle)
           end do
           r_coordTransform(iu,iv) = r_temp
           z_coordTransform(iu,iv) = z_temp
           !write(5,*) u_rootSolve_soln
        end do
     end do
     !close(unit=5)
     call system_clock(toc1)
     print *,"  Time for root solving:",real(toc1-tic1)/countrate

     ! Now that we have R and Z on a grid in the new coordinates, Fourier transform the results.

     ! The next bit of code is much like initFourierModesMod, but with 2 differences: 
     ! 1. We need to keep the m=n=0 mode.
     ! 2. We follow VMEC convention that n includes the factor of nfp.

     ! xm is nonnegative.
     ! xn can be negative, zero, or positive.
     ! When xm is 0, xn must be positive.
     mnmax = mpol*(ntor*2+1) + ntor + 1
    
     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmnc(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns(mnmax),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
    
     ! Handle the xm=0 modes:
     xm=0
     do jn=0,ntor
        xn(jn+1)=jn*nfp
     end do
    
     ! Handle the xm>0 modes:
     index = ntor + 1
     do jm = 1,mpol
        do jn = -ntor, ntor
           index = index + 1
           xn(index) = jn*nfp
           xm(index) = jm
        end do
     end do
     ! Initialization of xm and xn is now complete.

     call system_clock(tic1)
     do imn = 1, mnmax
        dnorm = (1.0_dp)/(nu_coordTransform*nv_coordTransform)
        if (xm(imn).ne.0 .or. xn(imn).ne.0) dnorm = 2*dnorm
        r_temp = 0
        z_temp = 0
        do iv = 1, nv_coordTransform
           v = (iv-1.0_dp)/nv_coordTransform
           do iu = 1, nu_coordTransform
              u = (iu-1.0_dp)/nu_coordTransform
              angle = twopi*(xm(imn)*u-xn(imn)*v/nfp)
              cosangle = cos(angle)
              sinangle = sin(angle)
              r_temp = r_temp + r_coordTransform(iu,iv) * cosangle
              z_temp = z_temp + z_coordTransform(iu,iv) * sinangle
           end do
        end do
        rmnc(imn) = r_temp*dnorm
        zmns(imn) = z_temp*dnorm
     end do
     call system_clock(toc1)
     print *,"  Time for Fourier transform:",real(toc1-tic1)/countrate

  case (5)
     ! EFIT

     lasym = .true.
     nfp = nfp_imposed
     mnmax = efit_num_modes
     allocate(xm(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(xn(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmnc(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmns(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmnc(efit_num_modes),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     call read_efit(efit_filename, efit_psiN, efit_num_modes, rmnc, zmns, rmns, zmnc)

     ! Set major radius equal to the zero-frequency component of R(theta)
     R0_plasma = rmnc(1)

     xn = 0
     do i=1,efit_num_modes
        xm(i) = i-1
     end do

  case default
     print *,"Error! Invalid setting for geometry_option_plasma:",geometry_option_plasma
     stop
  end select

!!$  ! Save plasma surface shape for NESCOIL
!!$  iunit = 15
!!$  call safe_open(iunit,ierr,'nescin_plasma_surface','replace','formatted')
!!$  if (ierr .ne. 0) stop "Unable to open nescin output file"
!!$  write (iunit, 10) '------ Plasma Surface ---- '
!!$  write (iunit, 10) 'Number of fourier modes in table'
!!$  write (iunit,*) mnmax
!!$  write (iunit, 10) 'Table of fourier coefficients'
!!$  write (iunit, 10) 'm,n,crc,czs,cls,crs,czc,clc'
!!$  do imn = 1, mnmax
!!$     write (iunit,'(x,2i6,1p6e20.12)') xm(imn), xn(imn), &
!!$          rmnc(imn), zmns(imn), cl(m,n), crs(m,n), czc(m,n), clc(m,n)
!!$  end do
!!$
!!$  close(iunit)
     


  nvl_plasma = nv_plasma * nfp
  nvl_middle = nv_middle * nfp
  nvl_outer = nv_outer * nfp

  allocate(u_plasma(nu_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(v_plasma(nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(vl_plasma(nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do i=1,nu_plasma
     u_plasma(i) = (i-1.0_dp)/nu_plasma
  end do

  do i=1,nv_plasma
     v_plasma(i) = (i-1.0_dp)/nv_plasma
  end do

  do i=1,nvl_plasma
     vl_plasma(i) = (i-1.0_dp)/nv_plasma
  end do

  ! First coordinate is the Cartesian component x, y, or z
  allocate(r_plasma(3,nu_plasma,nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdu_plasma(3,nu_plasma,nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdv_plasma(3,nu_plasma,nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(normal_plasma(3,nu_plasma,nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  r_plasma=0
  drdu_plasma=0
  drdv_plasma=0
 
  do iv = 1,nvl_plasma
     angle2 = twopi*vl_plasma(iv)/nfp
     sinangle2 = sin(angle2)
     cosangle2 = cos(angle2)
     dsinangle2dv = cosangle2*twopi/nfp
     dcosangle2dv = -sinangle2*twopi/nfp
     do iu = 1,nu_plasma
        do imn = 1,mnmax
           angle = twopi*(xm(imn)*u_plasma(iu) - xn(imn)*vl_plasma(iv)/nfp)
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledu = cosangle*twopi*xm(imn)
           dcosangledu = -sinangle*twopi*xm(imn)
           dsinangledv = -cosangle*twopi*xn(imn)/nfp
           dcosangledv = sinangle*twopi*xn(imn)/nfp

           r_plasma(1,iu,iv) = r_plasma(1,iu,iv) + rmnc(imn) * cosangle * cosangle2
           r_plasma(2,iu,iv) = r_plasma(2,iu,iv) + rmnc(imn) * cosangle * sinangle2
           r_plasma(3,iu,iv) = r_plasma(3,iu,iv) + zmns(imn) * sinangle

           drdu_plasma(1,iu,iv) = drdu_plasma(1,iu,iv) + rmnc(imn) * dcosangledu * cosangle2
           drdu_plasma(2,iu,iv) = drdu_plasma(2,iu,iv) + rmnc(imn) * dcosangledu * sinangle2
           drdu_plasma(3,iu,iv) = drdu_plasma(3,iu,iv) + zmns(imn) * dsinangledu

           drdv_plasma(1,iu,iv) = drdv_plasma(1,iu,iv) + rmnc(imn) * (dcosangledv * cosangle2 + cosangle * dcosangle2dv)
           drdv_plasma(2,iu,iv) = drdv_plasma(2,iu,iv) + rmnc(imn) * (dcosangledv * sinangle2 + cosangle * dsinangle2dv)
           drdv_plasma(3,iu,iv) = drdv_plasma(3,iu,iv) + zmns(imn) * dsinangledv

           if (lasym) then
              r_plasma(1,iu,iv) = r_plasma(1,iu,iv) + rmns(imn) * sinangle * cosangle2
              r_plasma(2,iu,iv) = r_plasma(2,iu,iv) + rmns(imn) * sinangle * sinangle2
              r_plasma(3,iu,iv) = r_plasma(3,iu,iv) + zmnc(imn) * cosangle

              drdu_plasma(1,iu,iv) = drdu_plasma(1,iu,iv) + rmns(imn) * dsinangledu * cosangle2
              drdu_plasma(2,iu,iv) = drdu_plasma(2,iu,iv) + rmns(imn) * dsinangledu * sinangle2
              drdu_plasma(3,iu,iv) = drdu_plasma(3,iu,iv) + zmnc(imn) * dcosangledu

              drdv_plasma(1,iu,iv) = drdv_plasma(1,iu,iv) + rmns(imn) * (dsinangledv * cosangle2 + sinangle * dcosangle2dv)
              drdv_plasma(2,iu,iv) = drdv_plasma(2,iu,iv) + rmns(imn) * (dsinangledv * sinangle2 + sinangle * dsinangle2dv)
              drdv_plasma(3,iu,iv) = drdv_plasma(3,iu,iv) + zmnc(imn) * dcosangledv
           end if
        end do
     end do
  end do

  ! Evaluate cross product
  normal_plasma(1,:,:) = drdv_plasma(2,:,:) * drdu_plasma(3,:,:) - drdu_plasma(2,:,:) * drdv_plasma(3,:,:)
  normal_plasma(2,:,:) = drdv_plasma(3,:,:) * drdu_plasma(1,:,:) - drdu_plasma(3,:,:) * drdv_plasma(1,:,:)
  normal_plasma(3,:,:) = drdv_plasma(1,:,:) * drdu_plasma(2,:,:) - drdu_plasma(1,:,:) * drdv_plasma(2,:,:)

  allocate(norm_normal_plasma(nu_plasma, nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  norm_normal_plasma = sqrt(normal_plasma(1,:,:)**2 + normal_plasma(2,:,:)**2 + normal_plasma(3,:,:)**2)

  du_plasma = u_plasma(2)-u_plasma(1)
  dv_plasma = v_plasma(2)-v_plasma(1)

  area_plasma = du_plasma * dv_plasma * sum(norm_normal_plasma)

  select case (geometry_option_plasma)
  case (2,3,4)
     ! A VMEC wout file is available
     print *,"Overriding net_poloidal_current_Amperes with value from the VMEC wout file."
     ! VMEC stores the toroidal Boozer component B_zeta as "bvco", using the HALF mesh
     net_poloidal_current_Amperes = 2*pi/mu0*(1.5_dp*bvco(ns)-0.5_dp*bvco(ns-1))
     ! curpol is a number which multiplies the data in the bnorm file.
     curpol = (2*pi/nfp)*(1.5_dp*bsubvmnc(1,ns) - 0.5_dp*bsubvmnc(1,ns-1))
  case default
     if (abs(net_poloidal_current_Amperes-1)<1e-12) then
        print *,"No VMEC file is available, and the default value of net_poloidal_current_Amperes (=1) will be used."
     else
        print *,"No VMEC file is available, so net_poloidal_current_Amperes will be taken from the bdistrib input file."
     end if
  end select

  call system_clock(toc)
  print *,"Done initializing plasma surface. Took ",real(toc-tic)/countrate," sec."

contains

  function fzero_residual(u_old)

    implicit none

    real(dp) :: u_old, fzero_residual
    integer :: imn

    ! residual = twopi*(u_new - u_new_target) = (twopi*u_old + lambda) - u_new_target*twopi
    fzero_residual = twopi*(u_old - u_rootSolve_target)

    do imn = 1, mnmax_vmec
       fzero_residual = fzero_residual + lmns(imn,ns)*sin(twopi*(xm_vmec(imn)*u_old - xn_vmec(imn)/nfp*v))
    end do

  end function fzero_residual

end subroutine init_plasma
