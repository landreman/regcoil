subroutine init_sensitivity()

  use global_variables
  use stel_constants
  use init_Fourier_modes_mod
  !use omp_lib

  implicit none

  integer :: iomega, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry, tic, toc, countrate,indexl_coil

  real(dp) :: angle, angle2, sinangle, cosangle
  real(dp) :: sinangle2, cosangle2
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
  integer :: izeta_plasma, itheta_plasma, izeta_coil, itheta_coil
  integer :: index_coil, l_coil, izetal_coil, index_plasma

  ! For volume calculation
  real(dp) :: major_R_squared_half_grid, dzdtheta_coil_half_grid
  integer :: index_coil_first, index_coil_last
  real(dp), dimension(:), allocatable :: dR_squared_domega_half_grid,domegadzdtheta_half_grid

  ! Initialize Fourier arrays - need to include m = 0, n = 0 mode
  call init_Fourier_modes_sensitivity(mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, nomega_coil, &
    xm_sensitivity, xn_sensitivity, omega_coil)

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

  call system_clock(tic,countrate)

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
  print *,"Loop over coil geom in init_sensitivity:",real(toc-tic)/countrate," sec."

  call system_clock(tic,countrate)

  ! Compute Volume Sensitivity - need arrays on half theta grid
  volume_coil_with_dzdtheta = 0
  allocate(dvolume_coildomega(nomega_coil))
  allocate(dR_squared_domega_half_grid(nomega_coil))
  allocate(domegadzdtheta_half_grid(nomega_coil))
  dvolume_coildomega = 0
  do itheta_coil=1,ntheta_coil-1
    do izeta_coil=1,nzeta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
      do l_coil = 0, (nfp-1)
        izetal_coil = izeta_coil + l_coil*nzeta_coil

        major_R_squared_half_grid = (r_coil(1,itheta_coil,izetal_coil) &
          *r_coil(1,itheta_coil,izetal_coil) &
          + r_coil(2,itheta_coil,izetal_coil)*r_coil(2,itheta_coil,izetal_coil) &
          + r_coil(1,itheta_coil+1,izetal_coil)*r_coil(1,itheta_coil+1,izetal_coil) &
          + r_coil(2,itheta_coil+1,izetal_coil)*r_coil(2,itheta_coil+1,izetal_coil))/2
        dzdtheta_coil_half_grid = (drdtheta_coil(3,itheta_coil,izetal_coil) &
          + drdtheta_coil(3,itheta_coil,izetal_coil))/2

        dR_squared_domega_half_grid = &
          drdomega(1,index_coil,l_coil+1,:)*r_coil(1,itheta_coil,izetal_coil) &
          + drdomega(2,index_coil,l_coil+1,:)*r_coil(2,itheta_coil,izetal_coil) &
          + drdomega(1,index_coil+1,l_coil+1,:)*r_coil(1,itheta_coil+1,izetal_coil) &
          + drdomega(2,index_coil+1,l_coil+1,:)*r_coil(2,itheta_coil+1,izetal_coil)
        domegadzdtheta_half_grid = &
          (domegadzdtheta(:,itheta_coil,izetal_coil) + domegadzdtheta(:,itheta_coil+1,izetal_coil))/2

        dvolume_coildomega = dvolume_coildomega + (dR_squared_domega_half_grid*dzdtheta_coil_half_grid &
          + major_R_squared_half_grid*domegadzdtheta_half_grid)*dtheta_coil*dzeta_coil/2
        volume_coil_with_dzdtheta = volume_coil_with_dzdtheta + major_R_squared_half_grid*dzdtheta_coil_half_grid &
          * dtheta_coil*dzeta_coil/2
      end do
    end do
  end do
  do izeta_coil=1,nzeta_coil
    index_coil_first = (izeta_coil-1)*ntheta_coil + 1
    index_coil_last = (izeta_coil-1)*ntheta_coil + ntheta_coil
    do l_coil = 0, (nfp-1)
      izetal_coil = izeta_coil + l_coil*nzeta_coil

      major_R_squared_half_grid = (r_coil(1,ntheta_coil,izetal_coil)*r_coil(1,ntheta_coil,izetal_coil) &
        + r_coil(2,ntheta_coil,izetal_coil)*r_coil(2,ntheta_coil,izetal_coil) &
        + r_coil(1,1,izetal_coil)*r_coil(1,1,izetal_coil) &
        + r_coil(2,1,izetal_coil)*r_coil(2,1,izetal_coil))/2
      dzdtheta_coil_half_grid = (drdtheta_coil(3,ntheta_coil,izetal_coil) + drdtheta_coil(3,1,izetal_coil))/2

      dR_squared_domega_half_grid = &
        drdomega(1,index_coil_first,l_coil+1,:)*r_coil(1,1,izetal_coil) &
        + drdomega(2,index_coil_first,l_coil+1,:)*r_coil(2,1,izetal_coil) &
        + drdomega(1,index_coil_last,l_coil+1,:)*r_coil(1,ntheta_coil,izetal_coil) &
        + drdomega(2,index_coil_last,l_coil+1,:)*r_coil(2,ntheta_coil,izetal_coil)
      domegadzdtheta_half_grid = &
        (domegadzdtheta(:,1,izetal_coil) + domegadzdtheta(1,ntheta_coil,izetal_coil))/2

      dvolume_coildomega = dvolume_coildomega + (dR_squared_domega_half_grid*dzdtheta_coil_half_grid &
        + major_R_squared_half_grid*domegadzdtheta_half_grid)*dtheta_coil*dzeta_coil/2
      volume_coil_with_dzdtheta = volume_coil_with_dzdtheta + major_R_squared_half_grid*dzdtheta_coil_half_grid &
        * dtheta_coil*dzeta_coil/2
    end do
  end do
  ! volume_coil_alt should be the same as volume_coil
  call system_clock(toc)
  print *,"Volume computed in init_sensitivity: ", real(volume_coil_with_dzdtheta), " m^3"
  call system_clock(toc)
  print *,"Coil volume computation:",real(toc-tic)/countrate," sec."
  deallocate(dR_squared_domega_half_grid)
  deallocate(domegadzdtheta_half_grid)

!  call system_clock(tic,countrate)
!
!  !$OMP PARALLEL
!  !$OMP MASTER
!  print *,"  Number of OpenMP threads:",omp_get_num_threads()
!  !$OMP END MASTER
!  !$OMP DO PRIVATE(itheta_plasma,index_plasma,x,y,z,izeta_coil,itheta_coil,index_coil,l_coil,izetal_coil,angle2,sinangle2,cosangle2,dx,dy,dz,dr2inv,dr32inv,dr52inv,dr_dot_norm_coil,dr_dot_norm_plasma,norm_plasma_dot_norm_coil,dinductancedr,dinductancednorm)
!  do izeta_plasma = 1, nzeta_plasma
!    do itheta_plasma = 1, ntheta_plasma
!      index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma
!      x = r_plasma(1,itheta_plasma,izeta_plasma)
!      y = r_plasma(2,itheta_plasma,izeta_plasma)
!      z = r_plasma(3,itheta_plasma,izeta_plasma)
!      do izeta_coil = 1,nzeta_coil
!        do itheta_coil = 1,ntheta_coil
!          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
!          do l_coil = 0, (nfp-1)
!            izetal_coil = izeta_coil + l_coil*nzeta_coil
!
!            angle2 = zetal_coil(izetal_coil)
!            sinangle2 = sin(angle2)
!            cosangle2 = cos(angle2)
!
!            dx = x - r_coil(1,itheta_coil,izetal_coil)
!            dy = y - r_coil(2,itheta_coil,izetal_coil)
!            dz = z - r_coil(3,itheta_coil,izetal_coil)
!
!            dr2inv = 1/(dx*dx + dy*dy + dz*dz)
!            dr32inv = dr2inv*sqrt(dr2inv)
!            dr52inv = dr2inv*dr32inv
!
!            dr_dot_norm_coil = dx*normal_coil(1,itheta_coil,izetal_coil) &
!              + dy*normal_coil(2,itheta_coil,izetal_coil) + dz*normal_coil(3,itheta_coil,izetal_coil)
!            dr_dot_norm_plasma = dx*normal_plasma(1,itheta_plasma,izeta_plasma) &
!              + dy*normal_plasma(2,itheta_plasma,izeta_plasma) + dz*normal_plasma(3, itheta_plasma,izeta_plasma)
!            norm_plasma_dot_norm_coil = normal_plasma(1,itheta_plasma,izeta_plasma) &
!              *normal_coil(1,itheta_coil,izetal_coil) + normal_plasma(2,itheta_plasma,izeta_plasma) &
!              *normal_coil(2,itheta_coil,izetal_coil) + normal_plasma(3,itheta_plasma,izeta_plasma) &
!              *normal_coil(3,itheta_coil,izetal_coil)
!
!            dinductancedr(1,l_coil+1) = (3*dx*norm_plasma_dot_norm_coil &
!              - 15*dx*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
!              + 3*(normal_plasma(1,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
!              + normal_coil(1,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
!            dinductancedr(2,l_coil+1) = (3*dy*norm_plasma_dot_norm_coil &
!              - 15*dy*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
!              + 3*(normal_plasma(2,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
!              + normal_coil(2,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
!            dinductancedr(3,l_coil+1) = (3*dz*norm_plasma_dot_norm_coil &
!              - 15*dz*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
!              + 3*(normal_plasma(3,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
!              + normal_coil(3,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
!
!            dinductancednorm(1,l_coil+1) = (normal_plasma(1,itheta_plasma,izeta_plasma) &
!              - 3*dr2inv*dr_dot_norm_plasma*dx)*(dr32inv*mu0/(4*pi))
!            dinductancednorm(2,l_coil+1) = (normal_plasma(2,itheta_plasma,izeta_plasma) &
!              - 3*dr2inv*dr_dot_norm_plasma*dy)*(dr32inv*mu0/(4*pi))
!            dinductancednorm(3,l_coil+1) = (normal_plasma(3,itheta_plasma,izeta_plasma) &
!              - 3*dr2inv*dr_dot_norm_plasma*dz)*(dr32inv*mu0/(4*pi))
!
!            dinductancedomega(:,index_plasma,index_coil) = &
!                dinductancedomega(:,index_plasma,index_coil) &
!              + dinductancednorm(1,l_coil+1)*dnormxdomega(:,index_coil,l_coil+1) &
!              + dinductancednorm(2,l_coil+1)*dnormydomega(:,index_coil,l_coil+1) &
!              + dinductancednorm(3,l_coil+1)*dnormzdomega(:,index_coil,l_coil+1) &
!              + dinductancedr(1,l_coil+1)*drdomega(1,index_coil,l_coil+1,:) &
!              + dinductancedr(2,l_coil+1)*drdomega(2,index_coil,l_coil+1,:) &
!              + dinductancedr(3,l_coil+1)*drdomega(3,index_coil,l_coil+1,:)
!          end do
!          ! Using matmuls seems to be about twice as slow
!          !dinductancedomega(:,index_plasma,index_coil) = &
!           !   matmul(dnormxdomega(:, index_coil,:),dinductancednorm(1,:)) &
!           ! + matmul(dnormydomega(:, index_coil,:),dinductancednorm(2,:)) &
!           ! + matmul(dnormzdomega(:, index_coil,:),dinductancednorm(3,:)) &
!           ! + matmul(dinductancedr(1,:),drdomega(1,index_coil,:,:)) &
!           ! + matmul(dinductancedr(2,:),drdomega(2,index_coil,:,:)) &
!           ! + matmul(dinductancedr(3,:),drdomega(3,index_coil,:,:))
!        end do
!      end do
!    end do
!  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

!  call system_clock(toc)
!  print *,"Inductance sensitivity loop in init_sensitivity:",real(toc-tic)/countrate," sec."
!
!  call system_clock(tic,countrate)

  ! Now need to multiply dinductancedomega by basis_functions to compute sensitivity of g_j
  ! dincutancedomega(nomega_coil, ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil)
  ! basis_functions(ntheta_coil*nzeta_coil, num_basis_functions)
  ! dgdomega(nomega_coil, ntheta_plasma*nzeta_plasma, num_basis_functions)

!  ! Now we need to compute dgd(...) = dinductanced(...)*basis_functions
!  M = ntheta_plasma*nzeta_plasma ! # rows of A
!  N = num_basis_functions ! # cols of B
!  K = ntheta_coil*nzeta_coil ! Common dimension of A and B
!  LDA = M
!  LDB = K
!  LDC = M
!  TRANSA = 'N' ! No transposes
!  TRANSB = 'N'
!  BLAS_ALPHA=dtheta_coil*dzeta_coil
!  BLAS_BETA=0
!  dgdomega = 0
!
!  !$OMP PARALLEL
!  !$OMP MASTER
!  print *,"  Number of OpenMP threads:",omp_get_num_threads()
!  !$OMP END MASTER
!  !$OMP DO
!  do iomega = 1, nomega_coil
!    ! With matmuls is about 10x slower so commenting out
!    !dgdomega(iomega,:,:) = dtheta_coil*dzeta_coil*matmul(dinductancedomega(iomega,:,:), basis_functions)
!    call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedomega(iomega,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdomega(iomega,:,:),LDC)
!  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL
!
!  call system_clock(toc)
!  print *,"loop with DGEMM for g sensitivity: ",real(toc-tic)/countrate," sec."

  call system_clock(tic,countrate)

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
  print *,"d sensitivity loop: ",real(toc-tic)/countrate," sec."

!  call system_clock(tic,countrate)

  ! Needed for computing sensitivity of f
!  select case (symmetry_option)
!  case (1)
!    minSymmetry = 1
!    maxSymmetry = 1
!  case (2)
!    minSymmetry = 2
!    maxSymmetry = 2
!  case (3)
!    minSymmetry = 1
!    maxSymmetry = 2
!  end select
!
!  do izeta_coil = 1, nzeta_coil
!    do itheta_coil = 1, ntheta_coil
!      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
!      do iomega = 1,nomega_coil
!
!        do whichSymmetry = minSymmetry, maxSymmetry
!
!          if (whichSymmetry==2 .and. symmetry_option==3) then
!            offset = mnmax_coil
!          else
!            offset = 0
!          end if
!
!          do imn_coil = 1, mnmax_coil
!
!            angle_coil = xm_coil(imn_coil)*theta_coil(itheta_coil) - xn_coil(imn_coil)*zeta_coil(izeta_coil)
!            cosangle_coil = cos(angle_coil)
!            sinangle_coil = sin(angle_coil)
!            cosangle_xn = cosangle_coil*xn_coil(imn_coil)
!            sinangle_xn = sinangle_coil*xn_coil(imn_coil)
!            cosangle_xm = cosangle_coil*xm_coil(imn_coil)
!            sinangle_xm = sinangle_coil*xm_coil(imn_coil)
!            if (whichSymmetry==1) then
!              dfxdomega(iomega, index_coil,imn_coil) = &
!                cosangle_xn*domegadxdtheta(iomega,itheta_coil,izeta_coil) &
!                + cosangle_xm*domegadxdzeta(iomega,itheta_coil,izeta_coil)
!              dfydomega(iomega, index_coil,imn_coil) = &
!                cosangle_xn*domegadydtheta(iomega,itheta_coil,izeta_coil) &
!                + cosangle_xm*domegadydzeta(iomega,itheta_coil,izeta_coil)
!              dfzdomega(iomega, index_coil,imn_coil) = &
!                cosangle_xn*domegadzdtheta(iomega,itheta_coil,izeta_coil) &
!                + cosangle_xm*domegadzdzeta(iomega,itheta_coil,izeta_coil)
!            else
!              dfxdomega(iomega, index_coil,imn_coil+offset) = &
!                -sinangle_xn*domegadxdtheta(iomega,itheta_coil,izeta_coil) &
!                -sinangle_xm*domegadxdzeta(iomega,itheta_coil,izeta_coil)
!              dfydomega(iomega, index_coil,imn_coil+offset) = &
!                -sinangle_xn*domegadydtheta(iomega,itheta_coil,izeta_coil) &
!                -sinangle_xm*domegadydzeta(iomega,itheta_coil,izeta_coil)
!              dfzdomega(iomega, index_coil,imn_coil+offset) = &
!                -sinangle_xn*domegadzdtheta(iomega,itheta_coil,izeta_coil) &
!                -sinangle_xm*domegadzdzeta(iomega,itheta_coil,izeta_coil)
!            end if
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo

!  call system_clock(toc)
!  print *,"f sensitivity loop: ",real(toc-tic)/countrate," sec."
  print *,"Init sensitivity complete."

end subroutine init_sensitivity





