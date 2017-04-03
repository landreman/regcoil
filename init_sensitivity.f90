module init_sensitivity

	implicit none

	contains 

		subroutine init_partials()

			use global_variables
			use stel_constants
      use init_Fourier_modes_mod

			implicit none
			
			integer :: imn, iflag, minSymmetry, maxSymmetry, offset, whichSymmetry
			real(dp), dimension(:,:,:,:), allocatable :: dnormdrmnc, dnormdrmns, dnormdzmnc, dnormdzmns
      real(dp), dimension(:), allocatable :: drdrmnc, drdrmns, drdzmnc, drdzmns
      real(dp), dimension(:), allocatable :: dinductancednorm, dinductancedr
      real(dp), dimension(:,:,:), allocatable :: dinductancedrmnc, dinductancedrmns, dinductancedzmnc, dinductancedzmns

			real(dp) :: angle, angle2, sinangle, cosangle
			real(dp) :: sinangle2, cosangle2, dr_dot_norm_coil, dr_dot_norm_plasma, norm_plasma_dot_norm_coil
      real(dp) :: drmncdxdtheta, drmncdxdzeta, drmncdydtheta, drmncdydzeta
      real(dp) :: drmncdzdtheta, drmncdzdzeta, drmnsdxdtheta, drmnsdxdzeta
			real(dp) :: drmnsdydtheta, drmnsdydzeta, drmnsdzdtheta, drmnsdzdzeta
			real(dp) :: dzmncdxdtheta, dzmncdxdzeta, dzmncdydtheta, dzmncdydzeta
      real(dp) :: dzmncdzdtheta, dzmncdzdzeta, dzmnsdydzeta, dzmnsdzdtheta
      real(dp) :: dzmnsdzdzeta, dzmnsdxdtheta, dzmnsdxdzeta, dzmnsdydtheta
			real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
			integer :: izeta_plasma, itheta_plasma, izeta_coil, itheta_coil
			integer :: index_coil, l_coil, izetal_coil, index_plasma, imn_coil
			real(dp) :: x, y, z, dx, dy, dz, dr2inv, dr32inv, dr52inv
      real(dp) :: angle_coil, sinangle_coil, cosangle_coil

      ! Variables needed by BLAS DGEMM:
      character :: TRANSA, TRANSB
      integer :: M, N, K, LDA, LDB, LDC
      real(dp) :: BLAS_ALPHA=1, BLAS_BETA=0


			allocate(dnormdrmnc(3, mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnormdrmns(3, mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnormdzmnc(3, mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnormdzmns(3, mnmax_sensitivity,ntheta_coil,nzetal_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'

      allocate(drdrmnc(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdrmns(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdzmnc(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdzmns(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      allocate(dinductancednorm(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dinductancedr(3),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      allocate(dinductancedrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
        ntheta_coil*nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dinductancedrmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
        ntheta_coil*nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dinductancedzmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
        ntheta_coil*nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dinductancedzmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, &
        ntheta_coil*nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

			allocate(dnorm_normaldrmnc(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnorm_normaldrmns(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnorm_normaldzmnc(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dnorm_normaldzmns(mnmax_sensitivity,ntheta_coil,nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'

			allocate(dddrmnc(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dddrmns(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dddzmnc(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'
			allocate(dddzmns(3, mnmax_sensitivity,ntheta_coil*nzeta_coil),stat=iflag)
			if (iflag .ne. 0) stop 'Allocation error!'

      allocate(dfdrmnc(3, mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dfdrmns(3, mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dfdzmnc(3, mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dfdzmns(3, mnmax_sensitivity, ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      allocate(dgdrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dgdrmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dgdzmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(dgdzmns(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      ! Initialize Fourier arrays
      call init_Fourier_modes(mmax_sensitivity, nmax_sensitivity, mnmax_sensitivity, &
        xm_sensitivity, xn_sensitivity)

      ! Computing chi2_B sensitivity

      ! These quantities used to compute dnormd(...)
			drmncdzdtheta = 0
			drmncdzdzeta = 0
			drmnsdzdtheta = 0
			drmnsdzdzeta = 0

			dzmncdxdtheta = 0
			dzmncdxdzeta = 0
			dzmncdydtheta = 0
			dzmncdydzeta = 0
			dzmncdzdtheta = 0

			dzmnsdxdtheta = 0
			dzmnsdxdzeta = 0
			dzmnsdydtheta = 0
			dzmnsdydzeta = 0
			dzmnsdzdtheta = 0

			do izeta_coil = 1,nzeta_coil

				do itheta_coil = 1,ntheta_coil

					index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil

					do l_coil = 0, (nfp-1)

						izetal_coil = izeta_coil + l_coil*nzeta_coil

            angle2 = zetal_coil(izetal_coil)
            sinangle2 = sin(angle2)
            cosangle2 = cos(angle2)

            dzmncdzdzeta = -sinangle2
            dzmnsdzdzeta = cosangle2

						dxdtheta = drdtheta_coil(1,itheta_coil,izetal_coil)
						dxdzeta = drdzeta_coil(1,itheta_coil,izetal_coil)
						dydzeta = drdzeta_coil(2,itheta_coil,izetal_coil)
						dydtheta = drdtheta_coil(2,itheta_coil,izetal_coil)
						dzdtheta = drdtheta_coil(3,itheta_coil,izetal_coil)
						dzdzeta = drdzeta_coil(3,itheta_coil,izetal_coil)

						do imn = 1,mnmax_sensitivity

							! For Fourier decomposition of surface, need to index using izetal_coil
							angle = xm_sensitivity(imn)*theta_coil(itheta_coil) - xn_sensitivity(imn)*zetal_coil(izetal_coil)
							sinangle = sin(angle)
							cosangle = cos(angle)

							drmncdxdtheta = -xm_sensitivity(imn)*sinangle*cosangle2
							drmncdxdzeta = xn_sensitivity(imn)*sinangle*cosangle2 - cosangle*sinangle2
							drmncdydtheta = -xm_sensitivity(imn)*sinangle*sinangle2
							drmncdydzeta = xn_sensitivity(imn)*sinangle*sinangle2 + cosangle*cosangle2

							drmnsdxdtheta = xm_sensitivity(imn)*cosangle*cosangle2
							drmnsdxdzeta = -xn_sensitivity(imn)*cosangle*cosangle2 - sinangle*sinangle2!
							drmnsdydtheta = xm_sensitivity(imn)*cosangle*sinangle2
							drmnsdydzeta = -xn_sensitivity(imn)*cosangle*sinangle2 + sinangle*cosangle2

							! dnormd(...) stored for use in chi_k^2 sensitivity
							dnormdrmnc(1, imn, itheta_coil, izetal_coil) = drmncdydzeta*dzdtheta + drmncdzdtheta*dydzeta &
								- drmncdydtheta*dzdzeta - drmncdzdzeta*dydtheta
							dnormdrmnc(2, imn, itheta_coil, izetal_coil) = drmncdzdzeta*dxdtheta + drmncdxdtheta*dzdzeta &
								- drmncdzdtheta*dxdzeta - drmncdxdzeta*dzdtheta
							dnormdrmnc(3, imn, itheta_coil, izetal_coil) = drmncdxdzeta*dydtheta + drmncdydtheta*dxdzeta &
								- drmncdxdtheta*dydzeta - drmncdydzeta*dxdtheta

							dnormdrmns(1, imn, itheta_coil, izetal_coil) = drmnsdydzeta*dzdtheta + drmnsdzdtheta*dydzeta &
                - drmnsdydtheta*dzdzeta - drmnsdzdzeta*dydtheta
							dnormdrmns(2, imn, itheta_coil, izetal_coil) = drmnsdzdzeta*dxdtheta + drmnsdxdtheta*dzdzeta &
								- drmnsdzdtheta*dxdzeta - drmnsdxdzeta*dzdtheta
							dnormdrmns(3, imn, itheta_coil, izetal_coil) = drmnsdxdzeta*dydtheta + drmnsdydtheta*dxdzeta &
								- drmnsdxdtheta*dydzeta - drmnsdydzeta*dxdtheta
							
							dnormdzmnc(1, imn, itheta_coil, izetal_coil) = dzmncdydzeta*dzdtheta + dzmncdzdtheta*dydzeta &
								- dzmncdydtheta*dzdzeta - dzmncdzdzeta*dydtheta
							dnormdzmnc(2, imn, itheta_coil, izetal_coil) = dzmncdzdzeta*dxdtheta + dzmncdxdtheta*dzdzeta &
								- dzmncdzdtheta*dxdzeta - dzmncdxdzeta*dzdtheta
							dnormdzmnc(3, imn, itheta_coil, izetal_coil) = dzmncdxdzeta*dydtheta + dzmncdydtheta*dxdzeta &
								- dzmncdxdtheta*dydzeta - dzmncdydzeta*dxdtheta

							dnormdzmns(1, imn, itheta_coil, izetal_coil) = dzmnsdydzeta*dzdtheta + dzmnsdzdtheta*dydzeta &
								- dzmnsdydtheta*dzdzeta - dzmnsdzdzeta*dydtheta
							dnormdzmns(2, imn, itheta_coil, izetal_coil) = dzmnsdzdzeta*dxdtheta + dzmnsdxdtheta*dzdzeta &
								- dzmnsdzdtheta*dxdzeta - dzmnsdxdzeta*dzdtheta
							dnormdzmns(3, imn, itheta_coil, izetal_coil) = dzmnsdxdzeta*dydtheta + dzmnsdydtheta*dxdzeta &
								- dzmnsdxdtheta*dydzeta - dzmnsdydzeta*dxdtheta

              ! dxdOmega
							drdrmnc(1) = cosangle*cosangle2
							drdrmns(1) = sinangle*cosangle2
							drdzmnc(1) = 0
							drdzmns(1) = 0
              ! dydOmega
							drdrmnc(2) = cosangle*sinangle2
							drdrmns(2) = sinangle*sinangle2
							drdzmnc(2) = 0
							drdzmns(2) = 0
              ! dzdOmega
							drdrmnc(3) = 0
							drdrmns(3) = 0
							drdzmnc(3) = cosangle2
							drdzmns(3) = sinangle2

							do izeta_plasma = 1, nzeta_plasma
								do itheta_plasma = 1, ntheta_plasma
									index_plasma = (izeta_plasma-1)*ntheta_plasma + itheta_plasma

									x = r_plasma(1,itheta_plasma,izeta_plasma)
									y = r_plasma(2,itheta_plasma,izeta_plasma)
									z = r_plasma(3,itheta_plasma,izeta_plasma)

                  dx = x - r_coil(1,itheta_coil,izetal_coil)
                  dy = y - r_coil(2,itheta_coil,izetal_coil)
                  dz = z - r_coil(3,itheta_coil,izetal_coil)

                  dr2inv = 1/(dx*dx + dy*dy + dz*dz)
									dr32inv = dr2inv*sqrt(dr2inv)
									dr52inv = dr2inv*dr32inv

									dr_dot_norm_coil = dx*normal_coil(1,itheta_coil,izetal_coil) &
                    + dy*normal_coil(2,itheta_coil,izetal_coil) + dz*normal_coil(3,itheta_coil,izetal_coil)
									dr_dot_norm_plasma = dx*normal_plasma(1,itheta_plasma,izeta_plasma) &
                    + dy*normal_plasma(2,itheta_plasma,izeta_plasma) + dz*normal_plasma(3, itheta_coil,izeta_coil)
									norm_plasma_dot_norm_coil = normal_plasma(1,itheta_plasma,izeta_plasma)*normal_coil(1,itheta_coil, izetal_coil) &
										+ normal_plasma(2,itheta_plasma,izeta_plasma)*normal_coil(2,itheta_coil,izetal_coil) &
										+ normal_plasma(3,itheta_plasma,izeta_plasma)*normal_coil(3,itheta_coil,izetal_coil)

                  ! dgdr before sum over l_coil
									dinductancedr(1) = (3*dx*norm_plasma_dot_norm_coil &
                    - 15*dx*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
										+ 3*(normal_plasma(1,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
										+ normal_coil(1,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
                  dinductancedr(2) = (3*dy*norm_plasma_dot_norm_coil &
                    - 15*dy*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                    + 3*(normal_plasma(2,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                    + normal_coil(2,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))
                  dinductancedr(3) = (3*dz*norm_plasma_dot_norm_coil &
                    - 15*dz*dr2inv*dr_dot_norm_coil*dr_dot_norm_plasma &
                    + 3*(normal_plasma(3,itheta_plasma,izeta_plasma)*dr_dot_norm_coil &
                    + normal_coil(3,itheta_coil,izetal_coil)*dr_dot_norm_plasma))*(dr52inv*mu0/(4*pi))

                  ! dgdN' before sum over l_coil
									dinductancednorm(1) = (normal_plasma(1,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dx)*(dr32inv*mu0/(4*pi))
									dinductancednorm(2) = (normal_plasma(2,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dy)*(dr32inv*mu0/(4*pi))
									dinductancednorm(3) = (normal_plasma(3,itheta_plasma,izeta_plasma) &
                    - 3*dr2inv*dr_dot_norm_plasma*dz)*(dr32inv*mu0/(4*pi))

                  ! These quantity are summed over l_coil
									dinductancedrmnc(imn, index_plasma,index_coil) = dinductancedrmnc(imn, index_plasma, index_coil) &
                    + dinductancednorm(1)*dnormdrmnc(1,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(2)*dnormdrmnc(2,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(3)*dnormdrmnc(3,imn,itheta_coil,izetal_coil) &
                    + dinductancedr(1)*drdrmnc(1) + dinductancedr(2)*drdrmnc(2) &
                    + dinductancedr(3)*drdrmnc(3)
									dinductancedrmns(imn, index_plasma,index_coil) = dinductancedrmns(imn, index_plasma,index_coil) &
                    + dinductancednorm(1)*dnormdrmns(1,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(2)*dnormdrmns(2,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(3)*dnormdrmns(3,imn,itheta_coil,izetal_coil) &
                    + dinductancedr(1)*drdrmns(1) + dinductancedr(2)*drdrmns(2) &
                    + dinductancedr(3)*drdrmns(3)
									dinductancedzmnc(imn, index_plasma,index_coil) = dinductancedzmnc(imn, index_plasma,index_coil) &
                    + dinductancednorm(1)*dnormdzmnc(1,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(2)*dnormdzmnc(2,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(3)*dnormdzmnc(3,imn,itheta_coil,izetal_coil) &
                    + dinductancedr(1)*drdzmnc(1) + dinductancedr(2)*drdzmnc(2) &
                    + dinductancedr(3)*drdzmnc(3)
									dinductancedzmns(imn, index_plasma,index_coil) = dinductancedzmns(imn, index_plasma,index_coil) &
                    + dinductancednorm(1)*dnormdzmns(1,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(2)*dnormdzmns(2,imn,itheta_coil,izetal_coil) &
                    + dinductancednorm(3)*dnormdzmns(3,imn,itheta_coil,izetal_coil) &
                    + dinductancedr(1)*drdzmns(1) + dinductancedr(2)*drdzmns(2) &
                    + dinductancedr(3)*drdzmns(3)
								end do
							end do
						end do
					end do
				end do
			end do

      ! Now need to multiply dinductanced(...) by basis_functions to compute sensitivity of g_j
      ! dincutanced(...)(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil)
      ! basis_functions(ntheta_coil*nzeta_coil, num_basis_functions)
      ! dgdrmnc(mnmax_sensitivity, ntheta_plasma*nzeta_plasma, num_basis_functions)

      ! Now we need to compute dgd(...) = dinductanced(...)*basis_functions
      M = ntheta_plasma*nzeta_plasma ! # rows of A
      N = num_basis_functions ! # cols of B
      K = ntheta_coil*nzeta_coil ! Common dimension of A and B
      LDA = M
      LDB = K
      LDC = M
      TRANSA = 'N' ! No transposes
      TRANSB = 'N'
      BLAS_ALPHA=dtheta_coil*dzeta_coil
      BLAS_BETA=0
      dgdrmnc = 0
      dgdrmns = 0
      dgdzmnc = 0
      dgdzmns = 0

      do imn = 1, mnmax_sensitivity
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedrmnc(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdrmnc(imn,:,:),LDC)
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedrmns(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdrmns(imn,:,:),LDC)
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedzmnc(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdzmnc(imn,:,:),LDC)
        call DGEMM(TRANSA,TRANSB,M,N,K,BLAS_ALPHA,dinductancedzmns(imn,:,:),LDA,basis_functions,LDB,BLAS_BETA,dgdzmns(imn,:,:),LDC)
      enddo

      ! debugging
      print *,"  Max dgdrmnc: ",maxval(dgdrmnc)

			! Computing quantities needed chi2_K sensitivity
      ! Partial derivatives of d, f, and N' computed

      ! Needed for computing sensitivity of f
      select case (symmetry_option)
      case (1)
        minSymmetry = 1
        maxSymmetry = 1
      case (2)
        minSymmetry = 2
        maxSymmetry = 2
      case (3)
        minSymmetry = 1
        maxSymmetry = 2
      end select

			do izeta_coil = 1, nzeta_coil

        angle2 = zeta_coil(izeta_coil)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)

        dzmncdzdzeta = -sinangle2
        dzmnsdzdzeta = cosangle2

				do itheta_coil = 1, ntheta_coil
          index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil
					do imn = 1, mnmax_sensitivity

						! We are computing a quantity symmetric in toroidal periods, but surface has Fourier components
            ! in izetal_coil
						angle = xm_sensitivity(imn)*theta_coil(itheta_coil) - xn_sensitivity(imn)*zetal_coil(izeta_coil)
						sinangle = sin(angle)
						cosangle = cos(angle)

						!! needs imn, itheta_coil, izeta_coil
						drmncdxdtheta = -xm_sensitivity(imn)*sinangle*cosangle2
						drmncdxdzeta = xn_sensitivity(imn)*sinangle*cosangle2 - cosangle*sinangle2
						drmncdydtheta = -xm_sensitivity(imn)*sinangle*sinangle2
						drmncdydzeta = xn_sensitivity(imn)*sinangle*sinangle2 + cosangle*cosangle2

						drmnsdxdtheta = xm_sensitivity(imn)*cosangle*cosangle2
						drmnsdxdzeta = -xn_sensitivity(imn)*cosangle*cosangle2 - sinangle*sinangle2
						drmnsdydtheta = xm_sensitivity(imn)*cosangle*sinangle2
						drmnsdydzeta = -xn_sensitivity(imn)*cosangle*sinangle2 + sinangle*cosangle2

						! Here the sensitivity is l_coil periodic - indexed by izeta_coil
						dnorm_normaldrmnc(imn,itheta_coil,izeta_coil) = &
               (normal_coil(1, itheta_coil, izeta_coil)*dnormdrmnc(1, imn, itheta_coil, izeta_coil) &
							+ normal_coil(2, itheta_coil, izeta_coil)*dnormdrmnc(2, imn, itheta_coil, izeta_coil) &
							+ normal_coil(3, itheta_coil, izeta_coil)*dnormdrmnc(3, imn, itheta_coil, izeta_coil)) &
                /norm_normal_coil(itheta_coil,izeta_coil)
						dnorm_normaldrmns(imn, itheta_coil,izeta_coil) = &
               (normal_coil(1, itheta_coil, izeta_coil)*dnormdrmns(1, imn, itheta_coil, izeta_coil) &
							+ normal_coil(2, itheta_coil, izeta_coil)*dnormdrmns(2, imn, itheta_coil, izeta_coil) &
							+ normal_coil(3, itheta_coil, izeta_coil)*dnormdrmns(3, imn, itheta_coil, izeta_coil)) &
                /norm_normal_coil(itheta_coil,izeta_coil)
						dnorm_normaldzmnc(imn, itheta_coil,izeta_coil) = &
               (normal_coil(1, itheta_coil, izeta_coil)*dnormdzmnc(1, imn, itheta_coil, izeta_coil) &
							+ normal_coil(2, itheta_coil, izeta_coil)*dnormdzmnc(2, imn, itheta_coil, izeta_coil) &
							+ normal_coil(3, itheta_coil, izeta_coil)*dnormdzmnc(3, imn, itheta_coil, izeta_coil)) &
               /norm_normal_coil(itheta_coil,izeta_coil)
						dnorm_normaldzmns(imn, itheta_coil,izeta_coil) = &
               (normal_coil(1, itheta_coil, izeta_coil)*dnormdzmns(1, imn, itheta_coil, izeta_coil) &
							+ normal_coil(2, itheta_coil, izeta_coil)*dnormdzmns(2, imn, itheta_coil, izeta_coil) &
							+ normal_coil(3, itheta_coil, izeta_coil)*dnormdzmns(3, imn, itheta_coil, izeta_coil)) &
               /norm_normal_coil(itheta_coil,izeta_coil)

						dddrmnc(1, imn, index_coil) = (net_poloidal_current_Amperes*drmncdxdtheta &
              - net_toroidal_current_Amperes*drmncdxdzeta)/twopi
						dddrmnc(2, imn, index_coil) = (net_poloidal_current_Amperes*drmncdydtheta &
              - net_toroidal_current_Amperes*drmncdydzeta)/twopi
						dddrmnc(3, imn, index_coil) = (net_poloidal_current_Amperes*drmncdzdtheta &
              - net_toroidal_current_Amperes*drmncdzdzeta)/twopi
						dddrmns(1, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdxdtheta &
              - net_toroidal_current_Amperes*drmnsdxdzeta)/twopi
						dddrmns(2, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdydtheta &
              - net_toroidal_current_Amperes*drmnsdydzeta)/twopi
						dddrmns(3, imn, index_coil) = (net_poloidal_current_Amperes*drmnsdzdtheta &
              - net_toroidal_current_Amperes*drmnsdzdzeta)/twopi
						dddzmnc(1, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdxdtheta &
              - net_toroidal_current_Amperes*dzmncdxdzeta)/twopi
						dddzmnc(2, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdydtheta &
              - net_toroidal_current_Amperes*dzmncdydzeta)/twopi
						dddzmnc(3, imn, index_coil) = (net_poloidal_current_Amperes*dzmncdzdtheta &
              - net_toroidal_current_Amperes*dzmncdzdzeta)/twopi
						dddzmns(1, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdxdtheta &
              - net_toroidal_current_Amperes*dzmnsdxdzeta)/twopi
						dddzmns(2, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdydtheta &
              - net_toroidal_current_Amperes*dzmnsdydzeta)/twopi
						dddzmns(3, imn, index_coil) = (net_poloidal_current_Amperes*dzmnsdzdtheta &
              - net_toroidal_current_Amperes*dzmnsdzdzeta)/twopi

            do whichSymmetry = minSymmetry, maxSymmetry

              if (whichSymmetry==2 .and. symmetry_option==3) then
                offset = mnmax_coil
                else
                  offset = 0
              end if

              do imn_coil = 1, mnmax_coil

                angle_coil = xm_coil(imn_coil)*theta_coil(itheta_coil) - xn_coil(imn_coil)*zetal_coil(izeta_coil)
                cosangle_coil = cos(angle_coil)
                sinangle_coil = sin(angle_coil)

                if (whichSymmetry==1) then
                  dfdrmnc(1, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdxdtheta &
                    + xm_coil(imn_coil)*drmncdxdzeta)
                  dfdrmns(1, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdxdtheta &
                    + xm_coil(imn_coil)*drmnsdxdzeta)
                  dfdzmnc(1, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdxdtheta &
                    + xm_coil(imn_coil)*dzmncdxdzeta)
                  dfdzmns(1, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdxdtheta &
                    + xm_coil(imn_coil)*dzmnsdxdzeta)

                  dfdrmnc(2, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdydtheta &
                    + xm_coil(imn_coil)*drmncdydzeta)
                  dfdrmns(2, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdydtheta &
                    + xm_coil(imn_coil)*drmnsdydzeta)
                  dfdzmnc(2, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdydtheta &
                    + xm_coil(imn_coil)*dzmncdydzeta)
                  dfdzmns(2, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdydtheta &
                    + xm_coil(imn_coil)*dzmnsdydzeta)

                  dfdrmnc(3, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmncdzdtheta &
                    + xm_coil(imn_coil)*drmncdzdzeta)
                  dfdrmns(3, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*drmnsdzdtheta &
                    + xm_coil(imn_coil)*drmnsdzdzeta)
                  dfdzmnc(3, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmncdzdtheta &
                    + xm_coil(imn_coil)*dzmncdzdzeta)
                  dfdzmns(3, imn, index_coil,imn_coil) = cosangle_coil*(xn_coil(imn_coil)*dzmnsdzdtheta &
                    + xm_coil(imn_coil)*dzmnsdzdzeta)

                else
                  dfdrmnc(1, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmncdxdtheta &
                    + xm_coil(imn_coil)*drmncdxdzeta)
                  dfdrmns(1, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdxdtheta &
                    + xm_coil(imn_coil)*drmnsdxdzeta)
                  dfdzmnc(1, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdxdtheta &
                    + xm_coil(imn_coil)*dzmncdxdzeta)
                  dfdzmns(1, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdxdtheta &
                    + xm_coil(imn_coil)*dzmnsdxdzeta)

                  dfdrmnc(2, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmncdydtheta &
                    + xm_coil(imn_coil)*drmncdydzeta)
                  dfdrmns(2, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdydtheta &
                    + xm_coil(imn_coil)*drmnsdydzeta)
                  dfdzmnc(2, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdydtheta &
                    + xm_coil(imn_coil)*dzmncdydzeta)
                  dfdzmns(2, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdydtheta &
                    + xm_coil(imn_coil)*dzmnsdydzeta)

                  dfdrmnc(3, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmncdzdtheta &
                    + xm_coil(imn_coil)*drmncdzdzeta)
                  dfdrmns(3, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*drmnsdzdtheta &
                    + xm_coil(imn_coil)*drmnsdzdzeta)
                  dfdzmnc(3, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmncdzdtheta &
                    + xm_coil(imn_coil)*dzmncdzdzeta)
                  dfdzmns(3, imn, index_coil,imn_coil) = -sinangle_coil*(xn_coil(imn_coil)*dzmnsdzdtheta &
                    + xm_coil(imn_coil)*dzmnsdzdzeta)
                end if
              enddo
            enddo
					enddo
				enddo
			enddo

		end subroutine init_partials

	end module init_sensitivity


