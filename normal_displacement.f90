subroutine normal_displacement()

  use global_variables
  use stel_constants

  real(dp), dimension(:,:), allocatable :: dchi2dx, dchi2dy, dchi2dz
  real(dp), dimension(:), allocatable :: this_dchi2dx, this_dchi2dy, this_dchi2dz
  integer :: iomega, iflag

  allocate(dchi2dx(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dy(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dz(ntheta_coil*nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dr_normal(nlambda,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_dchi2dx(nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_dchi2dy(nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_dchi2dz(nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  dchi2dx = 0
  dchi2dy = 0
  dchi2dz = 0
  dchi2dr_normal = 0
  this_dchi2dx = 0
  this_dchi2dy = 0
  this_dchi2dz = 0

  do izeta_coil = 1,nzeta_coil
    do itheta_coil = 1,ntheta_coil
      index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil

      angle2 = zeta_coil(izeta_coil)
      sinangle2 = sin(angle2)
      cosangle2 = cos(angle2)

      do iomega = 1,nomega_coil
        angle = xm_sensitivity(iomega)*theta_coil(itheta_coil) + nfp*xn_sensitivity(iomega)*zeta_coil(izeta_coil)
        sinangle = sin(angle)
        cosangle = cos(angle)

        if (omega_coil(iomega) == 1) then ! omega = rmnc
          this_dchi2dx = dchi2domega(iomega,:)*cosangle*cosangle2/(pi*pi)
          this_dchi2dy = dchi2domega(iomega,:)*cosangle*sinangle2/(pi*pi)
        else if (omega_coil(iomega) == 2) then ! omega = zmns
          this_dchi2dz = dchi2domega(iomega,:)*sinangle/(2*pi*pi)
        else if (omega_coil(iomega) == 3) then ! omega = rmns
          this_dchi2dx = dchi2domega(iomega,:)*sinangle*cosangle2/(pi*pi)
          this_dchi2dy = dchi2domega(iomega,:)*sinangle*sinangle2/(pi*pi)
        else if (omega_coil(iomega) == 4) then ! omega = zmnc
          this_dchi2dz = dchi2domega(iomega,:)*cosangle/(2*pi*pi)
        else
          print *,"omega_coil out of range!\n"
          call exit
        endif

        if (xm_sensitivity(iomega) == 0 .and. xn_sensitivity(iomega) == 0) then
          this_dchi2dx = this_dchi2dx/2
          this_dchi2dy = this_dchi2dy/2
          this_dchi2dz = this_dchi2dz/2
        endif

        ! Check for NaNs
        do ilambda = 1,nlambda
          if (isnan(this_dchi2dx(ilambda)) .or. isnan(this_dchi2dy(ilambda)) .or. isnan(this_dchi2dz(ilambda))) then
            print *,"lambda = %d\n", lambda(ilambda)
            print *,"NaN detected in normal displacement!\n"
            call exit
          endif
        enddo

        dchi2dx(index_coil,:) = dchi2dx(index_coil,:) + this_dchi2dx
        dchi2dy(index_coil,:) = dchi2dy(index_coil,:) + this_dchi2dy
        dchi2dz(index_coil,:) = dchi2dz(index_coil,:) + this_dchi2dz

      enddo

      ! Project onto normal direction
      if (norm_normal_coil(itheta_coil,izeta_coil) /= 0) then
        dchi2dr_normal(:,index_coil) = &
          (normal_coil(1,itheta_coil,izeta_coil)*dchi2dx(index_coil,:) &
          + normal_coil(2,itheta_coil,izeta_coil)*dchi2dy(index_coil,:) &
          + normal_coil(3,itheta_coil,izeta_coil)*dchi2dz(index_coil,:)) &
          /norm_normal_coil(itheta_coil,izeta_coil)
      endif

    enddo
  enddo

end subroutine normal_displacement