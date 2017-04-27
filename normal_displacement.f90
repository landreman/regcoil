subroutine normal_displacement()

  use global_variables

  implicit none

  integer :: iflag, iomega, ilambda

  real(dp), dimension(:,:,:), allocatable :: dchi2dx,dchi2dy,dchi2dz

  allocate(dchi2dr_normal(ntheta_coil,nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dx(ntheta_coil,nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dy(ntheta_coil,nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dchi2dz(ntheta_coil,nzeta_coil,nlambda),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do iomega = 1, nomega_coil
    do ilambda = 1, nlambda
      dchi2dx(:,:,ilambda) = dchi2domega(iomega,ilambda)/drdomega(1,:,:,iomega)
      dchi2dy(:,:,ilambda) = dchi2domega(iomega,ilambda)/drdomega(2,:,:,iomega)
      dchi2dz(:,:,ilambda) = dchi2domega(iomega,ilambda)/drdomega(3,:,:,iomega)
    enddo
  enddo

  do ilambda = 1, nlambda
    dchi2dr_normal(:,:,ilambda) = (normal_coil(1,:,:)*dchi2dx(:,:,ilambda) &
      + normal_coil(2,:,:)*dchi2dy(:,:,ilambda) &
      + normal_coil(3,:,:)*dchi2dz(:,:,ilambda))/norm_normal_coil
  enddo

  deallocate(dchi2dx)
  deallocate(dchi2dy)
  deallocate(dchi2dz)

end subroutine normal_displacement
