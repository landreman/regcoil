subroutine init_sensitivity_adjoint()

  use global_variables
  use stel_constants

  implicit none

  integer :: iomega, j_basis, tic, toc, countrate, iflag
  real(dp), dimension(:,:), allocatable :: Afactorx, Afactory, Afactorz, Afactor
  real(dp), dimension(:,:), allocatable :: bfactorx, bfactory, bfactorz

  real(dp), dimension(:), allocatable :: norm_normal_coil_inv1D
  real(dp), dimension(:), allocatable :: norm_normal_plasma_inv1D
  real(dp), dimension(:,:), allocatable :: dnorm_normaldomega2D

  allocate(dAKdomega(num_basis_functions,num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dABdomega(num_basis_functions,num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dbKdomega(num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dbBdomega(num_basis_functions,nomega_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactorx(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactory(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactorz(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Afactor(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactorx(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactory(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(bfactorz(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_normal_coil_inv1D(ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(norm_normal_plasma_inv1D(ntheta_plasma*nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(dnorm_normaldomega2D(nomega_coil,ntheta_coil*nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  if (sensitivity_option == 3) then
    call system_clock(tic,countrate)
  norm_normal_coil_inv1D   = reshape(1/norm_normal_coil,   (/ ntheta_coil*nzeta_coil /))
  do iomega = 1, nomega_coil
    ! dfxdomega(nomega_coil, ntheta_coil*nzeta_coil, num_basis_functions)
    ! f_x(ntheta_coil*nzeta_coil, num_basis_functions)
    dAKdomega(:,:,iomega) = dtheta_coil*dzeta_coil*(2*(matmul(transpose(f_x), dfxdomega(iomega,:,:)) &
    + matmul(transpose(f_y), dfydomega(iomega,:,:)) + matmul(transpose(f_z),dfzdomega(iomega,:,:))))
    dnorm_normaldomega2D(iomega,:) = reshape(dnorm_normaldomega(iomega,:,:),(/ ntheta_coil * nzeta_coil /))
    do j_basis = 1, num_basis_functions
      Afactorx(:, j_basis) = f_x(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      Afactory(:, j_basis) = f_y(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      Afactorz(:, j_basis) = f_z(:,j_basis)*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D
      bfactorx(:,j_basis) = f_x(:,j_basis)*norm_normal_coil_inv1D
      bfactory(:,j_basis) = f_y(:,j_basis)*norm_normal_coil_inv1D
      bfactorz(:,j_basis) = f_z(:,j_basis)*norm_normal_coil_inv1D
    enddo
    dAKdomega(:,:,iomega) = dAKdomega(:,:,iomega) - dtheta_coil*dzeta_coil*(matmul(transpose(f_x),Afactorx) &
      + matmul(transpose(f_y),Afactory) + matmul(transpose(f_z), Afactorz))
    dbKdomega(:,iomega) = dtheta_coil*dzeta_coil*(matmul(transpose(bfactorx),dddomega(1,iomega,:)) + &
      matmul(transpose(bfactory),dddomega(2,iomega,:)) + matmul(transpose(bfactorz),dddomega(3,iomega,:)))
    dbKdomega(:,iomega) = dbKdomega(:,iomega) + dtheta_coil*dzeta_coil*(matmul(transpose(dfxdomega(iomega,:,:)), &
      d_x*norm_normal_coil_inv1D) + matmul(transpose(dfydomega(iomega,:,:)),d_y*norm_normal_coil_inv1D) &
      + matmul(transpose(dfzdomega(iomega,:,:)),d_z*norm_normal_coil_inv1D))
      dbKdomega(:,iomega) = dbKdomega(:,iomega) - dtheta_coil*dzeta_coil*(matmul(transpose(f_x), &
      d_x*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D) + &
      matmul(transpose(f_y),d_y*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D) &
      + matmul(transpose(f_z),d_z*dnorm_normaldomega2D(iomega,:)*norm_normal_coil_inv1D*norm_normal_coil_inv1D))
  enddo

  call system_clock(toc)
  print *,"Computing matrices for dFKdomega: ",real(toc-tic)/countrate," sec."

  call system_clock(tic,countrate)

  norm_normal_plasma_inv1D   = reshape(1/norm_normal_plasma,   (/ ntheta_plasma*nzeta_plasma /))
  do iomega = 1, nomega_coil
    do j_basis = 1, num_basis_functions
      Afactor(:,j_basis) = dgdomega(iomega,:,j_basis)*norm_normal_plasma_inv1D
    enddo
    dABdomega(:,:,iomega) = 2*dtheta_plasma*dzeta_plasma*(matmul(transpose(g),Afactor))
    dbBdomega(:,iomega) = -dtheta_plasma*dzeta_plasma*matmul(transpose(dgdomega(iomega,:,:)), &
    reshape(Bnormal_from_plasma_current+Bnormal_from_net_coil_currents, (/ ntheta_plasma*nzeta_plasma /)))
  enddo

  call system_clock(toc)
  print *,"Computing matrices for dFBdomega: ",real(toc-tic)/countrate," sec."
  endif

  if (sensitivity_option == 3) then
    deallocate(Afactorx)
    deallocate(Afactory)
    deallocate(Afactorz)
    deallocate(bfactorx)
    deallocate(bfactory)
    deallocate(bfactorz)
    deallocate(Afactor)
    deallocate(norm_normal_coil_inv1D)
    deallocate(norm_normal_plasma_inv1D)
    deallocate(dnorm_normaldomega2D)
  endif

end subroutine init_sensitivity_adjoint
