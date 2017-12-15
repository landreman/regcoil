subroutine save_nescin()

  use global_variables
  use init_Fourier_modes_mod
  use stel_constants

  implicit none

  integer :: save_nescin_mnmax, imn, m, n, itheta, izeta
  integer :: iflag, mi, ni, ntheta_transform, nzetal_transform
  integer, dimension(:), allocatable :: save_nescin_xn, save_nescin_xm
  real(dp), dimension(:), allocatable :: save_nescin_rmnc, save_nescin_zmns
  real(dp), dimension(:,:), allocatable :: R_major_coil, z_coil
  real(dp) :: angle, dtransform
  real(dp), dimension(:), allocatable :: theta_transform, zetal_transform
  character(13) :: nescin_outname = "nescin.offset"
  integer, parameter :: out_unit=20
  character(1) tab
  tab = char(9)

  save_nescin_mnmax = (save_nescin_ntor+1) + (save_nescin_mpol)*(2*save_nescin_ntor+1)

  ! Eventually this should be changed - increase res. and interpolate w/ splines
  ntheta_transform = ntheta_coil
  nzetal_transform = nzetal_coil

  allocate(save_nescin_rmnc(save_nescin_mnmax),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(save_nescin_zmns(save_nescin_mnmax),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(save_nescin_xn(save_nescin_mnmax),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(save_nescin_xm(save_nescin_mnmax),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(R_major_coil(ntheta_transform,nzetal_transform),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(z_coil(ntheta_transform,nzetal_transform),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(theta_transform(ntheta_transform),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(zetal_transform(nzetal_transform),stat=iflag)

  do itheta=1,ntheta_transform
    theta_transform(itheta) = 2*pi*(itheta-1.0_dp)/ntheta_transform
  end do

  do izeta = 1,nzetal_transform
    zetal_transform(izeta) = 2*pi*(izeta-1.0_dp)/nzetal_transform
  end do

  ! Initialize m = 0 modes
  imn = 1
  do ni = 0,save_nescin_ntor
    save_nescin_xm(imn) = 0
    save_nescin_xn(imn) = ni
    imn = imn + 1
  end do
  ! Initialize m > 0 modes
  do mi = 1, save_nescin_mpol
    do ni = -save_nescin_ntor, save_nescin_ntor
      save_nescin_xm(imn) = mi
      save_nescin_xn(imn) = ni
      imn = imn + 1
    end do
  end do

  R_major_coil = sqrt(r_coil(1,:,:)**2 + r_coil(2,:,:)**2)
  z_coil = r_coil(3,:,:)

  dtransform = (1.0_dp)/(ntheta_transform*nzetal_transform)
  save_nescin_rmnc = 0
  save_nescin_zmns = 0
  do imn = 1,save_nescin_mnmax
    m = save_nescin_xm(imn)
    n = save_nescin_xn(imn)
    do itheta = 1,ntheta_transform
      do izeta = 1,nzetal_transform
        angle = m*theta_transform(itheta) + nfp*n*zetal_transform(izeta)
        save_nescin_rmnc(imn) = save_nescin_rmnc(imn) + R_major_coil(itheta,izeta)*cos(angle)
        save_nescin_zmns(imn) = save_nescin_zmns(imn) + z_coil(itheta,izeta)*sin(angle)
      end do
    end do
    if (m /= 0 .or. n /= 0) then
      save_nescin_rmnc(imn) = 2*save_nescin_rmnc(imn)
      save_nescin_zmns(imn) = 2*save_nescin_zmns(imn)
    end if
  end do
  save_nescin_rmnc = save_nescin_rmnc*dtransform
  save_nescin_zmns = save_nescin_zmns*dtransform

  ! Write to nescin file 
  open(unit=out_unit,file=nescin_outname,action="write",status="replace")
  write(out_unit,'(A,F8.3)') "------ Current Surface"
  write(out_unit,*) ""
  write(out_unit,*) save_nescin_mnmax
  write(out_unit,*) ""
  write(out_unit,*) ""
  do imn=1,save_nescin_mnmax
    write(out_unit,"(I4,I4,ES10.2,ES10.2,ES10.2,ES10.2,ES10.2)") save_nescin_xm(imn),save_nescin_xn(imn),save_nescin_rmnc(imn),save_nescin_zmns(imn),0.0,0.0
  end do
  close(out_unit)

end subroutine save_nescin
