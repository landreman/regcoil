subroutine read_nescin(nescin_filename, r, drdtheta, drdzeta, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, ntheta, nzetal, theta, zetal, compute_2nd_derivs)

  use regcoil_variables, only: nfp, xm, xn, mnmax, rmnc_global => rmnc, zmns_global => zmns, rmns_global => rmns, zmnc_global => zmnc
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none
  
  character(*) :: nescin_filename
  integer, intent(in) :: ntheta, nzetal
  real(dp), dimension(3,ntheta,nzetal) :: r, drdtheta, drdzeta
  real(dp), dimension(3,ntheta,nzetal) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2
  real(dp), dimension(ntheta)  :: theta
  real(dp), dimension(nzetal) :: zetal
  logical :: compute_2nd_derivs
  
  integer :: iunit = 7, itheta, izeta, iflag
  integer :: m, n, ntotal, k, mr, nr, istat
  real(dp) :: rmnc, zmns, rmns, zmnc
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
  real(dp) :: d2sinangle2dzeta2, d2cosangle2dzeta2
  real(dp) :: d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2
  real(dp) :: d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2

  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"

  r=0
  drdtheta=0
  drdzeta=0

  if (compute_2nd_derivs) then
     d2rdtheta2 = 0
     d2rdthetadzeta = 0
     d2rdzeta2 = 0
  end if

  call safe_open(iunit, istat, trim(nescin_filename), 'old', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file'
  endif


  ! Skip down to the line in the file that begins with matchString
  do
     read (iunit,"(a)") myline
     if (myline(:len(matchString)) == matchString) then
        exit
     end if
  end do
  read (iunit, *)

  read (iunit, *) ntotal
  print *,"  Reading",ntotal,"Fourier modes from nescin"


!!$  if (geometry_option_outer==4) then
!!$     ! Clear arrays associated with the plasma surface for offsetting,
!!$     ! and replace them with the nescin values.
!!$     deallocate(xm,xn,rmnc_global,zmns_global)
!!$     if (allocated(rmns_global)) then
!!$        deallocate(rmns_global)
!!$     end if
!!$     if (allocated(zmnc_global)) then
!!$        deallocate(zmnc_global)
!!$     end if
!!$     mnmax = ntotal
!!$     allocate(xm(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 1"
!!$     allocate(xn(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 2"
!!$     allocate(rmnc_global(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 3"
!!$     allocate(zmns_global(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 4"
!!$     allocate(rmns_global(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 5"
!!$     allocate(zmnc_global(mnmax),stat=iflag)
!!$     if (iflag .ne. 0) stop "Allocation error! read_nescin 6"
!!$  end if

  read (iunit, *)
  read (iunit, *)
  do k = 1, ntotal
     read (iunit, *) m, n, rmnc, zmns, rmns, zmnc

!!$     if (geometry_option_outer==4) then
!!$        ! Set arrays associated with offsetting surfaces
!!$        xm(k) = m
!!$        xn(k) = -n*nfp
!!$        rmnc_global(k) = rmnc
!!$        zmns_global(k) = zmns
!!$        rmns_global(k) = rmns
!!$        zmnc_global(k) = zmnc
!!$     end if

     do izeta = 1,nzetal
        angle2 = zetal(izeta)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dzeta = cosangle2
        dcosangle2dzeta = -sinangle2
        d2sinangle2dzeta2 = -sinangle2
        d2cosangle2dzeta2 = -cosangle2
        do itheta = 1,ntheta
           angle = m*theta(itheta) + n*nfp*zetal(izeta)
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledtheta = cosangle*m
           dcosangledtheta = -sinangle*m
           dsinangledzeta = cosangle*n*nfp
           dcosangledzeta = -sinangle*n*nfp
           d2sinangledtheta2  = -m*m*sinangle
           d2sinangledthetadzeta = -m*n*nfp*sinangle
           d2sinangledzeta2  = -n*n*nfp*nfp*sinangle
           d2cosangledtheta2  = -m*m*cosangle
           d2cosangledthetadzeta = -m*n*nfp*cosangle
           d2cosangledzeta2  = -n*n*nfp*nfp*cosangle
           
           r(1,itheta,izeta) = r(1,itheta,izeta) + rmnc * cosangle * cosangle2 + rmns * sinangle * cosangle2
           r(2,itheta,izeta) = r(2,itheta,izeta) + rmnc * cosangle * sinangle2 + rmns * sinangle * sinangle2
           r(3,itheta,izeta) = r(3,itheta,izeta) + zmns * sinangle             + zmnc * cosangle
           
           drdtheta(1,itheta,izeta) = drdtheta(1,itheta,izeta) + rmnc * dcosangledtheta * cosangle2 + rmns * dsinangledtheta * cosangle2
           drdtheta(2,itheta,izeta) = drdtheta(2,itheta,izeta) + rmnc * dcosangledtheta * sinangle2 + rmns * dsinangledtheta * sinangle2
           drdtheta(3,itheta,izeta) = drdtheta(3,itheta,izeta) + zmns * dsinangledtheta + zmnc * dcosangledtheta
           
           drdzeta(1,itheta,izeta) = drdzeta(1,itheta,izeta) + rmnc * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta) &
                + rmns * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
           drdzeta(2,itheta,izeta) = drdzeta(2,itheta,izeta) + rmnc * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta) &
                + rmns * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
           drdzeta(3,itheta,izeta) = drdzeta(3,itheta,izeta) + zmns * dsinangledzeta + zmnc * dcosangledzeta

           if (compute_2nd_derivs) then
              d2rdtheta2(1,itheta,izeta) = d2rdtheta2(1,itheta,izeta) + rmnc * d2cosangledtheta2 * cosangle2 + rmns * d2sinangledtheta2 * cosangle2
              d2rdtheta2(2,itheta,izeta) = d2rdtheta2(2,itheta,izeta) + rmnc * d2cosangledtheta2 * sinangle2 + rmns * d2sinangledtheta2 * sinangle2
              d2rdtheta2(3,itheta,izeta) = d2rdtheta2(3,itheta,izeta) + zmns * d2sinangledtheta2 + zmnc * d2cosangledtheta2

              d2rdthetadzeta(1,itheta,izeta) = d2rdthetadzeta(1,itheta,izeta) + rmnc * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta) &
                   + rmns * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
              d2rdthetadzeta(2,itheta,izeta) = d2rdthetadzeta(2,itheta,izeta) + rmnc * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta) &
                   + rmns * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
              d2rdthetadzeta(3,itheta,izeta) = d2rdthetadzeta(3,itheta,izeta) + zmns * d2sinangledthetadzeta + zmnc * d2cosangledthetadzeta
           
              d2rdzeta2(1,itheta,izeta) = d2rdzeta2(1,itheta,izeta) + rmnc * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
                   + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2) &
                   + rmns * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
                   + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
              d2rdzeta2(2,itheta,izeta) = d2rdzeta2(2,itheta,izeta) + rmnc * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
                   + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2) &
                   + rmns * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
                   + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
              d2rdzeta2(3,itheta,izeta) = d2rdzeta2(3,itheta,izeta) + zmns * d2sinangledzeta2 + zmnc * d2cosangledzeta2
           end if
        end do
     end do
  end do
  
  
  close(iunit)

end subroutine  read_nescin


subroutine write_nescin
  ! FFT the offset winding surface and write the harmonics to file; by czhu on 05/17/2018;
  use regcoil_variables, only: R0_plasma, nfp, volume_coil, &
       ntheta_coil, nzeta_coil, nzetal_coil, &
       theta_coil, zeta_coil, zetal_coil, &
       r_coil, drdtheta_coil, drdzeta_coil, &
       normal_coil, norm_normal_coil, area_coil, &
       geometry_option_coil, R0_coil, a_coil,  &
       separation, dtheta_coil, dzeta_coil, &
       nescin_filename, mpol_coil, ntor_coil, &
       nescin_filename, curpol

  use stel_kinds
  use stel_constants
  use safe_open_mod
  use omp_lib

  implicit none


  real(dp) :: angle
  integer :: i, itheta, izeta, im, in, mnmax_ws, istat = 0, iunit = 8
  real(dp) :: r_tmp, rc_tmp, rs_tmp, zc_tmp, zs_tmp, factor, tol

  !mnmax_ws = (mpol_coil+1)*(2*ntor_coil+1) ! (0:mf)*(-nf:nf)
  mnmax_ws = 0
  r_tmp = 0.0_dp
  !factor = dtheta_coil * dzeta_coil / (16*atan(1))**2
  factor = 1.0/(ntheta_coil*nzetal_coil)
  tol = separation * 1.0E-3_dp ! set the tolarence to ignore harmonics

  call safe_open(iunit, istat, trim(nescin_filename), 'replace', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file: file exsited or something wrong.'
  endif
  print *, "writing winding surface into ", nescin_filename

  ! to avoid repeated computation, we need to save the harmnics into arrays.
  do im = 0, mpol_coil
     do in = -ntor_coil, ntor_coil

        rc_tmp = 0.0_dp ; rs_tmp = 0.0_dp
        zc_tmp = 0.0_dp ; zs_tmp = 0.0_dp

        do itheta = 1,ntheta_coil
           do izeta = 1,nzetal_coil ! nzeta
              angle = im*theta_coil(itheta) + in*zetal_coil(izeta)  ! mu+nv
              r_tmp = sqrt(r_coil(1,itheta,izeta)*r_coil(1,itheta,izeta) + r_coil(2,itheta,izeta)*r_coil(2,itheta,izeta))
              rc_tmp = rc_tmp + r_tmp*cos(angle)
              rs_tmp = rs_tmp + r_tmp*sin(angle)
              zc_tmp = zc_tmp + r_coil(3,itheta,izeta)*cos(angle)
              zs_tmp = zs_tmp + r_coil(3,itheta,izeta)*sin(angle)
           enddo
        enddo

        if (im .eq. 0 .and. in .gt. 0) cycle
        if (im .ne. 0 .or. in .ne. 0) then
           rc_tmp = rc_tmp * 2.0_dp
           rs_tmp = rs_tmp * 2.0_dp
           zc_tmp = zc_tmp * 2.0_dp
           zs_tmp = zs_tmp * 2.0_dp
        endif
        rc_tmp = rc_tmp * factor
        rs_tmp = rs_tmp * factor
        zc_tmp = zc_tmp * factor
        zs_tmp = zs_tmp * factor

        if ( (abs(rc_tmp) + abs(rs_tmp) + abs(zc_tmp)+ abs(zs_tmp)) .ge. tol ) & ! count nonzero numbers
             mnmax_ws = mnmax_ws + 1
        
     enddo
  enddo


  write (iunit, '(a)') '------ Plasma information from VMEC ----'
  write (iunit, '(a)') 'np     iota_edge       phip_edge       curpol'
  write (iunit, '(I6, 3ES20.12)') nfp, 0.0, 0.0, curpol  ! write nfp and curpol information 

  write (iunit,*)
  write (iunit, '(a, 1pe20.12, a)') '------ Current Surface: Coil-Plasma separation = ', separation,' -----'
  write (iunit, '(a)') 'Number of fourier modes in table'
  write (iunit,*) mnmax_ws
  write (iunit, '(a)') 'Table of fourier coefficients'
  write (iunit, '(a)') 'm,n,crc2,czs2,crs2,czc2'

  do im = 0, mpol_coil
     do in = -ntor_coil, ntor_coil

        rc_tmp = 0.0_dp ; rs_tmp = 0.0_dp
        zc_tmp = 0.0_dp ; zs_tmp = 0.0_dp

        do itheta = 1,ntheta_coil
           do izeta = 1,nzetal_coil ! nzeta
              angle = im*theta_coil(itheta) + in*nfp*zetal_coil(izeta)  ! mu+nv
              r_tmp = sqrt(r_coil(1,itheta,izeta)*r_coil(1,itheta,izeta) + r_coil(2,itheta,izeta)*r_coil(2,itheta,izeta))
              rc_tmp = rc_tmp + r_tmp*cos(angle)
              rs_tmp = rs_tmp + r_tmp*sin(angle)
              zc_tmp = zc_tmp + r_coil(3,itheta,izeta)*cos(angle)
              zs_tmp = zs_tmp + r_coil(3,itheta,izeta)*sin(angle)
           enddo
        enddo

        if (im .eq. 0 .and. in .gt. 0) cycle
        if (im .ne. 0 .or. in .ne. 0) then
           rc_tmp = rc_tmp * 2.0_dp
           rs_tmp = rs_tmp * 2.0_dp
           zc_tmp = zc_tmp * 2.0_dp
           zs_tmp = zs_tmp * 2.0_dp
        endif
        rc_tmp = rc_tmp * factor
        rs_tmp = rs_tmp * factor
        zc_tmp = zc_tmp * factor
        zs_tmp = zs_tmp * factor

        if ( (abs(rc_tmp) + abs(rs_tmp) + abs(zc_tmp)+ abs(zs_tmp)) .ge. tol ) & ! write down only larger than tol
             write (iunit,'(x,2i6,1p4e20.12)') im, in, rc_tmp, zs_tmp, rs_tmp, zc_tmp
        
     enddo
  enddo

  print *, "Total harmonics number: ", mnmax_ws

  close(iunit)

end subroutine write_nescin
