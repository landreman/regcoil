subroutine regcoil_write_nescin
  ! FFT the offset winding surface and write the harmonics to file; by czhu on 05/17/2018;
  use regcoil_variables, only: R0_plasma, nfp, volume_coil, &
       ntheta_coil, nzeta_coil, nzetal_coil, &
       theta_coil, zeta_coil, zetal_coil, &
       r_coil, drdtheta_coil, drdzeta_coil, &
       normal_coil, norm_normal_coil, area_coil, &
       geometry_option_coil, R0_coil, a_coil,  &
       separation, dtheta_coil, dzeta_coil, &
       nescin_filename, curpol, max_mpol_coil, max_ntor_coil

  use stel_kinds
  use stel_constants
  use safe_open_mod
  use omp_lib

  implicit none


  real(dp) :: angle
  integer :: i, itheta, izeta, im, in, mnmax_ws, istat = 0, iunit = 8
  integer :: mpol_coil, ntor_coil
  real(dp) :: r_tmp, rc_tmp, rs_tmp, zc_tmp, zs_tmp, factor, tol

  mpol_coil = min(ntheta_coil / 2, max_mpol_coil)
  ntor_coil = min( nzeta_coil / 2, max_ntor_coil)
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

end subroutine regcoil_write_nescin
