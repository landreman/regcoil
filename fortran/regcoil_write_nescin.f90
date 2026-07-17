subroutine regcoil_write_nescin
  ! Write the harmonics of the coil winding surface to a file; by czhu on 05/17/2018;
  use regcoil_variables, only: mnmax_coil, xm_coil, xn_coil, &
       rmnc_coil, rmns_coil, zmnc_coil, zmns_coil, &
       separation, nescin_filename, nfp, curpol
  use stel_kinds
  use safe_open_mod

  implicit none

  integer :: imn, mnmax_ws, istat = 0, iunit = 8
  real(dp) :: tol

  tol = separation * 1.0E-3_dp ! set the tolarence to ignore harmonics
  tol = 0 ! MJL 20180617 Causes all harmonics to be written.

  call safe_open(iunit, istat, trim(nescin_filename), 'replace', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file: file exsited or something wrong.'
  end if
  print *, "Writing winding surface into ", trim(nescin_filename)

  ! Count number of non-"zero" mode amplitudes
  mnmax_ws = 0
  do imn = 1, mnmax_coil
     if ( (abs(rmnc_coil(imn)) + abs(rmns_coil(imn)) + abs(zmnc_coil(imn))+ abs(zmns_coil(imn))) .ge. tol ) then
        mnmax_ws = mnmax_ws + 1
     end if
  end do


  write (iunit, '(a)') '------ Plasma information from VMEC ----'
  write (iunit, '(a)') 'np     iota_edge       phip_edge       curpol'
  write (iunit, '(I6, 3ES20.12)') nfp, 0.0, 0.0, curpol  ! write nfp and curpol information 

  write (iunit,*)
  write (iunit, '(a, 1pe20.12, a)') '------ Current Surface: Coil-Plasma separation = ', separation,' -----'
  write (iunit, '(a)') 'Number of fourier modes in table'
  write (iunit,*) mnmax_ws
  write (iunit, '(a)') 'Table of fourier coefficients'
  write (iunit, '(a)') 'm,n,crc2,czs2,crs2,czc2'

  do imn = 1, mnmax_coil
     if ( (abs(rmnc_coil(imn)) + abs(rmns_coil(imn)) + abs(zmnc_coil(imn))+ abs(zmns_coil(imn))) .ge. tol ) then
        ! Note in the next line that we convert the toroidal mode number n from regcoil/vmec convention (angle = m*theta - n*zeta) to
        ! nescin convention (angle = m*theta + n*nfp*zeta).
        write (iunit,'(x,2i6,1p4e20.12)') xm_coil(imn), -xn_coil(imn)/nfp, rmnc_coil(imn), zmns_coil(imn), rmns_coil(imn), zmnc_coil(imn)
     end if
  end do

  print *, "Number of harmonics written: ", mnmax_ws

  close(iunit)

end subroutine regcoil_write_nescin
