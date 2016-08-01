subroutine read_bnorm()

  use global_variables, only: load_bnorm, bnorm_filename, ntheta_plasma, nzeta_plasma, &
       Bnormal_from_plasma_current, theta_plasma, zeta_plasma, nfp, curpol
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none

  integer :: iunit, i, mm, nn, itheta, izeta, index, iflag
  real(dp) :: bf
  integer :: tic, toc, countrate

  call system_clock(tic,countrate)

  allocate(Bnormal_from_plasma_current(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  Bnormal_from_plasma_current = 0

  if (.not. load_bnorm) then
     print *,"Not reading a bnorm file, so Bnormal_from_plasma_current arrays will all be 0."
     return
  end if

  print *,"Loading B_normal on the plasma surface due to plasma current from file ",trim(bnorm_filename)

  call safe_open(iunit, i, trim(bnorm_filename), 'old', 'formatted')
  if (i .ne. 0 ) then
     stop 'Unable to open bnorm_filename.'
  end if

  do
     read(iunit,*,iostat = i) mm, nn, bf
     if (i .ne. 0) exit
     print *,"  Adding a mode with m=",mm,",  n=",nn
     do izeta = 1,nzeta_plasma
        do itheta = 1,ntheta_plasma
           Bnormal_from_plasma_current(itheta,izeta) = Bnormal_from_plasma_current(itheta,izeta) + &
                bf*sin(mm*theta_plasma(itheta) + nn*nfp*zeta_plasma(izeta))

           ! To see that it should be (mu+nv) rather than (mu-nv) in the above line, you can examine
           ! BNORM/Sources/bn_fouri.f (where the arrays in the bnorm files are computed)
           ! or
           ! either NESCOIL/Sources/bnfld.f (where bnorm files are read)
        end do
     end do
  end do

  close(iunit)

  ! BNORM scales B_n by curpol=(2*pi/nfp)*bsubv(m=0,n=0)
  ! where bsubv is the extrapolation to the last full mesh point of
  ! bsubvmnc.  Let's undo this scaling now.
  Bnormal_from_plasma_current = Bnormal_from_plasma_current * curpol

!!$  ! I'll use an explicit loop here so there is no ambiguity about the order of theta vs zeta in the 1D vector:
!!$  do izeta = 1,nzeta_plasma
!!$     do itheta = 1,ntheta_plasma
!!$        Bnormal_from_plasma_current_1D((izeta-1)*ntheta_plasma+itheta) = Bnormal_from_plasma_current(itheta,izeta)
!!$     end do
!!$  end do
!!$  ! The syntax in the next line might be faster? But I should check the order.
!!$  !Bnormal_from_plasma_current_1D = reshape(Bnormal_from_plasma_current, (/ntheta_plasma*nzeta_plasma/))

  call system_clock(toc)
  print *,"Done reading B_normal on the plasma surface due to plasma current."
  print *,"Took ",real(toc-tic)/countrate," sec."

end subroutine  read_bnorm
