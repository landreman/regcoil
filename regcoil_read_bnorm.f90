subroutine regcoil_read_bnorm()

  use regcoil_variables, only: load_bnorm, bnorm_filename, ntheta_plasma, nzeta_plasma, &
       Bnormal_from_plasma_current, theta_plasma, zeta_plasma, nfp, curpol, verbose, &
       nbf, bfn, bfm, bfs, bfc
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none

  integer :: iunit, i, mm, nn, itheta, izeta, index, iflag, num_modes_added
  real(dp) :: bf
  integer :: tic, toc, countrate

  call system_clock(tic,countrate)

  if (allocated(Bnormal_from_plasma_current)) deallocate(Bnormal_from_plasma_current)
  allocate(Bnormal_from_plasma_current(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_read_bnorm Allocation error!'

  Bnormal_from_plasma_current = 0

  if (.not. load_bnorm) then
     if (verbose) print *,"Not reading a bnorm file, so Bnormal_from_plasma_current arrays will all be 0."
     return
  end if

  ! using data from FOCUS file ; by czhu
  if ( nbf > 0 ) then
     if (verbose) print *, "Using Bn coefficients in FOCUS file."

     do izeta = 1,nzeta_plasma
        do itheta = 1,ntheta_plasma
           do i = 1, nbf
              Bnormal_from_plasma_current(itheta,izeta) = Bnormal_from_plasma_current(itheta,izeta) + &
                   bfc(i)*cos(bfm(i)*theta_plasma(itheta) - bfn(i)*zeta_plasma(izeta)) + &
                   bfs(i)*sin(bfm(i)*theta_plasma(itheta) - bfn(i)*zeta_plasma(izeta))
              ! In FOCUS, the angle is using mu-nv. Both cosine and sin terms are retained.
              ! For more information, please check https://princetonuniversity.github.io/FOCUS/rdsurf.pdf.
           enddo
        end do
     end do
     return
  endif


  if (verbose) print *,"Loading B_normal on the plasma surface due to plasma current from file ",trim(bnorm_filename)

  call safe_open(iunit, i, trim(bnorm_filename), 'old', 'formatted')
  if (i .ne. 0 ) then
     stop 'Unable to open bnorm_filename.'
  end if

  num_modes_added = 0
  do
     read(iunit,*,iostat = i) mm, nn, bf
     if (i .ne. 0) exit
     num_modes_added = num_modes_added+1
     !print *,"  Adding a mode with m=",mm,",  n=",nn
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
  if (num_modes_added>0) then
     if (verbose) print *,"Number of modes read from bnorm file:",num_modes_added
  else
     print *,"WARNING!!! No modes found in the bnorm file."
  end if
  if (verbose) print *,"Done reading B_normal on the plasma surface due to plasma current."
  if (verbose) print *,"Took ",real(toc-tic)/countrate," sec."

end subroutine  regcoil_read_bnorm
