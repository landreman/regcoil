subroutine regcoil_read_bnorm(prob)

  use regcoil_variables, only: regcoil_t
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none


  type(regcoil_t), intent(inout) :: prob
  integer :: iunit, i, mm, nn, itheta, izeta, index, iflag, num_modes_added
  real(dp) :: bf
  integer :: tic, toc, countrate

    associate ( &
       ntheta_plasma => prob%plasma%ntheta_plasma, &
       nzeta_plasma => prob%plasma%nzeta_plasma, &
       nfp => prob%plasma%nfp, &
       nbf => prob%plasma%nbf, &
       verbose => prob%input%verbose, &
       load_bnorm => prob%input%load_bnorm, &
       bnorm_filename => prob%input%bnorm_filename, &
       curpol => prob%input%curpol &
       )
  call system_clock(tic,countrate)

  if (allocated(prob%plasma%Bnormal_from_plasma_current)) deallocate(prob%plasma%Bnormal_from_plasma_current)
  allocate(prob%plasma%Bnormal_from_plasma_current(ntheta_plasma,nzeta_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'regcoil_read_bnorm Allocation error!'

  prob%plasma%Bnormal_from_plasma_current = 0

  if (.not. load_bnorm) then
     if (verbose) print *,"Not reading a bnorm file, so prob%plasma%Bnormal_from_plasma_current arrays will all be 0."
     return
  end if

  ! using data from FOCUS file ; by czhu
  if ( nbf > 0 ) then
     if (verbose) print *, "Using Bn coefficients in FOCUS file."

     do izeta = 1,nzeta_plasma
        do itheta = 1,ntheta_plasma
           do i = 1, nbf
              prob%plasma%Bnormal_from_plasma_current(itheta,izeta) = prob%plasma%Bnormal_from_plasma_current(itheta,izeta) + &
                   prob%plasma%bfc(i)*cos(prob%plasma%bfm(i)*prob%plasma%theta_plasma(itheta) - prob%plasma%bfn(i)*prob%plasma%zeta_plasma(izeta)) + &
                   prob%plasma%bfs(i)*sin(prob%plasma%bfm(i)*prob%plasma%theta_plasma(itheta) - prob%plasma%bfn(i)*prob%plasma%zeta_plasma(izeta))
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
           prob%plasma%Bnormal_from_plasma_current(itheta,izeta) = prob%plasma%Bnormal_from_plasma_current(itheta,izeta) + &
                bf*sin(mm*prob%plasma%theta_plasma(itheta) + nn*nfp*prob%plasma%zeta_plasma(izeta))

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
  prob%plasma%Bnormal_from_plasma_current = prob%plasma%Bnormal_from_plasma_current * curpol

!!$  ! I'll use an explicit loop here so there is no ambiguity about the order of theta vs zeta in the 1D vector:
!!$  do izeta = 1,nzeta_plasma
!!$     do itheta = 1,ntheta_plasma
!!$        Bnormal_from_plasma_current_1D((izeta-1)*ntheta_plasma+itheta) = prob%plasma%Bnormal_from_plasma_current(itheta,izeta)
!!$     end do
!!$  end do
!!$  ! The syntax in the next line might be faster? But I should check the order.
!!$  !Bnormal_from_plasma_current_1D = reshape(prob%plasma%Bnormal_from_plasma_current, (/ntheta_plasma*nzeta_plasma/))

  call system_clock(toc)
  if (num_modes_added>0) then
     if (verbose) print *,"Number of modes read from bnorm file:",num_modes_added
  else
     print *,"WARNING!!! No modes found in the bnorm file."
  end if
  if (verbose) print *,"Done reading B_normal on the plasma surface due to plasma current."
  if (verbose) print *,"Took ",real(toc-tic)/countrate," sec."


  end associate
end subroutine  regcoil_read_bnorm
