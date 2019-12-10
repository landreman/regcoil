subroutine regcoil_write_single_Fourier
    ! Write the single Fourier series amplitudes to a file.
    use safe_open_mod
    use regcoil_variables
!    use regcoil_variables, only: nfp, B_0, lasym, nmax_axis, mnmax_plasma, xn_axis, xm_plasma, xn_plasma, &
!        raxis_cc, zaxis_cs, lmnc, lmns, curpol, R0_plasma, singleFourierFilename

    implicit none

    integer :: i, istat = 0, iunit = 8


    call safe_open(iunit, istat, trim(singleFourierFilename), 'replace', 'formatted')
    if (istat .ne. 0) then
        stop 'Error opening nescin file: file exsited or something wrong.'
    end if

    write (iunit, '(a)') 'nfp,    B_0,      net_poloidal_current_amperes, curpol,             R0_plasma,lasym'
    write (iunit, '(1I6,1p4e20.12,1L2)') nfp, B_0, net_poloidal_current_amperes, curpol, R0_plasma, lasym
    write (iunit, *)

    write (iunit, '(a)') 'nmax_axis'
    write (iunit, '(1I6)') nmax_axis   ! write the number of axis modes
    write (iunit, '(a)') 'xn_axis, raxis_cc,           raxis_cs'
    do i = 1,nmax_axis
        write (iunit,'(x,1I6,1p2e20.12)') xn_axis(i), raxis_cc(i), zaxis_cs(i)
    end do
    write (iunit, *)

    write (iunit, '(a)') 'mnmax_plasma'
    write (iunit, '(1I6)') mnmax_plasma   ! write the number of plasma modes
    write (iunit, '(a)') '      xm,  xn, lmnc,               lmns'
    do i = 1,mnmax_plasma
        if (lasym) then
            write (iunit,'(x,2I6,1p2e20.12)') xm_plasma(i), xn_plasma(i), lmnc(i), lmns(i)
        else
            write (iunit,'(x,2I6,1p1e20.12)') xm_plasma(i), xn_plasma(i), lmnc(i)
        end if

    end do

    close(iunit)

end subroutine regcoil_write_single_Fourier









