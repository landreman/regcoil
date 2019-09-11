subroutine regcoil_read_single_Fourier()

    use regcoil_variables
    use safe_open_mod
    use stel_kinds

    implicit none

    integer :: iunit = 7, k, istat

    call safe_open(iunit, istat, trim(singleFourierFilename), 'old', 'formatted')
    if (istat .ne. 0) then
     stop 'Error opening single Fourier file file'
    endif

    ! Skip the first line
    read (iunit, *)
    read (iunit, *) nfp, B_0, lasym

    ! Skip the next 2 lines
    read (iunit, *)
    read (iunit, *)
    read (iunit, *) nmax_axis
    if (verbose) print *,"  Reading",nmax_axis,"axis modes from single Fourier series file."

    allocate(xn_axis(nmax_axis))
    allocate(raxis_cc(nmax_axis))
    allocate(zaxis_cs(nmax_axis))

    read (iunit, *)
    do k = 1, nmax_axis
        read (iunit, *) xn_axis(k), raxis_cc(k), zaxis_cs(k)
    end do

    ! Skip the next 2 lines
    read (iunit, *)
    read (iunit, *)
    read (iunit, *) mnmax_plasma
    if (verbose) print *,"  Reading",mnmax_plasma,"Fourier modes from single Fourier series file."

    if (allocated(xm_plasma)) deallocate(xm_plasma)
    if (allocated(xn_plasma)) deallocate(xn_plasma)
!    if (allocated(lmnc)) deallocate(lmnc)
!    if (allocated(lmns)) deallocate(lmns)

    allocate(lmnc(mnmax_plasma))
    allocate(lmns(mnmax_plasma))
    allocate(xm_plasma(mnmax_plasma))
    allocate(xn_plasma(mnmax_plasma))

    read (iunit, *)
    do k = 1, mnmax_plasma
        if (lasym) then
            read (iunit, *) xm_plasma(k), xn_plasma(k), lmnc(k), lmns(k)
        else
            read (iunit, *) xm_plasma(k), xn_plasma(k), lmnc(k)
        end if
    end do

!    print *,"lmnc(1,n_max+1:2*n_max+1): ", lmnc(1,n_max+1:2*n_max+1)
!    print *,"lmns(1,n_max+1:2*n_max+1): ", lmns(1,n_max+1:2*n_max+1)
!    print *,"m_list: ", m_list
!    print *,"n_list: ", n_list

    close(iunit)

end subroutine  regcoil_read_single_Fourier
