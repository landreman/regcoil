module init_Fourier_modes_mod

  implicit none

contains

  subroutine init_Fourier_modes(mpol, ntor, mnmax, xm, xn)

    implicit none

    integer :: mpol, ntor, mnmax
    integer, dimension(:), allocatable :: xm, xn
    
    integer :: jn, jm, index, iflag

    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be positive.
    mnmax = mpol*(ntor*2+1) + ntor
    
    allocate(xm(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xn(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    ! Handle the xm=0 modes:
    xm=0
    do jn=1,ntor
       xn(jn)=jn
    end do
    
    ! Handle the xm>0 modes:
    index = ntor
    do jm = 1,mpol
       do jn = -ntor, ntor
          index = index + 1
          xn(index) = jn
          xm(index) = jm
       end do
    end do
    
    if (index .ne. mnmax) then
       print *,"Error!  index=",index," but mnmax=",mnmax
       stop
    end if

  end subroutine init_Fourier_modes

  subroutine init_Fourier_modes_sensitivity(mpol,ntor,mnmax,nomega,xm,xn,omega)

    implicit none

    integer :: mpol, ntor, mnmax, iomega, i, nomega
    integer, dimension(:), allocatable :: xm, xn, omega

    integer :: jn, jm, iflag

    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be non-negative
    mnmax = mpol*(ntor*2+1) + ntor+1

    ! nomega is the length of the number of fourier coefficients
    nomega = mnmax*4

    allocate(xm(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xn(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(omega(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    ! Handle the xm=0 modes:
    xm=0
    iomega = 0
    do jn=1,ntor+1
      do i=1,4
        iomega = iomega + 1
        omega(iomega) = i
        xn(iomega)=jn-1
      enddo
    end do

    ! Handle the xm>0 modes:
    do jm = 1,mpol
      do jn = -ntor, ntor
        do i=1,4
          iomega = iomega + 1
          xn(iomega) = jn
          xm(iomega) = jm
          omega(iomega) = i
        enddo
      end do
    end do

    if (iomega .ne. nomega) then
      print *,"Error!  iomega=",iomega," but nomega=",nomega
      stop
    end if

  end subroutine init_Fourier_modes_sensitivity

end module init_Fourier_modes_mod
