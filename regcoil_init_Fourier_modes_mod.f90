module regcoil_init_Fourier_modes_mod

  implicit none

contains

  subroutine regcoil_init_Fourier_modes(mpol, ntor, mnmax, xm, xn, include_00)

    implicit none

    integer :: mpol, ntor, mnmax
    integer, dimension(:), allocatable :: xm, xn
    logical, intent(in) :: include_00
    
    integer :: jn, jm, index, iflag
    integer, dimension(:), allocatable :: xm_temp, xn_temp
    
    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be 0 or positive.
    mnmax = mpol*(ntor*2+1) + ntor
    if (include_00) mnmax = mnmax + 1
   

    if (allocated(xm)) deallocate(xm)
    allocate(xm(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'init_Fourier Allocation error!'
 
    if (allocated(xn)) deallocate(xn)
    allocate(xn(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'init_Fourier Allocation error!'

    xm = 0
    xn = 0

    ! Handle the xm=0 modes:
    index = 0
    if (include_00) index = 1
    do jn = 1, ntor
       index = index + 1
       xn(index)=jn
    end do
    
    ! Handle the xm>0 modes:
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

  end subroutine regcoil_init_Fourier_modes
  
end module regcoil_init_Fourier_modes_mod
