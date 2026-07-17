subroutine regcoil_lambda_scan(prob)

  use regcoil_variables, only: regcoil_t

  implicit none

  type(regcoil_t), intent(inout) :: prob
  integer :: ilambda

  associate ( &
       nlambda => prob%input%nlambda, &
       lambda => prob%input%lambda &
       )

    do ilambda = 1,nlambda
       print "(a,es10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of ",nlambda,")"
       call regcoil_solve(prob, ilambda)
    end do

  end associate

end subroutine regcoil_lambda_scan
