subroutine regcoil_lambda_scan

  use regcoil_variables, only: nlambda, lambda

  implicit none

  integer :: ilambda


  do ilambda = 1,nlambda
     print "(a,es10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of ",nlambda,")"
     call regcoil_solve(ilambda)
  end do

end subroutine regcoil_lambda_scan
