module regcoil_compute_lambda

contains
 
subroutine compute_lambda(lscreen_optin)

  use regcoil_variables, only: nlambda, lambda_min, lambda_max, lambda, general_option
  use stel_kinds

  implicit none
  ! variables to handle printing to the screen
  logical, optional :: lscreen_optin
  logical :: lscreen

  integer :: j

  if(present(lscreen_optin)) then 
    lscreen = lscreen_optin
  else
    lscreen = .true.
  endif

  ! Adding a check to release previously allocated variable.
  ! This is because STELLOPT may call this function multiple times.
  if (allocated(lambda)) deallocate(lambda)
  allocate(lambda(nlambda))
  
  lambda(1) = 0
  do j = 1,nlambda-1
     lambda(j+1) = lambda_min * exp((log(lambda_max/lambda_min)*(j-1))/(nlambda-2))
  end do

  if (general_option==1) then
     print *,"We will use the following values of the regularization weight lambda:"
     print "(*(es10.3))",lambda
  end if

end subroutine compute_lambda

end module regcoil_compute_lambda
