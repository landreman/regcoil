subroutine regcoil_compute_lambda(prob)

  use regcoil_variables, only: regcoil_t
  use stel_kinds

  implicit none

  type(regcoil_t), intent(inout) :: prob
  integer :: j

  ! Allocate before associating — do not allocate/deallocate a selector while associated.
  if (allocated(prob%input%lambda)) deallocate(prob%input%lambda)
  allocate(prob%input%lambda(prob%input%nlambda))

  associate ( &
       nlambda => prob%input%nlambda, &
       lambda_min => prob%input%lambda_min, &
       lambda_max => prob%input%lambda_max, &
       lambda => prob%input%lambda, &
       general_option => prob%input%general_option, &
       verbose => prob%input%verbose &
       )

    lambda(1) = 0
    do j = 1,nlambda-1
       lambda(j+1) = lambda_min * exp((log(lambda_max/lambda_min)*(j-1))/(nlambda-2))
    end do

    if (general_option==1 .and. verbose) then
       print *,"We will use the following values of the regularization weight lambda:"
       print "(*(es10.3))",lambda
    end if

  end associate

end subroutine regcoil_compute_lambda
