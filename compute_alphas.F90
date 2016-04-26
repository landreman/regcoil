subroutine compute_alphas

  use global_variables, only: N_alphas, alpha_min, alpha_max, alphas
  use stel_kinds

  implicit none

  integer :: j

  allocate(alphas(N_alphas))
  
  alphas(1) = 0
  do j = 1,N_alphas-1
     alphas(j+1) = alpha_min * exp((log(alpha_max/alpha_min)*(j-1))/(N_alphas-1))
  end do

  print *,"We will use the following values of the regularization weight alpha:"
  print *,alphas

end subroutine compute_alphas


