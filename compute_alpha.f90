subroutine compute_alpha

  use global_variables, only: nalpha, alpha_min, alpha_max, alpha
  use stel_kinds

  implicit none

  integer :: j

  allocate(alpha(nalpha))
  
  alpha(1) = 0
  do j = 1,nalpha-1
     alpha(j+1) = alpha_min * exp((log(alpha_max/alpha_min)*(j-1))/(nalpha-2))
  end do

  print *,"We will use the following values of the regularization weight alpha:"
  print *,alpha

end subroutine compute_alpha


