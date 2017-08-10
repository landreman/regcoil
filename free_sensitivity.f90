subroutine free_sensitivity()

  use global_variables

  print *,"Deallocating sensitivity global variables."

  deallocate(omega_coil)
  deallocate(xm_sensitivity)
  deallocate(xn_sensitivity)
  deallocate(dddomega)
  deallocate(dfxdomega)
  deallocate(dfydomega)
  deallocate(dfzdomega)
  deallocate(dnorm_normaldomega)
  deallocate(dgdomega)
  deallocate(dinductancedomega)
  deallocate(dnormxdomega)
  deallocate(dnormydomega)
  deallocate(dnormzdomega)
  deallocate(dchi2Kdomega)
  deallocate(dchi2Bdomega)
  deallocate(dchi2domega)
  deallocate(drdomega)
  deallocate(domegadxdtheta)
  deallocate(domegadxdzeta)
  deallocate(domegadydtheta)
  deallocate(domegadydzeta)
  deallocate(domegadzdtheta)
  deallocate(domegadzdzeta)
  deallocate(dhdomega)
  !deallocate(dchi2dr_normal)
  if (normal_displacement_option == 2) then
    deallocate(domegadx)
    deallocate(domegady)
    deallocate(domegadz)
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 4) then
    deallocate(q_K)
  endif
  if (sensitivity_option == 3 .or. sensitivity_option == 5) then
    deallocate(q_B)
  endif
  if (sensitivity_option > 2) then
    deallocate(dmatrix_Kdomega)
    deallocate(dmatrix_Bdomega)
    deallocate(dmatrixdomega)
    deallocate(dRHS_Kdomega)
    deallocate(dRHS_Bdomega)
    deallocate(dRHSdomega)
  endif
  if (sensitivity_option == 3) then
    deallocate(adjoint_sum)
  endif

end subroutine free_sensitivity
