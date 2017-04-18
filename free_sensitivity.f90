subroutine free_sensitivity()

  use global_variables

  print *,"Deallocating sensitivity global variables."

  deallocate(xm_sensitivity)
  deallocate(xn_sensitivity)
  deallocate(dddomega)
  deallocate(dfxdomega)
  deallocate(dfydomega)
  deallocate(dfzdomega)
  deallocate(dnorm_normaldomega)
  deallocate(dchi2Kdomega)
  deallocate(dchi2Bdomega)
  deallocate(dgdomega)

end subroutine free_sensitivity
