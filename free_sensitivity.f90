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
  deallocate(dinductancedomega)
  deallocate(dnormxdomega)
  deallocate(dnormydomega)
  deallocate(dnormzdomega)
  if (sensitivity_option == 3) then
    deallocate(dAKdomega)
    deallocate(dABdomega)
    deallocate(dbKdomega)
    deallocate(dbBdomega)
  endif

end subroutine free_sensitivity
