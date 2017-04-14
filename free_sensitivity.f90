subroutine free_sensitivity()

  use global_variables

  print *,"Deallocating sensitivity global variables."

  deallocate(xm_sensitivity)
  deallocate(xn_sensitivity)
  deallocate(dddrmnc)
  deallocate(dddrmns)
  deallocate(dddzmnc)
  deallocate(dddzmns)
  !deallocate(dfdrmnc)
  !deallocate(dfdrmns)
  !deallocate(dfdzmnc)
  !deallocate(dfdzmns)
  deallocate(dnorm_normaldrmnc)
  deallocate(dnorm_normaldrmns)
  deallocate(dnorm_normaldzmnc)
  deallocate(dnorm_normaldzmns)
  deallocate(dchi2Kdrmnc)
  deallocate(dchi2Kdrmns)
  deallocate(dchi2Kdzmnc)
  deallocate(dchi2Kdzmns)
  deallocate(dchi2Bdrmnc)
  deallocate(dchi2Bdrmns)
  deallocate(dchi2Bdzmnc)
  deallocate(dchi2Bdzmns)
  deallocate(dgdrmnc)
  deallocate(dgdrmns)
  deallocate(dgdzmnc)
  deallocate(dgdzmns)
end subroutine free_sensitivity
