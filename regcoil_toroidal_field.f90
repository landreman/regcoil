subroutine regcoil_toroidal_field()
  use regcoil_variables, only: ntheta_plasma, nzeta_plasma, r_plasma, normal_plasma, verbose, &
       Bnormal_from_plasma_current, norm_normal_plasma, curpol, net_poloidal_current_Amperes
  use stel_kinds
  
  implicit none

  integer ::  itheta, izeta
  real(dp) :: rr, Bx, By, Bt
  real(dp), parameter :: two = 2, bsconstant =  1.0E-07

  do izeta = 1,nzeta_plasma
     do itheta = 1,ntheta_plasma
        ! Bt = u0*I/(2 pi R)
        rr = sqrt( r_plasma(1,itheta,izeta)**2 + r_plasma(2,itheta,izeta)**2 )
        Bt = two/rr * net_poloidal_current_Amperes * bsconstant
        Bx = - Bt * r_plasma(2,itheta,izeta)/rr
        By =   Bt * r_plasma(1,itheta,izeta)/rr
        ! Bz =   0
        Bnormal_from_plasma_current(itheta,izeta) = Bnormal_from_plasma_current(itheta,izeta) &
             & + Bx * normal_plasma(1,itheta,izeta) / norm_normal_plasma(itheta,izeta) &
             & + By * normal_plasma(2,itheta,izeta) / norm_normal_plasma(itheta,izeta) 
     enddo
  enddo

  if (verbose) print *,"Toroidal field contribution has been added to Bnormal_from_plasma_current. "

  return 
end subroutine regcoil_toroidal_field
