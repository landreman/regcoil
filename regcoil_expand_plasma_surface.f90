subroutine regcoil_expand_plasma_surface(theta, zeta, separation, x,y,z)
  
  use regcoil_variables, only: nfp, xm, xn, mnmax, rmnc, zmns, rmns, zmnc, lasym
  use stel_kinds
  use stel_constants

  implicit none

  real(dp), intent(in) :: theta,zeta,separation
  real(dp), intent(out) :: x,y,z

  integer :: imn
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
  real(dp) :: angle, sinangle, cosangle
  real(dp) :: dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
  real(dp) :: phi, sinphi, cosphi
  real(dp) :: dsinphidzeta, dcosphidzeta
  real(dp) :: normal_x, normal_y, normal_z, normal_abs

  x=0
  y=0
  z=0
  dxdtheta=0
  dydtheta=0
  dzdtheta=0
  dxdzeta=0
  dydzeta=0
  dzdzeta=0

  phi = zeta
  sinphi = sin(phi)
  cosphi = cos(phi)
  dsinphidzeta = cosphi
  dcosphidzeta = -sinphi

  do imn = 1,mnmax
     angle = xm(imn)*theta - xn(imn)*zeta
     cosangle = cos(angle)
     sinangle = sin(angle)
     dsinangledtheta = cosangle*xm(imn)
     dcosangledtheta = -sinangle*xm(imn)
     dsinangledzeta = -cosangle*xn(imn)
     dcosangledzeta = sinangle*xn(imn)

     x = x + rmnc(imn) * cosangle * cosphi
     y = y + rmnc(imn) * cosangle * sinphi
     z = z + zmns(imn) * sinangle
     
     dxdtheta = dxdtheta + rmnc(imn) * dcosangledtheta * cosphi
     dydtheta = dydtheta + rmnc(imn) * dcosangledtheta * sinphi
     dzdtheta = dzdtheta + zmns(imn) * dsinangledtheta
     
     dxdzeta = dxdzeta + rmnc(imn) * (dcosangledzeta * cosphi + cosangle * dcosphidzeta)
     dydzeta = dydzeta + rmnc(imn) * (dcosangledzeta * sinphi + cosangle * dsinphidzeta)
     dzdzeta = dzdzeta + zmns(imn) * dsinangledzeta

     if (lasym) then
        x = x + rmns(imn) * sinangle * cosphi
        y = y + rmns(imn) * sinangle * sinphi
        z = z + zmnc(imn) * cosangle
     
        dxdtheta = dxdtheta + rmns(imn) * dsinangledtheta * cosphi
        dydtheta = dydtheta + rmns(imn) * dsinangledtheta * sinphi
        dzdtheta = dzdtheta + zmnc(imn) * dcosangledtheta
     
        dxdzeta = dxdzeta + rmns(imn) * (dsinangledzeta * cosphi + sinangle * dcosphidzeta)
        dydzeta = dydzeta + rmns(imn) * (dsinangledzeta * sinphi + sinangle * dsinphidzeta)
        dzdzeta = dzdzeta + zmnc(imn) * dcosangledzeta
     end if
     
  end do

  ! Evaluate the (non-unit-magnitude) surface normal vector from N = (dr/dv) cross (dr/du):
  normal_x = dydzeta * dzdtheta - dydtheta * dzdzeta
  normal_y = dzdzeta * dxdtheta - dzdtheta * dxdzeta
  normal_z = dxdzeta * dydtheta - dxdtheta * dydzeta

  ! Now normalize the normal vector so it has unit magnitude:
  normal_abs = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
  normal_x = normal_x / normal_abs
  normal_y = normal_y / normal_abs
  normal_z = normal_z / normal_abs

  ! Move in the normal direction away from the plasma surface:
  x = x + normal_x * separation
  y = y + normal_y * separation
  z = z + normal_z * separation

end subroutine regcoil_expand_plasma_surface
