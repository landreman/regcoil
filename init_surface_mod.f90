module init_surface_mod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine init_surface(ntheta, nzeta, nzetal, theta, zeta, zetal, &
         r, drdtheta, drdzeta, &
         normal, norm_normal, area, &
         geometry_option, R_specified, a, separation, dtheta, dzeta, nescin_filename, which_surface)

      use compute_offset_surface_mod
      use global_variables, only: R0_plasma, nfp, volume_coil
      use stel_kinds
      use stel_constants
      use omp_lib

      implicit none

      character(*) :: nescin_filename
      integer :: ntheta, nzeta, nzetal, geometry_option, iflag, which_surface
      real(dp) :: R_specified, a, separation, dtheta, dzeta, area
      real(dp), dimension(:), allocatable :: theta, zeta, zetal
      real(dp), dimension(:,:,:), allocatable :: r, drdtheta, drdzeta, normal
      ! These next 2 arrays are not used now, but may be needed in the future for Merkel's regularization
      real(dp), dimension(:,:,:), allocatable :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2
      real(dp), dimension(:,:), allocatable :: norm_normal, major_R_squared
      real(dp) :: R0_to_use
      real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta
      real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
      real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
      integer :: i, itheta, izeta
      real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_theta, delta_zeta, temp

      allocate(theta(ntheta),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 1'
      allocate(zeta(nzeta),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 2'
      allocate(zetal(nzetal),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 3'

      do i = 1,ntheta
         theta(i) = twopi*(i-1.0_dp)/ntheta
      end do

      do i = 1,nzeta
         zeta(i) = twopi/nfp*(i-1.0_dp)/nzeta
      end do

      do i = 1,nzetal
         zetal(i) = twopi*(i-1.0_dp)/nzetal
      end do

      dtheta = theta(2)-theta(1)
      dzeta  = zeta(2)-zeta(1)


      ! First dimension is the Cartesian component x, y, or z.
      allocate(r(3,ntheta,nzetal),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 4'
      allocate(drdtheta(3,ntheta,nzetal),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 5'
      allocate(drdzeta(3,ntheta,nzetal),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 6'
      allocate(normal(3,ntheta,nzetal),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 7'

      r = 0
      drdtheta = 0
      drdzeta = 0

      ! We do not use the 2nd derivatives presently, but these placeholders are here in case we need them in the future:
!!$      allocate(d2rdtheta2(3,ntheta,nzetal),stat=iflag)
!!$      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 8'
!!$      allocate(d2rdthetadzeta(3,ntheta,nzetal),stat=iflag)
!!$      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 9'
!!$      allocate(d2rdzeta2(3,ntheta,nzetal),stat=iflag)
!!$      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 10'
!!$
!!$      d2rdtheta2 = 0
!!$      d2rdthetadzeta = 0
!!$      d2rdzeta2 = 0

      if (geometry_option==3 .or. geometry_option == 4) then
         print *,"  Reading coil surface from nescin file ",trim(nescin_filename)

         call read_nescin(nescin_filename, r, drdtheta, drdzeta, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, &
              ntheta, nzetal, theta, zetal, .false.)
      end if


      select case (geometry_option)
      case (0,1)
         ! Torus with circular cross-section

         print *,"  Building a plain circular torus."

         if (geometry_option==0) then
            R0_to_use = R0_plasma
         else
            R0_to_use = R_specified
         end if

         do itheta = 1,ntheta
            angle = theta(itheta)
            sinangle = sin(angle)
            cosangle = cos(angle)
            dsinangledtheta = cosangle
            dcosangledtheta = -sinangle
            d2sinangledtheta2 = -sinangle
            d2cosangledtheta2 = -cosangle
            do izeta = 1,nzetal
               angle2 = zetal(izeta)
               sinangle2 = sin(angle2)
               cosangle2 = cos(angle2)
               dsinangle2dzeta = cosangle2
               dcosangle2dzeta = -sinangle2
               d2sinangle2dzeta2 = -sinangle2
               d2cosangle2dzeta2 = -cosangle2

               r(1,itheta,izeta) = (R0_to_use + a * cosangle) * cosangle2
               r(2,itheta,izeta) = (R0_to_use + a * cosangle) * sinangle2
               r(3,itheta,izeta) = a * sinangle

               drdtheta(1,itheta,izeta) = (a * dcosangledtheta) * cosangle2
               drdtheta(2,itheta,izeta) = (a * dcosangledtheta) * sinangle2
               drdtheta(3,itheta,izeta) = a * dsinangledtheta

               drdzeta(1,itheta,izeta) = (R0_to_use + a * cosangle) * dcosangle2dzeta
               drdzeta(2,itheta,izeta) = (R0_to_use + a * cosangle) * dsinangle2dzeta
               !drdzeta(3,itheta,izeta) = 0, so no equation needed for it here.

! This next bit is remarked out since we don't need 2nd derivatives presently
!!$               if (transfer_matrix_option==2 .and. which_surface == 1) then
!!$                  d2rdu2(1,itheta,izeta) = a * d2cosangledu2 * cosangle2
!!$                  d2rdu2(2,itheta,izeta) = a * d2cosangledu2 * sinangle2
!!$                  d2rdu2(3,itheta,izeta) = a * d2sinangledu2
!!$
!!$                  d2rdudzeta(1,itheta,izeta) = a * dcosangledu * dcosangle2dv
!!$                  d2rdudzeta(2,itheta,izeta) = a * dcosangledu * dsinangle2dv
!!$                  !d2rdudzeta(3,itheta,izeta) = 0, so no equation needed for it here.
!!$
!!$                  d2rdv2(1,itheta,izeta) = (R0_to_use + a * cosangle) * d2cosangle2dv2
!!$                  d2rdv2(2,itheta,izeta) = (R0_to_use + a * cosangle) * d2sinangle2dv2
!!$                  !d2rdv2(3,itheta,izeta) = 0, so no equation needed for it here.
!!$               end if
            end do
         end do

      case (2,4)

!!$         if (transfer_matrix_option==2) then
!!$            stop "Error! This geometry_option_middle is not yet implemented for transfer_matrix_option=2."
!!$         end if

         if (geometry_option==2) then
            print "(a,f10.4,a)","   Constructing a surface offset from the plasma by ",separation," meters."
         else
            print "(a,f10.4,a)","   Constructing a surface offset from the nescin surface by ",separation," meters."
         end if

         ! Finite differences to use:
         ! (Numerical Recipes suggests (machine epsilon)^(1/3)
         delta_theta = 1e-5;
         delta_zeta = 1e-5;
         ! Trick from Numerical Recipes for improved accuracy:
         temp = 1.0_dp + delta_theta
         delta_theta = temp - 1.0_dp
         temp = 1.0_dp + delta_zeta
         delta_zeta = temp - 1.0_dp
 

         !$OMP PARALLEL

         !$OMP MASTER
         print *,"  Number of OpenMP threads:",omp_get_num_threads()
         !$OMP END MASTER

         !$OMP DO PRIVATE(x_new,y_new,z_new,x_old,y_old,z_old)
         do itheta = 1,ntheta
            do izeta = 1,nzetal

               ! Compute r:
               call compute_offset_surface_xyz_of_thetazeta(theta(itheta),zetal(izeta),x_new,y_new,z_new,separation)
               r(1,itheta,izeta) = x_new
               r(2,itheta,izeta) = y_new
               r(3,itheta,izeta) = z_new

               ! Compute dr/du:
               call compute_offset_surface_xyz_of_thetazeta(theta(itheta)-delta_theta,zetal(izeta),x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_thetazeta(theta(itheta)+delta_theta,zetal(izeta),x_new,y_new,z_new,separation)
               drdtheta(1,itheta,izeta) = (x_new-x_old)/(2*delta_theta)
               drdtheta(2,itheta,izeta) = (y_new-y_old)/(2*delta_theta)
               drdtheta(3,itheta,izeta) = (z_new-z_old)/(2*delta_theta)

               ! Compute dr/dv:
               call compute_offset_surface_xyz_of_thetazeta(theta(itheta),zetal(izeta)-delta_zeta,x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_thetazeta(theta(itheta),zetal(izeta)+delta_zeta,x_new,y_new,z_new,separation)
               drdzeta(1,itheta,izeta) = (x_new-x_old)/(2*delta_zeta)
               drdzeta(2,itheta,izeta) = (y_new-y_old)/(2*delta_zeta)
               drdzeta(3,itheta,izeta) = (z_new-z_old)/(2*delta_zeta)

            end do

         end do
         !$OMP END DO
         !$OMP END PARALLEL

      case (3)
         ! Nothing to do - we already read the nescin file.

      case default
         print *,"Invalid setting for geometry_option: ",geometry_option
         stop
      end select

      ! Evaluate cross product:
      normal(1,:,:) = drdzeta(2,:,:) * drdtheta(3,:,:) - drdtheta(2,:,:) * drdzeta(3,:,:)
      normal(2,:,:) = drdzeta(3,:,:) * drdtheta(1,:,:) - drdtheta(3,:,:) * drdzeta(1,:,:)
      normal(3,:,:) = drdzeta(1,:,:) * drdtheta(2,:,:) - drdtheta(1,:,:) * drdzeta(2,:,:)

      allocate(norm_normal(ntheta, nzeta),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 11'
      norm_normal = sqrt(normal(1,:,1:nzeta)**2 + normal(2,:,1:nzeta)**2 + normal(3,:,1:nzeta)**2)

      area = nfp * dtheta * dzeta * sum(norm_normal)

      ! Compute coil surface volume using \int (1/2) R^2 dZ dzeta.
      ! These quantities will be evaluated on the half theta grid, which is the natural grid for dZ,
      ! but we will need to interpolate R^2 from the full to half grid.
      allocate(major_R_squared(ntheta,nzeta))
      major_R_squared = r(1,:,:)*r(1,:,:) + r(2,:,:)*r(2,:,:)
      ! First handle the interior of the theta grid:
      volume_coil = sum((major_R_squared(1:ntheta-1,:) + major_R_squared(2:ntheta,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r(3,2:ntheta,:)-r(3,1:ntheta-1,:))) ! dZ
      ! Add the contribution from the ends of the theta grid:
      volume_coil = volume_coil + sum((major_R_squared(1,:) + major_R_squared(ntheta,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r(3,1,:)-r(3,ntheta,:))) ! dZ
      volume_coil = abs(volume_coil * dzeta / 2) ! r includes all nfp periods already, so no factor of nfp needed.
      deallocate(major_R_squared)
      print "(a,es10.3,a,es10.3,a)"," Coil surface area:",area," m^2. Volume:",volume_coil," m^3."

    end subroutine init_surface

  end module init_surface_mod
  

