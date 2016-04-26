module init_surface_mod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine init_surface(nu, nv, nvl, u, v, vl, &
         r, drdu, drdv, &
         d2rdu2, d2rdudv, d2rdv2, &
         normal, norm_normal, area, &
         geometry_option, R_specified, a, separation, du, dv, nescin_filename, which_surface)

      use global_variables, only: R0_plasma, nfp, transfer_matrix_option
      use stel_kinds
      use stel_constants
      use omp_lib

      implicit none

      character(*) :: nescin_filename
      integer :: nu, nv, nvl, geometry_option, iflag, which_surface
      real(dp) :: R_specified, a, separation, du, dv, area
      real(dp), dimension(:), allocatable :: u, v, vl
      real(dp), dimension(:,:,:), allocatable :: r, drdu, drdv, normal
      real(dp), dimension(:,:,:), allocatable :: d2rdu2, d2rdudv, d2rdv2
      real(dp), dimension(:,:), allocatable :: norm_normal
      real(dp) :: R0_to_use
      real(dp) :: angle, sinangle, cosangle, dsinangledu, dcosangledu
      real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv
      real(dp) :: d2sinangledu2, d2cosangledu2, d2sinangle2dv2, d2cosangle2dv2
      integer :: i, iu, iv
      real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_u, delta_v, temp

      allocate(u(nu),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 1'
      allocate(v(nv),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 2'
      allocate(vl(nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 3'

      do i = 1,nu
         u(i) = (i-1.0_dp)/nu
      end do

      do i = 1,nv
         v(i) = (i-1.0_dp)/nv
      end do

      do i = 1,nvl
         vl(i) = (i-1.0_dp)/nv
      end do

      du = u(2)-u(1)
      dv = v(2)-v(1)

!!$      if (transfer_matrix_option==2 .and. which_surface==2) then
!!$         u = u + du/2
!!$         v = v + dv/2
!!$         vl = vl + dv/2
!!$      end if

      ! First dimension is the Cartesian component x, y, or z.
      allocate(r(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 4'
      allocate(drdu(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 5'
      allocate(drdv(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 6'
      allocate(normal(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 7'

      r = 0
      drdu = 0
      drdv = 0

      if (transfer_matrix_option == 2 .and. which_surface == 1) then
         allocate(d2rdu2(3,nu,nvl),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 8'
         allocate(d2rdudv(3,nu,nvl),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 9'
         allocate(d2rdv2(3,nu,nvl),stat=iflag)
         if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 10'

         d2rdu2 = 0
         d2rdudv = 0
         d2rdv2 = 0
      end if

      if (geometry_option==3 .or. geometry_option == 4) then
         print *,"  Reading coil surface from nescin file ",trim(nescin_filename)

         call read_nescin(nescin_filename, r, drdu, drdv, d2rdu2, d2rdudv, d2rdv2, &
              nu, nvl, u, vl, transfer_matrix_option==2 .and. which_surface==1)
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

         do iu = 1,nu
            angle = twopi*u(iu)
            sinangle = sin(angle)
            cosangle = cos(angle)
            dsinangledu = cosangle*twopi
            dcosangledu = -sinangle*twopi
            d2sinangledu2 = -twopi*twopi*sinangle
            d2cosangledu2 = -twopi*twopi*cosangle
            do iv = 1,nvl
               angle2 = twopi * vl(iv) / nfp
               sinangle2 = sin(angle2)
               cosangle2 = cos(angle2)
               dsinangle2dv = cosangle2*twopi/nfp
               dcosangle2dv = -sinangle2*twopi/nfp
               d2sinangle2dv2 = -twopi*twopi/(nfp*nfp)*sinangle2
               d2cosangle2dv2 = -twopi*twopi/(nfp*nfp)*cosangle2

               r(1,iu,iv) = (R0_to_use + a * cosangle) * cosangle2
               r(2,iu,iv) = (R0_to_use + a * cosangle) * sinangle2
               r(3,iu,iv) = a * sinangle

               drdu(1,iu,iv) = (a * dcosangledu) * cosangle2
               drdu(2,iu,iv) = (a * dcosangledu) * sinangle2
               drdu(3,iu,iv) = a * dsinangledu

               drdv(1,iu,iv) = (R0_to_use + a * cosangle) * dcosangle2dv
               drdv(2,iu,iv) = (R0_to_use + a * cosangle) * dsinangle2dv
               !drdv(3,iu,iv) = 0, so no equation needed for it here.

               if (transfer_matrix_option==2 .and. which_surface == 1) then
                  d2rdu2(1,iu,iv) = a * d2cosangledu2 * cosangle2
                  d2rdu2(2,iu,iv) = a * d2cosangledu2 * sinangle2
                  d2rdu2(3,iu,iv) = a * d2sinangledu2

                  d2rdudv(1,iu,iv) = a * dcosangledu * dcosangle2dv
                  d2rdudv(2,iu,iv) = a * dcosangledu * dsinangle2dv
                  !d2rdudv(3,iu,iv) = 0, so no equation needed for it here.

                  d2rdv2(1,iu,iv) = (R0_to_use + a * cosangle) * d2cosangle2dv2
                  d2rdv2(2,iu,iv) = (R0_to_use + a * cosangle) * d2sinangle2dv2
                  !d2rdv2(3,iu,iv) = 0, so no equation needed for it here.
               end if
            end do
         end do

      case (2,4)

         if (transfer_matrix_option==2) then
            stop "Error! This geometry_option_middle is not yet implemented for transfer_matrix_option=2."
         end if

         if (geometry_option==2) then
            print *,"  Constructing a surface offset from the plasma by ",separation
         else
            print *,"  Constructing a surface offset from the nescin surface by ",separation
         end if

         ! Finite differences to use:
         ! (Numerical Recipes suggests (machine epsilon)^(1/3)
         delta_u = 1e-5;
         delta_v = 1e-5;
         ! Trick from Numerical Recipes for improved accuracy:
         temp = 1.0_dp + delta_u
         delta_u = temp - 1.0_dp
         temp = 1.0_dp + delta_v
         delta_v = temp - 1.0_dp
 

         !$OMP PARALLEL

         !$OMP MASTER
         print *,"  Number of OpenMP threads:",omp_get_num_threads()
         !$OMP END MASTER

         !$OMP DO PRIVATE(x_new,y_new,z_new,x_old,y_old,z_old)
         do iu = 1,nu
            do iv = 1,nvl

               ! Compute r:
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv),x_new,y_new,z_new,separation)
               r(1,iu,iv) = x_new
               r(2,iu,iv) = y_new
               r(3,iu,iv) = z_new

               ! Compute dr/du:
               call compute_offset_surface_xyz_of_uv(u(iu)-delta_u,vl(iv),x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_uv(u(iu)+delta_u,vl(iv),x_new,y_new,z_new,separation)
               drdu(1,iu,iv) = (x_new-x_old)/(2*delta_u)
               drdu(2,iu,iv) = (y_new-y_old)/(2*delta_u)
               drdu(3,iu,iv) = (z_new-z_old)/(2*delta_u)

               ! Compute dr/dv:
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv)-delta_v,x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv)+delta_v,x_new,y_new,z_new,separation)
               drdv(1,iu,iv) = (x_new-x_old)/(2*delta_v)
               drdv(2,iu,iv) = (y_new-y_old)/(2*delta_v)
               drdv(3,iu,iv) = (z_new-z_old)/(2*delta_v)

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
      normal(1,:,:) = drdv(2,:,:) * drdu(3,:,:) - drdu(2,:,:) * drdv(3,:,:)
      normal(2,:,:) = drdv(3,:,:) * drdu(1,:,:) - drdu(3,:,:) * drdv(1,:,:)
      normal(3,:,:) = drdv(1,:,:) * drdu(2,:,:) - drdu(1,:,:) * drdv(2,:,:)

      allocate(norm_normal(nu, nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 11'
      norm_normal = sqrt(normal(1,:,:)**2 + normal(2,:,:)**2 + normal(3,:,:)**2)

      area = du * dv * sum(norm_normal)

    end subroutine init_surface

  end module init_surface_mod
  

