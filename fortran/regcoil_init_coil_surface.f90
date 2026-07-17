subroutine regcoil_init_coil_surface(prob)

  use regcoil_compute_offset_surface_mod
  use regcoil_variables, only: regcoil_t
  use regcoil_init_Fourier_modes_mod
  use regcoil_splines
  
  use stel_kinds
  use stel_constants
  use omp_lib
  
  implicit none
  

  type(regcoil_t), intent(inout) :: prob
  integer :: iflag
  real(dp), dimension(:,:), allocatable :: major_R_coil, z_coil
  real(dp) :: R0_to_use, relative_variation_in_dl
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
  real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
  integer :: i, j, itheta, izeta, iteration, max_iterations
  real(dp) :: x_new, y_new, z_new, delta_theta, delta_zeta, temp, factor, factor2
  integer :: tic, toc, countrate, tic1, toc1
  integer :: mpol_coil, ntor_coil
  integer, allocatable, dimension(:) :: xm, xn
  real(dp), allocatable, dimension(:) :: rmnc, rmns, zmnc, zmns  
  real(dp), allocatable, dimension(:) :: theta_plasma_corresponding_to_theta_coil, dl, R_slice, z_slice
  real(dp), allocatable, dimension(:) :: constant_arclength_theta_on_old_theta_grid
  type (periodic_spline) :: theta_spline

    associate ( &
       R0_plasma => prob%plasma%R0_plasma, &
       nfp => prob%plasma%nfp, &
       lasym => prob%plasma%lasym, &
       constant_arclength_tolerance => prob%plasma%constant_arclength_tolerance, &
       ntheta_coil => prob%coil%ntheta_coil, &
       nzeta_coil => prob%coil%nzeta_coil, &
       nzetal_coil => prob%coil%nzetal_coil, &
       geometry_option_coil => prob%coil%geometry_option_coil, &
       R0_coil => prob%coil%R0_coil, &
       a_coil => prob%coil%a_coil, &
       separation => prob%coil%separation, &
       nescin_filename => prob%coil%nescin_filename, &
       dtheta_coil => prob%coil%dtheta_coil, &
       dzeta_coil => prob%coil%dzeta_coil, &
       mnmax_coil => prob%coil%mnmax_coil, &
       max_mpol_coil => prob%coil%max_mpol_coil, &
       max_ntor_coil => prob%coil%max_ntor_coil, &
       verbose => prob%input%verbose &
       )
  call system_clock(tic,countrate)
  if (verbose) print *,"Initializing coil surface."

  ! In case STELLOPT or some other code calls subroutines of regcoil multiple times, make sure arrays are deallocated:
  if(allocated(prob%coil%xm_coil)) deallocate(prob%coil%xm_coil)
  if(allocated(prob%coil%xn_coil)) deallocate(prob%coil%xn_coil)
  if(allocated(prob%coil%rmnc_coil)) deallocate(prob%coil%rmnc_coil)
  if(allocated(prob%coil%rmns_coil)) deallocate(prob%coil%rmns_coil)
  if(allocated(prob%coil%zmnc_coil)) deallocate(prob%coil%zmnc_coil)
  if(allocated(prob%coil%zmns_coil)) deallocate(prob%coil%zmns_coil)

  nzetal_coil   = nzeta_coil   * nfp

  if (allocated(prob%coil%theta_coil)) deallocate(prob%coil%theta_coil)
  allocate(prob%coil%theta_coil(ntheta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 1'

  if (allocated(prob%coil%zeta_coil)) deallocate(prob%coil%zeta_coil)
  allocate(prob%coil%zeta_coil(nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 2'

  if (allocated(prob%coil%zetal_coil)) deallocate(prob%coil%zetal_coil)
  allocate(prob%coil%zetal_coil(nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 3'

  do i = 1,ntheta_coil
     prob%coil%theta_coil(i) = twopi*(i-1.0_dp)/ntheta_coil
  end do

  do i = 1,nzeta_coil
     prob%coil%zeta_coil(i) = twopi/nfp*(i-1.0_dp)/nzeta_coil
  end do

  do i = 1,nzetal_coil
     prob%coil%zetal_coil(i) = twopi*(i-1.0_dp)/nzetal_coil
  end do

  dtheta_coil = prob%coil%theta_coil(2)-prob%coil%theta_coil(1)
  dzeta_coil  = prob%coil%zeta_coil(2)-prob%coil%zeta_coil(1)

  select case (geometry_option_coil)
  case (0,1)
     ! Torus with circular cross-section
     
     if (verbose) print *,"  Building a plain circular torus."
     
     if (geometry_option_coil==0) then
        R0_to_use = R0_plasma
     else
        R0_to_use = R0_coil
     end if
     
     mnmax_coil = 2
     allocate(prob%coil%xm_coil(mnmax_coil))
     allocate(prob%coil%xn_coil(mnmax_coil))
     allocate(prob%coil%rmnc_coil(mnmax_coil))
     allocate(prob%coil%rmns_coil(mnmax_coil))
     allocate(prob%coil%zmnc_coil(mnmax_coil))
     allocate(prob%coil%zmns_coil(mnmax_coil))
     prob%coil%xm_coil = (/ 0, 1 /)
     prob%coil%xn_coil = 0
     prob%coil%rmnc_coil = (/ R0_to_use, a_coil /)
     prob%coil%rmns_coil = 0
     prob%coil%zmnc_coil = 0
     prob%coil%zmns_coil = (/ 0.0d+0, a_coil /)

  case (2, 4)
     
     if (verbose) print "(a,f10.4,a)","   Constructing a surface offset from the plasma by ",separation," meters."
     if (verbose .and. geometry_option_coil==4) print *,"  Using the constant-arclength theta coordinate."

     allocate(major_R_coil(ntheta_coil, nzeta_coil))
     allocate(z_coil(ntheta_coil, nzeta_coil))
     major_R_coil = 0
     z_coil = 0
     max_iterations = 100

     call system_clock(tic1)
     !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iteration,x_new,y_new,z_new,theta_plasma_corresponding_to_theta_coil,dl,relative_variation_in_dl,R_slice,z_slice,itheta,constant_arclength_theta_on_old_theta_grid,theta_spline)
     
     allocate(theta_plasma_corresponding_to_theta_coil(Ntheta_coil))
     allocate(R_slice(Ntheta_coil))
     allocate(z_slice(Ntheta_coil))
     allocate(dl(Ntheta_coil))
     allocate(constant_arclength_theta_on_old_theta_grid(Ntheta_coil+1))
     constant_arclength_theta_on_old_theta_grid = 0

     !$OMP MASTER
     if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
     !$OMP END MASTER
     
     !$OMP DO
     do izeta = 1,nzeta_coil
        theta_plasma_corresponding_to_theta_coil = prob%coil%theta_coil
        do iteration = 1,max_iterations
           do itheta = 1,ntheta_coil           
              ! Compute r:
              call regcoil_compute_offset_surface_xyz_of_thetazeta(prob, theta_plasma_corresponding_to_theta_coil(itheta), &
                   prob%coil%zetal_coil(izeta),x_new,y_new,z_new,separation)
              R_slice(itheta) = sqrt(x_new * x_new + y_new * y_new)
              z_slice(itheta) = z_new
           end do

           if (geometry_option_coil==2) exit ! No need to iterate if we use the original theta coordinate.

           dl(1:Ntheta_coil-1) = sqrt((R_slice(2:Ntheta_coil) - R_slice(1:Ntheta_coil-1)) * (R_slice(2:Ntheta_coil) - R_slice(1:Ntheta_coil-1)) &
                + (z_slice(2:Ntheta_coil) - z_slice(1:Ntheta_coil-1)) * (z_slice(2:Ntheta_coil) - z_slice(1:Ntheta_coil-1)))
           ! Handle endpoint:
           dl(Ntheta_coil) = sqrt((R_slice(Ntheta_coil) - R_slice(1)) * (R_slice(Ntheta_coil) - R_slice(1)) &
                + (z_slice(Ntheta_coil) - z_slice(1)) * (z_slice(Ntheta_coil) - z_slice(1)))

           relative_variation_in_dl = (maxval(dl) - minval(dl)) / (sum(dl) / Ntheta_coil)
           if (relative_variation_in_dl < constant_arclength_tolerance) then
              !print *,omp_get_thread_num(),izeta," Tolerance achieved."
              exit
           end if

           constant_arclength_theta_on_old_theta_grid(1) = 0
           do itheta = 1, Ntheta_coil
              constant_arclength_theta_on_old_theta_grid(itheta+1) = constant_arclength_theta_on_old_theta_grid(itheta) + dl(itheta)
           end do
           constant_arclength_theta_on_old_theta_grid = constant_arclength_theta_on_old_theta_grid * (2 * pi / constant_arclength_theta_on_old_theta_grid(Ntheta_coil+1))

           call new_periodic_spline(Ntheta_coil, constant_arclength_theta_on_old_theta_grid(1:Ntheta_coil), &
                theta_plasma_corresponding_to_theta_coil - constant_arclength_theta_on_old_theta_grid(1:Ntheta_coil), 2*pi, theta_spline)
           ! In the line above, we subtract constant_arclength_theta_on_old_theta_grid so the function that is interpolated is periodic.

           do itheta = 1, Ntheta_coil
              theta_plasma_corresponding_to_theta_coil(itheta) = periodic_splint(prob%coil%theta_coil(itheta), theta_spline)
           end do
           theta_plasma_corresponding_to_theta_coil = theta_plasma_corresponding_to_theta_coil + prob%coil%theta_coil ! Here we add back the secular term.

           call delete_periodic_spline(theta_spline)
              
        end do

        if (iteration >= max_iterations) then
           print "(a,i3,a,i4,a,i5,a)","*** WARNING: constant-arclength conversion did not converge after",max_iterations," iterations. (proc=",omp_get_thread_num()," izeta=",izeta,")"
        end if

        major_R_coil(:,izeta) = R_slice
        z_coil(:,izeta) = z_slice
        
     end do
     !$OMP END DO

     deallocate(theta_plasma_corresponding_to_theta_coil,dl,R_slice,z_slice,constant_arclength_theta_on_old_theta_grid)

     !$OMP END PARALLEL


     call system_clock(toc1)
     if (verbose) print *,"  Computing offset points:",real(toc1-tic1)/countrate,"sec"

     ! Fourier transform the result.
     ! This is not a rate-limiting step, so for clarity of code, we don't bother with an FFT.
     call system_clock(tic1)
     mpol_coil = min(ntheta_coil / 2, max_mpol_coil)
     ntor_coil = min( nzeta_coil / 2, max_ntor_coil)
     if (verbose) print "(a,i3,a,i3)", "   Representing the offset surface with mpol=",mpol_coil,", ntor=",ntor_coil
     call regcoil_init_Fourier_modes(mpol_coil, ntor_coil, prob%coil%mnmax_coil, prob%coil%xm_coil, prob%coil%xn_coil, .true.)
     prob%coil%xn_coil = prob%coil%xn_coil * nfp
     allocate(prob%coil%rmnc_coil(mnmax_coil))
     allocate(prob%coil%rmns_coil(mnmax_coil))
     allocate(prob%coil%zmnc_coil(mnmax_coil))
     allocate(prob%coil%zmns_coil(mnmax_coil))
     prob%coil%rmnc_coil = 0
     prob%coil%rmns_coil = 0
     prob%coil%zmnc_coil = 0
     prob%coil%zmns_coil = 0
     factor = (2.0d+0) / (ntheta_coil * nzeta_coil)
     do izeta = 1, nzeta_coil
        do itheta = 1, ntheta_coil
           do j = 2, mnmax_coil
              angle = prob%coil%xm_coil(j) * prob%coil%theta_coil(itheta) - prob%coil%xn_coil(j) * prob%coil%zeta_coil(izeta)
              sinangle = sin(angle)
              cosangle = cos(angle)
              factor2 = factor
              ! The next 2 lines ensure inverse Fourier transform(Fourier transform) = identity
              if (mod(ntheta_coil,2) == 0 .and.     prob%coil%xm_coil(j)  ==    (ntheta_coil/2)) factor2 = factor2 / 2
              if (mod( nzeta_coil,2) == 0 .and. abs(prob%coil%xn_coil(j)) == nfp*(nzeta_coil/2)) factor2 = factor2 / 2
              prob%coil%rmnc_coil(j) = prob%coil%rmnc_coil(j) + major_R_coil(itheta, izeta) * cosangle * factor2
              prob%coil%rmns_coil(j) = prob%coil%rmns_coil(j) + major_R_coil(itheta, izeta) * sinangle * factor2
              prob%coil%zmnc_coil(j) = prob%coil%zmnc_coil(j) + z_coil(itheta, izeta) * cosangle * factor2
              prob%coil%zmns_coil(j) = prob%coil%zmns_coil(j) + z_coil(itheta, izeta) * sinangle * factor2
           end do
        end do
     end do
     prob%coil%rmnc_coil(1) = sum(major_R_coil) / (ntheta_coil * nzeta_coil)
     prob%coil%zmnc_coil(1) = sum(z_coil) / (ntheta_coil * nzeta_coil)

     if (.not. lasym) then
        prob%coil%rmns_coil = 0
        prob%coil%zmnc_coil = 0
     end if

     call regcoil_filter_coil_surface(prob)

     deallocate(major_R_coil, z_coil)
     call system_clock(toc1)
     if (verbose) print *,"  Fourier transform:",real(toc1-tic1)/countrate,"sec"
     
     if (verbose) call regcoil_write_nescin(prob)
     call system_clock(tic1)
     if (verbose) print *,"  Write nescin file:",real(tic1-toc1)/countrate,"sec"

  case (3)
     if (verbose) print *,"  Reading coil surface from nescin file ",trim(nescin_filename)
     call regcoil_read_nescin(prob)
     call regcoil_filter_coil_surface(prob)
     
  case default
     print *,"Invalid setting for geometry_option_coil: ",geometry_option_coil
     stop
  end select

  ! End of code that is specific to each value of geometry_option_coil.
  ! Now comes the code that is common to all settings of geometry_option_coil.

  ! First dimension is the Cartesian component x, y, or z.
  if (allocated(prob%coil%r_coil)) deallocate(prob%coil%r_coil)
  allocate(prob%coil%r_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 4'

  if (allocated(prob%coil%drdtheta_coil)) deallocate(prob%coil%drdtheta_coil)
  allocate(prob%coil%drdtheta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 5'

  if (allocated(prob%coil%drdzeta_coil)) deallocate(prob%coil%drdzeta_coil)
  allocate(prob%coil%drdzeta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 6'

  if (allocated(prob%coil%normal_coil)) deallocate(prob%coil%normal_coil)
  allocate(prob%coil%normal_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 7'

  if (allocated(prob%coil%d2rdtheta2_coil)) deallocate(prob%coil%d2rdtheta2_coil)
  allocate(prob%coil%d2rdtheta2_coil(3,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 8'

  if (allocated(prob%coil%d2rdthetadzeta_coil)) deallocate(prob%coil%d2rdthetadzeta_coil)
  allocate(prob%coil%d2rdthetadzeta_coil(3,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 9'

  if (allocated(prob%coil%d2rdzeta2_coil)) deallocate(prob%coil%d2rdzeta2_coil)
  allocate(prob%coil%d2rdzeta2_coil(3,ntheta_coil,nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 10'

  ! Convert Fourier representation of the coil winding surface to Cartesian coordinates:
  call regcoil_evaluate_coil_surface(prob)


  call system_clock(toc)
  if (verbose) print *,"Done initializing coil surface. Took ",real(toc-tic)/countrate," sec."
  

  end associate
end subroutine regcoil_init_coil_surface

subroutine regcoil_filter_coil_surface(prob)

  use regcoil_variables, only: regcoil_t

  implicit none


  type(regcoil_t), intent(inout) :: prob
  integer :: j

  ! Filter out high frequencies, if desired:
    associate ( &
       nfp => prob%plasma%nfp, &
       mnmax_coil => prob%coil%mnmax_coil, &
       mpol_coil_filter => prob%coil%mpol_coil_filter, &
       ntor_coil_filter => prob%coil%ntor_coil_filter &
       )
  do j = 1, mnmax_coil
     if (abs(prob%coil%xm_coil(j)) > mpol_coil_filter .or. abs(prob%coil%xn_coil(j)) > ntor_coil_filter*nfp) then
        prob%coil%rmnc_coil(j) = 0
        prob%coil%rmns_coil(j) = 0
        prob%coil%zmnc_coil(j) = 0
        prob%coil%zmns_coil(j) = 0
     end if
  end do


  end associate
end subroutine regcoil_filter_coil_surface
