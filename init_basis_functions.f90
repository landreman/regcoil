subroutine init_basis_functions()
  
  use global_variables, only: symmetry_option, nfp, mpol_coil, ntor_coil, xm_coil, xn_coil, mnmax_coil, num_basis_functions,&
       basis_functions, ntheta_coil, nzeta_coil, theta_coil, zeta_coil
  use init_Fourier_modes_mod
  use stel_kinds
  use stel_constants
  
  implicit none
  
  integer :: index, imn, iflag
  integer :: tic, toc, countrate
  integer :: whichSymmetry, minSymmetry, maxSymmetry, offset
  
  call system_clock(tic,countrate)
  print *,"Initializing basis functions"
  
  ! Initialize Fourier arrays
  call init_Fourier_modes(mpol_coil, ntor_coil, mnmax_coil, xm_coil, xn_coil)
  xn_coil = xn_coil * nfp
  
  select case (symmetry_option)
  case (1,2)
     num_basis_functions = mnmax_coil
  case (3)
     num_basis_functions = mnmax_coil * 2
  case default
     print *,"Error! Invalid setting for symmetry_option:",symmetry_option
     stop
  end select
  
  allocate(basis_functions(ntheta_coil*nzeta_coil, num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  select case (symmetry_option)
  case (1)
     minSymmetry = 1
     maxSymmetry = 1
  case (2)
     minSymmetry = 2
     maxSymmetry = 2
  case (3)
     minSymmetry = 1
     maxSymmetry = 2
  end select
  
  
  ! This loop could be made faster
  ! by using the sum-angle trig identities and pretabulating the trig functions.
  ! But these loops are not the rate-limiting step, so I'll use the more transparent direct method here.
  do whichSymmetry = minSymmetry, maxSymmetry
     
     if (whichSymmetry==2 .and. symmetry_option==3) then
        offset = mnmax_coil
     else
        offset = 0
     end if
     
     do itheta = 1, ntheta_coil
        do izeta = 1, nzeta_coil
           index = (izeta-1)*ntheta + itheta
           do imn = 1, mnmax_coil
              if (whichSymmetry==1) then
                 basis_functions(index, imn) = sin(xm_coil(imn)*theta_coil(itheta)-xn_coil(imn)*zeta_coil(izeta))
              else
                 basis_functions(index, imn+offset) = cos(xm_coil(imn)*theta_coil(itheta)-xn_coil(imn)*zeta_coil(izeta))
              end if
           end do
        end do
     end do
  end do
  
  call system_clock(toc)
  print *,"Done. Took",real(toc-tic)/countrate,"sec."
  
  
end subroutine init_basis_functions

  
