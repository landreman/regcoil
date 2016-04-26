subroutine constant_v_coil_field

  use global_variables
  use stel_constants
  use stel_kinds
  use omp_lib

  implicit none

  integer :: iv_1, iu_1, iv_2, iu_2, l_2, ivl_2, index_1
  real(DP) :: x, y, z, dx, dy, dz, dr2, dr32

  integer :: iflag
  real(dp) :: R2, temp, factors, prefactor, du, dv
  integer :: tic, toc, countrate, tic1, toc1
  real(dp), dimension(:), allocatable :: tempVec

  ! There are 2 differences between tempVec and Bnormal_from_const_v_coils_uv.
  ! 1) tempVec is 1D whereas Bnormal_from_const_v_coils_uv is 2D.
  ! 2) tempVec has an extra factor of |N| = norm_normal_plasma for the area integration.

  call system_clock(tic, countrate)
  print *,"Computing normal field from constant-v coils."
  
  allocate(Bnormal_from_const_v_coils(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_const_v_coils_inductance(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_const_v_coils_transfer(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_const_v_coils_uv(nu_plasma,nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(tempVec(nu_plasma*nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  du = u_plasma(2)-u_plasma(1)
  dv = v_plasma(2)-v_plasma(1)
  factors = nfp * du * dv

  Bnormal_from_const_v_coils_uv = 0
  call system_clock(tic1,countrate)

  !$OMP PARALLEL

  !$OMP MASTER
  print *,"  Number of OpenMP threads:",omp_get_num_threads()
  !$OMP END MASTER

  !$OMP DO PRIVATE(x,y,z,ivl_2,dx,dy,dz,dr2,dr32)
  do iv_1 = 1, nv_plasma
     do iu_1 = 1, nu_plasma
        x = r_plasma(1,iu_1,iv_1)
        y = r_plasma(2,iu_1,iv_1)
        z = r_plasma(3,iu_1,iv_1)
        do iv_2 = 1, nv_middle
           do iu_2 = 1, nu_middle
              do l_2 = 0, (nfp-1)
                 ivl_2 = iv_2 + l_2*nv_middle
                 dx = x - r_middle(1,iu_2,ivl_2)
                 dy = y - r_middle(2,iu_2,ivl_2)
                 dz = z - r_middle(3,iu_2,ivl_2)
                 
                 dr2 = dx*dx + dy*dy + dz*dz
                 dr32 = dr2*sqrt(dr2)

                 Bnormal_from_const_v_coils_uv(iu_1,iv_1) = Bnormal_from_const_v_coils_uv(iu_1,iv_1) +  &
                      (drdu_middle(1,iu_2,ivl_2) * dy * normal_plasma(3,iu_1,iv_1) + &
                      drdu_middle(2,iu_2,ivl_2) * dz * normal_plasma(1,iu_1,iv_1) + &
                      drdu_middle(3,iu_2,ivl_2) * dx * normal_plasma(2,iu_1,iv_1)  &
                      - drdu_middle(3,iu_2,ivl_2) * dy * normal_plasma(1,iu_1,iv_1) &
                      - drdu_middle(1,iu_2,ivl_2) * dz * normal_plasma(2,iu_1,iv_1) &
                      - drdu_middle(2,iu_2,ivl_2) * dx * normal_plasma(3,iu_1,iv_1) ) / dr32
                 ! We still need to divide by norm_normal_plasma.
                 ! This will be done 6 lines below.
              end do
           end do
        end do
        index_1 = (iv_1-1)*nu_plasma + iu_1
        tempVec(index_1) = factors*Bnormal_from_const_v_coils_uv(iu_1,iv_1)
        Bnormal_from_const_v_coils_uv(iu_1,iv_1) = Bnormal_from_const_v_coils_uv(iu_1,iv_1) / norm_normal_plasma(iu_1,iv_1)
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Factors in the next executable line:
  ! du*dv comes from the area integrals.
  ! mu0/(4pi) comes from the factor in front of the Biot-Savart law.
  ! -net_poloidal_current_Amperes/nfp comes from the factor multiplying v in eq (4) of the Nucl Fusion paper on NESCOIL.
  prefactor = (du * dv) * (mu0 / (4*pi)) * (-net_poloidal_current_Amperes / nfp)
  tempVec = tempVec * prefactor
  Bnormal_from_const_v_coils_uv = Bnormal_from_const_v_coils_uv * prefactor

  call system_clock(toc1)
  print *,"  assembly:",real(toc1-tic1)/countrate,"sec."

  call system_clock(tic1)
  ! This next line could be optimized using BLAS
  Bnormal_from_const_v_coils = matmul(tempVec,basis_functions_plasma)
  call system_clock(toc1)
  print *,"  matmul:",real(toc1-tic1)/countrate,"sec."

  deallocate(tempVec)
  call system_clock(toc)
  print *,"Done computing normal field from constant-v coils. Took ",real(toc-tic)/countrate," sec."

end subroutine constant_v_coil_field
