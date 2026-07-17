subroutine regcoil_solve(prob, ilambda)

  use regcoil_variables, only: regcoil_t

  implicit none

  type(regcoil_t), intent(inout) :: prob
  integer, intent(in) :: ilambda
  integer :: iflag, tic, toc, countrate

  associate ( &
       verbose => prob%input%verbose, &
       lambda => prob%input%lambda, &
       matrix_B => prob%work%matrix_B, &
       matrix_regularization => prob%work%matrix_regularization, &
       RHS_B => prob%work%RHS_B, &
       RHS_regularization => prob%work%RHS_regularization, &
       matrix => prob%work%matrix, &
       RHS => prob%work%RHS, &
       solution => prob%work%solution, &
       num_basis_functions => prob%work%num_basis_functions, &
       LAPACK_IPIV => prob%work%LAPACK_IPIV, &
       LAPACK_WORK => prob%work%LAPACK_WORK, &
       LAPACK_LWORK => prob%work%LAPACK_LWORK, &
       LAPACK_INFO => prob%work%LAPACK_INFO &
       )

    call system_clock(tic,countrate)

    ! The scaling of the terms below by 1/(1+lambda) ensures that matrix and RHS are O(1) regardless of whether lambda is >> 1 or << 1.
    matrix = (1 / (1 + lambda(ilambda))) * matrix_B + (lambda(ilambda) / (1 + lambda(ilambda))) * matrix_regularization
    RHS    = (1 / (1 + lambda(ilambda))) *    RHS_B + (lambda(ilambda) / (1 + lambda(ilambda))) *    RHS_regularization

    call system_clock(toc)
    if (verbose) print *,"  Additions: ",real(toc-tic)/countrate," sec."
    call system_clock(tic)

    ! Compute solution = matrix \ RHS.
    ! Use LAPACK's DSYSV since matrix is symmetric.
    ! Note: RHS will be over-written with the solution.
    call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, LAPACK_IPIV, RHS, num_basis_functions, LAPACK_WORK, LAPACK_LWORK, LAPACK_INFO)
    if (LAPACK_INFO /= 0) then
       print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", LAPACK_INFO
    end if
    solution = RHS

    call system_clock(toc)
    if (verbose) print *,"  DSYSV: ",real(toc-tic)/countrate," sec."
    call system_clock(tic)

    call regcoil_diagnostics(prob, ilambda)

  end associate

end subroutine regcoil_solve
