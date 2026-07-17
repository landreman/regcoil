! C-compatible API for the Python extension (Phase 4).
! Problem state still lives in regcoil_variables module globals (Phase 5 removes them).

module regcoil_c_api
  use iso_c_binding
  implicit none

contains

  subroutine regcoil_c_set_verbose(flag) bind(C, name="regcoil_c_set_verbose")
    use regcoil_variables, only: verbose
    integer(c_int), value, intent(in) :: flag
    verbose = (flag /= 0)
  end subroutine regcoil_c_set_verbose

  function regcoil_c_nlambda() result(n) bind(C, name="regcoil_c_nlambda")
    use regcoil_variables, only: nlambda
    integer(c_int) :: n
    n = int(nlambda, kind=c_int)
  end function regcoil_c_nlambda

  function regcoil_c_setup(c_path) result(ierr) bind(C, name="regcoil_c_setup")
    use regcoil_init_plasma_mod
    character(kind=c_char), intent(in) :: c_path(*)
    integer(c_int) :: ierr
    character(len=1024) :: path
    integer :: ios

    ierr = 0
    call c_to_f_string(c_path, path)

    call regcoil_read_input_file(trim(path), ios)
    if (ios /= 0) then
       ierr = int(ios, kind=c_int)
       return
    end if

    call regcoil_validate_input()
    call regcoil_compute_lambda()
    call regcoil_init_plasma()
    call regcoil_init_coil_surface()
    call regcoil_read_bnorm()
    call regcoil_build_matrices()
    call regcoil_prepare_solve()
  end function regcoil_c_setup

  function regcoil_c_build_matrices() result(ierr) bind(C, name="regcoil_c_build_matrices")
    integer(c_int) :: ierr
    ierr = 0
    call regcoil_build_matrices()
  end function regcoil_c_build_matrices

  function regcoil_c_prepare_solve() result(ierr) bind(C, name="regcoil_c_prepare_solve")
    integer(c_int) :: ierr
    ierr = 0
    call regcoil_prepare_solve()
  end function regcoil_c_prepare_solve

  function regcoil_c_solve_ilambda(ilambda_c, out_chi2_b, out_chi2_k, out_max_bn, out_max_k) result(ierr) &
       bind(C, name="regcoil_c_solve_ilambda")
    use regcoil_variables, only: nlambda, chi2_B, chi2_K, max_Bnormal, max_K, lambda
    integer(c_int), value, intent(in) :: ilambda_c
    real(c_double), intent(out) :: out_chi2_b, out_chi2_k, out_max_bn, out_max_k
    integer(c_int) :: ierr
    integer :: ilambda

    ierr = 0
    ilambda = int(ilambda_c)
    if (.not. allocated(lambda)) then
       ierr = 1
       return
    end if
    if (ilambda < 1 .or. ilambda > nlambda) then
       ierr = 2
       return
    end if

    call regcoil_solve(ilambda)
    out_chi2_b = real(chi2_B(ilambda), kind=c_double)
    out_chi2_k = real(chi2_K(ilambda), kind=c_double)
    out_max_bn = real(max_Bnormal(ilambda), kind=c_double)
    out_max_k = real(max_K(ilambda), kind=c_double)
  end function regcoil_c_solve_ilambda

  function regcoil_c_solve_lambda(lam, out_chi2_b, out_chi2_k, out_max_bn, out_max_k) result(ierr) &
       bind(C, name="regcoil_c_solve_lambda")
    use regcoil_variables, only: lambda, chi2_B, chi2_K, max_Bnormal, max_K
    real(c_double), value, intent(in) :: lam
    real(c_double), intent(out) :: out_chi2_b, out_chi2_k, out_max_bn, out_max_k
    integer(c_int) :: ierr

    ierr = 0
    if (.not. allocated(lambda)) then
       ierr = 1
       return
    end if

    ! One-λ solve uses the first diagnostic slot (globals still present).
    lambda(1) = real(lam, kind=kind(lambda))
    call regcoil_solve(1)
    out_chi2_b = real(chi2_B(1), kind=c_double)
    out_chi2_k = real(chi2_K(1), kind=c_double)
    out_max_bn = real(max_Bnormal(1), kind=c_double)
    out_max_k = real(max_K(1), kind=c_double)
  end function regcoil_c_solve_lambda

  subroutine c_to_f_string(c_path, f_path)
    character(kind=c_char), intent(in) :: c_path(*)
    character(len=*), intent(out) :: f_path
    integer :: i

    f_path = ''
    do i = 1, len(f_path)
       if (c_path(i) == c_null_char) exit
       f_path(i:i) = c_path(i)
    end do
  end subroutine c_to_f_string

end module regcoil_c_api
