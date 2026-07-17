! C-compatible API for the Python extension (Phase 5: opaque instance handles).

module regcoil_c_api
  use iso_c_binding
  use regcoil_variables, only: regcoil_t
  use regcoil_init_plasma_mod
  implicit none

contains

  function regcoil_c_create() result(handle) bind(C, name="regcoil_c_create")
    type(c_ptr) :: handle
    type(regcoil_t), pointer :: p
    allocate(p)
    handle = c_loc(p)
  end function regcoil_c_create

  subroutine regcoil_c_destroy(handle) bind(C, name="regcoil_c_destroy")
    type(c_ptr), value, intent(in) :: handle
    type(regcoil_t), pointer :: p
    if (.not. c_associated(handle)) return
    call c_f_pointer(handle, p)
    deallocate(p)
  end subroutine regcoil_c_destroy

  subroutine regcoil_c_set_verbose(handle, flag) bind(C, name="regcoil_c_set_verbose")
    type(c_ptr), value, intent(in) :: handle
    integer(c_int), value, intent(in) :: flag
    type(regcoil_t), pointer :: p
    call c_f_pointer(handle, p)
    p%input%verbose = (flag /= 0)
  end subroutine regcoil_c_set_verbose

  function regcoil_c_nlambda(handle) result(n) bind(C, name="regcoil_c_nlambda")
    type(c_ptr), value, intent(in) :: handle
    integer(c_int) :: n
    type(regcoil_t), pointer :: p
    call c_f_pointer(handle, p)
    n = int(p%input%nlambda, kind=c_int)
  end function regcoil_c_nlambda

  function regcoil_c_setup(handle, c_path) result(ierr) bind(C, name="regcoil_c_setup")
    type(c_ptr), value, intent(in) :: handle
    character(kind=c_char), intent(in) :: c_path(*)
    integer(c_int) :: ierr
    type(regcoil_t), pointer :: p
    character(len=1024) :: path
    integer :: ios

    ierr = 0
    call c_f_pointer(handle, p)
    call c_to_f_string(c_path, path)

    call regcoil_read_input_file(p, trim(path), ios)
    if (ios /= 0) then
       ierr = int(ios, kind=c_int)
       return
    end if

    call regcoil_validate_input(p)
    call regcoil_compute_lambda(p)
    call regcoil_init_plasma(p)
    call regcoil_init_coil_surface(p)
    call regcoil_read_bnorm(p)
    call regcoil_build_matrices(p)
    call regcoil_prepare_solve(p)
  end function regcoil_c_setup

  function regcoil_c_build_matrices(handle) result(ierr) bind(C, name="regcoil_c_build_matrices")
    type(c_ptr), value, intent(in) :: handle
    integer(c_int) :: ierr
    type(regcoil_t), pointer :: p
    ierr = 0
    call c_f_pointer(handle, p)
    call regcoil_build_matrices(p)
  end function regcoil_c_build_matrices

  function regcoil_c_prepare_solve(handle) result(ierr) bind(C, name="regcoil_c_prepare_solve")
    type(c_ptr), value, intent(in) :: handle
    integer(c_int) :: ierr
    type(regcoil_t), pointer :: p
    ierr = 0
    call c_f_pointer(handle, p)
    call regcoil_prepare_solve(p)
  end function regcoil_c_prepare_solve

  function regcoil_c_solve_ilambda(handle, ilambda_c, out_chi2_b, out_chi2_k, out_max_bn, out_max_k) result(ierr) &
       bind(C, name="regcoil_c_solve_ilambda")
    type(c_ptr), value, intent(in) :: handle
    integer(c_int), value, intent(in) :: ilambda_c
    real(c_double), intent(out) :: out_chi2_b, out_chi2_k, out_max_bn, out_max_k
    integer(c_int) :: ierr
    type(regcoil_t), pointer :: p
    integer :: ilambda

    ierr = 0
    call c_f_pointer(handle, p)
    ilambda = int(ilambda_c)
    if (.not. allocated(p%input%lambda)) then
       ierr = 1
       return
    end if
    if (ilambda < 1 .or. ilambda > p%input%nlambda) then
       ierr = 2
       return
    end if

    call regcoil_solve(p, ilambda)
    out_chi2_b = real(p%output%chi2_B(ilambda), kind=c_double)
    out_chi2_k = real(p%output%chi2_K(ilambda), kind=c_double)
    out_max_bn = real(p%output%max_Bnormal(ilambda), kind=c_double)
    out_max_k = real(p%output%max_K(ilambda), kind=c_double)
  end function regcoil_c_solve_ilambda

  function regcoil_c_solve_lambda(handle, lam, out_chi2_b, out_chi2_k, out_max_bn, out_max_k) result(ierr) &
       bind(C, name="regcoil_c_solve_lambda")
    type(c_ptr), value, intent(in) :: handle
    real(c_double), value, intent(in) :: lam
    real(c_double), intent(out) :: out_chi2_b, out_chi2_k, out_max_bn, out_max_k
    integer(c_int) :: ierr
    type(regcoil_t), pointer :: p

    ierr = 0
    call c_f_pointer(handle, p)
    if (.not. allocated(p%input%lambda)) then
       ierr = 1
       return
    end if

    p%input%lambda(1) = real(lam, kind=kind(p%input%lambda))
    call regcoil_solve(p, 1)
    out_chi2_b = real(p%output%chi2_B(1), kind=c_double)
    out_chi2_k = real(p%output%chi2_K(1), kind=c_double)
    out_max_bn = real(p%output%max_Bnormal(1), kind=c_double)
    out_max_k = real(p%output%max_K(1), kind=c_double)
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
