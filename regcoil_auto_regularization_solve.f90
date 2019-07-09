subroutine regcoil_auto_regularization_solve()

  use regcoil_variables
  use stel_kinds

  implicit none

  integer :: ilambda
  integer :: stage, next_stage
  logical :: initial_above_target, last_above_target, targeted_quantity_increases_with_lambda
  real(dp) :: Brendt_a, Brendt_b, Brendt_c, Brendt_fa, Brendt_fb, Brendt_fc, Brendt_d, Brendt_e
  real(dp) :: Brendt_p, Brendt_q, Brendt_r, Brendt_s, Brendt_tol1, Brendt_xm, Brendt_EPS, factor

  if (general_option==4) then
     stage = 1
  else
     stage = 10
  end if
  exit_code = -1
  do ilambda = 1,nlambda

     ! First, pick the next value of lambda, depending on what stage we are in in the search algorithm:
     ! ------------------------------------------------------------------------------------------------
     select case (stage)
     case (1)
        ! Initial guess for lambda:
        ! guess lambda = chi^2_B / chi^2_K, where the right hand side is evaluated taking the
        ! single-valued part of the current potential to be 0.

        KDifference_x = d_x !- matmul(f_x, solution)
        KDifference_y = d_y !- matmul(f_y, solution)
        KDifference_z = d_z !- matmul(f_z, solution)
        this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
             / norm_normal_coil
        chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)

        Bnormal_total(:,:,ilambda) = & ! (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) +
             Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

        chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
             * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)
        ! chi2_B, chi2_K, and Bnormal_total for this ilambda will be over-written with the real values below.

        lambda(ilambda) = chi2_B(ilambda) / chi2_K(ilambda) / 1000
        next_stage = 2

     case (2)
        ! Vary lambda by factors of 100 until we bracket the target current density.
        next_stage = 2

        if (initial_above_target) then
           factor = 100
        else
           factor = 0.01
        end if

        if (targeted_quantity_increases_with_lambda) factor = 1/factor

        lambda(ilambda) = lambda(ilambda-1) * factor

     case (3)
        ! Now that we've bracketed the target current density, search for the target using Brendt's algorithm.
        next_stage = 3

        if (abs(Brendt_e) >= Brendt_tol1 .and. abs(Brendt_fa)>abs(Brendt_fb)) then
           ! Attempt inverse quadratic interpolation
           Brendt_s = Brendt_fb / Brendt_fa
           if (Brendt_a == Brendt_c) then
              Brendt_p = 2.0*Brendt_xm*Brendt_s
              Brendt_q = 1.0-Brendt_s
           else
              Brendt_q = Brendt_fa / Brendt_fc
              Brendt_r = Brendt_fb / Brendt_fc
              Brendt_p = Brendt_s*(2.0*Brendt_xm*Brendt_q*(Brendt_q-Brendt_r)-(Brendt_b-Brendt_a)*(Brendt_r-1.0))
              Brendt_q = (Brendt_q-1.0)*(Brendt_r-1.0)*(Brendt_s-1.0)
           end if
           if (Brendt_p > 0) Brendt_q = -Brendt_q ! Check whether in bounds
           Brendt_p = abs(Brendt_p)
           if (2.0*Brendt_p < min(3.0*Brendt_xm*Brendt_q-abs(Brendt_tol1*Brendt_q),abs(Brendt_e*Brendt_q))) then
              ! Accept interpolation
              Brendt_e = Brendt_d
              Brendt_d = Brendt_p / Brendt_q
           else
              ! Interpolation failed, so use bisection
              Brendt_d = Brendt_xm
              Brendt_e = Brendt_d
           end if
        else
           ! Bounds are decreasing too slowly, so use bisection
           Brendt_d = Brendt_xm
           Brendt_e = Brendt_d
        end if
        ! Move last best guess to a.
        Brendt_a = Brendt_b
        Brendt_fa = Brendt_fb
        if (abs(Brendt_d) > Brendt_tol1) then
           ! Evaluate new trial root
           Brendt_b = Brendt_b + Brendt_d
        else
           Brendt_b = Brendt_b + sign(Brendt_tol1,Brendt_xm)
        end if

        lambda(ilambda) = exp(Brendt_b)

     case (10)
        ! Try lambda = infinity
        lambda(ilambda) = 1.0d200
        next_stage = 11

     case (11)
        lambda(ilambda) = 0
        next_stage = 1

     case default
        print *,"Invalid stage in auto_regularization_solve:",stage
        stop
     end select

     if (verbose) print "(a,es10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of at most ",nlambda,")"

     ! Done choosing the next lambda. Now comes the main solve.
     ! ------------------------------------------------------------------------------------------------

     call regcoil_solve(ilambda)

     last_above_target = (target_function(ilambda) > target_value)
     if (stage==1) initial_above_target = last_above_target
     if (stage==2 .and. (last_above_target .neqv. initial_above_target)) then
        ! If we've bracketed the target, move on to stage 3.
        next_stage = 3  
        if (verbose) print *,"Target_value has been bracketed."
     end if
     if (last_above_target .and. (((.not. targeted_quantity_increases_with_lambda) .and. stage==10) &
          .or. (targeted_quantity_increases_with_lambda .and. stage==11))) then
        if (verbose) then
           print *,"*******************************************************************************"
           print *,"*******************************************************************************"
           print *,"Error! The target_value you have set is not achievable because"
           print *,"it is too low."
           print *,"*******************************************************************************"
           print *,"*******************************************************************************"
        end if
        Nlambda = ilambda
        exit_code = -2
        ! Modified to change the behavior when current density target can not be
        ! reached.  Now, the chi2_B returned is the one that is calculated with
        ! infinite regularization 
        ! Put the worst achieved chi2_B into the chi2_B_target variable for STELLOPT.
        chi2_B_target = chi2_B(1)    ! 'Worst' achieved chi2_B
        exit
     end if
     if ((.not. last_above_target) .and. (((.not. targeted_quantity_increases_with_lambda) .and. stage==11) &
          .or. (targeted_quantity_increases_with_lambda .and. stage==10))) then
        if (verbose) then
           print *,"*******************************************************************************"
           print *,"*******************************************************************************"
           print *,"Error! The target_value you have set is not achievable because"
           print *,"it is too high."
           print *,"*******************************************************************************"
           print *,"*******************************************************************************"
        end if
        Nlambda = ilambda
        exit_code = -3
        ! Modified to change the behavior when current density target can not be
        ! reached.  Now, the chi2_B returned is the one that is calculated with
        ! infinite regularization 
        ! Put the worst achieved chi2_B into the chi2_B_target variable for STELLOPT.
        chi2_B_target = chi2_B(1)  ! 'Worst' achieved chi2_B
        exit
     end if
     if (stage==2 .and. next_stage == 3) then
        ! Initialize Brendt's algorithm for root-finding
        Brendt_a = log(lambda(ilambda-1))
        Brendt_b = log(lambda(ilambda))
        Brendt_fa = log(target_function(ilambda-1) / target_value)
        Brendt_fb = log(target_function(ilambda) / target_value)
        Brendt_c = Brendt_b
        Brendt_fc = Brendt_fb
        Brendt_d = Brendt_b - Brendt_a
        Brendt_e = Brendt_d
     end if
     if (stage==3) Brendt_fb = log(target_function(ilambda) / target_value)
     if (next_stage==3) then
        ! Analyze the most recent diagnostics for Brendt's algorithm.
        if ((Brendt_fb > 0 .and. Brendt_fc > 0) .or. (Brendt_fb < 0 .and. Brendt_fc < 0)) then
           Brendt_c = Brendt_a
           Brendt_fc = Brendt_fa
           Brendt_d = Brendt_b - Brendt_a
           Brendt_e = Brendt_d
        end if
        if (abs(Brendt_fc) < abs(Brendt_fb)) then
           Brendt_a = Brendt_b
           Brendt_b = Brendt_c
           Brendt_c = Brendt_a
           Brendt_fa = Brendt_fb
           Brendt_fb = Brendt_fc
           Brendt_fc = Brendt_fa
        end if
        Brendt_EPS = 1d-15
        Brendt_tol1 = 2.0*Brendt_EPS*abs(Brendt_b) + 0.5*lambda_search_tolerance
        Brendt_xm = 0.5*(Brendt_c - Brendt_b)
        if (abs(Brendt_xm) <= Brendt_tol1 .or. (Brendt_fb==0)) then
           ! We met the requested tolerance
           if (verbose) print *,"Requested tolerance has been met."
           exit_code=0
           Nlambda = ilambda
           chi2_B_target = chi2_B(Nlambda)
           exit
        end if
     end if
     stage = next_stage
  end do

 !   print *,"exit_code", exit_code, " chi2_B_target = ", chi2_B_target

  if (exit_code == -1) then
     print *,"*******************************************************************************"
     print *,"*******************************************************************************"
     print *,"The lambda search did not converge within Nlambda iterations!"
     print *,"*******************************************************************************"
     print *,"*******************************************************************************"
  end if

contains
 
  function target_function(jlambda)

    implicit none

    integer, intent(in) :: jlambda
    real(dp) :: target_function

    target_function = 0
    select case (trim(target_option))

    case (target_option_max_K)
       target_function = max_K(jlambda)
       targeted_quantity_increases_with_lambda = .false.

    case (target_option_rms_K)
       target_function = sqrt(chi2_K(jlambda) / area_coil)
       targeted_quantity_increases_with_lambda = .false.

    case (target_option_chi2_K)
       target_function = chi2_K(jlambda)
       targeted_quantity_increases_with_lambda = .false.

    case (target_option_max_Bnormal)
       target_function =  max_Bnormal(jlambda)
       targeted_quantity_increases_with_lambda = .true.

    case (target_option_rms_Bnormal)
       target_function = sqrt(chi2_B(jlambda) / area_plasma)
       targeted_quantity_increases_with_lambda = .true.

    case (target_option_chi2_B)
       target_function = chi2_B(jlambda)
       targeted_quantity_increases_with_lambda = .true.

    case (target_option_max_K_lse)
       target_function = max_K_lse(jlambda)
       targeted_quantity_increases_with_lambda = .false.

    case (target_option_lp_norm_K)
       target_function = lp_norm_K(jlambda)
       targeted_quantity_increases_with_lambda = .false.

    case default
       print *,"Invalid target_option: ",target_option
       stop
    end select

  end function target_function

end subroutine regcoil_auto_regularization_solve
