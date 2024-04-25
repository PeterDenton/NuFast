module Parameters
    implicit none
    real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
    real, parameter :: oneThird = 1. / 3
    real, parameter :: YerhoE2a = 1.52e-4
    real, parameter :: eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4
end module Parameters

! Compute the speed n times and calculate the mean over those n times
! N_Newton is used for NuFast only, ignored otherwise
real function Speed_Helper(Probability_Calculator, n, N_Newton)
    use Parameters
    implicit none
    interface
        subroutine Probability_Calculator(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
            real, intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye
            integer, intent(in) :: N_Newton
            real, intent(out) :: probs_returned(3, 3)
        end subroutine Probability_Calculator
    end interface
    integer, intent(in) :: n, N_Newton
    integer :: start_time, rate, end_time, i
    real :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, Emin, Emax, Estep, rho, Ye, probs_returned(3, 3)

    ! ------------------------------------- !
    ! Set the vacuum oscillation parameters !
    ! ------------------------------------- !
    s12sq = 0.31
    s13sq = 0.02
    s23sq = 0.55
    delta = -0.7 * PI
    Dmsq21 = 7.5e-5 ! eV^2
    Dmsq31 = 2.5e-3 ! eV^2

    ! ------------------------------- !
    ! Set the experimental parameters !
    ! ------------------------------- !
    L = 1300 ! km
    Emin = 0.5 ! GeV
    Emax = 5
    Estep = (Emax - Emin) / real(n)
    rho = 3 ! g/cc
    Ye = 0.5

    call system_clock(start_time, rate)
    do i = 1, n
        E = Emin + (i - 1) * Estep
        call Probability_Calculator(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
    end do ! i n
    call system_clock(end_time)
    Speed_Helper = real(end_time - start_time) / real(rate) / n
end function Speed_Helper

! Compute the speed many times and calculate the mean and standard deviation
! m is outer loop, n, is inner loop
! N_Newton is used for NuFast only, ignored otherwise
subroutine Speed(Probability_Calculator, m, n, N_Newton, mean, std)
    implicit none
    interface
        subroutine Probability_Calculator(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
            real, intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, Ye, rho
            integer, intent(in) :: N_Newton
            real, intent(out) :: probs_returned(3, 3)
        end subroutine Probability_Calculator
    end interface
    integer, intent(in) :: m, n, N_Newton
    real, intent(out) :: mean, std
    real, external :: Speed_Helper
    real :: s, speed_sum, speedsq_sum
    integer :: i

    speed_sum = 0
    speedsq_sum = 0
    do i = 1, m
        s = Speed_Helper(Probability_Calculator, n, N_Newton) * 1e9 ! in ns
        speed_sum = speed_sum + s
        speedsq_sum = speedsq_sum + s * s
    end do ! i, m
    mean = speed_sum / m
    std = sqrt(speedsq_sum / m - mean * mean)
end subroutine Speed

! The main program
! Computes the 9 oscillation probabilities for some
! oscillation parameters at some level of precision
! and prints them
program Benchmark
    use Parameters
    implicit none
    external :: Probability_Matter_LBL, Probability_Matter_LBL_Exact_Cubic, Probability_Vacuum_LBL, Speed
    real :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, Ye, rho, Emin, Emax, Estep
    real :: probs_returned(3, 3), probs_returned_Exact_Cubic(3, 3)
    real :: mean, std
    integer :: N_Newton, i, alpha, beta, m, n, dataf

    ! ------------------------------------- !
    ! Set the vacuum oscillation parameters !
    ! ------------------------------------- !
    s12sq = 0.31
    s13sq = 0.02
    s23sq = 0.55
    delta = -0.7 * PI
    Dmsq21 = 7.5e-5  ! eV^2
    Dmsq31 = +2.5e-3 ! eV^2

    ! ------------------------------- !
    ! Set the experimental parameters !
    ! ------------------------------- !
    L = 1300 ! km
    E = 2.5 ! GeV
    rho = 3 ! g/cc
    Ye = 0.5

    N_Newton = 0 ! <-- Sets the precision. 0 is close to the single precision limit and is better than DUNE/HK in the high statistics regime. Increasig N_Newton improves the precision at a modest computational cost

    ! ------------------------------------------ !
    ! Calculate all 9 oscillationa probabilities !
    ! ------------------------------------------ !
    call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
    call Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &
                                            probs_returned_Exact_Cubic)

    ! ------------------------------- !
    ! Print out the fractional errors !
    ! ------------------------------- !
    write (*,"(A,F8.1,A,F5.2,A,F5.2,A,F5.2,A,I1)") "L = ", L, " E = ", E, " rho = ", rho, " Ye = ", Ye, " N_Newton = ", N_Newton
    write (*,*) "Precision:"
    write (*,"(A6,A5,A10)") "alpha", "beta", "DeltaP/P"
    do alpha = 1, 3
        do beta = 1, 3
            write (*,"(I6,I5,E10.2)") alpha, beta, &
                ABS(probs_returned(alpha, beta) - probs_returned_Exact_Cubic(alpha, beta)) / probs_returned_Exact_Cubic(alpha, beta)
        end do ! beta, 1, 3
    end do ! alpha, 1, 3

    ! ----------------------------------- !
    ! Write the fractional errors to file !
    ! ----------------------------------- !
    n = 1e3
    ! Accelerator !
    Emin = 0.5
    Emax = 5
    Estep = (Emax - Emin) / n
    open (newunit = dataf, file = "data/Precision_accelerator.txt", status = "replace", action = "write")
    write (dataf, *) L, rho, Ye, N_Newton
    do i = 0, n
        E = Emin + i * Estep
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        call Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &
                                                probs_returned_Exact_Cubic)
        write (dataf, *) E, ABS(probs_returned(2, 2) - probs_returned_Exact_Cubic(2, 2)) / probs_returned_Exact_Cubic(2, 2), &
                 ABS(probs_returned(2, 1) - probs_returned_Exact_Cubic(2, 1)) / probs_returned_Exact_Cubic(2, 1)
    end do ! i
    close (dataf)

    ! Reactor !
    L = 50
    rho = 2.6
    Emin = 1e-3
    Emax = 5e-3
    Estep = (Emax - Emin) / n
    open (newunit = dataf, file = "data/Precision_reactor.txt", status = "replace", action = "write")
    write (dataf, *) L, rho, Ye, N_Newton
    do i = 0, n
        E = Emin + i * Estep
        call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
        call Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &
                                                probs_returned_Exact_Cubic)
        write (dataf, *) E, ABS(probs_returned(1, 1) - probs_returned_Exact_Cubic(1, 1)) / probs_returned_Exact_Cubic(1, 1)
    end do ! i
    close (dataf)

    ! ---------- !
    ! Speed test !
    ! ---------- !
    open (newunit = dataf, file = "data/Speed.txt", status = "replace", action = "write")
    m = int(1e3)
    n = int(1e6)

    write (*,*) ""
    write (*,*) "Speed (this takes a moment):"

    ! Vacuum probability
    write (*,"(A)") "Vacuum (1/6)"
    call Speed(Probability_Vacuum_LBL, m, n, 0, mean, std)
    write (*,"(F7.1,A,F7.1,A)") mean, " +- ", std, " ns"
    write (dataf,"(A,F7.1,A,F7.1,A)") "Vacuum", mean, " +- ", std, " ns"

    ! NuFast at various N_Newton levels
    do N_Newton = 0, 3
        write (*,"(A,I2,A,I1,A)") "N_Newton =", N_Newton, " (", N_Newton + 2, "/6)"
        call Speed(Probability_Matter_LBL, m, n, N_Newton, mean, std)
        write (*,"(F7.1,A,F7.1,A)") mean, " +- ", std, " ns"
        write (dataf,"(A,I1,F7.1,A,F7.1,A)") "N_Newton=", N_Newton, mean, " +- ", std, " ns"
    end do ! N_Newton, 0, 3

    ! Exact with cubic from Cardano/ZS
    write (*,"(A)") "Exact_Cubic (6/6)"
    call Speed(Probability_Matter_LBL_Exact_Cubic, m, n, 0, mean, std)
    write (*,"(F7.1,A,F7.1,A)") mean, " +- ", std, " ns"
    write (dataf,"(A,F7.1,A,F7.1,A)") "Exact_Cubic", mean, " +- ", std, " ns"
    close (dataf)
end program Benchmark
