subroutine Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, empty, probs_returned)
    use Parameters
    implicit none
    real, intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye
    integer, intent(in) :: empty
    real, intent(out) :: probs_returned(3, 3)

    real :: Amatter, c13sq, sind, cosd

    real :: A, B, C, S
    real :: See, Smm, Tee, Tmm, rootAsqm3B
    real :: lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32
    real :: Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq

    real :: sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin, PiDlambdaInv, Xp2, Xp3
    real :: D21, D32
    real :: sinD21, sinD31, sinD32
    real :: Lover4E, Jrr, Jmatter
    real :: Pee, Pme_CPC, Pme_CPV, Pmm


    ! --------------------------------------------------------------------- !
    ! First calculate useful simple functions of the oscillation parameters !
    ! --------------------------------------------------------------------- !
    c13sq = 1 - s13sq

    ! Ueisq's
    Ue2sq = c13sq * s12sq
    Ue3sq = s13sq

    ! Umisq's, Utisq's and Jvac     
    Um3sq = c13sq * s23sq
    ! Um2sq and Ut2sq are used here as temporary variables, will be properly defined later     
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)
      
    Jrr = sqrt(Um2sq * Ut2sq)
    sind = sin(delta)
    cosd = cos(delta)

    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jmatter = 8 * Jrr * c13sq * sind
    Amatter = Ye * rho * E * YerhoE2a ! assume Ye = 0.5

    ! calculate A, B, C, See, Tee, and part of Tmm
    A = Dmsq21 + Dmsq31 ! temporary variable
    See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq
    Tmm = Dmsq21 * Dmsq31 ! using Tmm as a temporary variable      
    Tee = Tmm * (1 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A + Amatter

    ! -------------------------------- !
    ! Get exact lambda3 from the cubic !
    ! -------------------------------- !
    B = Tmm + Amatter * See
    rootAsqm3B = sqrt(A * A - 3 * B)
    S = acos((A * A * A - 4.5 * A* B + 13.5 * C) / (rootAsqm3B * rootAsqm3B * rootAsqm3B))
    if (Dmsq31.LT.0) then
        S = S + 2 * PI ! adjust lambda3 for the IO
    endif
    lambda3 = (A + 2 * rootAsqm3B * cos(S * oneThird)) * oneThird

    ! ------------------- !
    ! Get  Delta lambda's !
    ! ------------------- !
    Dlambda21 = sqrt((A - lambda3) ** 2 - 4 * C / lambda3)
    lambda2 = 0.5 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21

    ! ----------------------- !
    ! Use Rosetta for Ueisq's !
    ! ----------------------- !
    ! denominators      
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31

    ! numerators
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2

    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq
    Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A)

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2

    ! ------------- !
    ! Use NHS for J !
    ! ------------- !
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv

    ! ----------------------- !
    ! Get all elements of Usq !
    ! ----------------------- !
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq

    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq

    ! ----------------------- !
    ! Get the kinematic terms !
    ! ----------------------- !
    Lover4E = eVsqkm_to_GeV_over4 * L / E

    D21 = Dlambda21 * Lover4E
    D32 = Dlambda32 * Lover4E
      
    sinD21 = sin(D21)
    sinD31 = sin(D32 + D21)
    sinD32 = sin(D32)

    triple_sin = sinD21 * sinD31 * sinD32

    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32

    ! ------------------------------------------------------------------- !
    ! Calculate the three necessary probabilities, separating CPC and CPV !
    ! ------------------------------------------------------------------- !
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 &
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 &
            + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    Pme_CPV = -Jmatter * triple_sin

    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 &
                 + Um3sq * Um1sq * sinsqD31_2 &
                 + Um3sq * Um2sq * sinsqD32_2)

    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 &
                 + Ue3sq * Ue1sq * sinsqD31_2 &
                 + Ue3sq * Ue2sq * sinsqD32_2)

    ! ---------------------------- !
    ! Assign all the probabilities !
    ! ---------------------------- !
    probs_returned(1, 1) = Pee                              ! Pee
    probs_returned(1, 2) = Pme_CPC - Pme_CPV                ! Pem
    probs_returned(1, 3) = 1 - Pee - probs_returned(1, 2)   ! Pet

    probs_returned(2, 1) = Pme_CPC + Pme_CPV                ! Pme
    probs_returned(2, 2) = Pmm                              ! Pmm
    probs_returned(2, 3) = 1 - probs_returned(2, 1) - Pmm   ! Pmt

    probs_returned(3, 1) = 1 - Pee - probs_returned(2, 1)
    probs_returned(3, 2) = 1 - probs_returned(1, 2) - Pmm
    probs_returned(3, 3) = 1 - probs_returned(1, 3) - probs_returned(2, 3)
end subroutine Probability_Matter_LBL_Exact_Cubic
