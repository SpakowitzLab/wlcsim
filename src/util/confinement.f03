#include "../defines.inc"
pure function in_confinement(RP, NT, IT1, IT2, wlc_p)
    use params, only : dp, wlcsim_params
    use inputparams, only : MAXPARAMLEN

    implicit none

    integer, intent(in) :: IT1, IT2, NT
    real(dp), intent(in) :: RP(3, NT)
    type(wlcsim_params), intent(in) :: wlc_p
    logical in_confinement
    integer i
    real(dp) rad, length, r2

    in_confinement = .True.
    ! if (WLC_P__CONFINETYPE == 'none') then
    !     in_confinement = .False.
    if (WLC_P__CONFINETYPE == 'platesInZperiodicXY') then
        ! Confinement only in the z-direction
        ! limits: 0 and LBox(3)
        do I = IT1,IT2
            if ((RP(3,I) < 0.0_dp) .or. (RP(3,I) > wlc_p%CONFINEMENTPARAMETER(1))) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'cube') then
        do I = IT1,IT2
            if ((RP(1,I) < 0.0_dp) &
                .or. (RP(1,I) > wlc_p%CONFINEMENTPARAMETER(1)) &
                .or. (RP(2,I) < 0.0_dp) &
                .or. (RP(2,I) > wlc_p%CONFINEMENTPARAMETER(1)) &
                .or. (RP(3,I) < 0.0_dp) &
                .or. (RP(3,I) > wlc_p%CONFINEMENTPARAMETER(1))) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'sphere') then
        ! sphere with given diameter, centered at (r,r,r)
        do I = IT1,IT2
            rad = wlc_p%CONFINEMENTPARAMETER(1)/2
            if ((RP(1,I) - rad)**2 + (RP(2,I) - rad)**2 + &
                (RP(3,I) - rad)**2 > rad) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'ecoli') then
        ! cylinder with hemispherical caps, one tip at origin
        ! full length - lbox(1)/confinementParameter(1)
        ! diameter - lbox(2:3)/confinementParameter(2)
        do I = IT1,IT2
            length = wlc_p%CONFINEMENTPARAMETER(1)
            rad = wlc_p%CONFINEMENTPARAMETER(2)/2
            r2 = RP(2,I)**2 + RP(3,I)**2
            if (r2 > rad &
                .or. RP(1,I) > length &
                .or. RP(1,I) < 0.0_dp) then
                in_confinement = .False.
                return
            elseif (RP(1,I) >= rad .and. RP(1,I) <= length - rad) then
                ! we're inside the main cylinder section, no need to check
                ! intersection with caps
                cycle
            ! if we are inside cap touching origin
            elseif (RP(1,I) < rad) then
                r2 = r2 + (RP(1,I) - rad)**2
            ! we are inside cap far from origin
            else
                r2 = r2 + (RP(1,I) - length + rad)**2
            endif
            if (r2 > rad*rad) then
                in_confinement = .False.
                return
            endif
        enddo
    ! always inside confinement otherwise
    ! elseif(WLC_P__CONFINETYPE == 'periodicUnequal') then
    ! else
    ! print*, "Undefined comfone Type"
    ! stop 1
    endif
    return

end function


