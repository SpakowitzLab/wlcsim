#include "../defines.inc"
pure function in_confinement(RP, NT, IT1, IT2)
    use params, only : dp, wlcsim_params
    use inputparams, only : MAXPARAMLEN

    implicit none

    integer, intent(in) :: IT1, IT2, NT
    real(dp), intent(in) :: RP(3, NT)
    logical in_confinement
    integer i
    real(dp) rad, length, r2
    real(dp), parameter :: center(3) = [WLC_P__LBOX_X/2.0_dp,&
                                        WLC_P__LBOX_Y/2.0_dp,&
                                        WLC_P__LBOX_Z/2.0_dp]

    in_confinement = .True.
    ! if (WLC_P__CONFINETYPE == 'none') then
    !     in_confinement = .False.
    if (WLC_P__CONFINETYPE == 'platesInZperiodicXY') then
        ! Confinement only in the z-direction
        ! limits: 0 and LBox(3)
        do I = IT1,IT2
            if ((RP(3,I) < 0.0_dp+WLC_P__DBIN) .or.&
                (RP(3,I) > WLC_P__CONFINEMENT_SLIT_WIDTH-WLC_P__DBIN)) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'cube') then
        do I = IT1,IT2
            if ((RP(1,I) < 0.0_dp+WLC_P__DBIN) &
                .or. (RP(1,I) > WLC_P__CONFINEMENT_CUBE_LENGTH-WLC_P__DBIN) &
                .or. (RP(2,I) < 0.0_dp+WLC_P__DBIN) &
                .or. (RP(2,I) > WLC_P__CONFINEMENT_CUBE_LENGTH-WLC_P__DBIN) &
                .or. (RP(3,I) < 0.0_dp+WLC_P__DBIN) &
                .or. (RP(3,I) > WLC_P__CONFINEMENT_CUBE_LENGTH-WLC_P__DBIN)) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'sphere') then
        ! sphere with given diameter
        do I = IT1,IT2
            rad = WLC_P__CONFINEMENT_SPHERE_DIAMETER/2
            if ((RP(1,I) - center(1))**2 + (RP(2,I) - center(2))**2 + &
                (RP(3,I) - center(3))**2 > rad) then
                in_confinement = .False.
                return
            endif
        enddo
    elseif (WLC_P__CONFINETYPE == 'ecoli') then
        ! cylinder with hemispherical caps, one tip at origin
        ! full length - lbox(1)/CONFINEMENT_ECOLI_LENGTH
        ! diameter - lbox(2:3)/CONFINEMENT_ECOLI_DIAMETER
        do I = IT1,IT2
            length = WLC_P__CONFINEMENT_ECOLI_LENGTH
            rad = WLC_P__CONFINEMENT_ECOLI_DIAMETER/2
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


