!
!       Code for MC move addaptations when you need to
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine MC_adapt(wlc_p,wlc_d,MCTYPE)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use params
    implicit none
    TYPE(wlcsim_params), intent(in) :: wlc_p
    TYPE(wlcsim_data), intent(inout) :: wlc_d
    integer, intent(in) :: MCTYPE   ! Type of move

!   Change the position if appropriate
    wlc_d%PHIT(MCTYPE) = real(wlc_d%SUCCESS(MCTYPE))/real(wlc_d%ATTEMPTS(MCTYPE))

    wlc_d%SUCCESS(MCTYPE) = 0
    wlc_d%ATTEMPTS(MCTYPE) = 0

    if ((MCTYPE.eq.8).or.(MCTYPE.eq.9)) then
        return
    endif

    ! If move type has no amplitude then only ajust window
    if (MCTYPE.eq.7) then
        if (wlc_d%PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
           wlc_d%WindoW(MCTYPE) = wlc_d%WindoW(MCTYPE)*1.05_dp
        else
           wlc_d%WindoW(MCTYPE) = wlc_d%WindoW(MCTYPE)*0.95_dp
        endif
        !window limits
        if (wlc_d%WindoW(MCTYPE) < wlc_p%MinWindoW(MCTYPE)) then
           wlc_d%WindoW(MCTYPE) = wlc_p%MinWindoW(MCTYPE)
        elseif (wlc_d%WindoW(MCTYPE) > wlc_p%MAXWindoW(MCTYPE)) then
           wlc_d%WindoW(MCTYPE) = wlc_p%MAXWindoW(MCTYPE)
        endif
        return
    endif

    ! If move has no window it doesn't need to be ajusted
    if ((MCTYPE.eq.4) .or. &
        (MCTYPE.eq.5) .or. &
        (MCTYPE.eq.6)) then
        ! Adjust Amplidtude
        if (wlc_d%PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
           wlc_d%MCAMP(MCTYPE) = wlc_d%MCAMP(MCTYPE)*1.05_dp
        else
           wlc_d%MCAMP(MCTYPE) = wlc_d%MCAMP(MCTYPE)*0.95_dp
        endif
        ! amplitude limits
        if (wlc_d%MCAMP(MCTYPE) > wlc_p%MAXAMP(MCTYPE)) then
           wlc_d%MCAMP(MCTYPE) = wlc_p%MAXAMP(MCTYPE)
        elseif (wlc_d%MCAMP(MCTYPE) < wlc_p%MinAMP(MCTYPE)) then
           wlc_d%MCAMP(MCTYPE) = wlc_p%MinAMP(MCTYPE)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (wlc_d%PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
       wlc_d%MCAMP(MCTYPE) = wlc_d%MCAMP(MCTYPE)*1.05_dp
       wlc_d%WindoW(MCTYPE) = wlc_d%WindoW(MCTYPE)*1.05_dp
    else
       wlc_d%MCAMP(MCTYPE) = wlc_d%MCAMP(MCTYPE)*0.95_dp
       wlc_d%WindoW(MCTYPE) = wlc_d%WindoW(MCTYPE)*0.95_dp
    endif

    ! Drift to target window
    wlc_d%WindoW(MCTYPE) = wlc_d%WindoW(MCTYPE)*0.98_dp + &
                      0.02_dp*wlc_p%winTarget(MCTYPE)


    ! amplitude limits
    if (wlc_d%MCAMP(MCTYPE) > wlc_p%MAXAMP(MCTYPE)) then
       wlc_d%MCAMP(MCTYPE) = wlc_p%MAXAMP(MCTYPE)
    elseif (wlc_d%MCAMP(MCTYPE) < wlc_p%MinAMP(MCTYPE)) then
       wlc_d%MCAMP(MCTYPE) = wlc_p%MinAMP(MCTYPE)
    endif

    !window limits
    if (wlc_d%WindoW(MCTYPE) < wlc_p%MinWindoW(MCTYPE)) then
       wlc_d%WindoW(MCTYPE) = wlc_p%MinWindoW(MCTYPE)
    elseif (wlc_d%WindoW(MCTYPE) > wlc_p%MAXWindoW(MCTYPE)) then
       wlc_d%WindoW(MCTYPE) = wlc_p%MAXWindoW(MCTYPE)
    endif

end subroutine
