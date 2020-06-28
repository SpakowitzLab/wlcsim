!
!       Code for MC move addaptations when you need to
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine mc_adapt(wlc_p,MCTYPE)
! values from wlcsim_data
use params, only: wlc_SUCCESS, wlc_ATTEMPTS, wlc_WindoW, wlc_PHIT, wlc_MCAMP
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use params
    implicit none
    TYPE(wlcsim_params), intent(in) :: wlc_p
    integer, intent(in) :: MCTYPE   ! Type of move

    if(wlc_ATTEMPTS(MCTYPE) < 10) then
        return
    endif

!   Change the position if appropriate
    wlc_PHIT(MCTYPE) = real(wlc_SUCCESS(MCTYPE),dp) &
                        /real(wlc_ATTEMPTS(MCTYPE),dp)

    wlc_SUCCESS(MCTYPE) = 0
    wlc_ATTEMPTS(MCTYPE) = 0

    if ((MCTYPE.eq.8).or.(MCTYPE.eq.9).or.(MCTYPE.eq.10).or.(MCTYPE.eq.11)) then
        return
    endif

    ! If move type has no amplitude then only ajust window
    if (MCTYPE.eq.7) then
        if (wlc_PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
           wlc_WindoW(MCTYPE) = wlc_WindoW(MCTYPE)*1.05_dp
        else
           wlc_WindoW(MCTYPE) = wlc_WindoW(MCTYPE)*0.95_dp
        endif
        !window limits
        if (wlc_WindoW(MCTYPE) < wlc_p%MINWINDOW(MCTYPE)) then
           wlc_WindoW(MCTYPE) = wlc_p%MINWINDOW(MCTYPE)
        elseif (wlc_WindoW(MCTYPE) > wlc_p%MAXWINDOW(MCTYPE)) then
           wlc_WindoW(MCTYPE) = wlc_p%MAXWINDOW(MCTYPE)
        endif
        return
    endif

    ! If move has no window it doesn't need to be ajusted
    if ((MCTYPE.eq.4) .or. &
        (MCTYPE.eq.5) .or. &
        (MCTYPE.eq.6) .or. &
        (MCTYPE.eq.12)) then
        ! Adjust Amplidtude
        if (wlc_PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
           wlc_MCAMP(MCTYPE) = wlc_MCAMP(MCTYPE)*1.05_dp
        else
           wlc_MCAMP(MCTYPE) = wlc_MCAMP(MCTYPE)*0.95_dp
        endif
        ! amplitude limits
        if (wlc_MCAMP(MCTYPE) > wlc_p%MAXAMP(MCTYPE)) then
           wlc_MCAMP(MCTYPE) = wlc_p%MAXAMP(MCTYPE)
        elseif (wlc_MCAMP(MCTYPE) < wlc_p%MINAMP(MCTYPE)) then
           wlc_MCAMP(MCTYPE) = wlc_p%MINAMP(MCTYPE)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (wlc_PHIT(MCTYPE) > wlc_p%PDESIRE(MCTYPE)) then
       wlc_MCAMP(MCTYPE) = wlc_MCAMP(MCTYPE)*1.05_dp
       wlc_WindoW(MCTYPE) = wlc_WindoW(MCTYPE)*1.05_dp
    else
       wlc_MCAMP(MCTYPE) = wlc_MCAMP(MCTYPE)*0.95_dp
       wlc_WindoW(MCTYPE) = wlc_WindoW(MCTYPE)*0.95_dp
    endif

    ! Drift to target window
    wlc_WindoW(MCTYPE) = wlc_WindoW(MCTYPE)*0.98_dp + &
                      0.02_dp*wlc_p%WINTARGET(MCTYPE)


    ! amplitude limits
    if (wlc_MCAMP(MCTYPE) > wlc_p%MAXAMP(MCTYPE)) then
       wlc_MCAMP(MCTYPE) = wlc_p%MAXAMP(MCTYPE)
    elseif (wlc_MCAMP(MCTYPE) < wlc_p%MINAMP(MCTYPE)) then
       wlc_MCAMP(MCTYPE) = wlc_p%MINAMP(MCTYPE)
    endif

    !window limits
    if (wlc_WindoW(MCTYPE) < wlc_p%MINWINDOW(MCTYPE)) then
       wlc_WindoW(MCTYPE) = wlc_p%MINWINDOW(MCTYPE)
    elseif (wlc_WindoW(MCTYPE) > wlc_p%MAXWINDOW(MCTYPE)) then
       wlc_WindoW(MCTYPE) = wlc_p%MAXWINDOW(MCTYPE)
    endif

end subroutine
