!
!       Code for MC move addaptations when you need to
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine MC_adapt(mc,md,MCTYPE)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use params
    IMPLICIT NONE
    TYPE(wlcsim_params), intent(in) :: mc
    TYPE(wlcsim_data), intent(inout) :: md
    INTEGER, intent(in) :: MCTYPE   ! Type of move

    ! Correct for turned down poor moves
    if ((md%PHit(MCTYPE).lt.mc%MIN_ACCEPT).and. &
        ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
        md%SUCCESS(MCTYPE)=md%SUCCESS(MCTYPE)*mc%reduce_move
    endif

!   Change the position if appropriate
    md%PHIT(MCTYPE)=real(md%SUCCESS(MCTYPE))/real(mc%NADAPT(MCTYPE))

    md%SUCCESS(MCTYPE)=0

    if ((MCTYPE.eq.8).or.(MCTYPE.eq.9)) then
        return
    endif

    ! If move type has no amplitude then only ajust window
    if (MCTYPE.eq.7) then
        if (md%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
           md%WINDOW(MCTYPE)=md%WINDOW(MCTYPE)*1.05_dp
        else
           md%WINDOW(MCTYPE)=md%WINDOW(MCTYPE)*0.95_dp
        endif
        !window limits
        if (md%WINDOW(MCTYPE).LT.mc%MINWINDOW(MCTYPE)) then
           md%WINDOW(MCTYPE)=mc%MINWINDOW(MCTYPE)
        elseif (md%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
           md%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
        endif
        return
    endif

    ! If move has no window it doesn't need to be ajusted
    if ((MCTYPE.eq.4) .or. &
        (MCTYPE.eq.5) .or. &
        (MCTYPE.eq.6)) then
        ! Adjust Amplidtude
        if (md%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
           md%MCAMP(MCTYPE)=md%MCAMP(MCTYPE)*1.05_dp
        else
           md%MCAMP(MCTYPE)=md%MCAMP(MCTYPE)*0.95_dp
        endif
        ! amplitude limits
        if (md%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
           md%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
        elseif (md%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
           md%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (md%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
       md%MCAMP(MCTYPE)=md%MCAMP(MCTYPE)*1.05_dp
       md%WINDOW(MCTYPE)=md%WINDOW(MCTYPE)*1.05_dp
    else
       md%MCAMP(MCTYPE)=md%MCAMP(MCTYPE)*0.95_dp
       md%WINDOW(MCTYPE)=md%WINDOW(MCTYPE)*0.95_dp
    endif

    ! Drift to target window
    md%WINDOW(MCTYPE)=md%WINDOW(MCTYPE)*0.98_dp+&
                      0.02_dp*mc%winTarget(MCTYPE)


    ! amplitude limits
    if (md%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
       md%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
    elseif (md%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
       md%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
    endif

    !window limits
    if (md%WINDOW(MCTYPE).LT.mc%MINWINDOW(MCTYPE)) then
       md%WINDOW(MCTYPE)=mc%MINWINDOW(MCTYPE)
    elseif (md%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
       md%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
    endif

end subroutine
