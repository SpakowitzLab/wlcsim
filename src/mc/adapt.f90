!
!       Code for MC move addaptations when you need to
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine mc_adapt(mc,md,MCTYPE)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use params
    IMPLICIT NONE
    TYPE(wlcsim_params), intent(inout) :: mc
    TYPE(wlcsim_data), intent(inout) :: md
    INTEGER, intent(in) :: mcTYPE   ! Type of move

    ! Correct for turned down poor moves
    if ((md%PHit(mctype).lt.md%MIN_ACCEPT).and. &
        ((mctype.eq.5).or.(mctype.eq.6))) then
        md%SUCCESS(mctype)=md%SUCCESS(mctype)*md%reduce_move
    endif

!   Change the position if appropriate
    md%PHIT(mctype)=real(md%SUCCESS(mctype))/real(md%NADAPT(mctype))
    md%SUCCESS(mctype)=0

    if ((mctype.eq.8).or.(mctype.eq.9)) then
        return
    endif

    ! If move type has no amplitude then only ajust window
    if (mctype.eq.7) then
        if (md%PHIT(mctype).GT.md%PDESIRE(mctype)) then
           md%WINDOW(mctype)=md%WINDOW(mctype)*1.05_dp
        else
           md%WINDOW(mctype)=md%WINDOW(mctype)*0.95_dp
        endif
        !window limits
        if (md%WINDOW(mctype).LT.md%MINWINDOW(mctype)) then
           md%WINDOW(mctype)=md%MINWINDOW(mctype)
        elseif (md%WINDOW(mctype).GT.md%MAXWINDOW(mctype)) then
           md%WINDOW(mctype)=md%MAXWINDOW(mctype)
        endif
        return
    endif

    ! If move has no window it doesn't need to be ajusted
    if ((mctype.eq.4) .or. &
        (mctype.eq.5) .or. &
        (mctype.eq.6)) then
        ! Adjust Amplidtude
        if (md%PHIT(mctype).GT.md%PDESIRE(mctype)) then
           md%mdAMP(mctype)=md%mdAMP(mctype)*1.05_dp
        else
           md%mdAMP(mctype)=md%mdAMP(mctype)*0.95_dp
        endif
        ! amplitude limits
        if (md%mdAMP(mctype).GT.md%MAXAMP(mctype)) then
           md%mdAMP(mctype)=md%MAXAMP(mctype)
        elseif (md%mdAMP(mctype).LT.md%MINAMP(mctype)) then
           md%mdAMP(mctype)=md%MINAMP(mctype)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (md%PHIT(mctype).GT.md%PDESIRE(mctype)) then
       md%mdAMP(mctype)=md%mdAMP(mctype)*1.05_dp
       md%WINDOW(mctype)=md%WINDOW(mctype)*1.05_dp
    else
       md%mdAMP(mctype)=md%mdAMP(mctype)*0.95_dp
       md%WINDOW(mctype)=md%WINDOW(mctype)*0.95_dp
    endif

    ! Drift to target window
    md%WINDOW(mctype)=md%WINDOW(mctype)*0.98_dp+&
                      0.02_dp*md%winTarget(mctype)


    ! amplitude limits
    if (md%mdAMP(mctype).GT.md%MAXAMP(mctype)) then
       md%mdAMP(mctype)=md%MAXAMP(mctype)
    elseif (md%mdAMP(mctype).LT.md%MINAMP(mctype)) then
       md%mdAMP(mctype)=md%MINAMP(mctype)
    endif

    !window limits
    if (md%WINDOW(mctype).LT.md%MINWINDOW(mctype)) then
       md%WINDOW(mctype)=md%MINWINDOW(mctype)
    elseif (md%WINDOW(mctype).GT.md%MAXWINDOW(mctype)) then
       md%WINDOW(mctype)=md%MAXWINDOW(mctype)
    endif

end subroutine
