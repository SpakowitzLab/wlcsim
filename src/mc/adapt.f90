!
!       Code for MC move addaptations when you need to 
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine MCvar_adapt(mc,MCTYPE)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use simMod
    IMPLICIT NONE
    TYPE(MCvar), intent(inout) :: mc
    INTEGER, intent(in) :: MCTYPE   ! Type of move

    ! Correct for turned down poor moves
    if ((mc%PHit(MCTYPE).lt.mc%MIN_ACCEPT).and. &
        ((MCTYPE.eq.5).or.(MCTYPE.eq.6))) then
        mc%SUCCESS(MCTYPE)=mc%SUCCESS(MCTYPE)*mc%reduce_move
    endif

!   Change the position if appropriate
    mc%PHIT(MCTYPE)=real(mc%SUCCESS(MCTYPE))/real(mc%NADAPT(MCTYPE))
    mc%SUCCESS(MCTYPE)=0

    if ((MCTYPE.eq.8).or.(MCTYPE.eq.9)) then
        return
    endif

    ! If move type has no amplitude then only ajust window
    if (MCTYPE.eq.7) then
        if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
           mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)*1.05_dp
        else
           mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)*0.95_dp
        endif
        !window limits
        if (mc%WINDOW(MCTYPE).LT.mc%MINWINDOW(MCTYPE)) then
           mc%WINDOW(MCTYPE)=mc%MINWINDOW(MCTYPE)
        elseif (mc%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
           mc%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
        endif
        return
    endif

    ! If move has no window it doesn't need to be ajusted
    if ((MCTYPE.eq.4) .or. & 
        (MCTYPE.eq.5) .or. &
        (MCTYPE.eq.6)) then
        ! Adjust Amplidtude
        if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
           mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*1.05_dp
        else
           mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*0.95_dp
        endif
        ! amplitude limits
        if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
           mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
        elseif (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
           mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*1.05_dp
       mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)*1.05_dp
    else
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*0.95_dp
       mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)*0.95_dp
    endif

    ! Drift to target window
    mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)*0.98_dp+&
                      0.02_dp*mc%winTarget(MCTYPE)
    

    ! amplitude limits
    if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
    elseif (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
    endif
    
    !window limits
    if (mc%WINDOW(MCTYPE).LT.mc%MINWINDOW(MCTYPE)) then
       mc%WINDOW(MCTYPE)=mc%MINWINDOW(MCTYPE)
    elseif (mc%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
       mc%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
    endif

end subroutine
