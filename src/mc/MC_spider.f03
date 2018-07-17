#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_spider(wlc_d,MCAMP,rand_stat,success,spider_id)

use mersenne_twister
use params, only: wlcsim_params, wlcsim_data, pi
use precision, only: dp, eps
use vector_utils, only: randomUnitVec, cross, distance, angle_of_triangle, &
                        round_into_pm1, rotateR, rotateU, axisAngle, rotateAintoB, random_perp
use windowTools, only: drawWindow

implicit none
type(wlcsim_data), intent(inout) :: wlc_d
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
logical, intent(out) :: success

real(dp), intent(in) :: MCAMP ! Amplitude of random change

integer leg_n, I, section_n
integer irnd(1) ! random intiger
real urnd(1) ! single random number
integer spider_id ! which spider
real(dp) dr(3) ! ranslation
integer hip, knee, toe
real(dp) thigh, shin, dold, dnew
real(dp) theta ! extra abount to rotate each leg
real(dp) dalpha ! amount to rotate shin
real(dp) dbeta ! amount to rotate thigh
real(dp), dimension(3) :: hipR,kneeR,toeR,direction, temp,rold,rnew,rinter
real(dp) ROT(3,4) ! Rotation matrix
logical swing_past
integer dib, IT1,IT2,IB1,IB2,IP,otherEnd

! choose spider
call random_number(urnd,rand_stat)
if (WLC_P__PROBABILITY_PRECALC_SPIDER < urnd(1)) then
    call random_index(wlc_d%numberOfSpiders,irnd,rand_stat)
    spider_id=irnd(1)
else
    call drawWindow(wlc_d,WLC_P__SPIDER_WINDOW,WLC_P__MAXWINDOW_CRANK_SHAFT,.true.,rand_stat,&
                    IT1,IT2,IB1,IB2,IP,DIB,success)
    if (success .eqv. .false.) return


    spider_id=wlc_d%numberOfSpiders+1
    wlc_d%spiders(spider_id)%nSections=1
    wlc_d%spiders(spider_id)%nLegs=2
    wlc_d%spiders(spider_id)%sections(:,1) = [IT1,IT2]

    if (WLC_P__EXPLICIT_BINDING) then
        !  Right Leg
        knee = IT2
        do I = 1,WLC_P__SPIDER_LEG_LENGTH
            knee = knee + 1
            otherEnd = wlc_d%explicitbindingpair(knee)
            if (otherEnd < 0) continue

            if (otherEnd < knee) then
                success = .False.
                return
            endif
            knee = otherEnd
        enddo
        toe=knee
        do I = 1,WLC_P__SPIDER_LEG_LENGTH
            toe = toe + 1
            otherEnd = wlc_d%explicitbindingpair(toe)
            if (otherEnd < 0) continue

            if (otherEnd < toe) then
                success = .False.
                return
            endif
            toe = otherEnd
        enddo
        wlc_d%spiders(spider_id)%legs(:,1) = [IT2,knee,toe]

        !  Left Leg
        knee = IT1
        do I = 1,WLC_P__SPIDER_LEG_LENGTH
            knee = knee - 1
            otherEnd = wlc_d%explicitbindingpair(knee)
            if (otherEnd < 0) continue

            if (otherEnd > knee) then
                success = .False.
                return
            endif
            knee = otherEnd
        enddo
        toe=knee
        do I = 1,WLC_P__SPIDER_LEG_LENGTH
            toe = toe - 1
            otherEnd = wlc_d%explicitbindingpair(knee)
            if (otherEnd < 0) continue

            if (otherEnd > toe) then
                success = .False.
                return
            endif
            toe = otherEnd
        enddo
        wlc_d%spiders(spider_id)%legs(:,2) = [IT1,knee,toe]
    else
        knee = IT2+WLC_P__SPIDER_LEG_LENGTH
        toe = knee + WLC_P__SPIDER_LEG_LENGTH
        wlc_d%spiders(spider_id)%legs(:,1) = [IT2,knee,toe]

        knee = IT1-WLC_P__SPIDER_LEG_LENGTH
        toe = knee - WLC_P__SPIDER_LEG_LENGTH
        wlc_d%spiders(spider_id)%legs(:,2) = [IT1,knee,toe]
    endif

    toe = wlc_d%spiders(spider_id)%legs(3,1)
    if ((toe-1)/WLC_P__NB .ne. (IT2-1)/WLC_P__NB) then ! is on a different polyme
        success = .False.
        return
    endif
    if (toe-IT2 > WLC_P__MAX_SPIDER_LEG_LENGTH) then
        success = .False.
        return
    endif
    wlc_d%spiders(spider_id)%moved_sections(2,1) = toe

    toe = wlc_d%spiders(spider_id)%legs(3,2)
    if ((toe-1)/WLC_P__NB .ne. (IT1-1)/WLC_P__NB) then
        success = .False.
        return
    endif
    if (IT1-toe > WLC_P__MAX_SPIDER_LEG_LENGTH) then
        success = .False.
        return
    endif
    wlc_d%spiders(spider_id)%moved_sections(1,1) = toe


endif


! choose random offset.  Here we use p~r^(-2) between 0 and MCAMP
call randomUnitVec(dr,rand_stat)
call random_number(urnd,rand_stat)
dr=dr*urnd(1)*MCAMP

! check to see if all legs will reach
do leg_n = 1,wlc_d%spiders(spider_id)%nLegs
    hip = wlc_d%spiders(spider_id)%legs(1,leg_n)
    knee= wlc_d%spiders(spider_id)%legs(2,leg_n)
    toe = wlc_d%spiders(spider_id)%legs(3,leg_n)
    hipR=wlc_d%R(:,hip)
    kneeR=wlc_d%R(:,knee)
    toeR=wlc_d%R(:,toe)
    ! check triangle inequlity
    if (distance(hipR,kneeR)+distance(kneeR,toeR) < distance(hipR+dr,toeR)) then
        success = .FALSE.
        return
    endif
    if (abs(distance(hipR,kneeR)-distance(kneeR,toeR)) > distance(hipR+dr,toeR)) then
        success = .FALSE.
        return
    endif
enddo
success= .TRUE.

! move legs
do leg_n = 1,wlc_d%spiders(spider_id)%nLegs
    hip = wlc_d%spiders(spider_id)%legs(1,leg_n)
    knee= wlc_d%spiders(spider_id)%legs(2,leg_n)
    toe = wlc_d%spiders(spider_id)%legs(3,leg_n)
    hipR=wlc_d%R(:,hip)
    kneeR=wlc_d%R(:,knee)
    toeR=wlc_d%R(:,toe)
    thigh = distance(kneeR,hipR)
    shin = distance(kneeR,toeR)

    rold = hipR-toeR
    dold = norm2(rold)
    rnew = hipR+dr-toeR
    dnew = norm2(rnew)

    if (dold<=0.001_dp*eps) then
        ! if hip and toe are in the same place, extend in random direction
        call random_perp((kneeR-toeR)/shin,direction,rinter,rand_stat)
        if (dot_product(rinter,rnew)<0.0_dp) then
            rinter=-1.0_dp*rinter
        endif
        direction = cross((kneeR-toeR),rinter)
        direction = direction/norm2(direction)
        rinter = rinter*dnew
        dalpha = asin(dnew/(2.0_dp*shin))
        dbeta = -dalpha
        swing_past = .False.
    else
        ! Calculate how much to rotate each leg to get new |distance|
        swing_past = dot_product(rold,rnew) < 0.0_dp  ! toe passes hip
        if (swing_past) then
            dalpha = angle_of_triangle(thigh,shin,dold) - (PI - angle_of_triangle(thigh,shin,dnew))
            dbeta = (PI - angle_of_triangle(shin,thigh,dnew)) - angle_of_triangle(shin,thigh,dold)
        else
            dalpha = angle_of_triangle(thigh,shin,dold) - angle_of_triangle(thigh,shin,dnew)
            dbeta = angle_of_triangle(shin,thigh,dnew) - angle_of_triangle(shin,thigh,dold)
        endif
        ! direction needs to be consistant with other angles
        direction = cross(kneeR-toeR,hipR-kneeR)
        if (norm2(direction)<0.01_dp*eps) then
            ! if fully exteded/contracted then contract/extend in random direction
            call random_perp(rnew/norm2(rnew),direction,rinter,rand_stat)
        else
            direction = direction/norm2(direction)
        endif
        if (swing_past) then
            rinter = -dnew*rold/dold
        else
            rinter = dnew*rold/dold
        endif
    endif

    ! angle to rotate shin
    call axisAngle(ROT,dalpha,direction,toeR)

    ! rotate shin to intermediate position
    do I = min(knee,toe),max(knee,toe)
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%R(:,I))
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%U(:,I))
    enddo
    temp = wlc_d%RP(:,toe)

    !Check toe stayed in the sampe place
    if ( distance(wlc_d%RP(:,toe),wlc_d%R(:,toe)) > eps ) then
        print*, "Broken toe in spider move"
        stop 1
    endif

    ! angle to rotate thigh
    call axisAngle(ROT,dbeta,direction,hipR)

    ! Rotation matrix for hip to intermediate position
    I = knee
    temp = rotateR(ROT,wlc_d%R(:,I)) + rinter - rold

    ! Check knee
    if ( distance(wlc_d%RP(:,knee),temp) > eps ) then
        print*, "Broken knee"
        stop 1
    endif

    ! rotate thigh
    do I = min(knee,hip),max(knee,hip)
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%R(:,I)) + rinter - rold
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%U(:,I))
    enddo

    ! rotate from intermediate position to final position
    call rotateAIntoB(rinter,rnew,toeR,ROT)
    do I = min(toe,hip),max(toe,hip)
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%RP(:,I))
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%UP(:,I))
    enddo

    ! don't rotatie hip and toe
    wlc_d%UP(:,hip)=wlc_d%U(:,hip)
    wlc_d%UP(:,toe)=wlc_d%U(:,toe)

    !Check toe stayed in the sampe place
    if ( distance(wlc_d%RP(:,toe),wlc_d%R(:,toe)) > eps ) then
        print*, "Broken toe in second half of spider move"
        stop 1
    endif
    ! chack to make sure hip moved the correct amount
    if ( distance(wlc_d%RP(:,hip),hipR+dr) > eps ) then
        print*, "Broken hip"
    endif

enddo

! translate sections
do section_n = 1,wlc_d%spiders(spider_id)%nLegs
    do I = wlc_d%spiders(spider_id)%sections(1,section_n), wlc_d%spiders(spider_id)%sections(2,section_n)
        wlc_d%RP(:,I) = wlc_d%R(:,I) + dr
        wlc_d%UP(:,I) = wlc_d%U(:,I)
    enddo
enddo
end subroutine
