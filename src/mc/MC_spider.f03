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
use vector_utils, only: randomUnitVec, cross, distance, angle_of_triangle, round_into_pm1, rotateR, rotateU

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
real(dp), dimension(3) :: hipR,kneeR,toeR,direction, temp,rold,rnew
real(dp) ROT(3,4) ! Rotation matrix
logical swing_past
! choose spider
call random_index(wlc_d%numberOfSpiders,irnd,rand_stat)
spider_id=irnd(1)

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
    if (distance(hipR,kneeR)+distance(kneeR,toeR) < distance(hipR+dr,toeR)) then
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

    ! Calculate how much to rotate each leg to get new |distance|
    swing_past = dot_product(rold,rnew) < 0.0_dp  ! toe passes hip
    if (swing_past) then
        dalpha = angle_of_triangle(thigh,shin,dold) - (PI - angle_of_triangle(thigh,shin,dnew))
        dbeta = (PI - angle_of_triangle(thigh,shin,dnew)) - angle_of_triangle(thigh,shin,dold)
    else
        dalpha = angle_of_triangle(thigh,shin,dold) - angle_of_triangle(thigh,shin,dnew)
        dbeta = angle_of_triangle(thigh,shin,dnew) - angle_of_triangle(thigh,shin,dold)
    endif

    ! direction needs to be consistant with other angles
    direction = cross(kneeR-toeR,hipR-kneeR)
    direction = direction/norm2(direction)

    ! claculate the angle the vector between toe and hip must turn
    if (swing_past) then
        theta = asin(round_into_pm1( dot_product(cross(-rold,rnew),direction)/(dold*dnew) ))
    else
        theta = asin(round_into_pm1( dot_product(cross( rold,rnew),direction)/(dold*dnew) ))
    endif

    ! angle to rotate shin
    call axisAngle(ROT,dalpha+theta,direction,toeR)

    ! rotate shin
    do I = min(knee,toe),max(knee,toe)
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%R(:,I))
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%U(:,I))
    enddo

    !Check toe stayed in the sampe place
    if ( distance(wlc_d%RP(:,toe),wlc_d%R(:,toe)) > eps ) then
        print*, "Broken toe in spider move"
        stop 1
    endif

    ! angle to rotate thigh
    call axisAngle(ROT,dbeta+theta,direction,hipR)

    ! Check knee
    I = knee
    temp = rotateR(ROT,wlc_d%R(:,I)) + dr
    if ( distance(wlc_d%RP(:,knee),temp) > eps ) then
        print*, "Broken knee"
        stop 1
    endif

    ! rotate thigh
    do I = min(knee,hip),max(knee,hip)
        wlc_d%RP(:,I) = rotateR(ROT,wlc_d%R(:,I)) + dr
        wlc_d%UP(:,I) = rotateU(ROT,wlc_d%U(:,I))
    enddo

    ! don't rotatie hip and toe
    wlc_d%UP(:,hip)=wlc_d%U(:,hip)
    wlc_d%UP(:,toe)=wlc_d%U(:,toe)

    ! chack to make sure hip moved the correct amount
    if ( distance(wlc_d%RP(:,hip),hipR+dr) > eps ) then
        print*, "Broken hip"
    endif

enddo

! translate sections
do section_n = 1,wlc_d%spiders(spider_id)%nLegs
    do I = wlc_d%spiders(spider_id)%sections(1,section_n), wlc_d%spiders(spider_id)%sections(2,section_n)
        wlc_d%RP(:,I) = wlc_d%R(:,I) + dr
    enddo
enddo
end subroutine
