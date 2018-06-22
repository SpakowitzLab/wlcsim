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
subroutine MC_spider(wlc_d,MCAMP,rand_stat)

use mersenne_twister
use params, only: wlcsim_params, wlcsim_data
use precision, only: dp, eps
use vector_utils, only: randomUnitVec, cross, distance, angle_of_triangle

implicit none
type(wlcsim_data), intent(inout) :: wlc_d
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator

real(dp), intent(in) :: MCAMP ! Amplitude of random change

integer leg_n, I, section_n
integer irnd(1) ! random intiger
real urnd(1) ! single random number
integer id ! which spider
real(dp) dr(3) ! ranslation
integer hip, knee, toe
real(dp) thigh, shin, dold, dnew
real(dp) dalpha ! amount to rotate shin
real(dp) dbeta ! amount to rotate thigh
real(dp), dimension(3) :: hipR,kneeR,toeR,direction, temp
logical success
real(dp) ROT(3,4) ! Rotation matrix

! choose spider
call random_index(wlc_d%numberOfSpiders,irnd,rand_stat)
id=irnd(1)

! choose random offset.  Here we use p~r^(-2) between 0 and MCAMP
call randomUnitVec(dr,rand_stat)
call random_number(urnd,rand_stat)
dr=dr*urnd(1)*MCAMP

! check to see if all legs will reach
do leg_n = 1,wlc_d%spiders(id)%nLegs
    hip = wlc_d%spiders(id)%legs(1,leg_n)
    knee= wlc_d%spiders(id)%legs(2,leg_n)
    toe = wlc_d%spiders(id)%legs(3,leg_n)
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
do leg_n = 1,wlc_d%spiders(id)%nLegs
    hip = wlc_d%spiders(id)%legs(1,leg_n)
    knee= wlc_d%spiders(id)%legs(2,leg_n)
    toe = wlc_d%spiders(id)%legs(3,leg_n)
    hipR=wlc_d%R(:,hip)
    kneeR=wlc_d%R(:,knee)
    toeR=wlc_d%R(:,toe)
    thigh = distance(kneeR,hipR)
    shin = distance(kneeR,toeR)
    
    dold = distance(hipR,toeR)
    dnew = distance(hipR+dr,toeR)

    dalpha = angle_of_triangle(thigh,shin,dold) - angle_of_triangle(thigh,shin,dnew)
    dbeta = angle_of_triangle(thigh,shin,dnew) - angle_of_triangle(thigh,shin,dold)
    
    direction = cross(kneeR-toeR,hipR-kneeR)

    ! angle to rotate shin
    call axisAngle(ROT,dalpha,direction,toeR)

    ! rotate shin
    do I = min(knee,toe),max(knee,toe)
        wlc_d%RP(1,I) = ROT(1,4) + ROT(1,1)*wlc_d%R(1,I) + ROT(1,2)*wlc_d%R(2,I) + ROT(1,3)*wlc_d%R(3,I)
        wlc_d%RP(2,I) = ROT(2,4) + ROT(2,1)*wlc_d%R(1,I) + ROT(2,2)*wlc_d%R(2,I) + ROT(2,3)*wlc_d%R(3,I)
        wlc_d%RP(3,I) = ROT(3,4) + ROT(3,1)*wlc_d%R(1,I) + ROT(3,2)*wlc_d%R(2,I) + ROT(3,3)*wlc_d%R(3,I)
        wlc_d%UP(1,I) = ROT(1,1)*wlc_d%U(1,I) + ROT(1,2)*wlc_d%U(2,I) + ROT(1,3)*wlc_d%U(3,I)
        wlc_d%UP(2,I) = ROT(2,1)*wlc_d%U(1,I) + ROT(2,2)*wlc_d%U(2,I) + ROT(2,3)*wlc_d%U(3,I)
        wlc_d%UP(3,I) = ROT(3,1)*wlc_d%U(1,I) + ROT(3,2)*wlc_d%U(2,I) + ROT(3,3)*wlc_d%U(3,I)
    enddo

    !Check toe stayed in the sampe place
    if ( distance(wlc_d%RP(:,toe),wlc_d%R(:,toe)) > eps ) then
        print*, "Broken toe in spider move"
        stop 1
    endif

    ! angle to rotate thigh
    call axisAngle(ROT,dbeta,direction,hipR)

    ! Check knee
    I = knee 
    temp(1) = ROT(1,4) + ROT(1,1)*wlc_d%R(1,I) + ROT(1,2)*wlc_d%R(2,I) + ROT(1,3)*wlc_d%R(3,I) + dr(1) 
    temp(2) = ROT(2,4) + ROT(2,1)*wlc_d%R(1,I) + ROT(2,2)*wlc_d%R(2,I) + ROT(2,3)*wlc_d%R(3,I) + dr(2)
    temp(3) = ROT(3,4) + ROT(3,1)*wlc_d%R(1,I) + ROT(3,2)*wlc_d%R(2,I) + ROT(3,3)*wlc_d%R(3,I) + dr(3)
    if ( distance(wlc_d%RP(:,knee),temp) > eps ) then 
        print*, "Broken knee"
        stop 1
    endif

    ! rotate thigh
    do I = min(knee,hip),max(knee,hip)
        wlc_d%RP(1,I) = ROT(1,4) + ROT(1,1)*wlc_d%R(1,I) + ROT(1,2)*wlc_d%R(2,I) + ROT(1,3)*wlc_d%R(3,I) + dr(1)
        wlc_d%RP(2,I) = ROT(2,4) + ROT(2,1)*wlc_d%R(1,I) + ROT(2,2)*wlc_d%R(2,I) + ROT(2,3)*wlc_d%R(3,I) + dr(2)
        wlc_d%RP(3,I) = ROT(3,4) + ROT(3,1)*wlc_d%R(1,I) + ROT(3,2)*wlc_d%R(2,I) + ROT(3,3)*wlc_d%R(3,I) + dr(3)
        wlc_d%UP(1,I) = ROT(1,1)*wlc_d%U(1,I) + ROT(1,2)*wlc_d%U(2,I) + ROT(1,3)*wlc_d%U(3,I)
        wlc_d%UP(2,I) = ROT(2,1)*wlc_d%U(1,I) + ROT(2,2)*wlc_d%U(2,I) + ROT(2,3)*wlc_d%U(3,I)
        wlc_d%UP(3,I) = ROT(3,1)*wlc_d%U(1,I) + ROT(3,2)*wlc_d%U(2,I) + ROT(3,3)*wlc_d%U(3,I)
    enddo

    ! chack to make sure hip moved the correct amount
    if ( distance(wlc_d%RP(:,hip),hipR+dr) > eps ) then 
        print*, "Broken hip"
    endif

enddo

! translate sections
do section_n = 1,wlc_d%spiders(id)%nLegs
    do I = wlc_d%spiders(id)%sections(1,section_n), wlc_d%spiders(id)%sections(2,section_n)
        wlc_d%RP(:,I) = wlc_d%R(:,I) + dr
    enddo
enddo
end subroutine
