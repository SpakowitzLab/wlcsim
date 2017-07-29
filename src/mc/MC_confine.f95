!-----------------------------------------------------------!
!
!     Make energy Huge when polymer goes outside boundary
!
!            Started by Quinn 2/17/15
!
!
!
! confineType  |  Discription
! _____________|_________________________________
!    0         |  No confinement
!    1         |  Betwene two plates in Z direction at 0 and LBox
!    2         |  Cube of size LBox**3,  range: 0-LBox
!    3         |  Circle of radius LBox, centered at LBox/2
!    4         |  Periodic, non-equal lengths

subroutine MC_confine(confineType, LBox, RP, NT, IT1, IT2, ECon)
use params, only: dp


implicit none

integer confineType  ! Specifier for type of confinement
real(dp) LBox(3) ! Side length of box
integer NT     ! Total number of beads in simulation
real(dp) RP(NT,3)  ! Bead positions
integer IT1    ! Start test bead
integer IT2    ! Final test bead
integer I      ! Index of bead being compared
real(dp) ECon

ECon = 0.0_dp


if (confineType == 0) then
    return
elseif(confineType == 1) then
    ! Confinement only in the z-direction
    ! limits: 0 and LBox
    do I = IT1,IT2
        if(RP(I,3)<0.0_dp) then
            ECon = 9990000.0_dp
        elseif (RP(I,3)>LBox(3)) then
            ECon = 9990000.0_dp
        endif
    ENDdo
elseif(confineType == 2) then
    do I = IT1,IT2
        if(RP(I,1)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(I,1)>LBox(1)) then
            ECon = 9990000.0_dp
        elseif(RP(I,2)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(I,2)>LBox(2)) then
            ECon = 9990000.0
        elseif(RP(I,3)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(I,3)>LBox(3)) then
            ECon = 9990000.0_dp
        endif
    ENDdo
elseif(confineType == 3) then
    do I = IT1,IT2
        if(((RP(I,1)-LBox(1)/2)**2 + (RP(I,2)-LBox(1)/2_dp)**2 + &
           (RP(I,3)-LBox(1)/2)**2) > dble(LBox(1)*LBox(1)*0.25_dp)) then
            ECon = 9990000.0
        endif
    Enddo
elseif(confineType == 4) then
    return
else
   print*, "Undefined comfone Type"
   stop 1
endif




END
