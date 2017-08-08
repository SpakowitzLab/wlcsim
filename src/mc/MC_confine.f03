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
use inputparams, only: MAXPARAMLEN
implicit none

character(MAXPARAMLEN), intent(in) :: confineType  ! Specifier for type of confinement
real(dp) LBox(3) ! Side length of box
integer NT     ! Total number of beads in simulation
real(dp) RP(3,NT)  ! Bead positions
integer IT1    ! Start test bead
integer IT2    ! Final test bead
integer I      ! Index of bead being compared
real(dp) ECon

ECon = 0.0_dp


if (confineType == 'none') then
    return
elseif(confineType == 'platesInZ') then
    ! Confinement only in the z-direction
    ! limits: 0 and LBox
    do I = IT1,IT2
        if(RP(3,I)<0.0_dp) then
            ECon = 9990000.0_dp
        elseif (RP(3,I)>LBox(3)) then
            ECon = 9990000.0_dp
        endif
    ENDdo
elseif(confineType == 'cube') then
    do I = IT1,IT2
        if(RP(1,I)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(1,I)>LBox(1)) then
            ECon = 9990000.0_dp
        elseif(RP(2,I)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(2,I)>LBox(2)) then
            ECon = 9990000.0
        elseif(RP(3,I)<0.0) then
            ECon = 9990000.0_dp
        elseif(RP(3,I)>LBox(3)) then
            ECon = 9990000.0_dp
        endif
    ENDdo
elseif(confineType == 'sphere') then
    do I = IT1,IT2
        if(((RP(1,I)-LBox(1)/2)**2 + (RP(2,I)-LBox(1)/2_dp)**2 + &
           (RP(3,I)-LBox(1)/2)**2) > dble(LBox(1)*LBox(1)*0.25_dp)) then
            ECon = 9990000.0
        endif
    Enddo
elseif(confineType == 'periodicUnequal') then
    return
else
   print*, "Undefined comfone Type"
   stop 1
endif




END
