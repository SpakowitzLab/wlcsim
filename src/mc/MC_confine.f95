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

SUBROUTINE MC_confine(confineType, LBox, RP, NT, IT1, IT2, ECon)
use params, only: dp


IMPLICIT NONE

INTEGER confineType  ! Specifier for type of confinement
DOUBLE PRECISION LBox(3) ! Side length of box
INTEGER NT     ! Total number of beads in simulation
DOUBLE PRECISION RP(NT,3)  ! Bead positions
INTEGER IT1    ! Start test bead
INTEGER IT2    ! Final test bead
INTEGER I      ! Index of bead being compared
DOUBLE PRECISION ECon

ECon=0.0_dp


if (confineType.EQ.0) then
    return
elseif(confineType.EQ.1) then
    ! Confinement only in the z-direction
    ! limits: 0 and LBox
    DO I=IT1,IT2
        if(RP(I,3)<0.0_dp) then
            ECon=9990000.0_dp
        elseif (RP(I,3)>LBox(3)) then
            ECon=9990000.0_dp
        endif
    ENDDO
elseif(confineType.EQ.2) then
    DO I=IT1,IT2
        if(RP(I,1)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,1)>LBox(1)) then
            ECon=9990000.0_dp
        elseif(RP(I,2)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,2)>LBox(2)) then
            ECon=9990000.0
        elseif(RP(I,3)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,3)>LBox(3)) then
            ECon=9990000.0_dp
        endif
    ENDDO
elseif(confineType.EQ.3) then
    DO I=IT1,IT2
        if(((RP(I,1)-LBox(1)/2)**2 + (RP(I,2)-LBox(1)/2_dp)**2 + &
           (RP(I,3)-LBox(1)/2)**2).GT.dble(LBox(1)*LBox(1)*0.25_dp)) then
            ECon=9990000.0
        endif
    Enddo
elseif(confineType.EQ.4) then
    return
else
   print*, "Undefined comfone Type"
   stop 1
endif




END
