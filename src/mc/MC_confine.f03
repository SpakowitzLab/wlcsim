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

subroutine MC_confine(RP, NT, IT1, IT2, ECon)
use params, only: dp, HUGE_ENERGY, wlcsim_params
implicit none

integer NT     ! Total number of beads in simulation
real(dp) RP(3,NT)  ! Bead positions
integer IT1    ! Start test bead
integer IT2    ! Final test bead
real(dp) ECon
logical in_confinement

ECon = 0.0_dp
if (.not. in_confinement(RP, NT, IT1, IT2)) then
    ECon = HUGE_ENERGY
endif

end
