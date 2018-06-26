#include "../defines.inc"
!---------------------------------------------------------------*
!
!      Chain backbone potential between two beads
!
! E_wlc(RM1,R,RP1,EB)
! E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
! E_GAUSS(R,RM1,EPAR)
!---------------------------------------------------------------
module MC_wlc
implicit none
contains

function E_wlc(RM1,R,RP1,EB)
    use params, only: dp
    implicit none
    real(dp), intent(in), dimension(3) :: RM1 ! R of bead i-1
    real(dp), intent(in), dimension(3) :: R ! R of bead i
    real(dp), intent(in), dimension(3) :: RP1 ! R of bead i+1
    real(dp), intent(in) :: EB
    real(dp) E_wlc
    real(dp), dimension(3) :: U0, U

    U0 = R - RM1
    U0 = U0/norm2(U0)
    U = RP1 - R
    U = U/norm2(U)
    E_wlc = EB*DOT_PRODUCT(U,U0)
end function E_wlc

function E_SSWLC(R,RM1,U,UM1,EB,EPAR,EPERP,ETA,GAM)
    use params, only: dp
    implicit none
    real(dp), intent(in), dimension(3) :: RM1 ! R of bead i-1
    real(dp), intent(in), dimension(3) :: R ! R of bead i
    real(dp), intent(in), dimension(3) :: UM1 ! U of bead i-1
    real(dp), intent(in), dimension(3) :: U ! U of bead i
    real(dp), intent(in) :: EB, EPAR, EPERP, ETA, GAM ! SSWLC constants
    real(dp), dimension(4) :: E_SSWLC
    real(dp) GI(3)
    real(dp) DR(3)
    real(dp) DRPAR
    real(dp) DRPERP(3)

    DR = R-RM1
    DRPAR = dot_product(DR,UM1)
    DRPERP = DR - DRPAR*UM1
    GI = U - UM1 - ETA*DRPERP
    E_SSWLC(1)=0.5_dp*EB*dot_product(GI,GI)
    E_SSWLC(2)=0.5_dp*EPAR*(DRPAR-GAM)**2
    E_SSWLC(3)=0.5_dp*EPERP*dot_product(DRPERP,DRPERP)
    E_SSWLC(4)=0.0_dp
end function E_SSWLC

function E_GAUSS(R,RM1,EPAR)
    use params, only: dp
    implicit none
    real(dp), intent(in), dimension(3) :: RM1 ! R of bead i-1
    real(dp), intent(in), dimension(3) :: R ! R of bead i
    real(dp), intent(in) :: EPAR
    real(dp) E_GAUSS 
    real(dp) DR(3)
    E_GAUSS = 0.5_dp*EPAR*dot_product(DR,DR)
    ! in gaussian chain, there's only parallel stretching energy. DEELAS init'd to zeros, so sum(DEELAS) == DEELAS(2) later
end function E_GAUSS

end module
