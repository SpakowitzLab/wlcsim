!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_reptation(R,U,RP,UP,NT,NB,NP,IP,IT1,IT2,IB1,IB2 &
                  ,rand_stat &
                  ,forward,ring,inTERP_BEAD_LENNARD_JONES)
use mersenne_twister
use params, only: dp

implicit none

integer, intent(in) :: NB     ! Number of beads on a polymer
integer, intent(in) :: NP     ! Number of polymers
integer, intent(in) :: NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,NT)  ! Bead positions
real(dp), intent(in) :: U(3,NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,NT)  ! Bead positions
real(dp), intent(out) :: UP(3,NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IB1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Index of test bead 2
logical, intent(in) :: ring
logical, intent(in) :: inTERP_BEAD_LENNARD_JONES

integer I  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
real(dp) MAG      ! Magnitude of vector
real(dp) DR(3)    ! Displacement for slide move
real(dp) Uvec(3) ! parallel component of triad
real(dp) pDir(3) ! perp component of triad
real(dp) tDir(3) ! twist component of triad
real(dp) r_relative(3) ! r in new coordinate system
real(dp) u_relative(3) ! u in new coordinate system
logical, intent(out) :: forward

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
endif

! single bead reptation
call random_number(urnd,rand_stat)
IP = ceiling(urnd(1)*NP)
IT1 = NB*(IP-1) + 1
IT2 = NB*(IP-1) + NB
IB1 = mod(IT1,NB)
IB2 = mod(IT2,NB)
! move forward or backward
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    forward = .true.
    dR(1) = R(1,IT1 + 1)-R(1,IT1)
    dR(2) = R(2,IT1 + 1)-R(2,IT1)
    dR(3) = R(3,IT1 + 1)-R(3,IT1)

    Uvec(1) = U(1,IT1); Uvec(2) = U(2,IT1); Uvec(3) = U(3,IT1)
    ! chose coordinate system
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! find next r and u in new coordinate system
    u_relative(1) = Uvec(1)*U(1,IT1 + 1) + &
                  Uvec(2)*U(2,IT1 + 1) + &
                  Uvec(3)*U(3,IT1 + 1)
    u_relative(2) = pDir(1)*U(1,IT1 + 1) + &
                  pDir(2)*U(2,IT1 + 1) + &
                  pDir(3)*U(3,IT1 + 1)
    u_relative(3) = tDir(1)*U(1,IT1 + 1) + &
                  tDir(2)*U(2,IT1 + 1) + &
                  tDir(3)*U(3,IT1 + 1)
    r_relative(1) = Uvec(1)*dR(1) + &
                  Uvec(2)*dR(2) + &
                  Uvec(3)*dR(3)
    r_relative(2) = pDir(1)*dR(1) + &
                  pDir(2)*dR(2) + &
                  pDir(3)*dR(3)
    r_relative(3) = tDir(1)*dR(1) + &
                  tDir(2)*dR(2) + &
                  tDir(3)*dR(3)


    ! orient coordinate system with end of chain
    Uvec(1) = U(1,IT2); Uvec(2) = U(2,IT2); Uvec(3) = U(3,IT2)
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! update UP and RP
    UP(1,IT2) = Uvec(1)*u_relative(1) + pDir(1)*u_relative(2) + tDir(1)*u_relative(3)
    UP(2,IT2) = Uvec(2)*u_relative(1) + pDir(2)*u_relative(2) + tDir(2)*u_relative(3)
    UP(3,IT2) = Uvec(3)*u_relative(1) + pDir(3)*u_relative(2) + tDir(3)*u_relative(3)
    mag = sqrt(UP(1,IT2)**2 + UP(2,IT2)**2 + UP(3,IT2)**2)
    UP(1,IT2) = UP(1,IT2)/mag
    UP(2,IT2) = UP(2,IT2)/mag
    UP(3,IT2) = UP(3,IT2)/mag
    RP(1,IT2) = R(1,IT2) + Uvec(1)*r_relative(1) + pDir(1)*r_relative(2) + tDir(1)*r_relative(3)
    RP(2,IT2) = R(2,IT2) + Uvec(2)*r_relative(1) + pDir(2)*r_relative(2) + tDir(2)*r_relative(3)
    RP(3,IT2) = R(3,IT2) + Uvec(3)*r_relative(1) + pDir(3)*r_relative(2) + tDir(3)*r_relative(3)

    do I = IT1,IT2-1
       RP(1,I) = R(1,I + 1)
       RP(2,I) = R(2,I + 1)
       RP(3,I) = R(3,I + 1)
       UP(1,I) = U(1,I + 1)
       UP(2,I) = U(2,I + 1)
       UP(3,I) = U(3,I + 1)
    enddo

   ! RperpMag = sqrt(r_relative(2)**2 + r_relative(3)**2)
   ! RparaMag = r_relative(1)
   ! call test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)

else
    forward = .false.
    dR(1) = R(1,IT2)-R(1,IT2-1)
    dR(2) = R(2,IT2)-R(2,IT2-1)
    dR(3) = R(3,IT2)-R(3,IT2-1)


    Uvec(1) = U(1,IT2); Uvec(2) = U(2,IT2); Uvec(3) = U(3,IT2)
    ! chose coordinate system
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! find next r and u in new coordinate system
    u_relative(1) = Uvec(1)*U(1,IT2-1) + &
                  Uvec(2)*U(2,IT2-1) + &
                  Uvec(3)*U(3,IT2-1)
    u_relative(2) = pDir(1)*U(1,IT2-1) + &
                  pDir(2)*U(2,IT2-1) + &
                  pDir(3)*U(3,IT2-1)
    u_relative(3) = tDir(1)*U(1,IT2-1) + &
                  tDir(2)*U(2,IT2-1) + &
                  tDir(3)*U(3,IT2-1)
    r_relative(1) = Uvec(1)*dR(1) + &
                  Uvec(2)*dR(2) + &
                  Uvec(3)*dR(3)
    r_relative(2) = pDir(1)*dR(1) + &
                  pDir(2)*dR(2) + &
                  pDir(3)*dR(3)
    r_relative(3) = tDir(1)*dR(1) + &
                  tDir(2)*dR(2) + &
                  tDir(3)*dR(3)

    ! orient coordinate system with end of chain
    Uvec(1) = U(1,IT1); Uvec(2) = U(2,IT1); Uvec(3) = U(3,IT1)
    call random_perp(Uvec,pDir,tDir,rand_stat)
    ! update UP and RP
    UP(1,IT1) = Uvec(1)*u_relative(1) + pDir(1)*u_relative(2) + tDir(1)*u_relative(3)
    UP(2,IT1) = Uvec(2)*u_relative(1) + pDir(2)*u_relative(2) + tDir(2)*u_relative(3)
    UP(3,IT1) = Uvec(3)*u_relative(1) + pDir(3)*u_relative(2) + tDir(3)*u_relative(3)
    mag = sqrt(UP(1,IT1)**2 + UP(2,IT1)**2 + UP(3,IT1)**2)
    UP(1,IT1) = UP(1,IT1)/mag
    UP(2,IT1) = UP(2,IT1)/mag
    UP(3,IT1) = UP(3,IT1)/mag
    RP(1,IT1) = R(1,IT1)-Uvec(1)*r_relative(1)-pDir(1)*r_relative(2)-tDir(1)*r_relative(3)
    RP(2,IT1) = R(2,IT1)-Uvec(2)*r_relative(1)-pDir(2)*r_relative(2)-tDir(2)*r_relative(3)
    RP(3,IT1) = R(3,IT1)-Uvec(3)*r_relative(1)-pDir(3)*r_relative(2)-tDir(3)*r_relative(3)

    do I = IT1 + 1,IT2
       RP(1,I) = R(1,I-1)
       RP(2,I) = R(2,I-1)
       RP(3,I) = R(3,I-1)
       UP(1,I) = U(1,I-1)
       UP(2,I) = U(2,I-1)
       UP(3,I) = U(3,I-1)
    enddo
endif

end subroutine
