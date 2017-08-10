!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn split out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_crank(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                  ,MCAMP,WindoW,rand_stat,winType &
                  ,dib,ring,inTERP_BEAD_LENNARD_JONES)

use mersenne_twister
use params, only: dp

implicit none

integer, intent(in) :: NB     ! Number of beads on a polymer
integer, intent(in) :: NP     ! Number of polymers
integer, intent(in) :: NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,NT)  ! Bead positions
real(dp), intent(in) :: U(3,NT)  ! Tangent vectorslet g:syntastic_fortran_compiler_options = '-Wall -pedantic -std=f2003'
real(dp), intent(out) :: RP(3,NT)  ! Bead positions
real(dp), intent(out) :: UP(3,NT)  ! Tangent vectors
integer, intent(out) :: IP    ! Test polymer
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: dib   ! number of beads moved by move
logical, intent(in) :: ring
logical, intent(in) :: inTERP_BEAD_LENNARD_JONES

integer I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
! Variables for the crank-shaft move

real(dp) TA(3)    ! Axis of rotation
real(dp) P1(3)    ! Point on rotation line
real(dp) MAG      ! Magnitude of vector
real(dp) ROT(4,4) ! Rotation matrix

real(dp) ALPHA    ! Angle of move

!     MC adaptation variables

integer, PARAMETER :: moveTypes = 10 ! Number of different move types
real(dp), intent(in) :: MCAMP(moveTypes) ! Amplitude of random change
integer, intent(in) :: MCTYPE            ! Type of MC move
integer, intent(in) :: winType
real(dp), intent(in) :: WindoW(moveTypes) ! Size of window for bead selection
integer TEMP

! Variables for change of binding state move
real(dp) d1,d2  !for testing

! variables for reptation move

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RinG .OR. inTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

!     Perform crank-shaft move (MCTYPE 1)


call random_number(urand,rand_stat)
IP = ceiling(urand(1)*NP)
IB1 = ceiling(urand(2)*NB)
! IB2 = ceiling(urand(2)*NB)
! instead of the above, we now choose only one random point and an
! exponentially-sized window after that point to move around, to ensure that
! small enough sections of chain are moved.
if (winType.eq.0) then
    IB2 = IB1 + nint((urand(3)-0.5_dp)*(2.0_dp*WindoW(MCTYPE) + 1.0))
elseif (winType.eq.1.and..not.RinG) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + (2*nint(urand(3))-1)* &
            nint(-1.0*log(urnd(1))*WindoW(MCTYPE))
elseif (winType.eq.1.and.RinG) then
    call random_number(urnd,rand_stat)
    IB2 = IB1 + nint(-1.0*log(urnd(1))*WindoW(MCTYPE))
else
    call stop_if_err(1, "Warning: winType not recognized")
endif

IT1 = NB*(IP-1) + IB1
IT2 = NB*(IP-1) + IB2

DIB = IB2-IB1
if (RinG) then                    !Polymer is a ring
   if (IB2 > NB) then
      IB2 = DIB-(NB-IB1)
   ENDif
   IT2 = NB*(IP-1) + IB2
   if (IB1 == IB2.AND.IB1 == 1) then
      TA(1) = R(1,IT1 + 1)-R(1,NB*IP)
      TA(2) = R(2,IT1 + 1)-R(2,NB*IP)
      TA(3) = R(3,IT1 + 1)-R(3,NB*IP)
   elseif (IB1 == IB2.AND.IB1 == NB) then
      TA(1) = R(1,NB*(IP-1) + 1)-R(1,IT1-1)
      TA(2) = R(2,NB*(IP-1) + 1)-R(2,IT1-1)
      TA(3) = R(3,NB*(IP-1) + 1)-R(3,IT1-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /= NB) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1-1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1-1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1-1)
   else
      TA(1) = R(1,IT2)-R(1,IT1)
      TA(2) = R(2,IT2)-R(2,IT1)
      TA(3) = R(3,IT2)-R(3,IT1)
   endif
else                                 !Polymer is not a ring
   if (IB2 > NB) then
      IB2 = NB
   elseif (IB2 < 1) then
      IB2 = 1
   endif
   IT2 = NB*(IP-1) + IB2

   if (IT1 > IT2) then
      TEMP = IT1
      IT1 = IT2
      IT2 = TEMP
      TEMP = IB1
      IB1 = IB2
      IB2 = TEMP
   endif
   DIB = IB2-IB1

   if (IB1 == IB2.AND.IB1 == 1) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1)
   elseif (IB1 == IB2.AND.IB1 == NB) then
      TA(1) = R(1,NB*IP)-R(1,NB*IP-1)
      TA(2) = R(2,NB*IP)-R(2,NB*IP-1)
      TA(3) = R(3,NB*IP)-R(3,NB*IP-1)
   elseif (IB1 == IB2.AND.IB1 /= 1.AND.IB2 /= NB) then
      TA(1) = R(1,IT1 + 1)-R(1,IT1-1)
      TA(2) = R(2,IT1 + 1)-R(2,IT1-1)
      TA(3) = R(3,IT1 + 1)-R(3,IT1-1)
   else
      TA(1) = R(1,IT2)-R(1,IT1)
      TA(2) = R(2,IT2)-R(2,IT1)
      TA(3) = R(3,IT2)-R(3,IT1)
   endif
endif


  MAG = sqrt(TA(1)**2. + TA(2)**2. + TA(3)**2.)
  TA(1) = TA(1)/MAG
  TA(2) = TA(2)/MAG
  TA(3) = TA(3)/MAG
  P1(1) = R(1,IT1)
  P1(2) = R(2,IT1)
  P1(3) = R(3,IT1)
  call random_number(urand,rand_stat)
  ALPHA = MCAMP(1)*(urand(1)-0.5)

  ROT(1,1) = TA(1)**2. + (TA(2)**2. + TA(3)**2.)*cos(ALPHA)
  ROT(1,2) = TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
  ROT(1,3) = TA(1)*TA(3)*(1.-cos(ALPHA)) + TA(2)*sin(ALPHA)
  ROT(1,4) = (P1(1)*(1.-TA(1)**2.) &
       -TA(1)*(P1(2)*TA(2) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

  ROT(2,1) = TA(1)*TA(2)*(1.-cos(ALPHA)) + TA(3)*sin(ALPHA)
  ROT(2,2) = TA(2)**2. + (TA(1)**2. + TA(3)**2.)*cos(ALPHA)
  ROT(2,3) = TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
  ROT(2,4) = (P1(2)*(1.-TA(2)**2.) &
       -TA(2)*(P1(1)*TA(1) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

  ROT(3,1) = TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
  ROT(3,2) = TA(2)*TA(3)*(1.-cos(ALPHA)) + TA(1)*sin(ALPHA)
  ROT(3,3) = TA(3)**2. + (TA(1)**2. + TA(2)**2.)*cos(ALPHA)
  ROT(3,4) = (P1(3)*(1.-TA(3)**2.) &
       -TA(3)*(P1(1)*TA(1) + P1(2)*TA(2)))*(1.-cos(ALPHA)) + (P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

  I = IT1

  do  J = 0,DIB
     if (I == (NB*IP + 1).AND.RinG) then
        I = NB*(IP-1) + 1
     endif
     RP(1,I) = ROT(1,4) + ROT(1,1)*R(1,I) + ROT(1,2)*R(2,I) + ROT(1,3)*R(3,I)
     RP(2,I) = ROT(2,4) + ROT(2,1)*R(1,I) + ROT(2,2)*R(2,I) + ROT(2,3)*R(3,I)
     RP(3,I) = ROT(3,4) + ROT(3,1)*R(1,I) + ROT(3,2)*R(2,I) + ROT(3,3)*R(3,I)
     UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
     UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
     UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
     I = I + 1

  ENDdo

!  ------begining testing---------
if(.false.) then
    ! This is a code block for testing
    if (abs(RP(1,IT1)-R(1,IT1)).gt.0.000001) then
        print*, "error in crank-shaft move"
        print*, RP(1,IT1), R(1,IT1)
        stop 1
    endif
    if (abs(RP(1,IT2)-R(1,IT2)).gt.0.000001) then
        print*, "error in crank-shaft move"
        print*, RP(1,IT1), R(1,IT1)
        stop 1
    endif
    if(IT1.ne.IT2) then
        d1 = (R(1,IT1 + 1)-R(1,IT1))**2 + &
           (R(2,IT1 + 1)-R(2,IT1))**2 + &
           (R(3,IT1 + 1)-R(3,IT1))**2
        d2 = (RP(1,IT1 + 1)-RP(1,IT1))**2 + &
           (RP(2,IT1 + 1)-RP(2,IT1))**2 + &
           (RP(3,IT1 + 1)-RP(3,IT1))**2
        if (abs(d1-d2).gt.0.000001) then
            print*, "error in crank-shaft move"
            print*, "distance change in 1"
            print*, "IT1",IT1," IT2",IT2
            print*, d1,d2
            stop 1
        endif
        d1 = (R(1,IT2-1)-R(1,IT2))**2 + &
           (R(2,IT2-1)-R(2,IT2))**2 + &
           (R(3,IT2-1)-R(3,IT2))**2
        d2 = (RP(1,IT2-1)-RP(1,IT2))**2 + &
           (RP(2,IT2-1)-RP(2,IT2))**2 + &
           (RP(3,IT2-1)-RP(3,IT2))**2
        if (abs(d1-d2).gt.0.000001) then
            print*, "error in crank-shaft move"
            print*, "distance change in 2"
            print*, d1,d2
            stop 1
        endif
    endif
endif
! --------end testing--------
end subroutine
