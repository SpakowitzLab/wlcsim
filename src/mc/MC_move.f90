!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn Made Changes to this file starting on 12/15/15
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_move(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                  ,MCAMP,WindoW,BPM,rand_stat,winType &
                  ,IT3,IT4,forward,dib,ring,inTERP_BEAD_LENNARD_JONES &
                  ,wlc_d)

use mersenne_twister
use params, only: dp, pi, wlcsim_data

!TODO: replace R,U,RP,UP .... with wlc_d

implicit none

!type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: NB     ! Number of beads on a polymer
integer, intent(in) :: NP     ! Number of polymers
integer, intent(in) :: NT     ! Total beads in simulation
real(dp), intent(in) :: R(3,NT)  ! Bead positions
real(dp), intent(in) :: U(3,NT)  ! Tangent vectors
real(dp), intent(out) :: RP(3,NT)  ! Bead positions
real(dp), intent(out) :: UP(3,NT)  ! Tangent vectors
integer, intent(in) :: BPM    ! Beads per monomer, aka G
integer, intent(out) :: IP    ! Test polymer
integer IP2   ! Second Test polymer if applicable
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IT3   ! Test bead position 3 if applicable
integer, intent(out) :: IT4   ! Test bead position 4 if applicable
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
real(dp) BETA     ! Angle of move

!     MC adaptation variables

integer, PARAMETER :: moveTypes = 10 ! Number of different move types
real(dp), intent(in) :: MCAMP(moveTypes) ! Amplitude of random change
integer, intent(in) :: MCTYPE            ! Type of MC move
integer, intent(in) :: winType
real(dp), intent(in) :: WindoW(moveTypes) ! Size of window for bead selection
real(dp) DR(3)    ! Displacement for slide move
integer TEMP

! Variables for change of binding state move
real(dp) d1,d2  !for testing

! variables for reptation move
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
    P1 = 0.0_dp
endif

!     Perform crank-shaft move (MCTYPE 1)

if (MCTYPE == 1) then

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

!     Perform slide move (MCTYPE 2)

elseif (MCTYPE == 2) then
   call random_number(urand,rand_stat)
   IP = ceiling(urand(1)*NP)
   IB1 = ceiling(urand(2)*NB)
   ! again, we use a window
   if (winType.eq.0) then
       IB2 = IB1 + nint((urand(3)-0.5_dp)*(2.0_dp*WindoW(MCTYPE) + 1.0))
   elseif (winType.eq.1.and..not.RinG) then
       call random_number(urnd,rand_stat)
       IB2 = IB1 + (2*nint(urand(3))-1)* &
               nint(-1.0*log(urnd(1))*WindoW(MCTYPE))
   elseif (winType.eq.1.and.RinG) then
       call random_number(urnd,rand_stat)
       IB2 = IB1 + nint(-1.0*log(urnd(1))*WindoW(MCTYPE))

   endif

   DIB = IB2-IB1

   if (RinG) then
    if (IB2 > NB) then
        IB2 = DIB-(NB-IB1)
    endif
   else
    if (IB2 > NB) then
        IB2 = NB
    endif
    if (IB2 < 1) then
       IB2 = 1
    endif
    if (IB2 < IB1) then
        TEMP = IB1
        IB1 = IB2
        IB2 = TEMP
    endif
    IT2 = NB*(IP-1) + IB2
    DIB = IB2-IB1
   endif

   IT1 = NB*(IP-1) + IB1
   IT2 = NB*(IP-1) + IB2

   call random_number(urand,rand_stat)
   DR(1) = MCAMP(2)*(urand(1)-0.5)
   DR(2) = MCAMP(2)*(urand(2)-0.5)
   DR(3) = MCAMP(2)*(urand(3)-0.5)

     I = IT1
     do  J = 0,DIB

        if (I == (NB*IP + 1).AND.RinG) then
           I = NB*(IP-1) + 1
        endif

        RP(1,I) = R(1,I) + DR(1)
        RP(2,I) = R(2,I) + DR(2)
        RP(3,I) = R(3,I) + DR(3)
        UP(1,I) = U(1,I)
        UP(2,I) = U(2,I)
        UP(3,I) = U(3,I)
        I = I + 1

     ENDdo
! We don't have to protect moves 4-10 with if ring because the code is identical in both cases
!     Perform pivot move (MCTYPE 3)

elseif (MCTYPE == 3) then

    call random_number(urnd,rand_stat)
    IP = ceiling(urnd(1)*NP)
    call random_number(urnd,rand_stat)
    if (urnd(1).gt.0.5_dp) then
        call random_number(urnd,rand_stat)
        IB2 = nint(-1.0_dp*log(urnd(1))*WindoW(MCTYPE)) + 1
        if (IB2 > NB) then
            IB2 = NB
        endif
        IB1 = 1
        IT1 = NB*(IP-1) + IB1
        IT2 = NB*(IP-1) + IB2
        P1(1) = R(1,IT2)
        P1(2) = R(2,IT2)
        P1(3) = R(3,IT2)
    else
        call random_number(urnd,rand_stat)
        IB1 = NB-nint(-1.0_dp*log(urnd(1))*WindoW(MCTYPE))
        if (IB1 < 1) then
            IB1 = 1
        endif
        IB2 = NB
        IT1 = NB*(IP-1) + IB1
        IT2 = NB*(IP-1) + IB2
        P1(1) = R(1,IT1)
        P1(2) = R(2,IT1)
        P1(3) = R(3,IT1)
    endif

   call random_number(urand,rand_stat)
   ALPHA = 2.*PI*urand(1)
   BETA = acos(2.*urand(2)-1.)
   TA(1) = sin(BETA)*cos(ALPHA)
   TA(2) = sin(BETA)*sin(ALPHA)
   TA(3) = cos(BETA)

   ALPHA = MCAMP(3)*(urand(3)-0.5)

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

   do I = IT1,IT2
      RP(1,I) = ROT(1,4) + ROT(1,1)*R(1,I) + ROT(1,2)*R(2,I) + ROT(1,3)*R(3,I)
      RP(2,I) = ROT(2,4) + ROT(2,1)*R(1,I) + ROT(2,2)*R(2,I) + ROT(2,3)*R(3,I)
      RP(3,I) = ROT(3,4) + ROT(3,1)*R(1,I) + ROT(3,2)*R(2,I) + ROT(3,3)*R(3,I)
      UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
      UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
      UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
   enddo

!     Perform rotate move (MCTYPE 4)
!     a.k.a. rotate a single bead
elseif (MCTYPE == 4) then

   call random_number(urand,rand_stat)
   IP = ceiling(urand(1)*NP)
   IB1 = ceiling(urand(2)*NB)
   IB2 = IB1
   IT1 = NB*(IP-1) + IB1
   IT2 = NB*(IP-1) + IB2

   call random_number(urand,rand_stat)
   ALPHA = 2.*PI*urand(1)
   BETA = acos(2.*urand(2)-1.)
   TA(1) = sin(BETA)*cos(ALPHA)
   TA(2) = sin(BETA)*sin(ALPHA)
   TA(3) = cos(BETA)

   ALPHA = MCAMP(4)*(urand(3)-0.5)

   ROT(1,1) = TA(1)**2. + (TA(2)**2. + TA(3)**2.)*cos(ALPHA)
   ROT(1,2) = TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
   ROT(1,3) = TA(1)*TA(3)*(1.-cos(ALPHA)) + TA(2)*sin(ALPHA)
   ROT(1,4) = 0.0

   ROT(2,1) = TA(1)*TA(2)*(1.-cos(ALPHA)) + TA(3)*sin(ALPHA)
   ROT(2,2) = TA(2)**2. + (TA(1)**2. + TA(3)**2.)*cos(ALPHA)
   ROT(2,3) = TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
   ROT(2,4) = 0.0

   ROT(3,1) = TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
   ROT(3,2) = TA(2)*TA(3)*(1.-cos(ALPHA)) + TA(1)*sin(ALPHA)
   ROT(3,3) = TA(3)**2. + (TA(1)**2. + TA(2)**2.)*cos(ALPHA)
   ROT(3,4) = 0.0

   I = IT1
   UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
   UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
   UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
   RP(1,I) = R(1,I)
   RP(2,I) = R(2,I)
   RP(3,I) = R(3,I)

!     Perform a full chain rotation

elseif (MCTYPE == 5) then

    call random_number(urand,rand_stat)
    IP = ceiling(urand(1)*NP)
    IB1 = 1
    IB2 = NB
    IT1 = NB*(IP-1) + IB1
    IT2 = NB*(IP-1) + IB2

    ALPHA = 2.0_dp*PI*urand(2)
    BETA = acos(2.0_dp*urand(3)-1.0_dp)
    TA(1) = sin(BETA)*cos(ALPHA)
    TA(2) = sin(BETA)*sin(ALPHA)
    TA(3) = cos(BETA)

    ! use ~central bead to put axes through
    ! you could also use center of mass if you wanted
    P1(1) = R(1,(IT1 + IT2)/2)
    P1(2) = R(2,(IT1 + IT2)/2)
    P1(3) = R(3,(IT1 + IT2)/2)

    call random_number(urnd,rand_stat)
    ALPHA = MCAMP(5)*(urnd(1)-0.5)

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

    do I = IT1,IT2
       RP(1,I) = ROT(1,4) + ROT(1,1)*R(1,I) + ROT(1,2)*R(2,I) + ROT(1,3)*R(3,I)
       RP(2,I) = ROT(2,4) + ROT(2,1)*R(1,I) + ROT(2,2)*R(2,I) + ROT(2,3)*R(3,I)
       RP(3,I) = ROT(3,4) + ROT(3,1)*R(1,I) + ROT(3,2)*R(2,I) + ROT(3,3)*R(3,I)
       UP(1,I) = ROT(1,1)*U(1,I) + ROT(1,2)*U(2,I) + ROT(1,3)*U(3,I)
       UP(2,I) = ROT(2,1)*U(1,I) + ROT(2,2)*U(2,I) + ROT(2,3)*U(3,I)
       UP(3,I) = ROT(3,1)*U(1,I) + ROT(3,2)*U(2,I) + ROT(3,3)*U(3,I)
    enddo

!     Perform full chain slide move (MCTYPE 6)
elseif (MCTYPE == 6) then

   call random_number(urnd,rand_stat)
   IP = ceiling(urnd(1)*NP)
   IB1 = 1
   IB2 = NB
   IT1 = NB*(IP-1) + IB1
   IT2 = NB*(IP-1) + IB2

   call random_number(urand,rand_stat)
   DR(1) = MCAMP(6)*(urand(1)-0.5_dp)
   DR(2) = MCAMP(6)*(urand(2)-0.5_dp)
   DR(3) = MCAMP(6)*(urand(3)-0.5_dp)

   do I = IT1,IT2
      RP(1,I) = R(1,I) + DR(1)
      RP(2,I) = R(2,I) + DR(2)
      RP(3,I) = R(3,I) + DR(3)
      UP(1,I) = U(1,I)
      UP(2,I) = U(2,I)
      UP(3,I) = U(3,I)
   enddo

elseif (MCTYPE == 7) then
   ! Change wlc_d%AB (a.k.a HP1 binding type fore section of polymer)
   ! Move amplitude is ignored for this move type
   call random_number(urand,rand_stat)
   IP = ceiling(urand(1)*NP)
   IB1 = ceiling(urand(2)*NB)
   call random_number(urnd,rand_stat)
   IB2 = IB1 + (2*nint(urand(3))-1)* &
           nint(-1.0*log(urnd(1))*WindoW(MCTYPE))

   if (IB2 < 1) then
      IB2 = 1
   endif
   if (IB2 > NB) then
      IB2 = NB
   endif

   if (IB2 < IB1) then
      TEMP = IB1
      IB1 = IB2
      IB2 = TEMP
   endif
   IT1 = NB*(IP-1) + IB1
   IT2 = NB*(IP-1) + IB2

   !keep binding constant within monomers
   IT1 = IT1-MOD(IT1-1,BPM)
   IT2 = IT2-MOD(IT2-1,BPM) + BPM-1

   do J = IT1,IT2
       wlc_d%ABP(J) = 1-wlc_d%AB(J)
   ENDdo

   !This loop may not be necessary
   do I = IT1,IT2
      RP(1,I) = R(1,I)
      RP(2,I) = R(2,I)
      RP(3,I) = R(3,I)
      UP(1,I) = U(1,I)
      UP(2,I) = U(2,I)
      UP(3,I) = U(3,I)
   ENDdo

! chain flip move
elseif (MCTYPE == 8) then
   call random_number(urand,rand_stat)
   IP = ceiling(urand(1)*NP)
   IB1 = 1
   IB2 = NB
   IT1 = NB*(IP-1) + IB1
   IT2 = NB*(IP-1) + IB2
   do I = 0,NB-1
      RP(1,IT1 + I) = R(1,IT2-I)
      RP(2,IT1 + I) = R(2,IT2-I)
      RP(3,IT1 + I) = R(3,IT2-I)
      UP(1,IT1 + I) = -U(1,IT2-I)
      UP(2,IT1 + I) = -U(2,IT2-I)
      UP(3,IT1 + I) = -U(3,IT2-I)
   ENDdo
! switch two chains
elseif(MCTYPE == 9) then
   call random_number(urnd,rand_stat)
   IP = ceiling(urnd(1)*NP)
   call random_number(urnd,rand_stat)
   IP2 = ceiling(urnd(1)*NP)
   ! Don't switch a chain with itself
   if (IP.eq.IP2) then
       IP2 = IP-1
       if (IP2.eq.0) then
           IP2 = 2
       endif
   endif
   IT1 = NB*(IP-1) + 1
   IT2 = NB*(IP-1) + NB
   IT3 = NB*(IP2-1) + 1
   IT4 = NB*(IP2-1) + NB
   do I = 0,NB-1
      RP(1,IT1 + I) = R(1,IT3 + I)
      RP(2,IT1 + I) = R(2,IT3 + I)
      RP(3,IT1 + I) = R(3,IT3 + I)
      UP(1,IT1 + I) = U(1,IT3 + I)
      UP(2,IT1 + I) = U(2,IT3 + I)
      UP(3,IT1 + I) = U(3,IT3 + I)
      RP(1,IT3 + I) = R(1,IT1 + I)
      RP(2,IT3 + I) = R(2,IT1 + I)
      RP(3,IT3 + I) = R(3,IT1 + I)
      UP(1,IT3 + I) = U(1,IT1 + I)
      UP(2,IT3 + I) = U(2,IT1 + I)
      UP(3,IT3 + I) = U(3,IT1 + I)
   ENDdo
   IB1 = -2000000
   IB2 = -2000000

! single bead reptation
elseif(MCTYPE == 10) then
    call random_number(urnd,rand_stat)
    IP = ceiling(urnd(1)*NP)
    IT1 = NB*(IP-1) + 1
    IT2 = NB*(IP-1) + NB
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
endif

RETURN
END
subroutine test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)
use params, only: dp
implicit none
! inputs
integer NT,IT1,IT2
real(dp) R(3,NT)  ! Bead positions
real(dp) U(3,NT)  ! Tangent vectors
real(dp) RP(3,NT)  ! Bead positions
real(dp) UP(3,NT)  ! Tangent vectors
real(dp) RparaMag, RperpMag

!defined
real(dp) drOld(3)
real(dp) drNew(3)
real(dp) drParOld, drParNew
real(dp) drPerpOld(3)
real(dp) drPerpNew(3)
real(dp) Eta
real(dp) GIOld(3)
real(dp) Ginew(3)
Eta = 1.89756278_dp

drOld(1) = R(1,IT1 + 1)-R(1,IT1)
drOld(2) = R(2,IT1 + 1)-R(2,IT1)
drOld(3) = R(3,IT1 + 1)-R(3,IT1)
DRPAROld = DROld(1)*U(1,IT1) + DROld(2)*U(2,IT1) + DROld(3)*U(3,IT1)
drNew(1) = RP(1,IT2)-RP(1,IT2-1)
drNew(2) = RP(2,IT2)-RP(2,IT2-1)
drNew(3) = RP(3,IT2)-RP(3,IT2-1)
DRPARNew = DRNew(1)*UP(1,IT2-1) + &
         DRNew(2)*UP(2,IT2-1) + &
         DRNew(3)*UP(3,IT2-1)
if (abs(drOld(1)**2 + drOld(2)**2 + drOld(3)**2&
      -(drNew(1)**2 + drNew(2)**2 + drNew(3)**2)).gt.0.000001) then
      print*, "drOld",drOld, " mag^2 = ",drOld(1)**2 + drOld(2)**2 + drOld(3)**2
      print*, "drNew",drNew, " mag^2 = ",drNew(1)**2 + drNew(2)**2 + drNew(3)**2
      print*, "Difference detected in test_equiv, 0"
      stop 1
endif

if (abs(drParOld-drParNew).gt.0.0000001_dp) then
    print*, "DRParOld",DRParOld,"DRParNew",DRParNew
    print*, "Difference detected in test_equiv, 1"
    stop 1
endif

drPerpOld(1) = drOld(1)-drParOld*U(1,IT1)
drPerpOld(2) = drOld(2)-drParOld*U(2,IT1)
drPerpOld(3) = drOld(3)-drParOld*U(3,IT1)
drPerpNew(1) = drNew(1)-drParNew*UP(1,IT2-1)
drPerpNew(2) = drNew(2)-drParNew*UP(2,IT2-1)
drPerpNew(3) = drNew(3)-drParNew*UP(3,IT2-1)

if (abs(drPerpOld(1)**2 + drPerpOld(2)**2 + drPerpOld(3)**2 &
      -(drPerpNew(1)**2 + drPerpNew(2)**2 + drPerpNew(3)**2)).gt.0.000001_dp) then
  print*, "drOld",sqrt(drOld(1)**2 + drOld(2)**2 + drOld(3)**2)
  print*, "drNew",sqrt(drNew(1)**2 + drNew(2)**2 + drNew(3)**2)
  print*, "dRparOld",dRparOld,"dRparNew",drParNew
  print*, "perp Old:", drPerpOld(1)**2 + drPerpOld(2)**2 + drPerpOld(3)**2
  print*, "perp New:", drPerpNew(1)**2 + drPerpNew(2)**2 + drPerpNew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  print*, "Difference detected in test_equiv, 2"
  stop 1
endif

GIOld(1) = U(1,IT1 + 1)-U(1,IT1)-Eta*dRperpOld(1)
GIOld(2) = U(2,IT1 + 1)-U(2,IT1)-Eta*dRperpOld(2)
GIOld(3) = U(3,IT1 + 1)-U(3,IT1)-Eta*dRperpOld(3)
Ginew(1) = UP(1,IT2)-UP(1,IT2-1)-Eta*dRperpNew(1)
Ginew(2) = UP(2,IT2)-UP(2,IT2-1)-Eta*dRperpNew(2)
Ginew(3) = UP(3,IT2)-UP(3,IT2-1)-Eta*dRperpNew(3)

if (abs(GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2&
      -(Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2)).gt.0.000001_dp) then
  print*, "Difference detected in test_equiv, 3"
  print*, "GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2", &
           GIOld(1)**2 + GIOld(2)**2 + GIOld(3)**2
  print*, "Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2", &
          Ginew(1)**2 + Ginew(2)**2 + Ginew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  stop 1
endif

return
end subroutine
subroutine random_perp(u,p,t,rand_stat)
! The subroutine generates the second two vectors in a unit triad
! The output vectors, p and t, are perpendicular to eachother and u
! The triad is randomly left or right handed
use mersenne_twister
use params, only: dp
implicit none
real(dp), PARAMETER :: PI = 3.141592654 ! Value of pi
type(random_stat) rand_stat  ! status of random number generator
real urnd(1) ! single random number

real(dp) v(2) ! random 2-vec
real(dp), intent(in) :: u(3) ! input
real(dp), intent(out) :: p(3) ! output: random perpendicular to u
real(dp), intent(out) :: t(3) ! orthogonal to p and u
real(dp) f

if (abs(u(1)**2 + u(2)**2 + u(3)**2-1.0_dp) .gt. 0.0000001_dp) then
    print*, u
    print*, "Error in random_perp, please give me a unit vector"
    stop 1
endif

call random_number(urnd,rand_stat)
v(1) = cos(2*PI*urnd(1))
v(2) = sin(2*PI*urnd(1))

if (u(3).gt.0.0) then
    f = 1.0_dp/(1 + u(3))
    p(1) = (u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2) = (u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3) = -1.0_dp*(u(2)*v(2) + u(1)*v(1))
else
    f = 1.0_dp/(1-u(3))
    p(1) = (-u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2) = (-u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3) = (u(2)*v(2) + u(1)*v(1))

endif

t(1) = u(2)*p(3)-u(3)*p(2)
t(2) = u(3)*p(1)-u(1)*p(3)
t(3) = u(1)*p(2)-u(2)*p(1)

! random sign
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    t(1) = -1.0_dp*t(1)
    t(2) = -1.0_dp*t(2)
    t(3) = -1.0_dp*t(3)
endif

! Testing
!if (abs(p(1)*u(1) + p(2)*u(2) + p(3)*u(3)).gt.0.000001_dp) then
!    print*, "Error in random_perp, 1"
!    stop 1
!endif
!if (abs(p(1)**2 + p(2)**2 + p(3)**2-1) .gt. 0.0000001_dp) then
!    print*, "Error in random_perp, 2"
!    stop 1
!endif
!if (abs(t(1)**2 + t(2)**2 + t(3)**2 -1).gt.0.000001_dp) then
!    print*, "Error in random_perp, 3"
!    stop 1
!endif
!if (abs(t(1)*p(1) + t(2)*p(2) + t(3)*p(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 4"
!    stop 1
!endif
!if (abs(t(1)*u(1) + t(2)*u(2) + t(3)*u(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 5"
!    stop 1
!endif
! END Testing

return
end subroutine
!---------------------------------------------------------------!
