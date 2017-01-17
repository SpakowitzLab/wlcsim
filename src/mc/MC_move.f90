!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn Made Changes to this file starting on 12/15/15
!
!---------------------------------------------------------------
SUBROUTINE MC_move(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                  ,MCAMP,WINDOW,AB,ABP,BPM,rand_stat,winType &
                  ,IT3,IT4,forward,dib,ring,INTERP_BEAD_LENNARD_JONES)

use mersenne_twister
use params, only: dp, pi

IMPLICIT NONE

INTEGER, intent(in) :: NB     ! Number of beads on a polymer
INTEGER, intent(in) :: NP     ! Number of polymers
INTEGER, intent(in) :: NT     ! Total beads in simulation
DOUBLE PRECISION, intent(in) :: R(NT,3)  ! Bead positions
DOUBLE PRECISION, intent(in) :: U(NT,3)  ! Tangent vectors
DOUBLE PRECISION, intent(out) :: RP(NT,3)  ! Bead positions
DOUBLE PRECISION, intent(out) :: UP(NT,3)  ! Tangent vectors
INTEGER, intent(in) :: BPM    ! Beads per monomer, aka G
INTEGER, intent(out) :: IP    ! Test polymer
INTEGER IP2   ! Second Test polymer if applicable
INTEGER, intent(out) :: IB1   ! Test bead position 1
INTEGER, intent(out) :: IT1   ! Index of test bead 1
INTEGER, intent(out) :: IB2   ! Test bead position 2
INTEGER, intent(out) :: IT2   ! Index of test bead 2
INTEGER, intent(out) :: IT3   ! Test bead position 3 if applicable
INTEGER, intent(out) :: IT4   ! Test bead position 4 if applicable
integer, intent(out) :: dib   ! number of beads moved by move
logical, intent(in) :: ring
logical, intent(in) :: INTERP_BEAD_LENNARD_JONES

INTEGER I,J  ! Test indices
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urand(3)  ! random vector
real urnd(1) ! single random number
! Variables for the crank-shaft move

DOUBLE PRECISION TA(3)    ! Axis of rotation
DOUBLE PRECISION P1(3)    ! Point on rotation line
DOUBLE PRECISION MAG      ! Magnitude of vector
DOUBLE PRECISION ROT(4,4) ! Rotation matrix

DOUBLE PRECISION ALPHA    ! Angle of move
DOUBLE PRECISION BETA     ! Angle of move

!     MC adaptation variables

INTEGER, PARAMETER :: moveTypes=10 ! Number of different move types
DOUBLE PRECISION, intent(in) :: MCAMP(moveTypes) ! Amplitude of random change
INTEGER, intent(in) :: MCTYPE            ! Type of MC move
INTEGER, intent(in) :: winType
Double precision, intent(in) :: WINDOW(moveTypes) ! Size of window for bead selection
DOUBLE PRECISION DR(3)    ! Displacement for slide move
INTEGER TEMP

! Variables for change of binding state move
INTEGER, intent(in) :: AB(NT)            ! Chemical (binding) state
INTEGER, intent(out) :: ABP(NT)          ! Underlying (methalation) state
Double precision d1,d2  !for testing

! variables for reptation move
double precision Uvec(3) ! parallel component of triad
double precision pDir(3) ! perp component of triad
double precision tDir(3) ! twist component of triad
double precision r_relative(3) ! r in new coordinate system
double precision u_relative(3) ! u in new coordinate system
logical, intent(out) :: forward

!TODO saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (RING .OR. INTERP_BEAD_LENNARD_JONES) then
    RP = R
    UP = U
    P1 = 0.0_dp
endif

!     Perform crank-shaft move (MCTYPE 1)

if (MCTYPE.EQ.1) then

   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
   ! IB2=ceiling(urand(2)*NB)
   ! instead of the above, we now choose only one random point and an
   ! exponentially-sized window after that point to move around, to ensure that
   ! small enough sections of chain are moved.
   if (winType.eq.0) then
       IB2=IB1+nint((urand(3)-0.5_dp)*(2.0_dp*WINDOW(MCTYPE)+1.0))
   elseif (winType.eq.1.and..not.RING) then
       call random_number(urnd,rand_stat)
       IB2=IB1+(2*nint(urand(3))-1)* &
               nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
   elseif (winType.eq.1.and.RING) then
       call random_number(urnd,rand_stat)
       IB2=IB1+nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
   endif

   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2

   DIB = IB2-IB1
   if (RING) then                    !Polymer is a ring
      IF (IB2.GT.NB) THEN
         IB2=DIB-(NB-IB1)
      ENDIF
      IT2=NB*(IP-1)+IB2
      if (IB1.EQ.IB2.AND.IB1.EQ.1) then
         TA(1)=R(IT1+1,1)-R(NB*IP,1)
         TA(2)=R(IT1+1,2)-R(NB*IP,2)
         TA(3)=R(IT1+1,3)-R(NB*IP,3)
      elseif (IB1.EQ.IB2.AND.IB1.EQ.NB) then
         TA(1)=R(NB*(IP-1)+1,1)-R(IT1-1,1)
         TA(2)=R(NB*(IP-1)+1,2)-R(IT1-1,2)
         TA(3)=R(NB*(IP-1)+1,3)-R(IT1-1,3)
      elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.NB) then
         TA(1)=R(IT1+1,1)-R(IT1-1,1)
         TA(2)=R(IT1+1,2)-R(IT1-1,2)
         TA(3)=R(IT1+1,3)-R(IT1-1,3)
      else
         TA(1)=R(IT2,1)-R(IT1,1)
         TA(2)=R(IT2,2)-R(IT1,2)
         TA(3)=R(IT2,3)-R(IT1,3)
      endif
   else                                 !Polymer is not a ring
      if (IB2.GT.NB) then
         IB2 = NB
      elseif (IB2.LT.1) then
         IB2 = 1
      endif
      IT2 = NB*(IP-1)+IB2

      if (IT1.GT.IT2) then
         TEMP=IT1
         IT1=IT2
         IT2=TEMP
         TEMP=IB1
         IB1=IB2
         IB2=TEMP
      endif
      DIB = IB2-IB1

      if (IB1.EQ.IB2.AND.IB1.EQ.1) then
         TA(1)=R(IT1+1,1)-R(IT1,1)
         TA(2)=R(IT1+1,2)-R(IT1,2)
         TA(3)=R(IT1+1,3)-R(IT1,3)
      elseif (IB1.EQ.IB2.AND.IB1.EQ.NB) then
         TA(1)=R(NB*IP,1)-R(NB*IP-1,1)
         TA(2)=R(NB*IP,2)-R(NB*IP-1,2)
         TA(3)=R(NB*IP,3)-R(NB*IP-1,3)
      elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.NB) then
         TA(1)=R(IT1+1,1)-R(IT1-1,1)
         TA(2)=R(IT1+1,2)-R(IT1-1,2)
         TA(3)=R(IT1+1,3)-R(IT1-1,3)
      else
         TA(1)=R(IT2,1)-R(IT1,1)
         TA(2)=R(IT2,2)-R(IT1,2)
         TA(3)=R(IT2,3)-R(IT1,3)
      endif
   endif


     MAG=sqrt(TA(1)**2.+TA(2)**2.+TA(3)**2.)
     TA(1)=TA(1)/MAG
     TA(2)=TA(2)/MAG
     TA(3)=TA(3)/MAG
     P1(1)=R(IT1,1)
     P1(2)=R(IT1,2)
     P1(3)=R(IT1,3)
     call random_number(urand,rand_stat)
     ALPHA=MCAMP(1)*(urand(1)-0.5)

     ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
     ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
     ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
          -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

     ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
     ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
     ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
          -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

     ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
     ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
     ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
     ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
          -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

     I=IT1

     DO  J=0,DIB
        if (I.EQ.(NB*IP+1).AND.RING) then
           I=NB*(IP-1)+1
        endif
        RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
        RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
        RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
        UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
        UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
        UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
        I=I+1


     ENDDO

  !  ------begining testing---------
  if(.false.) then
      ! This is a code block for testing
      if (abs(RP(IT1,1)-R(IT1,1)).gt.0.000001) then
          print*, "error in crank-shaft move"
          print*, RP(IT1,1), R(IT1,1)
          stop 1
      endif
      if (abs(RP(IT2,1)-R(IT2,1)).gt.0.000001) then
          print*, "error in crank-shaft move"
          print*, RP(IT1,1), R(IT1,1)
          stop 1
      endif
      if(IT1.ne.IT2) then
          d1=(R(IT1+1,1)-R(IT1,1))**2+&
             (R(IT1+1,2)-R(IT1,2))**2+&
             (R(IT1+1,3)-R(IT1,3))**2
          d2=(RP(IT1+1,1)-RP(IT1,1))**2+&
             (RP(IT1+1,2)-RP(IT1,2))**2+&
             (RP(IT1+1,3)-RP(IT1,3))**2
          if (abs(d1-d2).gt.0.000001) then
              print*, "error in crank-shaft move"
              print*, "distance change in 1"
              print*, "IT1",IT1," IT2",IT2
              print*, d1,d2
              stop 1
          endif
          d1=(R(IT2-1,1)-R(IT2,1))**2+&
             (R(IT2-1,2)-R(IT2,2))**2+&
             (R(IT2-1,3)-R(IT2,3))**2
          d2=(RP(IT2-1,1)-RP(IT2,1))**2+&
             (RP(IT2-1,2)-RP(IT2,2))**2+&
             (RP(IT2-1,3)-RP(IT2,3))**2
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

elseif (MCTYPE.EQ.2) then
   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
   ! again, we use a window
   if (winType.eq.0) then
       IB2=IB1+nint((urand(3)-0.5_dp)*(2.0_dp*WINDOW(MCTYPE)+1.0))
   elseif (winType.eq.1.and..not.RING) then
       call random_number(urnd,rand_stat)
       IB2=IB1+(2*nint(urand(3))-1)* &
               nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
   elseif (winType.eq.1.and.RING) then
       call random_number(urnd,rand_stat)
       IB2=IB1+nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))

   endif

   DIB=IB2-IB1

   if (RING) then
    if (IB2.GT.NB) then
        IB2=DIB-(NB-IB1)
    endif
   else
    if (IB2.GT.NB) then
        IB2=NB
    endif
    if (IB2.LT.1) then
       IB2=1
    endif
    if (IB2.LT.IB1) then
        TEMP=IB1
        IB1=IB2
        IB2=TEMP
    endif
    IT2=NB*(IP-1)+IB2
    DIB = IB2-IB1
   endif

   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2

   call random_number(urand,rand_stat)
   DR(1)=MCAMP(2)*(urand(1)-0.5)
   DR(2)=MCAMP(2)*(urand(2)-0.5)
   DR(3)=MCAMP(2)*(urand(3)-0.5)

     I=IT1
     DO  J=0,DIB

        if (I.EQ.(NB*IP+1).AND.RING) then
           I=NB*(IP-1)+1
        endif

        RP(I,1)=R(I,1)+DR(1)
        RP(I,2)=R(I,2)+DR(2)
        RP(I,3)=R(I,3)+DR(3)
        UP(I,1)=U(I,1)
        UP(I,2)=U(I,2)
        UP(I,3)=U(I,3)
        I=I+1

     ENDDO
! We don't have to protect moves 4-10 with if ring because the code is identical in both cases
!     Perform pivot move (MCTYPE 3)

elseif (MCTYPE.EQ.3) then

    call random_number(urnd,rand_stat)
    IP=ceiling(urnd(1)*NP)
    call random_number(urnd,rand_stat)
    if (urnd(1).gt.0.5_dp) then
        call random_number(urnd,rand_stat)
        IB2=nint(-1.0_dp*log(urnd(1))*WINDOW(MCTYPE))+1
        if (IB2.GT.NB) then
            IB2=NB
        endif
        IB1=1
        IT1=NB*(IP-1)+IB1
        IT2=NB*(IP-1)+IB2
        P1(1)=R(IT2,1)
        P1(2)=R(IT2,2)
        P1(3)=R(IT2,3)
    else
        call random_number(urnd,rand_stat)
        IB1=NB-nint(-1.0_dp*log(urnd(1))*WINDOW(MCTYPE))
        if (IB1.LT.1) then
            IB1=1
        endif
        IB2=NB
        IT1=NB*(IP-1)+IB1
        IT2=NB*(IP-1)+IB2
        P1(1)=R(IT1,1)
        P1(2)=R(IT1,2)
        P1(3)=R(IT1,3)
    endif

   call random_number(urand,rand_stat)
   ALPHA=2.*PI*urand(1)
   BETA=acos(2.*urand(2)-1.)
   TA(1)=sin(BETA)*cos(ALPHA)
   TA(2)=sin(BETA)*sin(ALPHA)
   TA(3)=cos(BETA)

   ALPHA=MCAMP(3)*(urand(3)-0.5)

   ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
   ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
   ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
   ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
   -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

   ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
   ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
   ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
   ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
   -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

   ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
   ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
   ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
   ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
   -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

   DO I=IT1,IT2
      RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
      RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
      RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
      UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
      UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
      UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
      ABP(I)=AB(I)
   enddo

!     Perform rotate move (MCTYPE 4)
!     a.k.a. rotate a single bead
elseif (MCTYPE.EQ.4) then

   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
   IB2=IB1
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2

   call random_number(urand,rand_stat)
   ALPHA=2.*PI*urand(1)
   BETA=acos(2.*urand(2)-1.)
   TA(1)=sin(BETA)*cos(ALPHA)
   TA(2)=sin(BETA)*sin(ALPHA)
   TA(3)=cos(BETA)

   ALPHA=MCAMP(4)*(urand(3)-0.5)

   ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
   ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
   ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
   ROT(1,4)=0.0

   ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
   ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
   ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
   ROT(2,4)=0.0

   ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
   ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
   ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
   ROT(3,4)=0.0

   I=IT1
   UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
   UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
   UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
   RP(I,1)=R(I,1)
   RP(I,2)=R(I,2)
   RP(I,3)=R(I,3)
   ABP(I)=AB(I)

!     Perform a full chain rotation

elseif (MCTYPE.EQ.5) then

    call random_number(urand,rand_stat)
    IP=ceiling(urand(1)*NP)
    IB1=1
    IB2=NB
    IT1=NB*(IP-1)+IB1
    IT2=NB*(IP-1)+IB2

    ALPHA=2.0_dp*PI*urand(2)
    BETA=acos(2.0_dp*urand(3)-1.0_dp)
    TA(1)=sin(BETA)*cos(ALPHA)
    TA(2)=sin(BETA)*sin(ALPHA)
    TA(3)=cos(BETA)

    ! use ~central bead to put axes through
    ! you could also use center of mass if you wanted
    P1(1)=R((IT1+IT2)/2,1)
    P1(2)=R((IT1+IT2)/2,2)
    P1(3)=R((IT1+IT2)/2,3)

    call random_number(urnd,rand_stat)
    ALPHA=MCAMP(5)*(urnd(1)-0.5)

    ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
    ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
    ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
    ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
    -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

    ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
    ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
    ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
    ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
    -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

    ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
    ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
    ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
    ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
    -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

    DO I=IT1,IT2
       RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
       RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
       RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
       UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
       UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
       UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
       ABP(I)=AB(I)
    enddo

!     Perform full chain slide move (MCTYPE 6)
elseif (MCTYPE.EQ.6) then

   call random_number(urnd,rand_stat)
   IP=ceiling(urnd(1)*NP)
   IB1=1
   IB2=NB
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2

   call random_number(urand,rand_stat)
   DR(1)=MCAMP(6)*(urand(1)-0.5_dp)
   DR(2)=MCAMP(6)*(urand(2)-0.5_dp)
   DR(3)=MCAMP(6)*(urand(3)-0.5_dp)

   DO I=IT1,IT2
      RP(I,1)=R(I,1)+DR(1)
      RP(I,2)=R(I,2)+DR(2)
      RP(I,3)=R(I,3)+DR(3)
      UP(I,1)=U(I,1)
      UP(I,2)=U(I,2)
      UP(I,3)=U(I,3)
      ABP(I)=AB(I)
   enddo

elseif (MCTYPE.EQ.7) then
   ! Change AB (a.k.a HP1 binding type fore section of polymer)
   ! Move amplitude is ignored for this move type

   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
   call random_number(urnd,rand_stat)
   IB2=IB1+(2*nint(urand(3))-1)* &
           nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))

   if (IB2.LT.1) then
      IB2=1
   endif
   if (IB2.GT.NB) then
      IB2=NB
   endif

   if (IB2.LT.IB1) then
      TEMP=IB1
      IB1=IB2
      IB2=TEMP
   endif
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2

   !keep binding constant within monomers
   IT1=IT1-MOD(IT1-1,BPM)
   IT2=IT2-MOD(IT2-1,BPM)+BPM-1

   DO J=IT1,IT2
       ABP(J)=1-AB(J)
   ENDDO

   !This loop may not be necessary
   DO I=IT1,IT2
      RP(I,1)=R(I,1)
      RP(I,2)=R(I,2)
      RP(I,3)=R(I,3)
      UP(I,1)=U(I,1)
      UP(I,2)=U(I,2)
      UP(I,3)=U(I,3)
   ENDDO

! chain flip move
elseif (MCTYPE.EQ.8) then
   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=1
   IB2=NB
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2
   DO I=0,NB-1
      RP(IT1+I,1)=R(IT2-I,1)
      RP(IT1+I,2)=R(IT2-I,2)
      RP(IT1+I,3)=R(IT2-I,3)
      UP(IT1+I,1)=-U(IT2-I,1)
      UP(IT1+I,2)=-U(IT2-I,2)
      UP(IT1+I,3)=-U(IT2-I,3)
      ABP(IT1+I)=AB(IT1+I)
   ENDDO
! switch two chains
elseif(MCTYPE.EQ.9) then
   call random_number(urnd,rand_stat)
   IP=ceiling(urnd(1)*NP)
   call random_number(urnd,rand_stat)
   IP2=ceiling(urnd(1)*NP)
   ! Don't switch a chain with itself
   if (IP.eq.IP2) then
       IP2=IP-1
       if (IP2.eq.0) then
           IP2=2
       endif
   endif
   IT1=NB*(IP-1)+1
   IT2=NB*(IP-1)+NB
   IT3=NB*(IP2-1)+1
   IT4=NB*(IP2-1)+NB
   DO I=0,NB-1
      RP(IT1+I,1)=R(IT3+I,1)
      RP(IT1+I,2)=R(IT3+I,2)
      RP(IT1+I,3)=R(IT3+I,3)
      UP(IT1+I,1)=U(IT3+I,1)
      UP(IT1+I,2)=U(IT3+I,2)
      UP(IT1+I,3)=U(IT3+I,3)
      ABP(IT1+I)=AB(IT1+I)
      RP(IT3+I,1)=R(IT1+I,1)
      RP(IT3+I,2)=R(IT1+I,2)
      RP(IT3+I,3)=R(IT1+I,3)
      UP(IT3+I,1)=U(IT1+I,1)
      UP(IT3+I,2)=U(IT1+I,2)
      UP(IT3+I,3)=U(IT1+I,3)
      ABP(IT3+I)=AB(IT3+I)
   ENDDO
   IB1=-2000000
   IB2=-2000000

! single bead reptation
elseif(MCTYPE.EQ.10) then
    call random_number(urnd,rand_stat)
    IP=ceiling(urnd(1)*NP)
    IT1=NB*(IP-1)+1
    IT2=NB*(IP-1)+NB
    ! move forward or backward
    call random_number(urnd,rand_stat)
    if (urnd(1).lt.0.5_dp) then
        forward=.true.
        dR(1)=R(IT1+1,1)-R(IT1,1)
        dR(2)=R(IT1+1,2)-R(IT1,2)
        dR(3)=R(IT1+1,3)-R(IT1,3)

        Uvec(1)=U(IT1,1); Uvec(2)=U(IT1,2); Uvec(3)=U(IT1,3)
        ! chose coordinate system
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! find next r and u in new coordinate system
        u_relative(1)=Uvec(1)*u(IT1+1,1)+&
                      Uvec(2)*u(IT1+1,2)+&
                      Uvec(3)*u(IT1+1,3)
        u_relative(2)=pDir(1)*u(IT1+1,1)+&
                      pDir(2)*u(IT1+1,2)+&
                      pDir(3)*u(IT1+1,3)
        u_relative(3)=tDir(1)*u(IT1+1,1)+&
                      tDir(2)*u(IT1+1,2)+&
                      tDir(3)*u(IT1+1,3)
        r_relative(1)=Uvec(1)*dR(1)+&
                      Uvec(2)*dR(2)+&
                      Uvec(3)*dR(3)
        r_relative(2)=pDir(1)*dR(1)+&
                      pDir(2)*dR(2)+&
                      pDir(3)*dR(3)
        r_relative(3)=tDir(1)*dR(1)+&
                      tDir(2)*dR(2)+&
                      tDir(3)*dR(3)


        ! orient coordinate system with end of chain
        Uvec(1)=U(IT2,1); Uvec(2)=U(IT2,2); Uvec(3)=U(IT2,3)
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! update UP and RP
        UP(IT2,1)=Uvec(1)*u_relative(1)+pDir(1)*u_relative(2)+tDir(1)*u_relative(3)
        UP(IT2,2)=Uvec(2)*u_relative(1)+pDir(2)*u_relative(2)+tDir(2)*u_relative(3)
        UP(IT2,3)=Uvec(3)*u_relative(1)+pDir(3)*u_relative(2)+tDir(3)*u_relative(3)
        mag=sqrt(UP(IT2,1)**2+UP(IT2,2)**2+UP(IT2,3)**2)
        UP(IT2,1)=UP(IT2,1)/mag
        UP(IT2,2)=UP(IT2,2)/mag
        UP(IT2,3)=UP(IT2,3)/mag
        RP(IT2,1)=R(IT2,1)+Uvec(1)*r_relative(1)+pDir(1)*r_relative(2)+tDir(1)*r_relative(3)
        RP(IT2,2)=R(IT2,2)+Uvec(2)*r_relative(1)+pDir(2)*r_relative(2)+tDir(2)*r_relative(3)
        RP(IT2,3)=R(IT2,3)+Uvec(3)*r_relative(1)+pDir(3)*r_relative(2)+tDir(3)*r_relative(3)

        DO I=IT1,IT2-1
           RP(I,1)=R(I+1,1)
           RP(I,2)=R(I+1,2)
           RP(I,3)=R(I+1,3)
           UP(I,1)=U(I+1,1)
           UP(I,2)=U(I+1,2)
           UP(I,3)=U(I+1,3)
        enddo

       ! RperpMag=sqrt(r_relative(2)**2+r_relative(3)**2)
       ! RparaMag=r_relative(1)
       ! call test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)

    else
        forward=.false.
        dR(1)=R(IT2,1)-R(IT2-1,1)
        dR(2)=R(IT2,2)-R(IT2-1,2)
        dR(3)=R(IT2,3)-R(IT2-1,3)


        Uvec(1)=U(IT2,1); Uvec(2)=U(IT2,2); Uvec(3)=U(IT2,3)
        ! chose coordinate system
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! find next r and u in new coordinate system
        u_relative(1)=Uvec(1)*u(IT2-1,1)+&
                      Uvec(2)*u(IT2-1,2)+&
                      Uvec(3)*u(IT2-1,3)
        u_relative(2)=pDir(1)*u(IT2-1,1)+&
                      pDir(2)*u(IT2-1,2)+&
                      pDir(3)*u(IT2-1,3)
        u_relative(3)=tDir(1)*u(IT2-1,1)+&
                      tDir(2)*u(IT2-1,2)+&
                      tDir(3)*u(IT2-1,3)
        r_relative(1)=Uvec(1)*dR(1)+&
                      Uvec(2)*dR(2)+&
                      Uvec(3)*dR(3)
        r_relative(2)=pDir(1)*dR(1)+&
                      pDir(2)*dR(2)+&
                      pDir(3)*dR(3)
        r_relative(3)=tDir(1)*dR(1)+&
                      tDir(2)*dR(2)+&
                      tDir(3)*dR(3)

        ! orient coordinate system with end of chain
        Uvec(1)=U(IT1,1); Uvec(2)=U(IT1,2); Uvec(3)=U(IT1,3)
        call random_perp(Uvec,pDir,tDir,rand_stat)
        ! update UP and RP
        UP(IT1,1)=Uvec(1)*u_relative(1)+pDir(1)*u_relative(2)+tDir(1)*u_relative(3)
        UP(IT1,2)=Uvec(2)*u_relative(1)+pDir(2)*u_relative(2)+tDir(2)*u_relative(3)
        UP(IT1,3)=Uvec(3)*u_relative(1)+pDir(3)*u_relative(2)+tDir(3)*u_relative(3)
        mag=sqrt(UP(IT1,1)**2+UP(IT1,2)**2+UP(IT1,3)**2)
        UP(IT1,1)=UP(IT1,1)/mag
        UP(IT1,2)=UP(IT1,2)/mag
        UP(IT1,3)=UP(IT1,3)/mag
        RP(IT1,1)=R(IT1,1)-Uvec(1)*r_relative(1)-pDir(1)*r_relative(2)-tDir(1)*r_relative(3)
        RP(IT1,2)=R(IT1,2)-Uvec(2)*r_relative(1)-pDir(2)*r_relative(2)-tDir(2)*r_relative(3)
        RP(IT1,3)=R(IT1,3)-Uvec(3)*r_relative(1)-pDir(3)*r_relative(2)-tDir(3)*r_relative(3)

        DO I=IT1+1,IT2
           RP(I,1)=R(I-1,1)
           RP(I,2)=R(I-1,2)
           RP(I,3)=R(I-1,3)
           UP(I,1)=U(I-1,1)
           UP(I,2)=U(I-1,2)
           UP(I,3)=U(I-1,3)
        enddo
    endif
    DO J=IT1,IT2
        ABP(J)=1-AB(J)
    ENDDO
endif

RETURN
END
subroutine test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)
use params, only: dp
IMPLICIT NONE
! inputs
INTEGER NT,IT1,IT2
DOUBLE PRECISION R(NT,3)  ! Bead positions
DOUBLE PRECISION U(NT,3)  ! Tangent vectors
DOUBLE PRECISION RP(NT,3)  ! Bead positions
DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
DOUBLE PRECISION RparaMag, RperpMag

!defined
double precision drOld(3)
double precision drNew(3)
double precision drParOld, drParNew
double precision drPerpOld(3)
double precision drPerpNew(3)
double precision Eta
double precision GIOld(3)
double precision GINew(3)
Eta=1.89756278_dp

drOld(1)=R(IT1+1,1)-R(IT1,1)
drOld(2)=R(IT1+1,2)-R(IT1,2)
drOld(3)=R(IT1+1,3)-R(IT1,3)
DRPAROld=DROld(1)*U(IT1,1)+DROld(2)*U(IT1,2)+DROld(3)*U(IT1,3)
drNew(1)=RP(IT2,1)-RP(IT2-1,1)
drNew(2)=RP(IT2,2)-RP(IT2-1,2)
drNew(3)=RP(IT2,3)-RP(IT2-1,3)
DRPARNew=DRNew(1)*UP(IT2-1,1)+&
         DRNew(2)*UP(IT2-1,2)+&
         DRNew(3)*UP(IT2-1,3)
if (abs(drOld(1)**2+drOld(2)**2+drOld(3)**2&
      -(drNew(1)**2+drNew(2)**2+drNew(3)**2)).gt.0.000001) then
      print*, "drOld",drOld, " mag^2=",drOld(1)**2+drOld(2)**2+drOld(3)**2
      print*, "drNew",drNew, " mag^2=",drNew(1)**2+drNew(2)**2+drNew(3)**2
      print*, "Difference detected in test_equiv, 0"
      stop 1
endif

if (abs(drParOld-drParNew).gt.0.0000001_dp) then
    print*, "DRParOld",DRParOld,"DRParNew",DRParNew
    print*, "Difference detected in test_equiv, 1"
    stop 1
endif

drPerpOld(1)=drOld(1)-drParOld*U(IT1,1)
drPerpOld(2)=drOld(2)-drParOld*U(IT1,2)
drPerpOld(3)=drOld(3)-drParOld*U(IT1,3)
drPerpNew(1)=drNew(1)-drParNew*UP(IT2-1,1)
drPerpNew(2)=drNew(2)-drParNew*UP(IT2-1,2)
drPerpNew(3)=drNew(3)-drParNew*UP(IT2-1,3)

if (abs(drPerpOld(1)**2+drPerpOld(2)**2+drPerpOld(3)**2 &
      -(drPerpNew(1)**2+drPerpNew(2)**2+drPerpNew(3)**2)).gt.0.000001_dp) then
  print*, "drOld",sqrt(drOld(1)**2+drOld(2)**2+drOld(3)**2)
  print*, "drNew",sqrt(drNew(1)**2+drNew(2)**2+drNew(3)**2)
  print*, "dRparOld",dRparOld,"dRparNew",drParNew
  print*, "perp Old:", drPerpOld(1)**2+drPerpOld(2)**2+drPerpOld(3)**2
  print*, "perp New:", drPerpNew(1)**2+drPerpNew(2)**2+drPerpNew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  print*, "Difference detected in test_equiv, 2"
  stop 1
endif

GIOld(1)=U(IT1+1,1)-U(IT1,1)-Eta*dRperpOld(1)
GIOld(2)=U(IT1+1,2)-U(IT1,2)-Eta*dRperpOld(2)
GIOld(3)=U(IT1+1,3)-U(IT1,3)-Eta*dRperpOld(3)
GINew(1)=UP(IT2,1)-UP(IT2-1,1)-Eta*dRperpNew(1)
GINew(2)=UP(IT2,2)-UP(IT2-1,2)-Eta*dRperpNew(2)
GINew(3)=UP(IT2,3)-UP(IT2-1,3)-Eta*dRperpNew(3)

if (abs(GIOld(1)**2+GIOld(2)**2+GIOld(3)**2&
      -(GINew(1)**2+GINew(2)**2+GINew(3)**2)).gt.0.000001_dp) then
  print*, "Difference detected in test_equiv, 3"
  print*, "GIOld(1)**2+GIOld(2)**2+GIOld(3)**2", &
           GIOld(1)**2+GIOld(2)**2+GIOld(3)**2
  print*, "GINew(1)**2+GINew(2)**2+GINew(3)**2", &
          GINew(1)**2+GINew(2)**2+GINew(3)**2
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
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi
type(random_stat) rand_stat  ! status of random number generator
real urnd(1) ! single random number

double precision v(2) ! random 2-vec
double precision, intent(in) :: u(3) ! input
double precision, intent(out) :: p(3) ! output: random perpendicular to u
double precision, intent(out) :: t(3) ! orthogonal to p and u
double precision f

if (abs(u(1)**2+u(2)**2+u(3)**2-1.0_dp) .gt. 0.0000001_dp) then
    print*, u
    print*, "Error in random_perp, please give me a unit vector"
    stop 1
endif

call random_number(urnd,rand_stat)
v(1)=cos(2*PI*urnd(1))
v(2)=sin(2*PI*urnd(1))

if (u(3).gt.0.0) then
    f=1.0_dp/(1+u(3))
    p(1)=(u(3)+f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2)=(u(3)+f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3)=-1.0_dp*(u(2)*v(2)+u(1)*v(1))
else
    f=1.0_dp/(1-u(3))
    p(1)=(-u(3)+f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2)=(-u(3)+f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3)=(u(2)*v(2)+u(1)*v(1))

endif

t(1)=u(2)*p(3)-u(3)*p(2)
t(2)=u(3)*p(1)-u(1)*p(3)
t(3)=u(1)*p(2)-u(2)*p(1)

! random sign
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    t(1)=-1.0_dp*t(1)
    t(2)=-1.0_dp*t(2)
    t(3)=-1.0_dp*t(3)
endif

! Testing
!if (abs(p(1)*u(1)+p(2)*u(2)+p(3)*u(3)).gt.0.000001_dp) then
!    print*, "Error in random_perp, 1"
!    stop 1
!endif
!if (abs(p(1)**2+p(2)**2+p(3)**2-1) .gt. 0.0000001_dp) then
!    print*, "Error in random_perp, 2"
!    stop 1
!endif
!if (abs(t(1)**2 + t(2)**2 + t(3)**2 -1).gt.0.000001_dp) then
!    print*, "Error in random_perp, 3"
!    stop 1
!endif
!if (abs(t(1)*p(1)+t(2)*p(2)+t(3)*p(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 4"
!    stop 1
!endif
!if (abs(t(1)*u(1)+t(2)*u(2)+t(3)*u(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 5"
!    stop 1
!endif
! END Testing

return
end subroutine
!---------------------------------------------------------------!
