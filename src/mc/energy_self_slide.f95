!---------------------------------------------------------------*

!
!     This subroutine calculates the self-interaction of a DNA-like
!     molecule.
!     The interaction is calculated as a repulsive Lennard Jones interaction
!     with an effective interaction diameter and interaction potential strength

!     This subroutine calcualtes the portion of the self-interaction that changes
!     during a segment slide move of the Monte Carlo simulation.
!
!     The interaction energy is determined used the distance of closest
!     approach between line segments. Energy of segments adjacent along the
!     chain is not calculated.

!     This routine currently only works when there is only one polymer





SUBROUTINE ENERGY_SELF_SLIDE(EPONP,R,NT,N,NP,PARA,RING,IB1,IB2)

  IMPLICIT NONE
  INTEGER N,NT,NP            ! Current number of beads
  DOUBLE PRECISION R(NT,3)   ! Bead positions
  DOUBLE PRECISION EPONP ! Self-interaction force
  DOUBLE PRECISION FMAG     ! Mag of force
  DOUBLE PRECISION RIJ      ! Interbead dist
  DOUBLE PRECISION EIJ(3)   ! Interbead unit vector
  INTEGER I, J              ! Index holders
  INTEGER SKIP              ! Bead skip index

  !     Variables for the calculation

  DOUBLE PRECISION U1(3),U2(3),U1U2
  DOUBLE PRECISION D1,D2
  DOUBLE PRECISION R12(3),D12,E12(3),R12T(3),R12C1(3),R12C2(3)
  DOUBLE PRECISION S1,S2
  DOUBLE PRECISION GI(3)
  INTEGER I1,J1,I2,J2
  INTEGER IB1,IB2,IB1P1,IB2P1
  INTEGER IO,II,IS1,IS2,IOP1,IIP1,IS1P1,IS2P1
  INTEGER DIO,DII,DIB

  !     Parameters in the simulation

  DOUBLE PRECISION PARA(10)
  DOUBLE PRECISION LHC      ! HC length
  DOUBLE PRECISION SIGP     ! HC diameter
  DOUBLE PRECISION VHC 	! Potential strengths
  DOUBLE PRECISION GAM
  DOUBLE PRECISION LBOX     ! Box edge length
  DOUBLE PRECISION SUM
  DOUBLE PRECISION DT
  DOUBLE PRECISION XIR
  DOUBLE PRECISION XIU
  DOUBLE PRECISION ETA
  DOUBLE PRECISION EPAR
  DOUBLE PRECISION EPERP
  DOUBLE PRECISION EB
  INTEGER RING              ! Is polymer a ring?
  INTEGER NMAX

  DOUBLE PRECISION D12MIN,FMAGMIN

  EB=PARA(1)
  EPAR=PARA(2)
  EPERP=PARA(3)
  GAM=PARA(4)
  ETA=PARA(5)
  XIR=PARA(6)
  XIU=PARA(7)
  LBOX=PARA(8)
  LHC=PARA(9)
  VHC=PARA(10)


  EPONP=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Setup. Determine the number of segments that lie
  !betwen the two terminal points of the segment slid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IF (IB2.GE.IB1) THEN
     DIB=IB2-IB1
  ELSE
     DIB=(N-IB1)+IB2
  ENDIF
  DIO=N-DIB-2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !The two segments immediately on the outside of the
  !beads that were slid are "stretched" or compressed.
  !The interaction of these segments with the remainder
  !of the outer segments is changed, and their interaction
  !with the segments inside the region slid is changed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  IS1P1=IB1                               ! Index of first stretched (or compressed) segment (not inside beads moved)
  IS2=IB2                                  ! Index of second stretched (or compressed) segment (not inside beads moved)

  IF (IS1P1.EQ.1) THEN
     IS1=N
  ELSE
     IS1=IS1P1-1
  ENDIF

  IF (IS2.EQ.N) THEN
     IS2P1=1
  ELSE
     IS2P1=IS2+1
  ENDIF

  II=IB1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Calculate the interaction between inner segments with all outer segments
  !The following loop is only performed if the two terminal beads are not equal
  !Note that special cases for end positions are currently only set-up for rings
  !This code can accomodate both linear and ring chains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (IB1.NE.IB2) then
     !Sum over segments inside segment slid (inner segments)
     DO  I=1,DIB

        IF (II.EQ.N.AND.RING.EQ.1) THEN
           IIP1=1
        ELSEIF (II.EQ.N+1) THEN
           II=1
           IIP1=II+1
        ELSE
           IIP1=II+1
        ENDIF
        IO=IB2+1
        !Sum over segments unchanged by slide (outer segments)
        DO J=1,DIO

           IF (IO.EQ.N.AND.RING.EQ.1) THEN
              IOP1=1
           ELSEIF (IO.EQ.N.AND.RING.EQ.0)THEN
              IO=0
              GOTO 110
           ELSEIF (IO.EQ.N+1) THEN
              IO=1
              IOP1=IO+1
           ELSE
              IOP1=IO+1
           ENDIF

           !Skip this pair if the segments are adjacent
           IF (IIP1.EQ.IO.OR.IOP1.EQ.II) THEN
              GOTO 70
           ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !Calculate the interaction energy between the inner segment
           !and outer segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           R12(1)=R(IO,1)-R(II,1)
           R12(2)=R(IO,2)-R(II,2)
           R12(3)=R(IO,3)-R(II,3)
           !Periodic Bounary conditions. Not used for single chains

           ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
           ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
           ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

           D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

           U1(1)=R(IIP1,1)-R(II,1)
           U1(2)=R(IIP1,2)-R(II,2)
           U1(3)=R(IIP1,3)-R(II,3)
           D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
           U1(1)=U1(1)/D1
           U1(2)=U1(2)/D1
           U1(3)=U1(3)/D1

           U2(1)=R(IOP1,1)-R(IO,1)
           U2(2)=R(IOP1,2)-R(IO,2)
           U2(3)=R(IOP1,3)-R(IO,3)
           D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
           U2(1)=U2(1)/D2
           U2(2)=U2(2)/D2
           U2(3)=U2(3)/D2

           U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
           if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
              D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
              GOTO 60
           endif

           GI(1)=U1(1)-U1U2*U2(1)
           GI(2)=U1(2)-U1U2*U2(2)
           GI(3)=U1(3)-U1U2*U2(3)

           S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

           if (S1.GT.D1.OR.S1.LT.0.) then
              R12T=R(IOP1,:)-R(IIP1,:)
              R12C1=R(IOP1,:)-R(II,:)
              R12C2=R(IIP1,:)-R(IO,:)
              D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                   & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
              GOTO 60
           endif

           GI(1)=U2(1)-U1U2*U1(1)
           GI(2)=U2(2)-U1U2*U1(2)
           GI(3)=U2(3)-U1U2*U1(3)

           S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

           if (S2.GT.D2.OR.S2.LT.0.) then
              R12T=R(IOP1,:)-R(IIP1,:)
              R12C1=R(IOP1,:)-R(II,:)
              R12C2=R(IIP1,:)-R(IO,:)
              D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                   & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
              GOTO 60
           endif

           R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
           R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
           R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

           D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

60         if (D12.GT.LHC) then
              goto 70
           endif

           FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

           EPONP=EPONP+FMAG

70         CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !Calculate the interaction between this outer segment and the
           !two stretched segments. Only do this once.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           IF (I.EQ.1) THEN

              !Skip this pair if the segments are adjacent
              IF (IS1.EQ.IOP1.OR.IO.EQ.IS1P1) THEN
                 GOTO 90
              ENDIF

              !If the chain is linear and IS1P1.EQ.1 then there is no first 'stretched' segment
              IF (IS1P1.EQ.1.AND.RING.EQ.0) THEN
                 GOTO 90
              ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !Calculate interaction between outer segment and first stretched
              !segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              R12(1)=R(IS1,1)-R(IO,1)
              R12(2)=R(IS1,2)-R(IO,2)
              R12(3)=R(IS1,3)-R(IO,3)
              !Periodic Bounary conditions. Not used for single chains

              ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
              ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
              ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

              D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

              U1(1)=R(IOP1,1)-R(IO,1)
              U1(2)=R(IOP1,2)-R(IO,2)
              U1(3)=R(IOP1,3)-R(IO,3)
              D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
              U1(1)=U1(1)/D1
              U1(2)=U1(2)/D1
              U1(3)=U1(3)/D1

              U2(1)=R(IS1P1,1)-R(IS1,1)
              U2(2)=R(IS1P1,2)-R(IS1,2)
              U2(3)=R(IS1P1,3)-R(IS1,3)
              D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
              U2(1)=U2(1)/D2
              U2(2)=U2(2)/D2
              U2(3)=U2(3)/D2

              U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
              if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
                 D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
                 GOTO 80
              endif

              GI(1)=U1(1)-U1U2*U2(1)
              GI(2)=U1(2)-U1U2*U2(2)
              GI(3)=U1(3)-U1U2*U2(3)

              S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

              if (S1.GT.D1.OR.S1.LT.0.) then

                 R12T=R(IOP1,:)-R(IS1P1,:)
                 R12C1=R(IOP1,:)-R(IS1,:)
                 R12C2=R(IS1P1,:)-R(IO,:)

                 D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                      & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
                 GOTO 80
              endif

              GI(1)=U2(1)-U1U2*U1(1)
              GI(2)=U2(2)-U1U2*U1(2)
              GI(3)=U2(3)-U1U2*U1(3)

              S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

              if (S2.GT.D2.OR.S2.LT.0.) then

                 R12T=R(IOP1,:)-R(IS1P1,:)
                 R12C1=R(IOP1,:)-R(IS1,:)
                 R12C2=R(IS1P1,:)-R(IO,:)

                 D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                      & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
                 GOTO 80
              endif

              R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
              R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
              R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

              D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

80            if (D12.GT.LHC) then
                 goto 90
              endif

              FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

              EPONP=EPONP+FMAG
90            CONTINUE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !Calculate the interaction between the outer segment and the second
              !stretched segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !Skip this pair if the segments are adjacent
              IF (IS2P1.EQ.IO.OR.IOP1.EQ.IS2) THEN
                 GOTO 110
              ENDIF

              !If the chain is linear and IS2=N, then there is no second 'stretched' segment
              IF (IS2.EQ.N.AND.RING.EQ.0) THEN
                 GOTO 110
              ENDIF

              R12(1)=R(IS2,1)-R(IO,1)
              R12(2)=R(IS2,2)-R(IO,2)
              R12(3)=R(IS2,3)-R(IO,3)
              !Periodic Bounary conditions. Not used for single chains

              ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
              ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
              ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

              D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

              U1(1)=R(IOP1,1)-R(IO,1)
              U1(2)=R(IOP1,2)-R(IO,2)
              U1(3)=R(IOP1,3)-R(IO,3)
              D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
              U1(1)=U1(1)/D1
              U1(2)=U1(2)/D1
              U1(3)=U1(3)/D1

              U2(1)=R(IS2P1,1)-R(IS2,1)
              U2(2)=R(IS2P1,2)-R(IS2,2)
              U2(3)=R(IS2P1,3)-R(IS2,3)
              D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
              U2(1)=U2(1)/D2
              U2(2)=U2(2)/D2
              U2(3)=U2(3)/D2

              U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
              if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
                 D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
                 GOTO 100
              endif

              GI(1)=U1(1)-U1U2*U2(1)
              GI(2)=U1(2)-U1U2*U2(2)
              GI(3)=U1(3)-U1U2*U2(3)

              S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

              if (S1.GT.D1.OR.S1.LT.0.) then

                 R12T=R(IOP1,:)-R(IS2P1,:)
                 R12C1=R(IOP1,:)-R(IS2,:)
                 R12C2=R(IS2P1,:)-R(IO,:)

                 D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                      & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
                 GOTO 100
              endif

              GI(1)=U2(1)-U1U2*U1(1)
              GI(2)=U2(2)-U1U2*U1(2)
              GI(3)=U2(3)-U1U2*U1(3)

              S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

              if (S2.GT.D2.OR.S2.LT.0.) then
                 R12T=R(IOP1,:)-R(IS2P1,:)
                 R12C1=R(IOP1,:)-R(IS2,:)
                 R12C2=R(IS2P1,:)-R(IO,:)
                 D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                      & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
                 GOTO 100
              endif

              R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
              R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
              R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

              D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

100           if (D12.GT.LHC) then
                 goto 110
              endif

              FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.



              EPONP=EPONP+FMAG


           ENDIF
           !step to next outer segment
110        CONTINUE
           IO=IO+1
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Calculate the interaction between the inner segment and the first stretched
        !segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Skip this pair if the segments are adjacent
        IF (IS1.EQ.IIP1.OR.II.EQ.IS1P1) THEN
           GOTO 130
        ENDIF

        !If the chain is linear and  IS1P1=1 there is no first "stretched" segment
        IF (IS1P1.EQ.1.AND.RING.EQ.0) THEN
           GOTO 130
        ENDIF

        R12(1)=R(IS1,1)-R(II,1)
        R12(2)=R(IS1,2)-R(II,2)
        R12(3)=R(IS1,3)-R(II,3)
        !Periodic Bounary conditions. Not used for single chains

        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

        U1(1)=R(IIP1,1)-R(II,1)
        U1(2)=R(IIP1,2)-R(II,2)
        U1(3)=R(IIP1,3)-R(II,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=R(IS1P1,1)-R(IS1,1)
        U2(2)=R(IS1P1,2)-R(IS1,2)
        U2(3)=R(IS1P1,3)-R(IS1,3)
        D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
        U2(1)=U2(1)/D2
        U2(2)=U2(2)/D2
        U2(3)=U2(3)/D2

        U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
        if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
           D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
           GOTO 120
        endif

        GI(1)=U1(1)-U1U2*U2(1)
        GI(2)=U1(2)-U1U2*U2(2)
        GI(3)=U1(3)-U1U2*U2(3)

        S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S1.GT.D1.OR.S1.LT.0.) then
           R12T=R(IS1P1,:)-R(IIP1,:)
           R12C1=R(IS1P1,:)-R(II,:)
           R12C2=R(IIP1,:)-R(IS1,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 120
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)

        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then

           R12T=R(IS1P1,:)-R(IIP1,:)
           R12C1=R(IS1P1,:)-R(II,:)
           R12C2=R(IIP1,:)-R(IS1,:)

           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 120
        endif

        R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
        R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
        R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

120     if (D12.GT.LHC) then
           goto 130
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

        EPONP=EPONP+FMAG
130     CONTINUE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Calculate interaction between inner segment and second stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF (IS2.EQ.IIP1.OR.II.EQ.IS2P1) THEN
           GOTO 150
        ENDIF

        !If the chain is linear and IS2=N, there is no second "stretched" segment
        IF (IS2.EQ.N.AND.RING.EQ.0) THEN
           GOTO 150
        ENDIF

        R12(1)=R(IS2,1)-R(II,1)
        R12(2)=R(IS2,2)-R(II,2)
        R12(3)=R(IS2,3)-R(II,3)
        !Periodic Bounary conditions. Not used for single chains

        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

        U1(1)=R(IIP1,1)-R(II,1)
        U1(2)=R(IIP1,2)-R(II,2)
        U1(3)=R(IIP1,3)-R(II,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=R(IS2P1,1)-R(IS2,1)
        U2(2)=R(IS2P1,2)-R(IS2,2)
        U2(3)=R(IS2P1,3)-R(IS2,3)
        D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
        U2(1)=U2(1)/D2
        U2(2)=U2(2)/D2
        U2(3)=U2(3)/D2

        U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
        if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
           D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
           GOTO 140
        endif

        GI(1)=U1(1)-U1U2*U2(1)
        GI(2)=U1(2)-U1U2*U2(2)
        GI(3)=U1(3)-U1U2*U2(3)

        S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S1.GT.D1.OR.S1.LT.0.) then
           R12T=R(IS2P1,:)-R(IIP1,:)
           R12C1=R(IS2P1,:)-R(II,:)
           R12C2=R(IIP1,:)-R(IS2,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 140
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)

        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then
           R12T=R(IS2P1,:)-R(IIP1,:)
           R12C1=R(IS2P1,:)-R(II,:)
           R12C2=R(IIP1,:)-R(IS2,:)

           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 140
        endif

        R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
        R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
        R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

140     if (D12.GT.LHC) then
           goto 150
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

        EPONP=EPONP+FMAG
150     CONTINUE



        II=II+1
     ENDDO
  ENDIF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Calculate the interaction between the two stretched segments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  R12(1)=R(IS2,1)-R(IS1,1)
  R12(2)=R(IS2,2)-R(IS1,2)
  R12(3)=R(IS2,3)-R(IS1,3)


  IF (IS1.EQ.IS2P1.OR.IS1P1.EQ.IS2) THEN
     GOTO 170
  ENDIF

  !If the chain is linear and IS1P1=1 or IS2=N, then there is no
  !first segment or second stretched segment

  IF (IS1P1.EQ.1.OR.IS2.EQ.N) THEN
     GOTO 170
  ENDIF

  !Periodic Bounary conditions. Not used for single chains

  ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
  ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
  ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

  D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

  U1(1)=R(IS1P1,1)-R(IS1,1)
  U1(2)=R(IS1P1,2)-R(IS1,2)
  U1(3)=R(IS1P1,3)-R(IS1,3)
  D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
  U1(1)=U1(1)/D1
  U1(2)=U1(2)/D1
  U1(3)=U1(3)/D1

  U2(1)=R(IS2P1,1)-R(IS2,1)
  U2(2)=R(IS2P1,2)-R(IS2,2)
  U2(3)=R(IS2P1,3)-R(IS2,3)
  D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
  U2(1)=U2(1)/D2
  U2(2)=U2(2)/D2
  U2(3)=U2(3)/D2

  U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
  if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
     D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
     GOTO 160
  endif

  GI(1)=U1(1)-U1U2*U2(1)
  GI(2)=U1(2)-U1U2*U2(2)
  GI(3)=U1(3)-U1U2*U2(3)

  S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

  if (S1.GT.D1.OR.S1.LT.0.) then
     R12T=R(IS2P1,:)-R(IS1P1,:)
     R12C1=R(IS2P1,:)-R(IS1,:)
     R12C2=R(IS1P1,:)-R(IS2,:)
     D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
          & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
     GOTO 160
  endif

  GI(1)=U2(1)-U1U2*U1(1)
  GI(2)=U2(2)-U1U2*U1(2)
  GI(3)=U2(3)-U1U2*U1(3)

  S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

  if (S2.GT.D2.OR.S2.LT.0.) then
     R12T=R(IS2P1,:)-R(IS1P1,:)
     R12C1=R(IS2P1,:)-R(IS1,:)
     R12C2=R(IS1P1,:)-R(IS2,:)
     D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
          & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
     GOTO 160
  endif

  R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
  R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
  R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

  D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

160 if (D12.GT.LHC) then
     goto 170
  endif

  FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

  EPONP=EPONP+FMAG


170 CONTINUE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !If only one bead was slid (IB1.EQ.IB2) then calculate the interaction between the
  !portion of the chain unaffected by the move and the two segments stretched/compressed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (IB1.EQ.IB2) THEN
     !If IB1.EQ.IB2, Get interaction between segment unaffected by move and two segments affected

     IO=IB1+1
     !Sum over segments unchanged by slide (outer segments)
     DO J=1,DIO

        IF (IO.EQ.N.AND.RING.EQ.1) THEN
           IOP1=1
        ELSEIF (IO.EQ.N.AND.RING.EQ.0) THEN
           IO=0
           GOTO 210
        ELSEIF (IO.EQ.N+1) THEN
           IO=1
           IOP1=IO+1
        ELSE
           IOP1=IO+1
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Calculate interaction between outer segment and first stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF (IS1.EQ.IOP1.OR.IO.EQ.IS1P1) THEN
           GOTO 190
        ENDIF

        !If the chain is linear and IP1P1=1, there is no first "stretched" segment
        IF (IS1P1.EQ.1.AND.RING.EQ.0) THEN
           GOTO 190
        ENDIF

        R12(1)=R(IS1,1)-R(IO,1)
        R12(2)=R(IS1,2)-R(IO,2)
        R12(3)=R(IS1,3)-R(IO,3)
        !Periodic Bounary conditions. Not used for single chains

        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX


        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

        U1(1)=R(IOP1,1)-R(IO,1)
        U1(2)=R(IOP1,2)-R(IO,2)
        U1(3)=R(IOP1,3)-R(IO,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=R(IS1P1,1)-R(IS1,1)
        U2(2)=R(IS1P1,2)-R(IS1,2)
        U2(3)=R(IS1P1,3)-R(IS1,3)
        D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
        U2(1)=U2(1)/D2
        U2(2)=U2(2)/D2
        U2(3)=U2(3)/D2

        U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
        if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
           D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
           GOTO 180
        endif

        GI(1)=U1(1)-U1U2*U2(1)
        GI(2)=U1(2)-U1U2*U2(2)
        GI(3)=U1(3)-U1U2*U2(3)

        S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)



        if (S1.GT.D1.OR.S1.LT.0.) then

           R12T=R(IOP1,:)-R(IS1P1,:)
           R12C1=R(IOP1,:)-R(IS1,:)
           R12C2=R(IS1P1,:)-R(IO,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 180
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)

        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then

           R12T=R(IOP1,:)-R(IS1P1,:)
           R12C1=R(IOP1,:)-R(IS1,:)
           R12C2=R(IS1P1,:)-R(IO,:)

           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 180
        endif



        R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
        R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
        R12(3)=R12(3)+S2*U2(3)-S1*U1(3)


        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)


180     if (D12.GT.LHC) then
           goto 190
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

        EPONP=EPONP+FMAG

190     CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Calculate the interaction between outer segment and second stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF (IS2.EQ.IOP1.OR.IO.EQ.IS2P1) THEN
           GOTO 210
        ENDIF

        !If the chain is linear and IS2=N, then there is no second "stretched" segment

        IF (IS2.EQ.N.AND.RING.EQ.0) THEN
           GOTO 210
        ENDIF

        R12(1)=R(IS2,1)-R(IO,1)
        R12(2)=R(IS2,2)-R(IO,2)
        R12(3)=R(IS2,3)-R(IO,3)


        !Periodic Bounary conditions. Not used for single chains

        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

        U1(1)=R(IOP1,1)-R(IO,1)
        U1(2)=R(IOP1,2)-R(IO,2)
        U1(3)=R(IOP1,3)-R(IO,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=R(IS2P1,1)-R(IS2,1)
        U2(2)=R(IS2P1,2)-R(IS2,2)
        U2(3)=R(IS2P1,3)-R(IS2,3)
        D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
        U2(1)=U2(1)/D2
        U2(2)=U2(2)/D2
        U2(3)=U2(3)/D2

        U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
        if (U1U2.EQ.1..OR.U1U2.EQ.-1.) then
           D12=SQRT(R12(1)**2+R12(2)**2+R12(3)**1)
           GOTO 200
        endif

        GI(1)=U1(1)-U1U2*U2(1)
        GI(2)=U1(2)-U1U2*U2(2)
        GI(3)=U1(3)-U1U2*U2(3)


        S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)


        if (S1.GT.D1.OR.S1.LT.0.) then
           R12T=R(IOP1,:)-R(IS2P1,:)
           R12C1=R(IOP1,:)-R(IS2,:)
           R12C2=R(IS2P1,:)-R(IO,:)


           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 200
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)


        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then
           R12T=R(IOP1,:)-R(IS2P1,:)
           R12C1=R(IOP1,:)-R(IS2,:)
           R12C2=R(IS2P1,:)-R(IO,:)


           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 200
        endif

        R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
        R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
        R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

200     if (D12.GT.LHC) then
           goto 210
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

        EPONP=EPONP+FMAG
210     CONTINUE


        IO=IO+1
     ENDDO
  ENDIF


  RETURN
ENDSUBROUTINE ENERGY_SELF_SLIDE

!---------------------------------------------------------------*

