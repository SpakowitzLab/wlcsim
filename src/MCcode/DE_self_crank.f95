!---------------------------------------------------------------*

!
!     This subroutine calculates the change in self-interaction of a DNA-like
!     molecule as a chain that interacts through a repulsive lennard-jones potential
!     
!     The change in self-interaction energy associated with a crankshaft move
!     is calculated.
!
!     The program takes the terminal bead indices of the segment that is 
!     rotated by the crankshaft move.
!
!
!     The interaction energy is determined used the distance of closest
!     approach between line segments. Energy of segments adjacent along the 
!     chain is not calculated.


SUBROUTINE DE_SELF_CRANK(DE,R,RP,NT,N,NP,PARA,RING,IB1,IB2)

  IMPLICIT NONE
  INTEGER N,NT,NP            ! Current number of beads
  DOUBLE PRECISION R(NT,3)   ! Bead positions
  DOUBLE PRECISION RP(NT,3)  ! Test bead positions
  DOUBLE PRECISION DE       ! Change in self-energy
  DOUBLE PRECISION E        ! Self-energy before move
  DOUBLE PRECISION EP       ! Self-energy of test polymer
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
  INTEGER IB1,IB2
  INTEGER IO,II,IOP1,IIP1
  INTEGER DIO,DII

  !     Parameters in the simulation

  DOUBLE PRECISION PARA(10)      
  DOUBLE PRECISION LHC      ! HC length
  DOUBLE PRECISION VHC 	! Potential strengths
  DOUBLE PRECISION GAM
  DOUBLE PRECISION LBOX     ! Box edge length
  DOUBLE PRECISION SUM
  DOUBLE PRECISION DT
  INTEGER RING              ! Is polymer a ring?
  INTEGER NMAX     

  DOUBLE PRECISION D12MIN,FMAGMIN


  GAM=PARA(4)
  LBOX=PARA(8)
  LHC=PARA(9)
  VHC=PARA(10)


  ! Determine number of segments inside segment moved and outside

  IF (RING.EQ.1) THEN
     IF (IB2.GE.IB1) THEN
        DII=IB2-IB1
     ELSE
        DII=(N-IB1)+IB2
     ENDIF
  ELSE
     DII=IB2-IB1
  ENDIF
  DIO=N-DII



  !///////////////////////////////////////////////////////////////////////
  !//////////////////////////////////////////////////////////////////////
  !Calculate the interaction energy for original polymer (before move)
  ! Calculate only the portion of the energy that changes during the move
  ! i.e. the interaction between the segment moved and the segment not moved
  !////////////////////////////////////////////////////////////////////////
  !///////////////////////////////////////////////////////////////////////

  E=0.
  II=IB1
  DO I=1,DII
     IF (II.EQ.N.AND.RING.EQ.1) THEN
        IIP1=1
     ELSEIF (II.EQ.N+1.AND.RING.EQ.1) THEN
        II=1
        IIP1=II+1
     ELSE
        IIP1=II+1
     ENDIF
     IO=IB2
     DO J=1,DIO
        IF (IO.EQ.N.AND.RING.EQ.1) THEN
           IOP1=1
        ELSEIF (IO.EQ.N+1.AND.RING.EQ.1) THEN
           IO=1
           IOP1=IO+1
        ELSEIF (IO.EQ.N.AND.RING.EQ.0) THEN
           IO=0
           GOTO 70
        ELSE
           IOP1=IO+1
        ENDIF

        IF (IIP1.EQ.IO.OR.IOP1.EQ.II) THEN
           GOTO 70
        ENDIF
        R12(1)=R(IO,1)-R(II,1)
        R12(2)=R(IO,2)-R(II,2)
        R12(3)=R(IO,3)-R(II,3)
        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
        ! if (D12.GT.(3.*GAM)) then
        !    goto 70
        ! endif

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

60      if (D12.GT.LHC) then
           goto 70
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.


        E=E+FMAG
70      CONTINUE
        IO=IO+1
     ENDDO

     II=II+1
  ENDDO

  !/////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  !Calculation the energy for the test polymer (polymer after the move)
  !Calculate only the portion of the energy that changes during the move.
  !////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  EP=0.
  II=IB1
  DO I=1,DII
     IF (II.EQ.N.AND.RING.EQ.1) THEN
        IIP1=1
     ELSEIF (II.EQ.N+1.AND.RING.EQ.1) THEN
        II=1
        IIP1=II+1
     ELSE
        IIP1=II+1
     ENDIF
     IO=IB2
     DO J=1,DIO
        IF (IO.EQ.N.AND.RING.EQ.1) THEN
           IOP1=1
        ELSEIF (IO.EQ.N+1.AND.RING.EQ.1) THEN
           IO=1
           IOP1=IO+1
        ELSEIF (IO.EQ.N.AND.RING.EQ.0) THEN
           IO=0
           GOTO 90
        ELSE
           IOP1=IO+1
        ENDIF

        IF (IIP1.EQ.IO.OR.IOP1.EQ.II) THEN
           GOTO 90
        ENDIF
        R12(1)=RP(IO,1)-RP(II,1)
        R12(2)=RP(IO,2)-RP(II,2)
        R12(3)=RP(IO,3)-RP(II,3)
        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
        ! if (D12.GT.(3.*GAM)) then
        !    goto 90
        ! endif

        U1(1)=RP(IIP1,1)-RP(II,1)
        U1(2)=RP(IIP1,2)-RP(II,2)
        U1(3)=RP(IIP1,3)-RP(II,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=RP(IOP1,1)-RP(IO,1)
        U2(2)=RP(IOP1,2)-RP(IO,2)
        U2(3)=RP(IOP1,3)-RP(IO,3)
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
           R12T=RP(IOP1,:)-RP(IIP1,:)
           R12C1=RP(IOP1,:)-RP(II,:)
           R12C2=RP(IIP1,:)-RP(IO,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 80
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)

        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then
           R12T=RP(IOP1,:)-RP(IIP1,:)
           R12C1=RP(IOP1,:)-RP(II,:)
           R12C2=RP(IIP1,:)-RP(IO,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 80
        endif

        R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
        R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
        R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

80      if (D12.GT.LHC) then
           goto 90
        endif

        FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.


        EP=EP+FMAG
90      CONTINUE
        IO=IO+1
     ENDDO

     II=II+1
  ENDDO


  !Get the energy difference
 
  DE=EP-E


  RETURN
ENDSUBROUTINE DE_SELF_CRANK
      
!---------------------------------------------------------------*
      
