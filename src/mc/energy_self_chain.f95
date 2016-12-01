!---------------------------------------------------------------*

!
!     This subroutine calcualtes the self-interaction of a single
!     chain with itself. Chain may be linear or a ring.
!
!     The interaction energy is determined used the distance of closest
!     approach between line segments. Energy of segments adjacent along the
!     chain is not calculated.

!     Interaction is calculated for a repulsive lennard jones potential.
!     Arbitrary pair-wise potentials may be substituted

SUBROUTINE ENERGY_SELF_CHAIN(EPONP,R,NT,N,NP,PARA,RING)
    !TODO change so works with multiple chains

    use params, only: dp

    implicit none

  INTEGER N,NT,NP            ! Current number of beads
  DOUBLE PRECISION R(NT,3)   ! Bead positions
  DOUBLE PRECISION EPONP ! Self-interaction force
  DOUBLE PRECISION FMAG     ! Mag of force

  !     Variables for the calculation

  DOUBLE PRECISION U1(3),U2(3),U1U2
  DOUBLE PRECISION D1,D2
  DOUBLE PRECISION R12(3),D12,R12T(3),R12C1(3),R12C2(3)
  DOUBLE PRECISION S1,S2
  DOUBLE PRECISION GI(3)
  INTEGER IT1,IT2,IT1P1,IT2P1

  !     Parameters in the simulation

  real(dp) PARA(10)
  DOUBLE PRECISION LHC      ! HC length
  DOUBLE PRECISION VHC      ! Potential strengths
  DOUBLE PRECISION GAM
  DOUBLE PRECISION LBOX     ! Box edge length
  DOUBLE PRECISION XIR
  INTEGER RING              ! Is polymer a ring?
  INTEGER NMAX


  ! EB=PARA(1)
  ! EPAR=PARA(2)
  ! EPERP=PARA(3)
  GAM=PARA(4)
  ! ETA=PARA(5)
  XIR=PARA(6)
  ! XIU=PARA(7)
  LBOX=PARA(8)
  LHC=PARA(9)
  VHC=PARA(10)


  !     Calculate the self-interaction forces


  IF (RING.EQ.1) THEN
     NMAX=N
  ELSE
     NMAX=N-1
  ENDIF


  EPONP=0.
  DO IT1=3,NMAX
     IF (IT1.EQ.N.AND.RING.EQ.1) THEN
        IT1P1=1
     ELSE
        IT1P1=IT1+1
     ENDIF
     DO IT2=1,IT1-2
        IF (IT2.EQ.N.AND.RING.EQ.1) THEN
           IT2P1=1
        ELSE
           IT2P1=IT2+1
        ENDIF

        IF (IT1P1.EQ.IT2.OR.IT2P1.EQ.IT1) THEN
           GOTO 70
        ENDIF
        R12(1)=R(IT2,1)-R(IT1,1)
        R12(2)=R(IT2,2)-R(IT1,2)
        R12(3)=R(IT2,3)-R(IT1,3)
        ! R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

        D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
        U1(1)=R(IT1P1,1)-R(IT1,1)
        U1(2)=R(IT1P1,2)-R(IT1,2)
        U1(3)=R(IT1P1,3)-R(IT1,3)
        D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
        U1(1)=U1(1)/D1
        U1(2)=U1(2)/D1
        U1(3)=U1(3)/D1

        U2(1)=R(IT2P1,1)-R(IT2,1)
        U2(2)=R(IT2P1,2)-R(IT2,2)
        U2(3)=R(IT2P1,3)-R(IT2,3)
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
           R12T=R(IT2P1,:)-R(IT1P1,:)
           R12C1=R(IT2P1,:)-R(IT1,:)
           R12C2=R(IT1P1,:)-R(IT2,:)
           D12=SQRT(MIN(R12(1)**2+R12(2)**2+R12(3)**2,R12T(1)**2+R12T(2)**2+R12T(3)**2,&
                & R12C1(1)**2+R12C1(2)**2+R12C1(3)**2,R12C2(1)**2+R12C2(2)**2+R12C2(3)**2))
           GOTO 60
        endif

        GI(1)=U2(1)-U1U2*U1(1)
        GI(2)=U2(2)-U1U2*U1(2)
        GI(3)=U2(3)-U1U2*U1(3)

        S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2.GT.D2.OR.S2.LT.0.) then
           R12T=R(IT2P1,:)-R(IT1P1,:)
           R12C1=R(IT2P1,:)-R(IT1,:)
           R12C2=R(IT1P1,:)-R(IT2,:)
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


        EPONP=EPONP+FMAG
70      CONTINUE
     ENDDO
  ENDDO



  RETURN
END SUBROUTINE ENERGY_SELF_CHAIN

!---------------------------------------------------------------*

