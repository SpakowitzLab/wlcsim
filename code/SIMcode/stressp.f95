!! ---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      SUBROUTINE stressp(COR,R,U,R0,U0,NT,N,NP,PARA,INTON,SIMTYPE)

      DOUBLE PRECISION FELAS(NT,3) ! Elastic force
      DOUBLE PRECISION FPONP(NT,3) ! Self-interaction force
      DOUBLE PRECISION TELAS(NT,3) ! Elastic force
      DOUBLE PRECISION FELAS0(NT,3) ! Elastic force
      DOUBLE PRECISION FPONP0(NT,3) ! Self-interaction force
      INTEGER INTON             ! Include polymer interactions
      DOUBLE PRECISION TELAS0(NT,3) ! Elastic force
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION R0(NT,3)  ! Bead positions
      DOUBLE PRECISION U0(NT,3)  ! Tangent vectors
      DOUBLE PRECISION COR

      DOUBLE PRECISION L0       ! Bead separation
      DOUBLE PRECISION FTOT(3)  ! Compress force
      INTEGER I,J,IB                 ! Index holders
      INTEGER N,NT,NP           ! Number of bead
      INTEGER SIMTYPE

!     Variables in the simulation

      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION DT

!     Variables for force and torque calculations

      DOUBLE PRECISION RCOM(3)  ! Center of mass
      DOUBLE PRECISION SIG(3,3)
      DOUBLE PRECISION SIG0(3,3)

!     Load the input parameters

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

      COR=0.

      DT=0.0001
      call force_elas(FELAS0,TELAS0,R0,U0,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)
      if (INTON.EQ.1) then
         call force_ponp(FPONP0,R0,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR)
      endif

      call force_elas(FELAS,TELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)
      if (INTON.EQ.1) then
         call force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR)
      endif

      DO 10 I=1,NP
         RCOM(1)=0.
         RCOM(2)=0.
         RCOM(3)=0.
         DO 20 J=1,N
            IB=J+N*(I-1.)
            RCOM(1)=RCOM(1)+R(IB,1)/N
            RCOM(2)=RCOM(2)+R(IB,2)/N
            RCOM(3)=RCOM(3)+R(IB,3)/N
 20      CONTINUE

         SIG(1,1)=0.
         SIG(1,2)=0.
         SIG(1,3)=0.
         SIG(2,1)=0.
         SIG(2,2)=0.
         SIG(2,3)=0.
         SIG(3,1)=0.
         SIG(3,2)=0.
         SIG(3,3)=0.
         DO 30 J=1,N
            IB=J+N*(I-1.)
            FTOT(1)=FELAS(IB,1)+INTON*FPONP(IB,1)
            FTOT(2)=FELAS(IB,2)+INTON*FPONP(IB,2)
            FTOT(3)=FELAS(IB,3)+INTON*FPONP(IB,3)
            SIG(1,1)=SIG(1,1)-(R(IB,1)-RCOM(1))*FTOT(1)
            SIG(1,2)=SIG(1,2)-(R(IB,1)-RCOM(1))*FTOT(2)
            SIG(1,3)=SIG(1,3)-(R(IB,1)-RCOM(1))*FTOT(3)
            SIG(2,1)=SIG(2,1)-(R(IB,2)-RCOM(2))*FTOT(1)
            SIG(2,2)=SIG(2,2)-(R(IB,2)-RCOM(2))*FTOT(2)
            SIG(2,3)=SIG(2,3)-(R(IB,2)-RCOM(2))*FTOT(3)
            SIG(3,1)=SIG(3,1)-(R(IB,3)-RCOM(3))*FTOT(1)
            SIG(3,2)=SIG(3,2)-(R(IB,3)-RCOM(3))*FTOT(2)
            SIG(3,3)=SIG(3,3)-(R(IB,3)-RCOM(3))*FTOT(3)
 30      CONTINUE

         RCOM(1)=0.
         RCOM(2)=0.
         RCOM(3)=0.
         DO 40 J=1,N
            IB=J+N*(I-1.)
            RCOM(1)=RCOM(1)+R0(IB,1)/N
            RCOM(2)=RCOM(2)+R0(IB,2)/N
            RCOM(3)=RCOM(3)+R0(IB,3)/N
 40      CONTINUE

         SIG0(1,1)=0.
         SIG0(1,2)=0.
         SIG0(1,3)=0.
         SIG0(2,1)=0.
         SIG0(2,2)=0.
         SIG0(2,3)=0.
         SIG0(3,1)=0.
         SIG0(3,2)=0.
         SIG0(3,3)=0.
         DO 50 J=1,N
            IB=J+N*(I-1.)
            FTOT(1)=FELAS0(IB,1)+INTON*FPONP0(IB,1)
            FTOT(2)=FELAS0(IB,2)+INTON*FPONP0(IB,2)
            FTOT(3)=FELAS0(IB,3)+INTON*FPONP0(IB,3)
            SIG0(1,1)=SIG0(1,1)-(R0(IB,1)-RCOM(1))*FTOT(1)
            SIG0(1,2)=SIG0(1,2)-(R0(IB,1)-RCOM(1))*FTOT(2)
            SIG0(1,3)=SIG0(1,3)-(R0(IB,1)-RCOM(1))*FTOT(3)
            SIG0(2,1)=SIG0(2,1)-(R0(IB,2)-RCOM(2))*FTOT(1)
            SIG0(2,2)=SIG0(2,2)-(R0(IB,2)-RCOM(2))*FTOT(2)
            SIG0(2,3)=SIG0(2,3)-(R0(IB,2)-RCOM(2))*FTOT(3)
            SIG0(3,1)=SIG0(3,1)-(R0(IB,3)-RCOM(3))*FTOT(1)
            SIG0(3,2)=SIG0(3,2)-(R0(IB,3)-RCOM(3))*FTOT(2)
            SIG0(3,3)=SIG0(3,3)-(R0(IB,3)-RCOM(3))*FTOT(3)
 50      CONTINUE

         COR=COR+SIG(1,2)*SIG0(1,2)+SIG(1,3)*SIG0(1,3)+SIG(2,1)*SIG0(2,1)+SIG(2,3)*SIG0(2,3)+SIG(3,1)*SIG0(3,1)+SIG(3,2)*SIG0(3,2)

 10   CONTINUE

      COR=COR/6.

      RETURN
      END

!---------------------------------------------------------------*
