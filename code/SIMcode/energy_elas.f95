!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      SUBROUTINE energy_elas(EELAS,R,U,NT,N,NP,PARA)

      DOUBLE PRECISION EELAS(3) ! Elastic force
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION B(NT)  ! Tangent vectors
      DOUBLE PRECISION UR(NT,3)  ! Tangent vectors
      DOUBLE PRECISION KAP,EPS  ! Elastic props
      DOUBLE PRECISION L0       ! Bead separation
      DOUBLE PRECISION FCOM(3)  ! Compress force
      DOUBLE PRECISION FBEND(3) ! Bend force
      INTEGER I,J,IB            ! Index holders
      INTEGER N,NT,NP           ! Number of bead

!     Polymer properties

      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA

!     Variables for force and torque calculations

      DOUBLE PRECISION DR(3),DRPAR,DRPERP(3)
      DOUBLE PRECISION FI(3),TI(3)
      DOUBLE PRECISION U1U2,GI(3),DOTGU,HI(3)

!     Calculate the forces and torques

      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)

      EELAS(1)=0.
      EELAS(2)=0.
      EELAS(3)=0.
      IB=1
      DO 30 I=1,NP
         DO 40 J=1,(N-1)
            DR(1)=R(IB+1,1)-R(IB,1)
            DR(2)=R(IB+1,2)-R(IB,2)
            DR(3)=R(IB+1,3)-R(IB,3)
            DRPAR=DR(1)*U(IB,1)+DR(2)*U(IB,2)+DR(3)*U(IB,3)

            DRPERP(1)=DR(1)-DRPAR*U(IB,1)
            DRPERP(2)=DR(2)-DRPAR*U(IB,2)
            DRPERP(3)=DR(3)-DRPAR*U(IB,3)
            U1U2=U(IB,1)*U(IB+1,1)+U(IB,2)*U(IB+1,2)+U(IB,3)*U(IB+1,3)

            GI(1)=(U(IB+1,1)-U(IB,1)-ETA*DRPERP(1))
            GI(2)=(U(IB+1,2)-U(IB,2)-ETA*DRPERP(2))
            GI(3)=(U(IB+1,3)-U(IB,3)-ETA*DRPERP(3))

            EELAS(1)=EELAS(1)+0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.)
            EELAS(2)=EELAS(2)+0.5*EPAR*(DRPAR-GAM)**2.
            EELAS(3)=EELAS(3)+0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

            IB=IB+1
 40      CONTINUE
         IB=IB+1
 30   CONTINUE


      RETURN
      END

!---------------------------------------------------------------*
