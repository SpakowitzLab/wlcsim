!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      SUBROUTINE energy_elas(EELAS,R,U,NT,NB,NP,PARA)
      use setPrecision
      IMPLICIT NONE 
      INTEGER, intent(in) :: NB           ! Number of beads in a polymer
      INTEGER, intent(in) :: NT           ! Number of beads total
      INTEGER, intent(in) :: NP           ! Number of polymers
      DOUBLE PRECISION, intent(in) :: R(NT,3)  ! Bead positions
      DOUBLE PRECISION, intent(in) :: U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION, intent(out):: EELAS(3) ! Elastic force
      INTEGER I,J,IB            ! Index holders

!     Polymer properties

      DOUBLE PRECISION, intent(in) :: PARA(10)
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA

      DOUBLE PRECISION DR(3),DRPAR,DRPERP(3)
      DOUBLE PRECISION GI(3) !,U1U2


      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)

      EELAS(1)=0.0_dp
      EELAS(2)=0.0_dp
      EELAS(3)=0.0_dp
      IB=1
      DO I=1,NP
         DO J=1,(NB-1)
            DR(1)=R(IB+1,1)-R(IB,1)
            DR(2)=R(IB+1,2)-R(IB,2)
            DR(3)=R(IB+1,3)-R(IB,3)
            DRPAR=DR(1)*U(IB,1)+DR(2)*U(IB,2)+DR(3)*U(IB,3)

            DRPERP(1)=DR(1)-DRPAR*U(IB,1)
            DRPERP(2)=DR(2)-DRPAR*U(IB,2)
            DRPERP(3)=DR(3)-DRPAR*U(IB,3)
            ! U1U2=U(IB,1)*U(IB+1,1)+U(IB,2)*U(IB+1,2)+U(IB,3)*U(IB+1,3)

            GI(1)=(U(IB+1,1)-U(IB,1)-ETA*DRPERP(1))
            GI(2)=(U(IB+1,2)-U(IB,2)-ETA*DRPERP(2))
            GI(3)=(U(IB+1,3)-U(IB,3)-ETA*DRPERP(3))

            EELAS(1) = EELAS(1) + 0.5_dp*EB &
                *(GI(1)**2.0_dp + GI(2)**2.0_dp + GI(3)**2)
            EELAS(2) = EELAS(2) + 0.5_dp*EPAR*(DRPAR - GAM)**2
            EELAS(3) = EELAS(3) + 0.5_dp*EPERP &
                *(DRPERP(1)**2 + DRPERP(2)**2 + DRPERP(3)**2)

            IB=IB+1
         ENDDO
         IB=IB+1
      ENDDO


      RETURN
      END

!---------------------------------------------------------------*
