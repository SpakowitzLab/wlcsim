!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      subroutine stress(SIG,R,U,NT,N,NP,PARA,inTON,SIMTYPE)

      use params, only : dp

      real(dp) FELAS(NT,3) ! Elastic force
      real(dp) FPONP(NT,3) ! Self-interaction force
      real(dp) TELAS(NT,3) ! Elastic force
      real(dp) R(NT,3)  ! Bead positions
      real(dp) U(NT,3)  ! Tangent vectors
      real(dp) L0       ! Bead separation
      real(dp) FTOT(3)  ! Compress force
      real(dp) FBEND(3) ! Bend force
      integer inTON             ! Include polymer interactions
      integer I,J,IB                 ! Index holders
      integer N,NT,NP           ! Number of bead
      integer SIMTYPE

!     Variables in the simulation

      real(dp) EB,EPAR,EPERP
      real(dp) GAM,ETA
      real(dp) XIR,XIU
      real(dp) LBOX     ! Box edge length
      real(dp) LHC      ! Length of HC int
      real(dp) VHC      ! HC strength
      real(dp) PARA(10)
      real(dp) DT

!     Variables for force and torque calculations

      real(dp) RCOM(3)  ! Center of mass
      real(dp) SIG(3,3)

!     Load the input parameters

      EB = PARA(1)
      EPAR = PARA(2)
      EPERP = PARA(3)
      GAM = PARA(4)
      ETA = PARA(5)
      XIR = PARA(6)
      XIU = PARA(7)
      LBOX = PARA(8)
      LHC = PARA(9)
      VHC = PARA(10)

      DT = 0.0001
      call force_elas(FELAS,TELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)

      if (inTON == 1) then
         call force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR)
      endif

      SIG(1,1) = 0.
      SIG(1,2) = 0.
      SIG(1,3) = 0.
      SIG(2,1) = 0.
      SIG(2,2) = 0.
      SIG(2,3) = 0.
      SIG(3,1) = 0.
      SIG(3,2) = 0.
      SIG(3,3) = 0.
      do 10 I = 1,NP
         RCOM(1) = 0.
         RCOM(2) = 0.
         RCOM(3) = 0.
         do 20 J = 1,N
            IB = J + N*(I-1)
            RCOM(1) = RCOM(1) + R(IB,1)/N
            RCOM(2) = RCOM(2) + R(IB,2)/N
            RCOM(3) = RCOM(3) + R(IB,3)/N
 20      continue

         do 30 J = 1,N
            IB = J + N*(I-1)
            FTOT(1) = FELAS(IB,1) + inTON*FPONP(IB,1)
            FTOT(2) = FELAS(IB,2) + inTON*FPONP(IB,2)
            FTOT(3) = FELAS(IB,3) + inTON*FPONP(IB,3)
            SIG(1,1) = SIG(1,1)-(R(IB,1)-RCOM(1))*FTOT(1)
            SIG(1,2) = SIG(1,2)-(R(IB,1)-RCOM(1))*FTOT(2)
            SIG(1,3) = SIG(1,3)-(R(IB,1)-RCOM(1))*FTOT(3)
            SIG(2,1) = SIG(2,1)-(R(IB,2)-RCOM(2))*FTOT(1)
            SIG(2,2) = SIG(2,2)-(R(IB,2)-RCOM(2))*FTOT(2)
            SIG(2,3) = SIG(2,3)-(R(IB,2)-RCOM(2))*FTOT(3)
            SIG(3,1) = SIG(3,1)-(R(IB,3)-RCOM(3))*FTOT(1)
            SIG(3,2) = SIG(3,2)-(R(IB,3)-RCOM(3))*FTOT(2)
            SIG(3,3) = SIG(3,3)-(R(IB,3)-RCOM(3))*FTOT(3)
 30      continue
 10   continue

      RETURN
      END

!---------------------------------------------------------------*
