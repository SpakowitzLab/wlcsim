!! ---------------------------------------------------------------*

!
!     This subroutine calculates the elastic forces for a wormlike
!     chain with a stretching potential.  The stretch and bend
!     moduli are fed along with the bead positions.
!
!     Andrew Spakowitz
!     Written 9-1-04

      subroutine stressp(COR,R,U,R0,U0,NT,N,NP,PARA,inTON,SIMTYPE)
          use params, only : dp
      implicit none
      real(dp) FELAS(3,NT) ! Elastic force
      real(dp) FPONP(3,NT) ! Self-interaction force
      real(dp) TELAS(3,NT) ! Elastic force
      real(dp) FELAS0(NT,3) ! Elastic force
      real(dp) FPONP0(NT,3) ! Self-interaction force
      integer inTON             ! Include polymer interactions
      real(dp) TELAS0(NT,3) ! Elastic force
      real(dp) R(3,NT)  ! Bead positions
      real(dp) U(3,NT)  ! Tangent vectors
      real(dp) R0(3,NT)  ! Bead positions
      real(dp) U0(3,NT)  ! Tangent vectors
      real(dp) COR

      real(dp) FTOT(3)  ! Compress force
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
      real(dp) SIG0(3,3)

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

      COR = 0.

      DT = 0.0001
      call force_elas(FELAS0,TELAS0,R0,U0,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)
      if (inTON == 1) then
         call force_ponp(FPONP0,R0,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR)
      endif

      call force_elas(FELAS,TELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)
      if (inTON == 1) then
         call force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR)
      endif

      do 10 I = 1,NP
         RCOM(1) = 0.
         RCOM(2) = 0.
         RCOM(3) = 0.
         do 20 J = 1,N
            IB = J + N*(I-1)
            RCOM(1) = RCOM(1) + R(1,IB)/N
            RCOM(2) = RCOM(2) + R(2,IB)/N
            RCOM(3) = RCOM(3) + R(3,IB)/N
 20      continue

         SIG(1,1) = 0.
         SIG(1,2) = 0.
         SIG(1,3) = 0.
         SIG(2,1) = 0.
         SIG(2,2) = 0.
         SIG(2,3) = 0.
         SIG(3,1) = 0.
         SIG(3,2) = 0.
         SIG(3,3) = 0.
         do 30 J = 1,N
            IB = J + N*(I-1)
            FTOT(1) = FELAS(1,IB) + inTON*FPONP(1,IB)
            FTOT(2) = FELAS(2,IB) + inTON*FPONP(2,IB)
            FTOT(3) = FELAS(3,IB) + inTON*FPONP(3,IB)
            SIG(1,1) = SIG(1,1)-(R(1,IB)-RCOM(1))*FTOT(1)
            SIG(1,2) = SIG(1,2)-(R(1,IB)-RCOM(1))*FTOT(2)
            SIG(1,3) = SIG(1,3)-(R(1,IB)-RCOM(1))*FTOT(3)
            SIG(2,1) = SIG(2,1)-(R(2,IB)-RCOM(2))*FTOT(1)
            SIG(2,2) = SIG(2,2)-(R(2,IB)-RCOM(2))*FTOT(2)
            SIG(2,3) = SIG(2,3)-(R(2,IB)-RCOM(2))*FTOT(3)
            SIG(3,1) = SIG(3,1)-(R(3,IB)-RCOM(3))*FTOT(1)
            SIG(3,2) = SIG(3,2)-(R(3,IB)-RCOM(3))*FTOT(2)
            SIG(3,3) = SIG(3,3)-(R(3,IB)-RCOM(3))*FTOT(3)
 30      continue

         RCOM(1) = 0.
         RCOM(2) = 0.
         RCOM(3) = 0.
         do 40 J = 1,N
            IB = J + N*(I-1)
            RCOM(1) = RCOM(1) + R0(1,IB)/N
            RCOM(2) = RCOM(2) + R0(2,IB)/N
            RCOM(3) = RCOM(3) + R0(3,IB)/N
 40      continue

         SIG0(1,1) = 0.
         SIG0(1,2) = 0.
         SIG0(1,3) = 0.
         SIG0(2,1) = 0.
         SIG0(2,2) = 0.
         SIG0(2,3) = 0.
         SIG0(3,1) = 0.
         SIG0(3,2) = 0.
         SIG0(3,3) = 0.
         do 50 J = 1,N
            IB = J + N*(I-1)
            FTOT(1) = FELAS0(IB,1) + inTON*FPONP0(IB,1)
            FTOT(2) = FELAS0(IB,2) + inTON*FPONP0(IB,2)
            FTOT(3) = FELAS0(IB,3) + inTON*FPONP0(IB,3)
            SIG0(1,1) = SIG0(1,1)-(R0(1,IB)-RCOM(1))*FTOT(1)
            SIG0(1,2) = SIG0(1,2)-(R0(1,IB)-RCOM(1))*FTOT(2)
            SIG0(1,3) = SIG0(1,3)-(R0(1,IB)-RCOM(1))*FTOT(3)
            SIG0(2,1) = SIG0(2,1)-(R0(2,IB)-RCOM(2))*FTOT(1)
            SIG0(2,2) = SIG0(2,2)-(R0(2,IB)-RCOM(2))*FTOT(2)
            SIG0(2,3) = SIG0(2,3)-(R0(2,IB)-RCOM(2))*FTOT(3)
            SIG0(3,1) = SIG0(3,1)-(R0(3,IB)-RCOM(3))*FTOT(1)
            SIG0(3,2) = SIG0(3,2)-(R0(3,IB)-RCOM(3))*FTOT(2)
            SIG0(3,3) = SIG0(3,3)-(R0(3,IB)-RCOM(3))*FTOT(3)
 50      continue

         COR = COR + SIG(1,2)*SIG0(1,2) + SIG(1,3)*SIG0(1,3) + SIG(2,1)*SIG0(2,1) &
             + SIG(2,3)*SIG0(2,3) + SIG(3,1)*SIG0(3,1) + SIG(3,2)*SIG0(3,2)

 10   continue

      COR = COR/6.

      RETURN
      END

!---------------------------------------------------------------*
