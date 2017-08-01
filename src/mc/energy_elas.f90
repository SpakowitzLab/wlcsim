!---------------------------------------------------------------*

!
!     This subroutine calculates the elastic energies for a wormlike
!     chain with a stretching potential. The stretch and bend
!     moduli are fed along with the bead positions.


      subroutine energy_elas(EELAS,R,U,NT,NB,NP,PARA,RinG,TWIST,Lk,lt,L)
      use params, only: dp, pi
      implicit none
      integer, intent(in) :: NB           ! Number of beads in a polymer
      integer, intent(in) :: NT           ! Number of beads total
      integer, intent(in) :: NP           ! Number of polymers
      integer, intent(in) :: lk
      real(dp), intent(in) :: l, lt
      real(dp), intent(in) :: R(3,NT)  ! Bead positions
      real(dp), intent(in) :: U(3,NT)  ! Tangent vectors
      real(dp), intent(out):: EELAS(6) ! Elastic force
      integer WR,TW ! writhe, twist
      integer I,J,IB,ibp1            ! Index holders
      LOGICAL, intent(in) :: RinG, TWIST

!     Polymer properties

      real(dp), intent(in) :: PARA(10)
      real(dp) EB,EPAR,EPERP
      real(dp) GAM,ETA

      real(dp) DR(3),DRPAR,DRPERP(3)
      real(dp) GI(3) !,U1U2


      EB = PARA(1)
      EPAR = PARA(2)
      EPERP = PARA(3)
      GAM = PARA(4)
      ETA = PARA(5)

      EELAS(1) = 0.0_dp
      EELAS(2) = 0.0_dp
      EELAS(3) = 0.0_dp
      IB = 1
      do I = 1,NP
         do J = 1,NB
            if (RinG) then
                if (J == NB) then
                    IBP1 = 1 + (I-1)*NB
                else
                    IBP1 = IB + 1
                ENDif
            elseif (J == NB) then
                CYCLE
            else
                IBP1 = IB + 1
            ENDif
            DR(1) = R(1,IBP1)-R(1,IB)
            DR(2) = R(2,IBP1)-R(2,IB)
            DR(3) = R(3,IBP1)-R(3,IB)
            DRPAR = DR(1)*U(1,IB) + DR(2)*U(2,IB) + DR(3)*U(3,IB)

            DRPERP(1) = DR(1)-DRPAR*U(1,IB)
            DRPERP(2) = DR(2)-DRPAR*U(2,IB)
            DRPERP(3) = DR(3)-DRPAR*U(3,IB)
            ! energy does not depend on u\cdot{}u
            !U1U2 = U(1,IB)*U(1,IBP1) + U(2,IB)*U(2,IBP1) + U(3,IB)*U(3,IBP1)

            GI(1) = (U(1,IBP1)-U(1,IB)-ETA*DRPERP(1))
            GI(2) = (U(2,IBP1)-U(2,IB)-ETA*DRPERP(2))
            GI(3) = (U(3,IBP1)-U(3,IB)-ETA*DRPERP(3))

            EELAS(1) = EELAS(1) + 0.5*EB*(GI(1)**2. + GI(2)**2. + GI(3)**2.)
            EELAS(2) = EELAS(2) + 0.5*EPAR*(DRPAR-GAM)**2.
            EELAS(3) = EELAS(3) + 0.5*EPERP*(DRPERP(1)**2. + DRPERP(2)**2. + DRPERP(3)**2.)

            IB = IB + 1
         ENDdo
         IB = IB + 1
      ENDdo

      ! Get Twist Energy
      if (TWIST) then
          call WRITHE(R,NB,Wr)
          Tw = Lk-Wr
          EELAS(4) = ((2*PI*Tw)**2)*LT/(2*L)
      ENDif

      RETURN
      END

!---------------------------------------------------------------*
