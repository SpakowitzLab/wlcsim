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

subroutine energy_self_slide(EPONP, R, NT, N, NP, PARA, RinG, IB1, IB2)
   use params, only: dp

   implicit none
   integer N, NT, NP            ! Current number of beads
   real(dp) R(3, NT)   ! Bead positions
   real(dp) EPONP ! Self-interaction force
   real(dp) FMAG     ! Mag of force
   real(dp) RIJ      ! Interbead dist
   real(dp) EIJ(3)   ! Interbead unit vector
   integer I, J              ! Index holders
   integer SKIP              ! Bead skip index

   !     Variables for the calculation

   real(dp) U1(3), U2(3), U1U2
   real(dp) D1, D2
   real(dp) R12(3), D12, E12(3), R12T(3), R12C1(3), R12C2(3)
   real(dp) S1, S2
   real(dp) GI(3)
   integer I1, J1, I2, J2
   integer IB1, IB2, IB1P1, IB2P1
   integer IO, II, IS1, IS2, IOP1, IIP1, IS1P1, IS2P1
   integer DIO, DII, DIB

   !     Parameters in the simulation

   real(dp) PARA(10)
   real(dp) LHC      ! HC length
   real(dp) SIGP     ! HC diameter
   real(dp) VHC         ! Potential strengths
   real(dp) GAM
   real(dp) LBOX     ! Box edge length
   real(dp) SUM
   real(dp) DT
   real(dp) XIR
   real(dp) XIU
   real(dp) ETA
   real(dp) EPAR
   real(dp) EPERP
   real(dp) EB
   logical RinG              ! Is polymer a ring?
   integer NMAX

   real(dp) D12Min, FMAGMin

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

   EPONP = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Setup. Determine the number of segments that lie
   !betwen the two terminal points of the segment slid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (IB2 >= IB1) then
      DIB = IB2 - IB1
   else
      DIB = (N - IB1) + IB2
   ENDif
   DIO = N - DIB - 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !The two segments immediately on the outside of the
   !beads that were slid are "stretched" or compressed.
   !The interaction of these segments with the remainder
   !of the outer segments is changed, and their interaction
   !with the segments inside the region slid is changed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IS1P1 = IB1                               ! Index of first stretched (or compressed) segment (not inside beads moved)
   IS2 = IB2                                  ! Index of second stretched (or compressed) segment (not inside beads moved)

   if (IS1P1 == 1) then
      IS1 = N
   else
      IS1 = IS1P1 - 1
   ENDif

   if (IS2 == N) then
      IS2P1 = 1
   else
      IS2P1 = IS2 + 1
   ENDif

   II = IB1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Calculate the interaction between inner segments with all outer segments
   !The following loop is only performed if the two terminal beads are not equal
   !Note that special cases for end positions are currently only set-up for rings
   !This code can accomodate both linear and ring chains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (IB1 /= IB2) then
      !Sum over segments inside segment slid (inner segments)
      do I = 1, DIB

         if (II == N .AND. RinG) then
            IIP1 = 1
         elseif (II == N + 1) then
            II = 1
            IIP1 = II + 1
         else
            IIP1 = II + 1
         ENDif
         IO = IB2 + 1
         !Sum over segments unchanged by slide (outer segments)
         do J = 1, DIO

            if (IO == N .AND. RinG) then
               IOP1 = 1
            elseif (IO == N .AND. (.not. RinG)) then
               IO = 0
               GOTO 110
            elseif (IO == N + 1) then
               IO = 1
               IOP1 = IO + 1
            else
               IOP1 = IO + 1
            ENDif

            !Skip this pair if the segments are adjacent
            if (IIP1 == IO .OR. IOP1 == II) then
               GOTO 70
            ENDif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Calculate the interaction energy between the inner segment
            !and outer segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            R12(1) = R(1, IO) - R(1, II)
            R12(2) = R(2, IO) - R(2, II)
            R12(3) = R(3, IO) - R(3, II)
            !Periodic Bounary conditions. Not used for single chains

            ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
            ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
            ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

            D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

            U1(1) = R(1, IIP1) - R(1, II)
            U1(2) = R(2, IIP1) - R(2, II)
            U1(3) = R(3, IIP1) - R(3, II)
            D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
            U1(1) = U1(1)/D1
            U1(2) = U1(2)/D1
            U1(3) = U1(3)/D1

            U2(1) = R(1, IOP1) - R(1, IO)
            U2(2) = R(2, IOP1) - R(2, IO)
            U2(3) = R(3, IOP1) - R(3, IO)
            D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
            U2(1) = U2(1)/D2
            U2(2) = U2(2)/D2
            U2(3) = U2(3)/D2

            U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
            if (U1U2 == 1. .OR. U1U2 == -1.) then
               D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
               GOTO 60
            endif

            GI(1) = U1(1) - U1U2*U2(1)
            GI(2) = U1(2) - U1U2*U2(2)
            GI(3) = U1(3) - U1U2*U2(3)

            S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

            if (S1 > D1 .OR. S1 < 0.) then
               R12T = R(:, IOP1) - R(:, IIP1)
               R12C1 = R(:, IOP1) - R(:, II)
               R12C2 = R(:, IIP1) - R(:, IO)
               D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                    & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
               GOTO 60
            endif

            GI(1) = U2(1) - U1U2*U1(1)
            GI(2) = U2(2) - U1U2*U1(2)
            GI(3) = U2(3) - U1U2*U1(3)

            S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

            if (S2 > D2 .OR. S2 < 0.) then
               R12T = R(:, IOP1) - R(:, IIP1)
               R12C1 = R(:, IOP1) - R(:, II)
               R12C2 = R(:, IIP1) - R(:, IO)
               D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                    & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
               GOTO 60
            endif

            R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
            R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
            R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

            D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

60          if (D12 > LHC) then
               goto 70
            endif

            FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

            EPONP = EPONP + FMAG

70          continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Calculate the interaction between this outer segment and the
            !two stretched segments. Only do this once.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (I == 1) then

               !Skip this pair if the segments are adjacent
               if (IS1 == IOP1 .OR. IO == IS1P1) then
                  GOTO 90
               ENDif

               !If the chain is linear and IS1P1 == 1 then there is no first 'stretched' segment
               if (IS1P1 == 1 .AND. (.not. RinG)) then
                  GOTO 90
               ENDif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Calculate interaction between outer segment and first stretched
               !segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               R12(1) = R(1, IS1) - R(1, IO)
               R12(2) = R(2, IS1) - R(2, IO)
               R12(3) = R(3, IS1) - R(3, IO)
               !Periodic Bounary conditions. Not used for single chains

               ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
               ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
               ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

               D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

               U1(1) = R(1, IOP1) - R(1, IO)
               U1(2) = R(2, IOP1) - R(2, IO)
               U1(3) = R(3, IOP1) - R(3, IO)
               D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
               U1(1) = U1(1)/D1
               U1(2) = U1(2)/D1
               U1(3) = U1(3)/D1

               U2(1) = R(1, IS1P1) - R(1, IS1)
               U2(2) = R(2, IS1P1) - R(2, IS1)
               U2(3) = R(3, IS1P1) - R(3, IS1)
               D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
               U2(1) = U2(1)/D2
               U2(2) = U2(2)/D2
               U2(3) = U2(3)/D2

               U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
               if (U1U2 == 1. .OR. U1U2 == -1.) then
                  D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
                  GOTO 80
               endif

               GI(1) = U1(1) - U1U2*U2(1)
               GI(2) = U1(2) - U1U2*U2(2)
               GI(3) = U1(3) - U1U2*U2(3)

               S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S1 > D1 .OR. S1 < 0.) then

                  R12T = R(:, IOP1) - R(:, IS1P1)
                  R12C1 = R(:, IOP1) - R(:, IS1)
                  R12C2 = R(:, IS1P1) - R(:, IO)

                  D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                       & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
                  GOTO 80
               endif

               GI(1) = U2(1) - U1U2*U1(1)
               GI(2) = U2(2) - U1U2*U1(2)
               GI(3) = U2(3) - U1U2*U1(3)

               S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S2 > D2 .OR. S2 < 0.) then

                  R12T = R(:, IOP1) - R(:, IS1P1)
                  R12C1 = R(:, IOP1) - R(:, IS1)
                  R12C2 = R(:, IS1P1) - R(:, IO)

                  D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                       & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
                  GOTO 80
               endif

               R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
               R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
               R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

               D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

80             if (D12 > LHC) then
                  goto 90
               endif

               FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

               EPONP = EPONP + FMAG
90             continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Calculate the interaction between the outer segment and the second
               !stretched segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               !Skip this pair if the segments are adjacent
               if (IS2P1 == IO .OR. IOP1 == IS2) then
                  GOTO 110
               ENDif

               !If the chain is linear and IS2 = N, then there is no second 'stretched' segment
               if (IS2 == N .AND. (.not. RinG)) then
                  GOTO 110
               ENDif

               R12(1) = R(1, IS2) - R(1, IO)
               R12(2) = R(2, IS2) - R(2, IO)
               R12(3) = R(3, IS2) - R(3, IO)
               !Periodic Bounary conditions. Not used for single chains

               ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
               ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
               ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

               D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

               U1(1) = R(1, IOP1) - R(1, IO)
               U1(2) = R(2, IOP1) - R(2, IO)
               U1(3) = R(3, IOP1) - R(3, IO)
               D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
               U1(1) = U1(1)/D1
               U1(2) = U1(2)/D1
               U1(3) = U1(3)/D1

               U2(1) = R(1, IS2P1) - R(1, IS2)
               U2(2) = R(2, IS2P1) - R(2, IS2)
               U2(3) = R(3, IS2P1) - R(3, IS2)
               D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
               U2(1) = U2(1)/D2
               U2(2) = U2(2)/D2
               U2(3) = U2(3)/D2

               U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
               if (U1U2 == 1. .OR. U1U2 == -1.) then
                  D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
                  GOTO 100
               endif

               GI(1) = U1(1) - U1U2*U2(1)
               GI(2) = U1(2) - U1U2*U2(2)
               GI(3) = U1(3) - U1U2*U2(3)

               S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S1 > D1 .OR. S1 < 0.) then

                  R12T = R(:, IOP1) - R(:, IS2P1)
                  R12C1 = R(:, IOP1) - R(:, IS2)
                  R12C2 = R(:, IS2P1) - R(:, IO)

                  D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                       & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
                  GOTO 100
               endif

               GI(1) = U2(1) - U1U2*U1(1)
               GI(2) = U2(2) - U1U2*U1(2)
               GI(3) = U2(3) - U1U2*U1(3)

               S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S2 > D2 .OR. S2 < 0.) then
                  R12T = R(:, IOP1) - R(:, IS2P1)
                  R12C1 = R(:, IOP1) - R(:, IS2)
                  R12C2 = R(:, IS2P1) - R(:, IO)
                  D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                       & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
                  GOTO 100
               endif

               R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
               R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
               R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

               D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

100            if (D12 > LHC) then
                  goto 110
               endif

               FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

               EPONP = EPONP + FMAG

            ENDif
            !step to next outer segment
110         continue
            IO = IO + 1
         ENDdo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Calculate the interaction between the inner segment and the first stretched
         !segment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !Skip this pair if the segments are adjacent
         if (IS1 == IIP1 .OR. II == IS1P1) then
            GOTO 130
         ENDif

         !If the chain is linear and  IS1P1 = 1 there is no first "stretched" segment
         if (IS1P1 == 1 .AND. (.not. RinG)) then
            GOTO 130
         ENDif

         R12(1) = R(1, IS1) - R(1, II)
         R12(2) = R(2, IS1) - R(2, II)
         R12(3) = R(3, IS1) - R(3, II)
         !Periodic Bounary conditions. Not used for single chains

         ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
         ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
         ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

         U1(1) = R(1, IIP1) - R(1, II)
         U1(2) = R(2, IIP1) - R(2, II)
         U1(3) = R(3, IIP1) - R(3, II)
         D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
         U1(1) = U1(1)/D1
         U1(2) = U1(2)/D1
         U1(3) = U1(3)/D1

         U2(1) = R(1, IS1P1) - R(1, IS1)
         U2(2) = R(2, IS1P1) - R(2, IS1)
         U2(3) = R(3, IS1P1) - R(3, IS1)
         D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
         U2(1) = U2(1)/D2
         U2(2) = U2(2)/D2
         U2(3) = U2(3)/D2

         U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
         if (U1U2 == 1. .OR. U1U2 == -1.) then
            D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
            GOTO 120
         endif

         GI(1) = U1(1) - U1U2*U2(1)
         GI(2) = U1(2) - U1U2*U2(2)
         GI(3) = U1(3) - U1U2*U2(3)

         S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S1 > D1 .OR. S1 < 0.) then
            R12T = R(:, IS1P1) - R(:, IIP1)
            R12C1 = R(:, IS1P1) - R(:, II)
            R12C2 = R(:, IIP1) - R(:, IS1)
            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 120
         endif

         GI(1) = U2(1) - U1U2*U1(1)
         GI(2) = U2(2) - U1U2*U1(2)
         GI(3) = U2(3) - U1U2*U1(3)

         S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S2 > D2 .OR. S2 < 0.) then

            R12T = R(:, IS1P1) - R(:, IIP1)
            R12C1 = R(:, IS1P1) - R(:, II)
            R12C2 = R(:, IIP1) - R(:, IS1)

            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 120
         endif

         R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
         R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
         R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

120      if (D12 > LHC) then
            goto 130
         endif

         FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

         EPONP = EPONP + FMAG
130      continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Calculate interaction between inner segment and second stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (IS2 == IIP1 .OR. II == IS2P1) then
            GOTO 150
         ENDif

         !If the chain is linear and IS2 = N, there is no second "stretched" segment
         if (IS2 == N .AND. (.not. RinG)) then
            GOTO 150
         ENDif

         R12(1) = R(1, IS2) - R(1, II)
         R12(2) = R(2, IS2) - R(2, II)
         R12(3) = R(3, IS2) - R(3, II)
         !Periodic Bounary conditions. Not used for single chains

         ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
         ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
         ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

         U1(1) = R(1, IIP1) - R(1, II)
         U1(2) = R(2, IIP1) - R(2, II)
         U1(3) = R(3, IIP1) - R(3, II)
         D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
         U1(1) = U1(1)/D1
         U1(2) = U1(2)/D1
         U1(3) = U1(3)/D1

         U2(1) = R(1, IS2P1) - R(1, IS2)
         U2(2) = R(2, IS2P1) - R(2, IS2)
         U2(3) = R(3, IS2P1) - R(3, IS2)
         D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
         U2(1) = U2(1)/D2
         U2(2) = U2(2)/D2
         U2(3) = U2(3)/D2

         U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
         if (U1U2 == 1. .OR. U1U2 == -1.) then
            D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
            GOTO 140
         endif

         GI(1) = U1(1) - U1U2*U2(1)
         GI(2) = U1(2) - U1U2*U2(2)
         GI(3) = U1(3) - U1U2*U2(3)

         S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S1 > D1 .OR. S1 < 0.) then
            R12T = R(:, IS2P1) - R(:, IIP1)
            R12C1 = R(:, IS2P1) - R(:, II)
            R12C2 = R(:, IIP1) - R(:, IS2)
            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 140
         endif

         GI(1) = U2(1) - U1U2*U1(1)
         GI(2) = U2(2) - U1U2*U1(2)
         GI(3) = U2(3) - U1U2*U1(3)

         S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S2 > D2 .OR. S2 < 0.) then
            R12T = R(:, IS2P1) - R(:, IIP1)
            R12C1 = R(:, IS2P1) - R(:, II)
            R12C2 = R(:, IIP1) - R(:, IS2)

            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 140
         endif

         R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
         R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
         R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

140      if (D12 > LHC) then
            goto 150
         endif

         FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

         EPONP = EPONP + FMAG
150      continue

         II = II + 1
      ENDdo
   ENDif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Calculate the interaction between the two stretched segments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   R12(1) = R(1, IS2) - R(1, IS1)
   R12(2) = R(2, IS2) - R(2, IS1)
   R12(3) = R(3, IS2) - R(3, IS1)

   if (IS1 == IS2P1 .OR. IS1P1 == IS2) then
      GOTO 170
   ENDif

   !If the chain is linear and IS1P1 = 1 or IS2 = N, then there is no
   !first segment or second stretched segment

   if (IS1P1 == 1 .OR. IS2 == N) then
      GOTO 170
   ENDif

   !Periodic Bounary conditions. Not used for single chains

   ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
   ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
   ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

   D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

   U1(1) = R(1, IS1P1) - R(1, IS1)
   U1(2) = R(2, IS1P1) - R(2, IS1)
   U1(3) = R(3, IS1P1) - R(3, IS1)
   D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
   U1(1) = U1(1)/D1
   U1(2) = U1(2)/D1
   U1(3) = U1(3)/D1

   U2(1) = R(1, IS2P1) - R(1, IS2)
   U2(2) = R(2, IS2P1) - R(2, IS2)
   U2(3) = R(3, IS2P1) - R(3, IS2)
   D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
   U2(1) = U2(1)/D2
   U2(2) = U2(2)/D2
   U2(3) = U2(3)/D2

   U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
   if (U1U2 == 1. .OR. U1U2 == -1.) then
      D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
      GOTO 160
   endif

   GI(1) = U1(1) - U1U2*U2(1)
   GI(2) = U1(2) - U1U2*U2(2)
   GI(3) = U1(3) - U1U2*U2(3)

   S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

   if (S1 > D1 .OR. S1 < 0.) then
      R12T = R(:, IS2P1) - R(:, IS1P1)
      R12C1 = R(:, IS2P1) - R(:, IS1)
      R12C2 = R(:, IS1P1) - R(:, IS2)
      D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
           & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
      GOTO 160
   endif

   GI(1) = U2(1) - U1U2*U1(1)
   GI(2) = U2(2) - U1U2*U1(2)
   GI(3) = U2(3) - U1U2*U1(3)

   S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

   if (S2 > D2 .OR. S2 < 0.) then
      R12T = R(:, IS2P1) - R(:, IS1P1)
      R12C1 = R(:, IS2P1) - R(:, IS1)
      R12C2 = R(:, IS1P1) - R(:, IS2)
      D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
           & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
      GOTO 160
   endif

   R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
   R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
   R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

   D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

160 if (D12 > LHC) then
      goto 170
   endif

   FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

   EPONP = EPONP + FMAG

170 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !If only one bead was slid (IB1 == IB2) then calculate the interaction between the
   !portion of the chain unaffected by the move and the two segments stretched/compressed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (IB1 == IB2) then
      !If IB1 == IB2, Get interaction between segment unaffected by move and two segments affected

      IO = IB1 + 1
      !Sum over segments unchanged by slide (outer segments)
      do J = 1, DIO

         if (IO == N .AND. RinG) then
            IOP1 = 1
         elseif (IO == N .AND. (.not. RinG)) then
            IO = 0
            GOTO 210
         elseif (IO == N + 1) then
            IO = 1
            IOP1 = IO + 1
         else
            IOP1 = IO + 1
         ENDif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Calculate interaction between outer segment and first stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (IS1 == IOP1 .OR. IO == IS1P1) then
            GOTO 190
         ENDif

         !If the chain is linear and IP1P1 = 1, there is no first "stretched" segment
         if (IS1P1 == 1 .AND. (.not. RinG)) then
            GOTO 190
         ENDif

         R12(1) = R(1, IS1) - R(1, IO)
         R12(2) = R(2, IS1) - R(2, IO)
         R12(3) = R(3, IS1) - R(3, IO)
         !Periodic Bounary conditions. Not used for single chains

         ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
         ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
         ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

         U1(1) = R(1, IOP1) - R(1, IO)
         U1(2) = R(2, IOP1) - R(2, IO)
         U1(3) = R(3, IOP1) - R(3, IO)
         D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
         U1(1) = U1(1)/D1
         U1(2) = U1(2)/D1
         U1(3) = U1(3)/D1

         U2(1) = R(1, IS1P1) - R(1, IS1)
         U2(2) = R(2, IS1P1) - R(2, IS1)
         U2(3) = R(3, IS1P1) - R(3, IS1)
         D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
         U2(1) = U2(1)/D2
         U2(2) = U2(2)/D2
         U2(3) = U2(3)/D2

         U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
         if (U1U2 == 1. .OR. U1U2 == -1.) then
            D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
            GOTO 180
         endif

         GI(1) = U1(1) - U1U2*U2(1)
         GI(2) = U1(2) - U1U2*U2(2)
         GI(3) = U1(3) - U1U2*U2(3)

         S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S1 > D1 .OR. S1 < 0.) then

            R12T = R(:, IOP1) - R(:, IS1P1)
            R12C1 = R(:, IOP1) - R(:, IS1)
            R12C2 = R(:, IS1P1) - R(:, IO)
            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 180
         endif

         GI(1) = U2(1) - U1U2*U1(1)
         GI(2) = U2(2) - U1U2*U1(2)
         GI(3) = U2(3) - U1U2*U1(3)

         S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S2 > D2 .OR. S2 < 0.) then

            R12T = R(:, IOP1) - R(:, IS1P1)
            R12C1 = R(:, IOP1) - R(:, IS1)
            R12C2 = R(:, IS1P1) - R(:, IO)

            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 180
         endif

         R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
         R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
         R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

180      if (D12 > LHC) then
            goto 190
         endif

         FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

         EPONP = EPONP + FMAG

190      continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Calculate the interaction between outer segment and second stretched segment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (IS2 == IOP1 .OR. IO == IS2P1) then
            GOTO 210
         ENDif

         !If the chain is linear and IS2 = N, then there is no second "stretched" segment

         if (IS2 == N .AND. (.not. RinG)) then
            GOTO 210
         ENDif

         R12(1) = R(1, IS2) - R(1, IO)
         R12(2) = R(2, IS2) - R(2, IO)
         R12(3) = R(3, IS2) - R(3, IO)

         !Periodic Bounary conditions. Not used for single chains

         ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
         ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
         ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

         U1(1) = R(1, IOP1) - R(1, IO)
         U1(2) = R(2, IOP1) - R(2, IO)
         U1(3) = R(3, IOP1) - R(3, IO)
         D1 = sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
         U1(1) = U1(1)/D1
         U1(2) = U1(2)/D1
         U1(3) = U1(3)/D1

         U2(1) = R(1, IS2P1) - R(1, IS2)
         U2(2) = R(2, IS2P1) - R(2, IS2)
         U2(3) = R(3, IS2P1) - R(3, IS2)
         D2 = sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
         U2(1) = U2(1)/D2
         U2(2) = U2(2)/D2
         U2(3) = U2(3)/D2

         U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
         if (U1U2 == 1. .OR. U1U2 == -1.) then
            D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
            GOTO 200
         endif

         GI(1) = U1(1) - U1U2*U2(1)
         GI(2) = U1(2) - U1U2*U2(2)
         GI(3) = U1(3) - U1U2*U2(3)

         S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S1 > D1 .OR. S1 < 0.) then
            R12T = R(:, IOP1) - R(:, IS2P1)
            R12C1 = R(:, IOP1) - R(:, IS2)
            R12C2 = R(:, IS2P1) - R(:, IO)

            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 200
         endif

         GI(1) = U2(1) - U1U2*U1(1)
         GI(2) = U2(2) - U1U2*U1(2)
         GI(3) = U2(3) - U1U2*U1(3)

         S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

         if (S2 > D2 .OR. S2 < 0.) then
            R12T = R(:, IOP1) - R(:, IS2P1)
            R12C1 = R(:, IOP1) - R(:, IS2)
            R12C2 = R(:, IS2P1) - R(:, IO)

            D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2, R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                 & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2, R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
            GOTO 200
         endif

         R12(1) = R12(1) + S2*U2(1) - S1*U1(1)
         R12(2) = R12(2) + S2*U2(2) - S1*U1(2)
         R12(3) = R12(3) + S2*U2(3) - S1*U1(3)

         D12 = sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

200      if (D12 > LHC) then
            goto 210
         endif

         FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

         EPONP = EPONP + FMAG
210      continue

         IO = IO + 1
      ENDdo
   ENDif

   RETURN
ENDsubroutine energy_self_slide

!---------------------------------------------------------------*

