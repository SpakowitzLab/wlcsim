!---------------------------------------------------------------*

!
! This subroutine calculates the change in the self energy for
! a small Monte Carlo move in the position.
!
! Corrections to force magnitude made 6-3-04.
!
! Andrew Spakowitz
! Written 6-29-04 + +
!
! Should probably be update so that s1/s2 are outside of e.g. [0,d1], we need to check that the endpoints (the beads themselves)
! are not actually colliding.

      subroutine mc_self(DESELF,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
      use params, only : dp
      implicit none
      real(dp) R(3,NT)  ! Bead positions
      real(dp) U(3,NT)  ! Tangent vectors
      real(dp) RP(3,NT)  ! Bead positions
      real(dp) UP(3,NT)  ! Tangent vectors
      integer N,NP,NT           ! Number of beads

!     Variables for the calculation

      real(dp) U1(3),U2(3),U1U2
      real(dp) D1,D2
      real(dp) R12(3),D12,E12(3)
      real(dp) S1,S2
      real(dp) GI(3)
      integer I1,J1,I2,J2,IP,IT1,IT2
      integer IMin,IMAX
      integer IB1,IB2
      integer inD1,inD2

      integer I                ! Current test index
      integer J                ! Index holder
      integer SKIP             ! Bead skip index
      real(dp) DESELF
      real(dp) EMAG

!     Parameters in the simulation

      real(dp) LHC      ! HC length
      real(dp) SIGP     ! HC diameter
      real(dp) VHC     ! Potential strengths
      real(dp) GAM
      real(dp) LBOX     ! Box edge length
      real(dp) SUM
      real(dp) DT
      real(dp) XIR

!     Calculate the self-interaction forces

      DESELF = 0.
      if (IB1 == 1) then
         IMin = 1
      else
         IMin = IB1-1
      endif
      if (IB2 == N) then
         IMAX = (N-1)
      else
         IMAX = IB2
      endif

      do 30 I1 = 1,NP
         if (I1 == IP) then
            goto 100
         endif

         do 40 J1 = 1,(N-1)
            inD1 = J1 + N*(I1-1)
            I2 = IP

            do 50 J2 = IMin,IMAX
               inD2 = J2 + N*(I2-1)
               R12(1) = R(1,inD2)-R(1,inD1)
               R12(2) = R(2,inD2)-R(2,inD1)
               R12(3) = R(3,inD2)-R(3,inD1)
               R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
               R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
               R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

               D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
               if (D12 > (3.*GAM)) then
                  goto 70
               endif

               U1(1) = R(1,inD1 + 1)-R(1,inD1)
               U1(2) = R(2,inD1 + 1)-R(2,inD1)
               U1(3) = R(3,inD1 + 1)-R(3,inD1)
               D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
               U1(1) = U1(1)/D1
               U1(2) = U1(2)/D1
               U1(3) = U1(3)/D1

               U2(1) = R(1,inD2 + 1)-R(1,inD2)
               U2(2) = R(2,inD2 + 1)-R(2,inD2)
               U2(3) = R(3,inD2 + 1)-R(3,inD2)
               D2 = sqrt(U2(1)**2. + U2(2)**2. + U2(3)**2.)
               U2(1) = U2(1)/D2
               U2(2) = U2(2)/D2
               U2(3) = U2(3)/D2

               U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
               if (U1U2 == 1.) then
                  goto 70
               endif

               GI(1) = U1(1)-U1U2*U2(1)
               GI(2) = U1(2)-U1U2*U2(2)
               GI(3) = U1(3)-U1U2*U2(3)

               S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S1 > D1.OR.S1 < 0.) then
                  goto 70
               endif

               GI(1) = U2(1)-U1U2*U1(1)
               GI(2) = U2(2)-U1U2*U1(2)
               GI(3) = U2(3)-U1U2*U1(3)

               S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S2 > D2.OR.S2 < 0.) then
                  goto 70
               endif

               R12(1) = R12(1) + S2*U2(1)-S1*U1(1)
               R12(2) = R12(2) + S2*U2(2)-S1*U1(2)
               R12(3) = R12(3) + S2*U2(3)-S1*U1(3)

               D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)

               if (D12 > LHC) then
                  goto 70
               endif

               EMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6. + 1.)/12.

               DESELF = DESELF-EMAG

 70            continue

 50         continue

            do 80 J2 = IMin,IMAX
               inD2 = J2 + N*(I2-1)
               R12(1) = RP(1,inD2)-R(1,inD1)
               R12(2) = RP(2,inD2)-R(2,inD1)
               R12(3) = RP(3,inD2)-R(3,inD1)
               R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
               R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
               R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

               D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
               if (D12 > (3.*GAM)) then
                  goto 90
               endif

               U1(1) = R(1,inD1 + 1)-R(1,inD1)
               U1(2) = R(2,inD1 + 1)-R(2,inD1)
               U1(3) = R(3,inD1 + 1)-R(3,inD1)
               D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
               U1(1) = U1(1)/D1
               U1(2) = U1(2)/D1
               U1(3) = U1(3)/D1

               U2(1) = RP(1,inD2 + 1)-RP(1,inD2)
               U2(2) = RP(2,inD2 + 1)-RP(2,inD2)
               U2(3) = RP(3,inD2 + 1)-RP(3,inD2)
               D2 = sqrt(U2(1)**2. + U2(2)**2. + U2(3)**2.)
               U2(1) = U2(1)/D2
               U2(2) = U2(2)/D2
               U2(3) = U2(3)/D2

               U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
               if (U1U2 == 1.) then
                  goto 90
               endif

               GI(1) = U1(1)-U1U2*U2(1)
               GI(2) = U1(2)-U1U2*U2(2)
               GI(3) = U1(3)-U1U2*U2(3)

               S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S1 > D1.OR.S1 < 0.) then
                  goto 90
               endif

               GI(1) = U2(1)-U1U2*U1(1)
               GI(2) = U2(2)-U1U2*U1(2)
               GI(3) = U2(3)-U1U2*U1(3)

               S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

               if (S2 > D2.OR.S2 < 0.) then
                  goto 90
               endif

               R12(1) = R12(1) + S2*U2(1)-S1*U1(1)
               R12(2) = R12(2) + S2*U2(2)-S1*U1(2)
               R12(3) = R12(3) + S2*U2(3)-S1*U1(3)

               D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)

               if (D12 > LHC) then
                  goto 90
               endif

               EMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6. + 1.)/12.

               DESELF = DESELF + EMAG

 90            continue
 80         continue


 40      continue
 100     continue
 30   continue

      RETURN
      END

!---------------------------------------------------------------*

