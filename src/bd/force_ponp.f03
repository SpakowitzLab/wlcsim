!---------------------------------------------------------------*

!
!     This subroutine calculates the self-interaction polymer
!     using the closest-point interpolation strategy
!
!     Andrew Spakowitz
!     Written 11-12-13

      subroutine force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR,SWDT)

      use params, only : dp

      real(dp) R(NT,3)   ! Bead positions
      integer N,NT,NP            ! Current number of beads
      real(dp) FPONP(NT,3) ! Self-interaction force
      real(dp) FMAG     ! Mag of force
      real(dp) RIJ      ! Interbead dist
      real(dp) EIJ(3)   ! Interbead unit vector
      integer I, J              ! Index holders
      integer SKIP              ! Bead skip index

!     Variables for the calculation

      real(dp) U1(3),U2(3),U1U2
      real(dp) D1,D2
      real(dp) R12(3),D12,E12(3)
      real(dp) S1,S2
      real(dp) GI(3)
      integer I1,J1,I2,J2
      integer IB1,IB2

!     Parameters in the simulation

      real(dp) LHC      ! HC length
      real(dp) SIGP     ! HC diameter
      real(dp) VHC     ! Potential strengths
      real(dp) GAM
      real(dp) LBOX     ! Box edge length
      real(dp) SUM
      real(dp) DT
      real(dp) XIR

!     Variables for the timestep switch

      integer SWDT

!     Setup the parameters

      IB1 = 1
      do 10 I1 = 1,NP
         do 20 J1 = 1,N
            FPONP(IB1,1) = 0.
            FPONP(IB1,2) = 0.
            FPONP(IB1,3) = 0.
            IB1 = IB1 + 1
 20      continue
 10   continue

!     Calculate the self-interaction forces

      do 30 I1 = 1,(NP-1)
         do 40 J1 = 1,(N-1)
            IB1 = J1 + N*(I1-1)
            do 50 I2 = (I1 + 1),NP
               do 60 J2 = 1,(N-1)
                  IB2 = J2 + N*(I2-1)
                  R12(1) = R(IB2,1)-R(IB1,1)
                  R12(2) = R(IB2,2)-R(IB1,2)
                  R12(3) = R(IB2,3)-R(IB1,3)

                  D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
                  if (D12 > (3.*GAM)) then
                     goto 70
                  endif

                  U1(1) = R(IB1 + 1,1)-R(IB1,1)
                  U1(2) = R(IB1 + 1,2)-R(IB1,2)
                  U1(3) = R(IB1 + 1,3)-R(IB1,3)
                  D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
                  U1(1) = U1(1)/D1
                  U1(2) = U1(2)/D1
                  U1(3) = U1(3)/D1

                  U2(1) = R(IB2 + 1,1)-R(IB2,1)
                  U2(2) = R(IB2 + 1,2)-R(IB2,2)
                  U2(3) = R(IB2 + 1,3)-R(IB2,3)
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

                  E12(1) = R12(1)/D12
                  E12(2) = R12(2)/D12
                  E12(3) = R12(3)/D12
                  FMAG = VHC*((LHC/D12)**13.-(LHC/D12)**7.)/LHC

                  if ((FMAG/XIR) > (0.001/DT)) then
                     DT = XIR*0.0005/FMAG
                     SWDT = 1
                  endif

                  FPONP(IB1,1) = FPONP(IB1,1) + FMAG*E12(1)*(-1. + S1/D1)
                  FPONP(IB1,2) = FPONP(IB1,2) + FMAG*E12(2)*(-1. + S1/D1)
                  FPONP(IB1,3) = FPONP(IB1,3) + FMAG*E12(3)*(-1. + S1/D1)
                  FPONP(IB1 + 1,1) = FPONP(IB1 + 1,1) + FMAG*E12(1)*(-S1/D1)
                  FPONP(IB1 + 1,2) = FPONP(IB1 + 1,2) + FMAG*E12(2)*(-S1/D1)
                  FPONP(IB1 + 1,3) = FPONP(IB1 + 1,3) + FMAG*E12(3)*(-S1/D1)
                  FPONP(IB2,1) = FPONP(IB2,1) + FMAG*E12(1)*(1.-S2/D2)
                  FPONP(IB2,2) = FPONP(IB2,2) + FMAG*E12(2)*(1.-S2/D2)
                  FPONP(IB2,3) = FPONP(IB2,3) + FMAG*E12(3)*(1.-S2/D2)
                  FPONP(IB2 + 1,1) = FPONP(IB2 + 1,1) + FMAG*E12(1)*(S2/D2)
                  FPONP(IB2 + 1,2) = FPONP(IB2 + 1,2) + FMAG*E12(2)*(S2/D2)
                  FPONP(IB2 + 1,3) = FPONP(IB2 + 1,3) + FMAG*E12(3)*(S2/D2)

 70               continue
 60            continue
 50         continue
 40      continue
 30   continue

      RETURN
      END

!---------------------------------------------------------------*

