!---------------------------------------------------------------*

!
!     This subroutine calculates the self-interaction of a DNA-like
!     molecule as a linear chain with charges along the centerline.
!     Within the program is specified the Bjerrum length, the
!     Debye length, and the length of hard-core repulsion.  Exactly
!     as given in AJSclamp4-16-04 and elsewhere.
!
!     Corrections to force magnitude made 6-3-04.
!
!     Andrew Spakowitz
!     Written 1-31-05

      subroutine energy_ponp(wlc_p, wlc_d)

      

      real(dp) R(NT,3)   ! Bead positions
      integer N,NT,NP            ! Current number of beads
      real(dp) EPONP ! Self-interaction force
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

      real(dp) PARA(10)
      real(dp) LHC      ! HC length
      real(dp) SIGP     ! HC diameter
      real(dp) VHC     ! Potential strengths
      real(dp) GAM
      real(dp) LBOX     ! Box edge length
      real(dp) SUM
      real(dp) DT
      real(dp) XIR


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


!     Calculate the self-interaction forces

      EPONP = 0.
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

                  FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6. + 1.)/12.

                  EPONP = EPONP + FMAG

 70               continue
 60            continue
 50         continue
 40      continue
 30   continue

      RETURN
      END

!---------------------------------------------------------------*

