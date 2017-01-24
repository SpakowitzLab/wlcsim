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

      SUBROUTINE energy_ponp(wlc_p, wlc_d)

      

      DOUBLE PRECISION R(NT,3)   ! Bead positions
      INTEGER N,NT,NP            ! Current number of beads
      DOUBLE PRECISION EPONP ! Self-interaction force
      DOUBLE PRECISION FMAG     ! Mag of force
      DOUBLE PRECISION RIJ      ! Interbead dist
      DOUBLE PRECISION EIJ(3)   ! Interbead unit vector
      INTEGER I, J              ! Index holders
      INTEGER SKIP              ! Bead skip index

!     Variables for the calculation

      DOUBLE PRECISION U1(3),U2(3),U1U2
      DOUBLE PRECISION D1,D2
      DOUBLE PRECISION R12(3),D12,E12(3)
      DOUBLE PRECISION S1,S2
      DOUBLE PRECISION GI(3)
      INTEGER I1,J1,I2,J2
      INTEGER IB1,IB2

!     Parameters in the simulation

      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION LHC      ! HC length
      DOUBLE PRECISION SIGP     ! HC diameter
      DOUBLE PRECISION VHC     ! Potential strengths
      DOUBLE PRECISION GAM
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION SUM
      DOUBLE PRECISION DT
      DOUBLE PRECISION XIR


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


!     Calculate the self-interaction forces

      EPONP=0.
      DO 30 I1=1,(NP-1)
         DO 40 J1=1,(N-1)
            IB1=J1+N*(I1-1)
            DO 50 I2=(I1+1),NP
               DO 60 J2=1,(N-1)
                  IB2=J2+N*(I2-1)
                  R12(1)=R(IB2,1)-R(IB1,1)
                  R12(2)=R(IB2,2)-R(IB1,2)
                  R12(3)=R(IB2,3)-R(IB1,3)

                  D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
                  if (D12.GT.(3.*GAM)) then
                     goto 70
                  endif

                  U1(1)=R(IB1+1,1)-R(IB1,1)
                  U1(2)=R(IB1+1,2)-R(IB1,2)
                  U1(3)=R(IB1+1,3)-R(IB1,3)
                  D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
                  U1(1)=U1(1)/D1
                  U1(2)=U1(2)/D1
                  U1(3)=U1(3)/D1

                  U2(1)=R(IB2+1,1)-R(IB2,1)
                  U2(2)=R(IB2+1,2)-R(IB2,2)
                  U2(3)=R(IB2+1,3)-R(IB2,3)
                  D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
                  U2(1)=U2(1)/D2
                  U2(2)=U2(2)/D2
                  U2(3)=U2(3)/D2

                  U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
                  if (U1U2.EQ.1.) then
                     goto 70
                  endif

                  GI(1)=U1(1)-U1U2*U2(1)
                  GI(2)=U1(2)-U1U2*U2(2)
                  GI(3)=U1(3)-U1U2*U2(3)

                  S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

                  if (S1.GT.D1.OR.S1.LT.0.) then
                     goto 70
                  endif

                  GI(1)=U2(1)-U1U2*U1(1)
                  GI(2)=U2(2)-U1U2*U1(2)
                  GI(3)=U2(3)-U1U2*U1(3)

                  S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

                  if (S2.GT.D2.OR.S2.LT.0.) then
                     goto 70
                  endif

                  R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
                  R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
                  R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

                  D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

                  if (D12.GT.LHC) then
                     goto 70
                  endif

                  FMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

                  EPONP=EPONP+FMAG

 70               CONTINUE
 60            CONTINUE
 50         CONTINUE
 40      CONTINUE
 30   CONTINUE

      RETURN
      END

!---------------------------------------------------------------*

