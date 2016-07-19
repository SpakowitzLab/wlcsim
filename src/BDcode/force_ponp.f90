!---------------------------------------------------------------*

!
!     This subroutine calculates the self-interaction polymer
!     using the closest-point interpolation strategy
!
!     Andrew Spakowitz
!     Written 11-12-13

      SUBROUTINE force_ponp(FPONP,R,NT,N,NP,LHC,VHC,LBOX,GAM,DT,XIR,SWDT)

      DOUBLE PRECISION R(NT,3)   ! Bead positions
      INTEGER N,NT,NP            ! Current number of beads
      DOUBLE PRECISION FPONP(NT,3) ! Self-interaction force
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

      DOUBLE PRECISION LHC      ! HC length
      DOUBLE PRECISION SIGP     ! HC diameter
      DOUBLE PRECISION VHC     ! Potential strengths
      DOUBLE PRECISION GAM
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION SUM
      DOUBLE PRECISION DT
      DOUBLE PRECISION XIR

!     Variables for the timestep switch

      INTEGER SWDT

!     Setup the parameters

      IB1=1
      DO 10 I1=1,NP
         DO 20 J1=1,N
            FPONP(IB1,1)=0.
            FPONP(IB1,2)=0.
            FPONP(IB1,3)=0.
            IB1=IB1+1
 20      CONTINUE
 10   CONTINUE

!     Calculate the self-interaction forces

      DO 30 I1=1,(NP-1)
         DO 40 J1=1,(N-1)
            IB1=J1+N*(I1-1)
            DO 50 I2=(I1+1),NP
               DO 60 J2=1,(N-1)
                  IB2=J2+N*(I2-1)
                  R12(1)=R(IB2,1)-R(IB1,1)
                  R12(2)=R(IB2,2)-R(IB1,2)
                  R12(3)=R(IB2,3)-R(IB1,3)
                  R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
                  R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
                  R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

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

                  E12(1)=R12(1)/D12
                  E12(2)=R12(2)/D12
                  E12(3)=R12(3)/D12
                  FMAG=VHC*((LHC/D12)**13.-(LHC/D12)**7.)/LHC

                  if ((FMAG/XIR).GT.(0.001/DT)) then
                     DT=XIR*0.0005/FMAG
                     SWDT=1
                  endif

                  FPONP(IB1,1)=FPONP(IB1,1)+FMAG*E12(1)*(-1.+S1/D1)
                  FPONP(IB1,2)=FPONP(IB1,2)+FMAG*E12(2)*(-1.+S1/D1)
                  FPONP(IB1,3)=FPONP(IB1,3)+FMAG*E12(3)*(-1.+S1/D1)
                  FPONP(IB1+1,1)=FPONP(IB1+1,1)+FMAG*E12(1)*(-S1/D1)
                  FPONP(IB1+1,2)=FPONP(IB1+1,2)+FMAG*E12(2)*(-S1/D1)
                  FPONP(IB1+1,3)=FPONP(IB1+1,3)+FMAG*E12(3)*(-S1/D1)
                  FPONP(IB2,1)=FPONP(IB2,1)+FMAG*E12(1)*(1.-S2/D2)
                  FPONP(IB2,2)=FPONP(IB2,2)+FMAG*E12(2)*(1.-S2/D2)
                  FPONP(IB2,3)=FPONP(IB2,3)+FMAG*E12(3)*(1.-S2/D2)
                  FPONP(IB2+1,1)=FPONP(IB2+1,1)+FMAG*E12(1)*(S2/D2)
                  FPONP(IB2+1,2)=FPONP(IB2+1,2)+FMAG*E12(2)*(S2/D2)
                  FPONP(IB2+1,3)=FPONP(IB2+1,3)+FMAG*E12(3)*(S2/D2)

 70               CONTINUE
 60            CONTINUE
 50         CONTINUE
 40      CONTINUE
 30   CONTINUE

      RETURN
      END

!---------------------------------------------------------------*

