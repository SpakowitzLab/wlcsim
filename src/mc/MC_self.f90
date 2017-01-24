!---------------------------------------------------------------*

!
! This subroutine calculates the change in the self energy for
! a small Monte Carlo move in the position.
!
! Corrections to force magnitude made 6-3-04.
!
! Andrew Spakowitz
! Written 6-29-04++
!
! Should probably be update so that s1/s2 are outside of e.g. [0,d1], we need to check that the endpoints (the beads themselves)
! are not actually colliding.

      SUBROUTINE MC_self(DESELF,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)

      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      INTEGER N,NP,NT           ! Number of beads

!     Variables for the calculation

      DOUBLE PRECISION U1(3),U2(3),U1U2
      DOUBLE PRECISION D1,D2
      DOUBLE PRECISION R12(3),D12,E12(3)
      DOUBLE PRECISION S1,S2
      DOUBLE PRECISION GI(3)
      INTEGER I1,J1,I2,J2
      INTEGER IMIN,IMAX
      INTEGER IB1,IB2
      INTEGER IND1,IND2

      INTEGER I                ! Current test index
      INTEGER J                ! Index holder
      INTEGER SKIP             ! Bead skip index
      DOUBLE PRECISION DESELF
      DOUBLE PRECISION EMAG

!     Parameters in the simulation

      DOUBLE PRECISION LHC      ! HC length
      DOUBLE PRECISION SIGP     ! HC diameter
      DOUBLE PRECISION VHC     ! Potential strengths
      DOUBLE PRECISION GAM
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION SUM
      DOUBLE PRECISION DT
      DOUBLE PRECISION XIR

!     Calculate the self-interaction forces

      DESELF=0.
      if (IB1.EQ.1) then
         IMIN=1
      else
         IMIN=IB1-1
      endif
      if (IB2.EQ.N) then
         IMAX=(N-1)
      else
         IMAX=IB2
      endif

      DO 30 I1=1,NP
         if (I1.EQ.IP) then
            goto 100
         endif

         DO 40 J1=1,(N-1)
            IND1=J1+N*(I1-1)
            I2=IP

            DO 50 J2=IMIN,IMAX
               IND2=J2+N*(I2-1)
               R12(1)=R(IND2,1)-R(IND1,1)
               R12(2)=R(IND2,2)-R(IND1,2)
               R12(3)=R(IND2,3)-R(IND1,3)
               R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
               R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
               R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

               D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
               if (D12.GT.(3.*GAM)) then
                  goto 70
               endif

               U1(1)=R(IND1+1,1)-R(IND1,1)
               U1(2)=R(IND1+1,2)-R(IND1,2)
               U1(3)=R(IND1+1,3)-R(IND1,3)
               D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
               U1(1)=U1(1)/D1
               U1(2)=U1(2)/D1
               U1(3)=U1(3)/D1

               U2(1)=R(IND2+1,1)-R(IND2,1)
               U2(2)=R(IND2+1,2)-R(IND2,2)
               U2(3)=R(IND2+1,3)-R(IND2,3)
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

               EMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

               DESELF=DESELF-EMAG

 70            CONTINUE

 50         CONTINUE

            DO 80 J2=IMIN,IMAX
               IND2=J2+N*(I2-1)
               R12(1)=RP(IND2,1)-R(IND1,1)
               R12(2)=RP(IND2,2)-R(IND1,2)
               R12(3)=RP(IND2,3)-R(IND1,3)
               R12(1)=R12(1)-nint(R12(1)/LBOX)*LBOX
               R12(2)=R12(2)-nint(R12(2)/LBOX)*LBOX
               R12(3)=R12(3)-nint(R12(3)/LBOX)*LBOX

               D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)
               if (D12.GT.(3.*GAM)) then
                  goto 90
               endif

               U1(1)=R(IND1+1,1)-R(IND1,1)
               U1(2)=R(IND1+1,2)-R(IND1,2)
               U1(3)=R(IND1+1,3)-R(IND1,3)
               D1=sqrt(U1(1)**2.+U1(2)**2.+U1(3)**2.)
               U1(1)=U1(1)/D1
               U1(2)=U1(2)/D1
               U1(3)=U1(3)/D1

               U2(1)=RP(IND2+1,1)-RP(IND2,1)
               U2(2)=RP(IND2+1,2)-RP(IND2,2)
               U2(3)=RP(IND2+1,3)-RP(IND2,3)
               D2=sqrt(U2(1)**2.+U2(2)**2.+U2(3)**2.)
               U2(1)=U2(1)/D2
               U2(2)=U2(2)/D2
               U2(3)=U2(3)/D2

               U1U2=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
               if (U1U2.EQ.1.) then
                  goto 90
               endif

               GI(1)=U1(1)-U1U2*U2(1)
               GI(2)=U1(2)-U1U2*U2(2)
               GI(3)=U1(3)-U1U2*U2(3)

               S1=(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

               if (S1.GT.D1.OR.S1.LT.0.) then
                  goto 90
               endif

               GI(1)=U2(1)-U1U2*U1(1)
               GI(2)=U2(2)-U1U2*U1(2)
               GI(3)=U2(3)-U1U2*U1(3)

               S2=-(R12(1)*GI(1)+R12(2)*GI(2)+R12(3)*GI(3))/(1.-U1U2**2.)

               if (S2.GT.D2.OR.S2.LT.0.) then
                  goto 90
               endif

               R12(1)=R12(1)+S2*U2(1)-S1*U1(1)
               R12(2)=R12(2)+S2*U2(2)-S1*U1(2)
               R12(3)=R12(3)+S2*U2(3)-S1*U1(3)

               D12=sqrt(R12(1)**2.+R12(2)**2.+R12(3)**2.)

               if (D12.GT.LHC) then
                  goto 90
               endif

               EMAG=VHC*((LHC/D12)**12.-2.*(LHC/D12)**6.+1.)/12.

               DESELF=DESELF+EMAG

 90            CONTINUE
 80         CONTINUE


 40      CONTINUE
 100     CONTINUE
 30   CONTINUE

      RETURN
      END

!---------------------------------------------------------------*

