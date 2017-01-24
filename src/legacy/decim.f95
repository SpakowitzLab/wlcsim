!---------------------------------------------------------------*

      SUBROUTINE decim(R,U,NT,N,NP,PARA,DT)

      PARAMETER (PI=3.141593)

      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION RT(NT,3)  ! Bead positions
      DOUBLE PRECISION UT(NT,3)  ! Tangent vectors
      DOUBLE PRECISION TTOT     ! Time of BD simulation
      INTEGER N,NP,NT           ! Number of beads
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION PVEC(60,8)
      INTEGER IND,CRS
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION M
      DOUBLE PRECISION DT
      INTEGER I,J,IP
      DOUBLE PRECISION L,LP

!     Reset the positions and orientations

      IND=1
      do 10 IP=1,NP
         I=1
         J=1
         do while (I.LE.N)
            RT(IND,1)=R(I+N*(IP-1),1)
            RT(IND,2)=R(I+N*(IP-1),2)
            RT(IND,3)=R(I+N*(IP-1),3)
            UT(IND,1)=U(I+N*(IP-1),1)
            UT(IND,2)=U(I+N*(IP-1),2)
            UT(IND,3)=U(I+N*(IP-1),3)
            I=I+2
            J=J+1
            IND=IND+1
         enddo
 10   CONTINUE
      N=J-1

      do 20 I=1,NP
         do 30 J=1,N
            IND=J+N*(I-1)
            R(IND,1)=RT(IND,1)
            R(IND,2)=RT(IND,2)
            R(IND,3)=RT(IND,3)
            U(IND,1)=UT(IND,1)
            U(IND,2)=UT(IND,2)
            U(IND,3)=UT(IND,3)
 30      continue
 20   continue

!     Load in the parameters for the simulation

      open (unit=5, file='input/input')
      read (unit=5, fmt='(4(/))')
      read (unit=5, fmt=*) LP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) L
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LBOX
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LHC
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) VHC
      close(5)
      L=L/LP
      DEL=L/(N-1.)

!     Load the tabulated parameters

      OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
      DO 40 I=1,60
         READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
 40   CONTINUE
      CLOSE(5)

      if (DEL.LT.PVEC(1,1)) then
         DEL=PVEC(1,1)
      endif
      if (DEL.GT.PVEC(60,1)) then
         DEL=PVEC(60,1)
      endif

      CRS=0
      IND=1
      do while (CRS.EQ.0)
         if (DEL.LE.PVEC(IND,1)) then
            CRS=1
         else
            IND=IND+1
         endif
      enddo

      I=2
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EB=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

      I=3
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      GAM=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

      I=4
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EPAR=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

      I=5
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EPERP=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

      I=6
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      ETA=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

      I=7
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      XIU=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

!      I=8
!      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
!      DT=XIU*(M*(DEL-PVEC(IND,1))+PVEC(IND,I))

      EB=EB/DEL
      EPAR=EPAR/DEL
      EPERP=EPERP/DEL
      GAM=DEL*GAM

      XIU=XIU*L/N
      XIR=L/N
      DT=0.5*XIU/(EPERP*GAM**2.)

      PARA(1)=EB
      PARA(2)=EPAR
      PARA(3)=EPERP
      PARA(4)=GAM
      PARA(5)=ETA
      PARA(6)=XIR
      PARA(7)=XIU
      PARA(8)=LBOX
      PARA(9)=LHC
      PARA(10)=VHC

      RETURN
      END

!---------------------------------------------------------------*
