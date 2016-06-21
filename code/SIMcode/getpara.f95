! *---------------------------------------------------------------*
!
!     subroutine getpara.f95
!     Setup the parameters for the simulation
!
!     1. Determine the simulation type
!     2. Evaluate the polymer elastic parameters
!     3. Determine the parameters for Brownian dynamics simulation
!
!     Andrew Spakowitz
!     8/17/15
!      

      SUBROUTINE getpara(PARA,DT,SIMTYPE)
      
      PARAMETER (PI=3.141593)
	  
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION PVEC(679,8)
      INTEGER IND,CRS
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION REND     ! Fixed end-to-end distance (dimensionless)
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION M
      DOUBLE PRECISION DT
      INTEGER I,N
      DOUBLE PRECISION L,LP
      INTEGER SIMTYPE           ! Simulation method (WLC=1,SSWLC=2,GC=3)
      
!     Load in the parameters for the simulation

      open (unit=5, file='input/input')
      read (unit=5, fmt='(4(/))')
      read (unit=5, fmt=*) LP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) L
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LBOX
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) REND
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) VHC
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) N
      close(5)
      L=L/LP
      LP=1.
      DEL=L/(N-1.)
      REND=REND*L

!     Load the tabulated parameters
	  
      OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
      DO 10 I=1,679
         READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
 10   CONTINUE 
      CLOSE(5)


!     Setup the parameters for WLC simulation
      
      if (DEL.LT.PVEC(1,1)) then
         EB=LP/DEL
         GAM=DEL
         XIR=L/N
         SIMTYPE=1
      endif

!    Setup the parameters for GC simulation

      if (DEL.GT.PVEC(679,1)) then
         EPAR=1.5/(DEL*LP**2.)
         GAM=0.
         SIMTYPE=3
      endif

!    Setup the parameters for ssWLC simulation

      if (DEL.GE.PVEC(1,1).AND.DEL.LE.PVEC(679,1)) then
         SIMTYPE=2
      
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

!         I=8
!         M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
!         DT=XIU*(M*(DEL-PVEC(IND,1))+PVEC(IND,I))
      
         EB=EB/DEL
         EPAR=EPAR/DEL
         EPERP=EPERP/DEL
         GAM=DEL*GAM
      
         XIU=XIU*L/N
         XIR=L/N
         DT=0.5*XIU/(EPERP*GAM**2.)

      endif

      PARA(1)=EB
      PARA(2)=EPAR
      PARA(3)=EPERP
      PARA(4)=GAM
      PARA(5)=ETA
      PARA(6)=XIR
      PARA(7)=XIU
      PARA(8)=LBOX
      PARA(9)=REND
      PARA(10)=DEL

      RETURN     
      END
      
!---------------------------------------------------------------*
