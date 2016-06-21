! ---------------------------------------------------------------*
      
!     
!     subroutine initcond.f95
!     Set the initial conformation of the polymer chains
!     from file or from random initialization.
!     
!     Andrew Spakowitz
!     8/17/15
      
      SUBROUTINE initcond(R,U,NT,N,NP,IDUM,FRMFILE,PARA)

      use mt19937, only : grnd, init_genrand, rnorm, mt, mti
      
      PARAMETER (PI=3.141593)
      
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      INTEGER N,NP,NT           ! Number of beads
      DOUBLE PRECISION GAM       ! Equil bead separation
      DOUBLE PRECISION LBOX     ! Box edge length
      INTEGER I,J,IB            ! Index Holders
      REAL ran1                 ! Random number generator
      INTEGER FRMFILE           ! Is conformation in file?
      INTEGER INPUT             ! Is input file set?
      DOUBLE PRECISION RMIN
      DOUBLE PRECISION R0(3)
      DOUBLE PRECISION PARA(10)

!     Parameters for end constraint

      INTEGER IMID
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION REND
      DOUBLE PRECISION DEL
      
!     Variables in the simulation
      
      DOUBLE PRECISION KAP,EPS  ! Elastic properties
      DOUBLE PRECISION XI       ! Drag coefficients

!     Random number generator initiation

      integer IDUM
      character*8 datedum
      character*10 timedum
      character*5 zonedum
      integer seedvalues(8)
      
!     Setup the choice parameters
      
      INPUT=1

!     Seed the random number generator off the computer clock

      call date_and_time(datedum,timedum,zonedum,seedvalues)
	
! concatenate filename, time within mins, secs, millisecs to seed random number generator	

      IDUM=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
      call init_genrand(IDUM)   
      
!     Input the conformation if FRMFILE=1
      
      if(FRMFILE.EQ.1)then
         OPEN (UNIT = 5, FILE = 'initial/r', STATUS = 'OLD')
         IB=1
         DO 10 I=1,NP
            DO 20 J=1,N
               READ(5,*) R(IB,1),R(IB,2),R(IB,3)
               IB=IB+1
 20         CONTINUE
 10      CONTINUE 
         CLOSE(5)

         OPEN (UNIT = 5, FILE = 'initial/u', STATUS = 'OLD')
         IB=1
         DO 30 I=1,NP
            DO 40 J=1,N
               READ(5,*) U(IB,1),U(IB,2),U(IB,3)
               IB=IB+1
 40         CONTINUE
 30      CONTINUE 
         CLOSE(5)
         
      endif
      
!     Set the initial conformation to a straight chain if CHOICE=1
      
      if(FRMFILE.EQ.0)then
         
!     Fix the initial condition
         
         if (INPUT.EQ.0) then
            LBOX=10.
            GAM=1.
         else
            GAM=PARA(4)
            LBOX=PARA(8)
         endif

         REND=PARA(9)
         DEL=PARA(10)
         IMID=nint((N+1.)/2.)
         ALPHA=acos(((REND/DEL)**2.-((N-IMID)**2.+(IMID-1)**2.))/(2.*(IMID-1.)*(N-IMID)))

         IB=1
         DO 50 I=1,NP
!            R0(1)=grnd()*LBOX
!            R0(2)=grnd()*LBOX
!            R0(3)=grnd()*LBOX
            DO 50 J=1,N
               if (J.LE.IMID) then
                  R(IB,1)=DEL*(J-1.)
                  R(IB,2)=0
                  R(IB,3)=0
                  U(IB,1)=1.
                  U(IB,2)=0.
                  U(IB,3)=0.			   
               else
                  R(IB,1)=DEL*(IMID-1.)+DEL*(J-IMID)*cos(ALPHA)
                  R(IB,2)=DEL*(J-IMID)*sin(ALPHA)
                  R(IB,3)=0
                  U(IB,1)=cos(ALPHA)
                  U(IB,2)=sin(ALPHA)
                  U(IB,3)=0.			   
               endif
               IB=IB+1
 60         CONTINUE
 50      CONTINUE
         
      endif      
      
      RETURN     
      END
      
!---------------------------------------------------------------*
