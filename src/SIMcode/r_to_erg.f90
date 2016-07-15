!! ---------------------------------------------------------------*

!     Find the energetic terms from R

!     Revised 6-22-04

      SUBROUTINE r_to_erg(R,NT,N,ECOM,EBEND)

      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(N-1,3) ! Tangent
      DOUBLE PRECISION B(N-1)   ! Bond length
      INTEGER N,NT              ! Number of beads
      INTEGER I
      DOUBLE PRECISION KAP,EPS  ! DNA mat props
      DOUBLE PRECISION L0       ! Equilibrium
                                ! bead separation
      DOUBLE PRECISION DOT      ! Dot product

!     Variables in the tube interaction

      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION RAD      ! Radius of int
      DOUBLE PRECISION SIGP     ! Poly diameter of int
      DOUBLE PRECISION FORCE

!     Energy variables

      DOUBLE PRECISION ECOM     ! Compression energy
      DOUBLE PRECISION EBEND    ! Bending energy
      DOUBLE PRECISION EEX      ! External energy
      DOUBLE PRECISION EPONP    ! Poly energy

!     Setup the properties

      open (unit=5, file='input/input')
      read (unit=5, fmt='(4(/))')
      read (unit=5, fmt=*) KAP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) EPS
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) L0
      close(5)

!     Setup the necessary quantities

      DO 10 I=1,(N-1)
         B(I)=sqrt((R(I+1,1)-R(I,1))**2.+(R(I+1,2)-R(I,2))**2.+(R(I+1,3)-R(I,3))**2.)
         U(I,1)=(R(I+1,1)-R(I,1))/B(I)
         U(I,2)=(R(I+1,2)-R(I,2))/B(I)
         U(I,3)=(R(I+1,3)-R(I,3))/B(I)
 10   CONTINUE


!     Calculate the energetic contributions

      ECOM=0.
      DO 30 I=1,(N-1)
         ECOM=ECOM+(B(I)-L0)**2.
 30   CONTINUE
      ECOM=0.5*KAP*ECOM

      EBEND=0.
      DO 40 I=1,(N-2)
         DOT=U(I+1,1)*U(I,1)+U(I+1,2)*U(I,2)+U(I+1,3)*U(I,3)
         EBEND=EBEND+1.-DOT
 40   CONTINUE
      EBEND=EPS*EBEND

!     call energy_ex(EEX,R,N,LHC,VHC,RAD)
!     call energy_ponp(EPONP,R,N,LB,LHC,VHC,SIGP)

      RETURN
      END

!---------------------------------------------------------------*
