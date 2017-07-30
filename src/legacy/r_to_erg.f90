!! ---------------------------------------------------------------*

!     Find the energetic terms from R

!     Revised 6-22-04

      subroutine r_to_erg(R,NT,N,ECOM,EBEND)

      real(dp) R(NT,3)  ! Bead positions
      real(dp) U(N-1,3) ! Tangent
      real(dp) B(N-1)   ! Bond length
      integer N,NT              ! Number of beads
      integer I
      real(dp) KAP,EPS  ! DNA mat props
      real(dp) L0       ! Equilibrium
                                ! bead separation
      real(dp) doT      ! Dot product

!     Variables in the tube interaction

      real(dp) LHC      ! Length of HC int
      real(dp) VHC      ! HC strength
      real(dp) RAD      ! Radius of int
      real(dp) SIGP     ! Poly diameter of int
      real(dp) forCE

!     Energy variables

      real(dp) ECOM     ! Compression energy
      real(dp) EBEND    ! Bending energy
      real(dp) EEX      ! External energy
      real(dp) EPONP    ! Poly energy

!     Setup the properties

      open (unit = 5, file = 'input/input')
      read (unit = 5, fmt = '(4(/))')
      read (unit = 5, fmt = *) KAP
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) EPS
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) L0
      close(5)

!     Setup the necessary quantities

      do 10 I = 1,(N-1)
         B(I) = sqrt((R(I + 1,1)-R(I,1))**2. + (R(I + 1,2)-R(I,2))**2. + (R(I + 1,3)-R(I,3))**2.)
         U(I,1) = (R(I + 1,1)-R(I,1))/B(I)
         U(I,2) = (R(I + 1,2)-R(I,2))/B(I)
         U(I,3) = (R(I + 1,3)-R(I,3))/B(I)
 10   continue


!     Calculate the energetic contributions

      ECOM = 0.
      do 30 I = 1,(N-1)
         ECOM = ECOM + (B(I)-L0)**2.
 30   continue
      ECOM = 0.5*KAP*ECOM

      EBEND = 0.
      do 40 I = 1,(N-2)
         doT = U(I + 1,1)*U(I,1) + U(I + 1,2)*U(I,2) + U(I + 1,3)*U(I,3)
         EBEND = EBEND + 1.-doT
 40   continue
      EBEND = EPS*EBEND

!     call energy_ex(EEX,R,N,LHC,VHC,RAD)
!     call energy_ponp(EPONP,R,N,LB,LHC,VHC,SIGP)

      RETURN
      END

!---------------------------------------------------------------*
