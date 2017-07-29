!---------------------------------------------------------------*

!
!     This subroutine analyses the structure of the expanding
!     filament.
!
!     Andrew Spakowitz
!     Written 9-1-04

!     Check the calculator of eigenvalue/eigenvector

      subroutine find_struc(R,NT,N,RCOM,DELR)

      real(dp) R(NT,3)  ! Bead positions
      real(dp) U(N-1,3) ! Tangent
      real(dp) B(N-1)   ! Bond length
      integer N,NT              ! Number of beads
      integer I

!     Structure analysis

      real(dp) RCOM(3)  ! Center of mass
      real(dp) DELR(3)  ! Mag of gyration tensor
      real(dp) TEMP(3)
      real(dp) DELRR(3) ! Real part
      real(dp) DELRI(3) ! Imaginary part
      real(dp) EVEC(3,3)
      real(dp) T(3,3)   ! Radius of gyration tensor
      real(dp) fv1(3)
      integer iv1(3)
      integer ERR

!     Find the center of mass

      RCOM(1) = 0.
      RCOM(2) = 0.
      RCOM(3) = 0.
      do 10 I = 1,N
         RCOM(1) = RCOM(1) + R(I,1)
         RCOM(2) = RCOM(2) + R(I,2)
         RCOM(3) = RCOM(3) + R(I,3)
 10   continue
      RCOM(1) = RCOM(1)/N
      RCOM(2) = RCOM(2)/N
      RCOM(3) = RCOM(3)/N

!     Find the principle radii of gyration

      T(1,1) = 0.
      T(1,2) = 0.
      T(1,3) = 0.
      T(2,1) = 0.
      T(2,2) = 0.
      T(2,3) = 0.
      T(3,1) = 0.
      T(3,2) = 0.
      T(3,3) = 0.
      do 20 I = 1,N
         T(1,1) = T(1,1) + (R(I,1)-RCOM(1))*(R(I,1)-RCOM(1))
         T(1,2) = T(1,2) + (R(I,1)-RCOM(1))*(R(I,2)-RCOM(2))
         T(1,3) = T(1,3) + (R(I,1)-RCOM(1))*(R(I,3)-RCOM(3))
         T(2,1) = T(2,1) + (R(I,2)-RCOM(2))*(R(I,1)-RCOM(1))
         T(2,2) = T(2,2) + (R(I,2)-RCOM(2))*(R(I,2)-RCOM(2))
         T(2,3) = T(2,3) + (R(I,2)-RCOM(2))*(R(I,3)-RCOM(3))
         T(3,1) = T(3,1) + (R(I,3)-RCOM(3))*(R(I,1)-RCOM(1))
         T(3,2) = T(3,2) + (R(I,3)-RCOM(3))*(R(I,2)-RCOM(2))
         T(3,3) = T(3,3) + (R(I,3)-RCOM(3))*(R(I,3)-RCOM(3))
 20   continue
      T(1,1) = T(1,1)/N
      T(1,2) = T(1,2)/N
      T(1,3) = T(1,3)/N
      T(2,1) = T(2,1)/N
      T(2,2) = T(2,2)/N
      T(2,3) = T(2,3)/N
      T(3,1) = T(3,1)/N
      T(3,2) = T(3,2)/N
      T(3,3) = T(3,3)/N

!      call rg(3,3,T,DELRR,DELRI,1,EVEC,iv1,fv1,ERR)

      DELR(1) = sqrt(DELRR(1))
      DELR(2) = sqrt(DELRR(2))
      DELR(3) = sqrt(DELRR(3))

      if (DELR(1) > DELR(2).AND.DELR(2) > DELR(3)) then
         TEMP(1) = DELR(1)
         TEMP(2) = DELR(2)
         TEMP(3) = DELR(3)
      elseif (DELR(2) > DELR(1).AND.DELR(1) > DELR(3)) then
         TEMP(1) = DELR(2)
         TEMP(2) = DELR(1)
         TEMP(3) = DELR(3)
      elseif (DELR(1) > DELR(3).AND.DELR(3) > DELR(2)) then
         TEMP(1) = DELR(1)
         TEMP(2) = DELR(3)
         TEMP(3) = DELR(2)
      elseif (DELR(2) > DELR(3).AND.DELR(3) > DELR(1)) then
         TEMP(1) = DELR(2)
         TEMP(2) = DELR(3)
         TEMP(3) = DELR(1)
      elseif (DELR(3) > DELR(1).AND.DELR(1) > DELR(2)) then
         TEMP(1) = DELR(3)
         TEMP(2) = DELR(1)
         TEMP(3) = DELR(2)
      elseif (DELR(3) > DELR(2).AND.DELR(2) > DELR(1)) then
         TEMP(1) = DELR(3)
         TEMP(2) = DELR(2)
         TEMP(3) = DELR(1)
      endif

      DELR(1) = TEMP(1)
      DELR(2) = TEMP(2)
      DELR(3) = TEMP(3)

      RETURN
      END

!---------------------------------------------------------------*
