!---------------------------------------------------------------*

!
! This subroutine calculates the change in the energy due to
! interaction with a spherical capsid.
!
! Corrections to force magnitude made 6-3-04.
!
! Andrew Spakowitz
! Written 6-29-04

      SUBROUTINE MC_ex(DEEX,R,NT,N,DR,I,RAD,VCAP)

      DOUBLE PRECISION R(NT,3)        ! Bead positions
      DOUBLE PRECISION DR(3)        ! Change in bead position
      DOUBLE PRECISION RMAG,RMAGD    ! Dist from origin
      INTEGER I                ! Current test index
      DOUBLE PRECISION RAD        ! Radius of capsid
      DOUBLE PRECISION VCAP        ! Capsid interaction
      INTEGER N,NT            ! Number of beads
      DOUBLE PRECISION DEEX        ! Change in external energy

! Calculate the change in energy

      RMAG=sqrt(R(I,1)**2+R(I,2)**2+R(I,3)**2)
      RMAGD=sqrt((R(I,1)+DR(1))**2+(R(I,2)+DR(2))**2+(R(I,3)+DR(3))**2)
      if (RMAG.GT.RAD.AND.RMAGD.GT.RAD) then
         DEEX=VCAP*((RMAGD-RAD)**4.-(RMAG-RAD)**4.)
      elseif (RMAG.GT.RAD.AND.RMAGD.LE.RAD) then
         DEEX=-VCAP*((RMAG-RAD)**4.)
      elseif (RMAG.LE.RAD.AND.RMAGD.GT.RAD) then
         DEEX=VCAP*((RMAGD-RAD)**4.)
      endif

      RETURN
      END

!---------------------------------------------------------------*