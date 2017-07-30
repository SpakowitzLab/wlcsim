!---------------------------------------------------------------*

!
! This subroutine calculates the change in the energy due to
! interaction with a spherical capsid.
!
! Corrections to force magnitude made 6-3-04.
!
! Andrew Spakowitz
! Written 6-29-04

      subroutine MC_ex(DEEX,R,NT,N,DR,I,RAD,VCAP)
      use params, only : dp

      real(dp) R(NT,3)        ! Bead positions
      real(dp) DR(3)        ! Change in bead position
      real(dp) RMAG,RMAGD    ! Dist from origin
      integer I                ! Current test index
      real(dp) RAD        ! Radius of capsid
      real(dp) VCAP        ! Capsid interaction
      integer N,NT            ! Number of beads
      real(dp) DEEX        ! Change in external energy

! Calculate the change in energy

      RMAG = sqrt(R(I,1)**2 + R(I,2)**2 + R(I,3)**2)
      RMAGD = sqrt((R(I,1) + DR(1))**2 + (R(I,2) + DR(2))**2 + (R(I,3) + DR(3))**2)
      if (RMAG > RAD.AND.RMAGD > RAD) then
         DEEX = VCAP*((RMAGD-RAD)**4.-(RMAG-RAD)**4.)
      elseif (RMAG > RAD.AND.RMAGD <= RAD) then
         DEEX = -VCAP*((RMAG-RAD)**4.)
      elseif (RMAG <= RAD.AND.RMAGD > RAD) then
         DEEX = VCAP*((RMAGD-RAD)**4.)
      endif

      RETURN
      END

!---------------------------------------------------------------*
