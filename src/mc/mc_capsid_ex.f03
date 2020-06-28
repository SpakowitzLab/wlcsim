!---------------------------------------------------------------*

!
! This subroutine calculates the change in the energy due to
! interaction with a spherical capsid.
!
! Corrections to force magnitude made 6-3-04.
!
! Andrew Spakowitz
! Written 6-29-04

      subroutine mc_capsid_ex(DEEX,R,NT,DR,I,RAD,VCAP)
      use params, only : dp
      implicit none
      integer NT            ! Number of beads
      real(dp) R(3,NT)        ! Bead positions
      real(dp) DR(3)        ! Change in bead position
      real(dp) RMAG,RMAGD    ! Dist from origin
      integer I                ! Current test index
      real(dp) RAD        ! Radius of capsid
      real(dp) VCAP        ! Capsid interaction
      real(dp) DEEX        ! Change in external energy

! Calculate the change in energy

      RMAG = sqrt(R(1,I)**2 + R(2,I)**2 + R(3,I)**2)
      RMAGD = sqrt((R(1,I) + DR(1))**2 + (R(2,I) + DR(2))**2 + (R(3,I) + DR(3))**2)
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
