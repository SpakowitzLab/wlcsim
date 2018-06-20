!--------------------------------------------------------------*
!
!    Calculates axis angle rotation matrix
!    The forth column is the offset
!
!    Quinn split out this file on 6/20/18
!
!---------------------------------------------------------------

subroutine axisAngle(ROT,alpha,TAin,P1)

use params, only: dp

implicit none
real(dp), intent(in) :: TAin(3)    ! Axis of rotation  (Need not be unit vector!!)
real(dp), intent(in) :: P1(3)    ! Point on rotation line
real(dp), intent(in) :: ALPHA    ! Angle of move
real(dp), intent(out) :: ROT(3,4) ! Rotation matrix
real(dp) MAG      ! Magnitude of vector
real(dp) TA(3)    ! Normalized vector

MAG = sqrt(TA(1)**2 + TA(2)**2 + TA(3)**2)
TA(1) = TAin(1)/MAG
TA(2) = TAin(2)/MAG
TA(3) = TAin(3)/MAG

ROT(1,1) = TA(1)**2 + (TA(2)**2 + TA(3)**2)*cos(ALPHA)
ROT(1,2) = TA(1)*TA(2)*(1.0_dp-cos(ALPHA))-TA(3)*sin(ALPHA)
ROT(1,3) = TA(1)*TA(3)*(1.0_dp-cos(ALPHA)) + TA(2)*sin(ALPHA)
ROT(1,4) = (P1(1)*(1.0_dp-TA(1)**2) &
     -TA(1)*(P1(2)*TA(2) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

ROT(2,1) = TA(1)*TA(2)*(1.0_dp-cos(ALPHA)) + TA(3)*sin(ALPHA)
ROT(2,2) = TA(2)**2 + (TA(1)**2 + TA(3)**2)*cos(ALPHA)
ROT(2,3) = TA(2)*TA(3)*(1.0_dp-cos(ALPHA))-TA(1)*sin(ALPHA)
ROT(2,4) = (P1(2)*(1.0_dp-TA(2)**2) &
     -TA(2)*(P1(1)*TA(1) + P1(3)*TA(3)))*(1.-cos(ALPHA)) + (P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

ROT(3,1) = TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
ROT(3,2) = TA(2)*TA(3)*(1.-cos(ALPHA)) + TA(1)*sin(ALPHA)
ROT(3,3) = TA(3)**2 + (TA(1)**2 + TA(2)**2)*cos(ALPHA)
ROT(3,4) = (P1(3)*(1.0_dp-TA(3)**2) &
     -TA(3)*(P1(1)*TA(1) + P1(2)*TA(2)))*(1.0_dp-cos(ALPHA)) + (P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
end subroutine
