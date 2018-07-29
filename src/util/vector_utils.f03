module vector_utils
use precision, only: dp, eps, pi
implicit none

contains

function distance(a,b)
    implicit none
    real(dp), dimension(3), intent(in) :: a,b
    real(dp) distance
    distance = sqrt( (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2 )
end function distance

function angle(a,b,c)
    implicit none
    real(dp) angle
    real(dp), dimension(3), intent(in) :: a,b,c
    angle = acos(dot_product(a-b,c-b)/(norm2(a-b)*norm2(c-b)))

end function angle

function cross(a,b)
    implicit none
    real(dp), dimension(3) :: cross
    real(dp), dimension(3), intent(in) :: a,b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

function rotateR(ROT,R)
    implicit none
    real(dp), intent(in) :: ROT(3,4) ! Rotation matrix
    real(dp), intent(in) :: R(3)
    real(dp) rotateR(3)
    rotateR(1) = ROT(1,4) + ROT(1,1)*R(1) + ROT(1,2)*R(2) + ROT(1,3)*R(3)
    rotateR(2) = ROT(2,4) + ROT(2,1)*R(1) + ROT(2,2)*R(2) + ROT(2,3)*R(3)
    rotateR(3) = ROT(3,4) + ROT(3,1)*R(1) + ROT(3,2)*R(2) + ROT(3,3)*R(3)
end function rotateR

function rotateU(ROT,U)
    implicit none
    real(dp), intent(in) :: ROT(3,4) ! Rotation matrix
    real(dp), intent(in) :: U(3)
    real(dp) rotateU(3)
    rotateU(1) = ROT(1,1)*U(1) + ROT(1,2)*U(2) + ROT(1,3)*U(3)
    rotateU(2) = ROT(2,1)*U(1) + ROT(2,2)*U(2) + ROT(2,3)*U(3)
    rotateU(3) = ROT(3,1)*U(1) + ROT(3,2)*U(2) + ROT(3,3)*U(3)
end function rotateU

function angle_of_triangle(a,b,c)
    implicit none
    real(dp), intent(in) :: a,b,c
    real(dp) angle_of_triangle
    angle_of_triangle=acos((a**2 - b**2 -c**2)/(-2.0_dp*b*c))
end function angle_of_triangle

function round_into_pm1(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp) round_into_pm1
    round_into_pm1 = min(x,1.0_dp)
    round_into_pm1 = max(x,-1.0_dp)
end function round_into_pm1

!! ---------------------------------------------------------------
!
!    Generates random unit vector
!
!------------------------------------------------------------
subroutine randomUnitVec(U,rand_stat)
    use mersenne_twister

    implicit none
    type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
    real(dp), intent(out) :: U(3)
    real(dp) z
    real(dp) theta
    real(dp) temp
    real(dp) urand(2)

    call random_number(urand,rand_stat)
    theta = urand(1)*2.0_dp*PI
    z = urand(2)*2.0_dp-1.0_dp
    temp = sqrt(1.0_dp-z**2)
    U(1) = temp*cos(theta)
    U(2) = temp*sin(theta)
    U(3) = z

    ! Alturnative algorithem
    !call random_number(urand,rand_stat)
    !ALPHA = 2.0_dp*PI*urand(1)
    !BETA = acos(2.0_dp*urand(2)-1.0_dp)
    !U(1) = sin(BETA)*cos(ALPHA)
    !U(2) = sin(BETA)*sin(ALPHA)
    !U(3) = cos(BETA)

    ! Alturnative algorithem
    !real(dp) urand(3)
    !call random_gauss(urand, rand_stat)
    !U = urand
    !U = U/norm2(U)
end subroutine


!--------------------------------------------------------------*
!
!           Calculates random perpendicular vectors
!
!    Quinn Made Changes to this file starting on 12/15/15
!    Quinn made this an independent file on 6/20/18
!
!---------------------------------------------------------------

subroutine random_perp(u,p,t,rand_stat)
    ! The subroutine generates the second two vectors in a unit triad
    ! The output vectors, p and t, are perpendicular to eachother and u
    ! The triad is randomly left or right handed
    use mersenne_twister
    implicit none
    !real(dp), PARAMETER :: PI = 3.141592654 ! Value of pi
    type(random_stat) rand_stat  ! status of random number generator
    real(dp) urnd(1) ! single random number

    real(dp) v(2) ! random 2-vec
    real(dp), intent(in) :: u(3) ! input
    real(dp), intent(out) :: p(3) ! output: random perpendicular to u
    real(dp), intent(out) :: t(3) ! orthogonal to p and u
    real(dp) f

    if (abs(u(1)**2 + u(2)**2 + u(3)**2-1.0_dp) .gt.eps) then
        print*, u
        print*, "Error in random_perp, please give me a unit vector"
        stop 1
    endif

    call random_number(urnd,rand_stat)
    v(1) = cos(2.0_dp*PI*urnd(1))
    v(2) = sin(2.0_dp*PI*urnd(1))

    if (u(3).gt.0.0) then
        f = 1.0_dp/(1.0_dp + u(3))
        p(1) = (u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
        p(2) = (u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
        p(3) = -1.0_dp*(u(2)*v(2) + u(1)*v(1))
    else
        f = 1.0_dp/(1.0_dp-u(3))
        p(1) = (-u(3) + f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
        p(2) = (-u(3) + f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
        p(3) = (u(2)*v(2) + u(1)*v(1))

    endif

    t(1) = u(2)*p(3)-u(3)*p(2)
    t(2) = u(3)*p(1)-u(1)*p(3)
    t(3) = u(1)*p(2)-u(2)*p(1)

    ! random sign
    call random_number(urnd,rand_stat)
    if (urnd(1).lt.0.5_dp) then
        t(1) = -1.0_dp*t(1)
        t(2) = -1.0_dp*t(2)
        t(3) = -1.0_dp*t(3)
    endif

    ! Testing
    !if (abs(dot_product(p,u)) > 0.000001_dp) then
    !    print*, "Error in random_perp, 1"
    !    stop 1
    !endif
    !if (abs(p(1)**2 + p(2)**2 + p(3)**2-1) .gt. 0.0000001_dp) then
    !    print*, "Error in random_perp, 2"
    !    stop 1
    !endif
    !if (abs(t(1)**2 + t(2)**2 + t(3)**2 -1).gt.0.000001_dp) then
    !    print*, "Error in random_perp, 3"
    !    stop 1
    !endif
    !if (abs(dot_product(t,p)) > 0.0000001_dp) then
    !    print*, "Error in random_perp, 4"
    !    stop 1
    !endif
    !if (abs(dot_product(t,u)) > 0.0000001_dp) then
    !    print*, "Error in random_perp, 5"
    !    stop 1
    !endif
    ! END Testing

    return
end subroutine

!--------------------------------------------------------------*
!
!    Calculates axis angle rotation matrix
!    Calculates the rotation in the right handed sense by angle alpha
!    in the TA direction.
!    The magnitude of TA is disreguarded.
!    The forth column is the offset
!
!    Quinn split out this file on 6/20/18
!
!---------------------------------------------------------------

subroutine axisAngle(ROT,alpha,TAin,P1)
implicit none
real(dp), intent(in) :: TAin(3)    ! Axis of rotation  (Need not be unit vector!!)
real(dp), intent(in) :: P1(3)    ! Point on rotation line
real(dp), intent(in) :: ALPHA    ! Angle of move
real(dp), intent(out) :: ROT(3,4) ! Rotation matrix
real(dp) MAG      ! Magnitude of vector
real(dp) TA(3)    ! Normalized vector

!print*, "PA", P1
MAG = sqrt(TAin(1)**2 + TAin(2)**2 + TAin(3)**2)
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

subroutine rotateAIntoB(A,B,P1,ROT)
implicit none
real(dp), intent(in) :: A(3)
real(dp), intent(in) :: B(3)
real(dp), intent(in) :: P1(3)    ! Point on rotation line
real(dp), intent(out) :: ROT(3,4) ! Rotation matrix
real(dp) ALPHA    ! Angle of move
real(dp) TA(3)    ! Normalized vector

TA = cross(A,B)
alpha = asin(norm2(TA)/(norm2(A)*norm2(B)))
call axisAngle(ROT,alpha,TA,P1)
end subroutine rotateAIntoB
end module
