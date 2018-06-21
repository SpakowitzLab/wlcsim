module vector_utils
implicit none

contains 


function cross(a,b)
    use params, only: dp
    real(dp), dimension(3) :: cross
    real(dp), dimension(3), intent(in) :: a,b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

!! ---------------------------------------------------------------
!
!    Generates random unit vector
!
!------------------------------------------------------------
subroutine randomUnitVec(U,rand_stat)

    use mersenne_twister
    use params, only: dp, pi
    
    implicit none
    type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
    real(dp), intent(out) :: U(3)
    real(dp) z
    real(dp) theta
    real(dp) temp
    real urand(2)
    
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
    !real urand(3)
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
    use params, only: pi, dp, eps
    implicit none
    !real(dp), PARAMETER :: PI = 3.141592654 ! Value of pi
    type(random_stat) rand_stat  ! status of random number generator
    real urnd(1) ! single random number
    
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
    !if (abs(p(1)*u(1) + p(2)*u(2) + p(3)*u(3)).gt.0.000001_dp) then
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
    !if (abs(t(1)*p(1) + t(2)*p(2) + t(3)*p(3)).gt.0.0000001_dp) then
    !    print*, "Error in random_perp, 4"
    !    stop 1
    !endif
    !if (abs(t(1)*u(1) + t(2)*u(2) + t(3)*u(3)).gt.0.0000001_dp) then
    !    print*, "Error in random_perp, 5"
    !    stop 1
    !endif
    ! END Testing
    
    return
end subroutine

end module