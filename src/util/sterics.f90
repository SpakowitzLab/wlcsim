! npagane | risca lab | dec 2019 | sterics calculation file

! define line-line intersection module that contains calculation and test
! inspired from http://paulbourke.net/geometry/pointlineplane/
MODULE LineLineIntersection
    use params, only: dp
    implicit none
    contains

    ! calculation
    FUNCTION LineLineIntersectionCalculation(A1,A2,B1,B2)
        implicit none
        ! A1 is the point: (ax1, ay1, az1)
        ! A2 is the point: (ax2, ay2, az2)
        ! A1->A2 makes the line segment A (likewise for B)
        real(dp), intent(in), dimension(3) :: A1
        real(dp), intent(in), dimension(3) :: A2
        real(dp), intent(in), dimension(3) :: B1
        real(dp), intent(in), dimension(3) :: B2
        real(dp), parameter ::  tol = 1.0e-8! tolerance for cooccupancy (should be small to disallow overlap)
        real(dp), dimension(3) :: pA ! closest point on A to B
        real(dp), dimension(3) :: pB ! closest point on B to A
        real(dp) dotA1B1B2B1, dotB2B1A2A1, dotA1B1A2A1, dotB2B1B2B1, dotA2A1A2A1
        real(dp) muA, muB
        real(dp), dimension(3) :: vecA, vecB, tA2, tB1, tB2
        integer LineLineIntersectionCalculation

        ! default value set to 0
        LineLineIntersectionCalculation = 0

        ! check for overlap of points
        if (ALL(ABS(A1-B1) <= tol) .OR. ALL(ABS(A1-B2) <= tol) .OR. ALL(ABS(A2-B1) <= tol) .OR. ALL(ABS(A2-B2) <= tol)) then
            print*, "collision, point overlap"
            LineLineIntersectionCalculation = 1
            return
        endif

        ! check if lines are parallel
        vecA = (A2-A1)/NORM2(A2-A1)
        vecB = (B2-B1)/NORM2(B2-B1)
        if ( ALL(ABS(vecA) - ABS(vecB) <= tol) ) then
            ! try to find overlap
            tA2 = (A2-A1)/vecA
            tB1 = (B1-A1)/vecA
            tB2 = (B2-A1)/vecA
            ! ensure all the components are the same
            if ( tA2(1)==tA2(2) .AND. tA2(2)==tA2(3) .AND. tB1(1)==tB1(2) .AND. tB1(2)==tB1(3) &
              .AND. tB2(1)==tB2(2) .AND. tB2(2)==tB2(3) ) then
                ! take just the first component
                if (( (tA2(1)*tB1(1) > 0) .OR. (tA2(1)*tB2(1) > 0) ) &
                  .AND. ( (abs(tA2(1)) >= abs(tB1(1))) .OR. (abs(tA2(1)) >= abs(tB2(1))) )) then
                    !print*, "collision, parallel overlap"
                    !print*, A1(1), ',',A1(2), ',',A1(3)
                    !print*, A2(1), ',',A2(2), ',',A2(3)
                    !print*, B1(1),',', B1(2), ',',B1(3)
                    !print*, B2(1), ',',B2(2), ',',B2(3)
                    LineLineIntersectionCalculation = 1
                    return ! quit early, coincident lines
                else
                    return ! quit early, non-coincident lines
                endif
            else
                !print*, "no collision, parallel but not overlap"
                return ! quit early, parallel and separated lines
            endif
        endif

        ! find shortest line between A and B
        dotA1B1B2B1 = dot_product(A1-B1, B2-B1)
        dotB2B1A2A1 = dot_product(B2-B1, A2-A1)
        dotA1B1A2A1 = dot_product(A1-B1, A2-A1)
        dotB2B1B2B1 = dot_product(B2-B1, B2-B1)
        dotA2A1A2A1 = dot_product(A2-A1, A2-A1)
        muA = (dotA1B1B2B1*dotB2B1A2A1 - dotA1B1A2A1*dotB2B1B2B1) / (dotA2A1A2A1*dotB2B1B2B1 - dotB2B1A2A1*dotB2B1A2A1)
        muB = (dotA1B1B2B1 + mua*dotB2B1A2A1) / dotB2B1B2B1
        pA = A1 + muA * (A2-A1)
        pB = B1 + muB * (B2-B1)
        ! check if dist <= tol
        if ((sqrt(dot_product(pA-pB, pA-pB)) <= tol)  .AND. (abs(muA) <= 1) .AND. (abs(muB) <= 1) &
          .AND. .NOT. (abs(-1 - muA*muB) <= tol)) then ! check for -1 and 1 pairs  
            print*, "collision, intersect"
            !print*, muA, muB, sqrt(dot_product(pA-pB, pA-pB)), abs(-1 - muA*muB)
            !print*, A1(1), ',',A1(2),',', A1(3)
            !print*, A2(1), ',',A2(2), ',',A2(3)
            !print*, B1(1),',', B1(2),',', B1(3)
            !print*, B2(1), ',',B2(2), ',',B2(3)
            LineLineIntersectionCalculation = 1
            return 
        endif
        return 

    END FUNCTION LineLineIntersectionCalculation

    ! testing functions used to live here !
    ! have since migrated to nicole_util/3D_collisions !

END MODULE LineLineIntersection

! define sphere-line intersection module that contains calculation and test
! inspired from http://paulbourke.net/geometry/circlesphere/
MODULE SphereLineIntersection
    use params, only: dp
    implicit none
    contains

    ! calculation
    FUNCTION SphereLineIntersectionCalculation(A1,A2,B1,r)
        implicit none
        ! A1 is the point: (ax1, ay1, az1)
        ! A2 is the point: (ax2, ay2, az2)
        ! B1 is the center of the sphere (x,y,z) of radius r
        real(dp), intent(in), dimension(3) :: A1
        real(dp), intent(in), dimension(3) :: A2
        real(dp), intent(in), dimension(3) :: B1
        real(dp), intent(in) :: r
        real(dp) :: a, b, c, discr, um, up ! quadratic variables
        integer SphereLineIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereLineIntersectionCalculation = 0

        ! calculate discriminant
        a = dot_product(A2-A1, A2-A1)
        b = 2*dot_product(A2-A1, A1-B1)
        c = dot_product(B1, B1) + dot_product(A1, A1) - 2*dot_product(B1, A1) - r*r
        discr = b*b - 4*a*c
        if ( discr < 0 ) then 
            return
        endif
        up = (-b + sqrt(discr))/(2*a)
        um = (-b - sqrt(discr))/(2*a)
        if ((um > 1 .AND. up < 0) .OR. (um < 0 .AND. up > 1)) then 
            SphereLineIntersectionCalculation = 1
            !print*, "line in sphere"
            !print*, um, up
            !print*, A1(1), ',',A1(2), ',',A1(3)
            !print*, A2(1), ',',A2(2), ',',A2(3)
            !print*, B1(1),',', B1(2), ',',B1(3), ',',r
            return
        else if ((um >= 0 .AND. um <= 1) .OR. (up >= 0 .AND. up <= 1)) then 
            SphereLineIntersectionCalculation = 1
            !print*, "one or more line-sphere intersections"
            !print*, um, up
            !print*, A1(1), ',', A1(2), ',', A1(3)
            !print*, A2(1), ',', A2(2), ',', A2(3)
            !print*, B1(1), ',', B1(2), ',', B1(3), ',', r
            return
        endif
        return 

    END FUNCTION SphereLineIntersectionCalculation

    ! testing functions used to live here !
    ! have since migrated to nicole_util/3D_collisions !

END MODULE SphereLineIntersection

! define sphere-sphere intersection module that contains calculation and test
MODULE SphereSphereIntersection
    use params, only: dp
    implicit none
    contains

    ! calculation
    FUNCTION SphereSphereIntersectionCalculation(A1,ra,B1,rb)
        implicit none
        ! A1 is the center of the sphere (x,y,z) of radius ra
        ! B1 is the center of the sphere (x,y,z) of radius rb
        real(dp), intent(in), dimension(3) :: A1
        real(dp), intent(in) :: ra
        real(dp), intent(in), dimension(3) :: B1
        real(dp), intent(in) :: rb
        real(dp), parameter :: dist = 1.0e-5 ! tolerance for collision (should be small)
        integer SphereSphereIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereSphereIntersectionCalculation = 0

        ! see if radii overalp
        if ( sqrt(dot_product(A1-B1, A1-B1)) - (ra+rb) < dist ) then
            SphereSphereIntersectionCalculation = 1
            !print*, "sphere collision"
            !print*, A1(1), ',', A1(2), ',', A1(3), ',', ra
            !print*, B1(1), ',', B1(2), ',', B1(3), ',', rb
            return
        endif
        return 

    END FUNCTION SphereSphereIntersectionCalculation

    ! testing functions used to live here !
    ! have since migrated to nicole_util/3D_collisions !

END MODULE SphereSphereIntersection