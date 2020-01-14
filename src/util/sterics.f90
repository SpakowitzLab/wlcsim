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
        real(dp), parameter ::  tol = 1.0e-3 ! tolerance for cooccupancy (should be small to disallow overlap)
        real(dp), parameter :: dist = 1.0e-5 ! tolerance for collision (pseudo thickness of line)
        real(dp), dimension(3) :: pA ! closest point on A to B
        real(dp), dimension(3) :: pB ! closest point on B to A
        real(dp) dotA1B1B2B1, dotB2B1A2A1, dotA1B1A2A1, dotB2B1B2B1, dotA2A1A2A1
        real(dp) muA, muB
        real(dp), dimension(3) :: vecA, vecB, tA2, tB1, tB2
        logical LineLineIntersectionCalculation 

        ! default value (i.e. no intersection)
        LineLineIntersectionCalculation = .FALSE.

        ! check for overlap of points
        if (ALL(ABS(A1-B1) <= tol) .OR. ALL(ABS(A1-B2) <= tol) .OR. ALL(ABS(A2-B1) <= tol) .OR. ALL(ABS(A2-B2) <= tol)) then
            !print*, "collision, point overlap"
            LineLineIntersectionCalculation = .TRUE.
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
            ! ensure all the components are the same (HMMM... should i use a tol instead of equiv)
            if ( tA2(1)==tA2(2) .AND. tA2(2)==tA2(3) .AND. tB1(1)==tB1(2) .AND. tB1(2)==tB1(3) &
              .AND. tB2(1)==tB2(2) .AND. tB2(2)==tB2(3) ) then
                ! take just the first component
                if (( (tA2(1)*tB1(1) > 0) .OR. (tA2(1)*tB2(1) > 0) ) &
                  .AND. ( (abs(tA2(1)) >= abs(tB1(1))) .OR. (abs(tA2(1)) >= abs(tB2(1))) )) then
                    print*, "collision, parallel overlap"
                    LineLineIntersectionCalculation = .TRUE.
                    return ! quit early
                else
                    return ! quit early
                endif
            else
                !print*, "no collision, parallel but not overlap"
                return ! quit early
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
        if ((sqrt(dot_product(pA-pB, pA-pB)) <= dist)  .AND. (abs(muA) <= 1) .AND. (abs(muB) <= 1)) then 
            print*, "collision, intersect"
            print*, sqrt(dot_product(pA-pB, pA-pB)), muA, muB
            print*, A1
            print*, A2
            print*, B1
            print*, B2
            LineLineIntersectionCalculation = .TRUE.
            return 
        endif
        return 

    END FUNCTION LineLineIntersectionCalculation

    ! testing routines !
    SUBROUTINE LineLineIntersectionTestOverlapA1B1
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/) 
        real(dp), dimension(3) :: A2 = (/1,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp), dimension(3) :: B2 = (/-1,-1,-1/)
        logical :: val 
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then 
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA1B1
  
    SUBROUTINE LineLineIntersectionTestOverlapA1B2
        implicit none
        real(dp), dimension(3) :: A1 = (/-11,-11,-11/)
        real(dp), dimension(3) :: A2 = (/-10,-10,-10/)
        real(dp), dimension(3) :: B1 = (/1,1,1/)
        real(dp), dimension(3) :: B2 = (/-10,-10,-10/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then 
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA1B2

    SUBROUTINE LineLineIntersectionTestOverlapA2B1
        implicit none
        real(dp), dimension(3) :: A1 = (/1,0,0/)
        real(dp), dimension(3) :: A2 = (/0,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp), dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA2B1

    SUBROUTINE LineLineIntersectionTestOverlapA2B2
        implicit none
        real(dp), dimension(3) :: A1 = (/1,0,0/)
        real(dp), dimension(3) :: A2 = (/0.0001,0.00002,0.00004/)
        real(dp), dimension(3) :: B1 = (/0,0,1/)
        real(dp), dimension(3) :: B2 = (/0.0001,0.00002,0.00004/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA2B2

    SUBROUTINE LineLineIntersectionTestSameLine
        implicit none
        real(dp), dimension(3) :: A1 = (/1,1,1/)
        real(dp), dimension(3) :: A2 = (/0,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp), dimension(3) :: B2 = (/1,1,1/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestSameLine"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestSameLine

    SUBROUTINE LineLineIntersectionTestParallelOverlapA1B1
        implicit none
        real(dp), dimension(3) :: A1 = (/2,2,2/)
        real(dp), dimension(3) :: A2 = (/0,0,0/)
        real(dp), dimension(3) :: B1 = (/1,1,1/)
        real(dp), dimension(3) :: B2 = (/3,3,3/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA1B1

    SUBROUTINE LineLineIntersectionTestParallelOverlapA1B2
        implicit none
        real(dp), dimension(3) :: A1 = (/2,2,2/)
        real(dp), dimension(3) :: A2 = (/-100,-100,-100/)
        real(dp), dimension(3) :: B1 = (/100,100,100/)
        real(dp), dimension(3) :: B2 = (/1.5,1.5,1.5/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA1B2

    SUBROUTINE LineLineIntersectionTestParallelOverlapA2B1
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/-100,-100,-100/)
        real(dp), dimension(3) :: B1 = (/-99,-99,-99/)
        real(dp), dimension(3) :: B2 = (/-1001,-1001,-1001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA2B1

    SUBROUTINE LineLineIntersectionTestParallelOverlapA2B2
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/1,1,1/)
        real(dp), dimension(3) :: B1 = (/2,2,2/)
        real(dp), dimension(3) :: B2 = (/0.999, 0.999, 0.999/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA2B2

    SUBROUTINE LineLineIntersectionTestParallelA1B1
        implicit none
        real(dp), dimension(3) :: A1 = (/2,2,2/)
        real(dp), dimension(3) :: A2 = (/0,0,0/)
        real(dp), dimension(3) :: B1 = (/2.1,2.1,2.1/)
        real(dp), dimension(3) :: B2 = (/3,3,3/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA1B1

    SUBROUTINE LineLineIntersectionTestParallelA1B2
        implicit none
        real(dp), dimension(3) :: A1 = (/1,1,1/)
        real(dp), dimension(3) :: A2 = (/-100,-100,-100/)
        real(dp), dimension(3) :: B1 = (/100,100,100/)
        real(dp), dimension(3) :: B2 = (/1.5,1.5,1.5/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA1B2

    SUBROUTINE LineLineIntersectionTestParallelA2B1
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/-100,-100,-100/)
        real(dp), dimension(3) :: B1 = (/-100.0001,-100.0001,-100.0001/)
        real(dp), dimension(3) :: B2 = (/-1001,-1001,-1001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA2B1

    SUBROUTINE LineLineIntersectionTestParallelA2B2
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/0.999, 0.999, 0.999/)
        real(dp), dimension(3) :: B1 = (/2,2,2/)
        real(dp), dimension(3) :: B2 = (/1.001,1.001,1.001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA2B2

    SUBROUTINE LineLineIntersectionTestIntersectA1
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/1,0,0/)
        real(dp), dimension(3) :: B1 = (/0,-1,0/)
        real(dp), dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectA1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectA1

    SUBROUTINE LineLineIntersectionTestIntersectA2
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/0,0,100/)
        real(dp), dimension(3) :: B1 = (/-200,0,100/)
        real(dp), dimension(3) :: B2 = (/200,0,100/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectA2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectA2

    SUBROUTINE LineLineIntersectionTestIntersectB1
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,0,1/)
        real(dp), dimension(3) :: A2 = (/1,0,1/)
        real(dp), dimension(3) :: B1 = (/0,0,1/)
        real(dp), dimension(3) :: B2 = (/0,0,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectB1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectB1

    SUBROUTINE LineLineIntersectionTestIntersectB2
        implicit none
        real(dp), dimension(3) :: A1 = (/0,-1,0/)
        real(dp), dimension(3) :: A2 = (/0,1,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp), dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectB2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectB2

    SUBROUTINE LineLineIntersectionTestIntersectMiddle
        implicit none
        real(dp), dimension(3) :: A1 = (/0,-2,0/)
        real(dp), dimension(3) :: A2 = (/0,2,0/)
        real(dp), dimension(3) :: B1 = (/-2,0,0/)
        real(dp), dimension(3) :: B2 = (/2,0,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectMiddle"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectMiddle

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideZ
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/0,0,5/)
        real(dp), dimension(3) :: B1 = (/-1,0,2/)
        real(dp), dimension(3) :: B2 = (/1,0,2/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideZ"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideZ

    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideZ
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,0,5/)
        real(dp), dimension(3) :: A2 = (/1,0,5/)
        real(dp), dimension(3) :: B1 = (/-1,0,2/)
        real(dp), dimension(3) :: B2 = (/1,0,2/)
        logical :: val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideZ"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideZ

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideY
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/0,5,0/)
        real(dp), dimension(3) :: B1 = (/-1,2,0/)
        real(dp), dimension(3) :: B2 = (/1,2,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideY"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideY

    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideY
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,5,0/)
        real(dp), dimension(3) :: A2 = (/1,5,0/)
        real(dp), dimension(3) :: B1 = (/-1,2,0/)
        real(dp), dimension(3) :: B2 = (/1,2,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideY"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideY

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideX
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/5,0,0/)
        real(dp), dimension(3) :: B1 = (/2,0,-1/)
        real(dp), dimension(3) :: B2 = (/2,0,1/)
        logical :: val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideX"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideX
        
    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideX
        implicit none
        real(dp), dimension(3) :: A1 = (/5,0,-1/)
        real(dp), dimension(3) :: A2 = (/5,0,1/)
        real(dp), dimension(3) :: B1 = (/2,0,-1/)
        real(dp), dimension(3) :: B2 = (/2,0,1/)
        logical :: val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideX"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideX

    SUBROUTINE LineLineIntersectionTestIntersectProjection
        implicit none
        real(dp), dimension(3) :: A1 = (/0,0,0/)
        real(dp), dimension(3) :: A2 = (/5,5,0/)
        real(dp), dimension(3) :: B1 = (/0,2,2/)
        real(dp), dimension(3) :: B2 = (/2,0,1/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjection"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjection

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
        logical SphereLineIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereLineIntersectionCalculation = .FALSE.

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
            SphereLineIntersectionCalculation = .TRUE.
            !print*, "line in sphere"
            return
        else if ((um >= 0 .AND. um <= 1) .OR. (up >= 0 .AND. up <= 1)) then 
            SphereLineIntersectionCalculation = .TRUE.
            !print*, "one or more intersections"
            return
        endif
        return 

    END FUNCTION SphereLineIntersectionCalculation

    ! testing routines !
    SUBROUTINE SphereLineIntersectionTestLineInside
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,0,0/)
        real(dp), dimension(3) :: A2 = (/-2,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 20.0
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInside"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInside

    SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA1
        implicit none
        real(dp), dimension(3) :: A1 = (/-2,0,0/)
        real(dp), dimension(3) :: A2 = (/-1,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInsideEdgeA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA1

    SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA2
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,0,0/)
        real(dp), dimension(3) :: A2 = (/-2,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInsideEdgeA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA2

    SUBROUTINE SphereLineIntersectionTestLineTangent
        implicit none
        real(dp), dimension(3) :: A1 = (/-1,0,0/)
        real(dp), dimension(3) :: A2 = (/1,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineTangent"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineTangent

    SUBROUTINE SphereLineIntersectionTestLineOutsideA1
        implicit none
        real(dp), dimension(3) :: A1 = (/1,0,0/)
        real(dp), dimension(3) :: A2 = (/2,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideA1

    SUBROUTINE SphereLineIntersectionTestLineOutsideA2
        implicit none
        real(dp), dimension(3) :: A1 = (/4,0,0/)
        real(dp), dimension(3) :: A2 = (/3,0,0/)
        real(dp), dimension(3) :: B1 = (/1,0,0/)
        real(dp) :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideA2

    SUBROUTINE SphereLineIntersectionTestLineOutsideBoth
        implicit none
        real(dp), dimension(3) :: A1 = (/5,0,0/)
        real(dp), dimension(3) :: A2 = (/-5,0,0/)
        real(dp), dimension(3) :: B1 = (/0,0,0/)
        real(dp) :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideBoth"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideBoth


    SUBROUTINE SphereLineIntersectionTestLineCloseA1
        implicit none
        real(dp), dimension(3) :: A1 = (/12,0,0/)
        real(dp), dimension(3) :: A2 = (/11.001,0.0,0.0/)
        real(dp), dimension(3) :: B1 = (/10,10,10/)
        real(dp) :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineCloseA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineCloseA1

    SUBROUTINE SphereLineIntersectionTestLineCloseA2
        implicit none
        real(dp), dimension(3) :: A1 = (/-55,0,0/)
        real(dp), dimension(3) :: A2 = (/-50.003,0.0,0.0/)
        real(dp), dimension(3) :: B1 = (/-10,-10,-10/)
        real(dp) :: r = 40
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineCloseA2

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
        logical SphereSphereIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereSphereIntersectionCalculation = .FALSE.

        ! see if radii overalp
        if ( sqrt(dot_product(A1-B1, A1-B1)) - (ra+rb) < dist ) then
            SphereSphereIntersectionCalculation = .TRUE.
            !print*, "sphere collision"
            return
        endif
        return 

    END FUNCTION SphereSphereIntersectionCalculation

    ! testing routines !
    SUBROUTINE SphereSphereIntersectionTestAinB
        implicit none
        real(dp), dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real(dp) :: ra = 2
        real(dp), dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real(dp) :: rb = 20.0
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestAinB"
            stop
        endif
   
    END SUBROUTINE SphereSphereIntersectionTestAinB

    SUBROUTINE SphereSphereIntersectionTestBinA
        implicit none
        real(dp), dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real(dp) :: ra = 20
        real(dp), dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real(dp) :: rb = 2
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestBinA"
            stop
        endif
   
    END SUBROUTINE SphereSphereIntersectionTestBinA

    SUBROUTINE SphereSphereIntersectionTestTangent
        implicit none
        real(dp), dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real(dp) :: ra = 1
        real(dp), dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real(dp) :: rb = 1
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestTangent"
            stop
        endif

    END SUBROUTINE SphereSphereIntersectionTestTangent

    SUBROUTINE SphereSphereIntersectionTestOverlap
        implicit none
        real(dp), dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real(dp) :: ra = 1
        real(dp), dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real(dp) :: rb = 1
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlap"
            stop 
        endif

    END SUBROUTINE SphereSphereIntersectionTestOverlap

    SUBROUTINE SphereSphereIntersectionTestNoOverlap
        implicit none
        real(dp), dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real(dp) :: ra = 1
        real(dp), dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real(dp) :: rb = 0.99
        logical val
   
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestNoOverlap"
            stop
        endif

    END SUBROUTINE SphereSphereIntersectionTestNoOverlap

END MODULE SphereSphereIntersection