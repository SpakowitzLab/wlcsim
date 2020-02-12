! npagane | risca lab | dec 2019 | sterics calculation file

! define line-line intersection module that contains calculation and test
! inspired from http://paulbourke.net/geometry/pointlineplane/
MODULE LineLineIntersection
    use precision, only: dp
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
        real(dp), parameter ::  neighbor = 1.0e-3 ! tolerance for nearest collision, set decently high
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
            !print*, "collision, point overlap"
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
                    print*, "collision, parallel overlap"
                    !print*, A1(1), ',',A1(2), ',',A1(3)
                    !print*, A2(1), ',',A2(2), ',',A2(3)
                    !print*, B1(1),',', B1(2), ',',B1(3)
                    !print*, B2(1), ',',B2(2), ',',B2(3)
                    LineLineIntersectionCalculation = 10
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
        ! check if dist == 0
        if ((sqrt(dot_product(pA-pB, pA-pB)) == 0)  .AND. (abs(muA) <= 1) .AND. (abs(muB) <= 1) &
          .AND. (abs(-1 - muA*muB) > neighbor)) then ! check for -1 and 1 pairs  
            print*, "collision, intersect"
            print*, muA, muB, sqrt(dot_product(pA-pB, pA-pB)), abs(-1 - muA*muB)
            print*, A1(1), ',',A1(2),',', A1(3)
            print*, A2(1), ',',A2(2), ',',A2(3)
            print*, B1(1),',', B1(2),',', B1(3)
            print*, B2(1), ',',B2(2), ',',B2(3)
            LineLineIntersectionCalculation = 3
            !print*, "stopping for now bc this should not happen"
            !stop
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
        integer val 
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then 
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then 
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        print*, 'here'
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
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
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideX"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideX

    SUBROUTINE LineLineIntersectionTestIntersectProjection
        implicit none
        real(dp), dimension(3) :: A1 = (/-8.8910375546638831E-002 ,  0.45626156257315231      ,   6.6668324441179081/)!(/0,0,0/)
        real(dp), dimension(3) :: A2 = (/5.5511151231257827E-017 ,  -4.4408920985006262E-016 ,   8.2500003278255480/)!(/5,5,0))
        real(dp), dimension(3) :: B1 = (/8.5784016210477236E-018 ,  -6.3001451707972334E-017 ,   9.9000003933906555/)!(/0,2,2/)
        real(dp), dimension(3) :: B2 = (/9.0729117231340004E-019 ,  -7.8930676690238801E-017 ,   11.550000458955765 /)!(/2,0,1/)
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val /= 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjection"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjection

END MODULE LineLineIntersection


! define sphere-line intersection module that contains calculation and test
! inspired from http://paulbourke.net/geometry/circlesphere/
MODULE SphereLineIntersection
    use precision, only: dp
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
            SphereLineIntersectionCalculation = 4
            !print*, "line in sphere"
            !print*, um, up
            !print*, A1(1), ',',A1(2), ',',A1(3)
            !print*, A2(1), ',',A2(2), ',',A2(3)
            !print*, B1(1),',', B1(2), ',',B1(3), ',',r
            return
        else if ((um >= 0 .AND. um <= 1) .OR. (up >= 0 .AND. up <= 1)) then 
            SphereLineIntersectionCalculation = 2
            !print*, "one or more line-sphere intersections"
            !print*, um, up
            !print*, A1(1), ',', A1(2), ',', A1(3)
            !print*, A2(1), ',', A2(2), ',', A2(3)
            !print*, B1(1), ',', B1(2), ',', B1(3), ',', r
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val == 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val /= 0) then
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
        integer val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val /= 0) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineCloseA2

END MODULE SphereLineIntersection

! define sphere-sphere intersection module that contains calculation and test
MODULE SphereSphereIntersection
    use precision, only: dp
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
            SphereSphereIntersectionCalculation = 30
            !print*, "sphere collision"
            !print*, A1(1), ',', A1(2), ',', A1(3), ',', ra
            !print*, B1(1), ',', B1(2), ',', B1(3), ',', rb
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
        integer val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val == 0) then
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
        integer val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val == 0) then
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
        integer val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val == 0) then
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
        integer val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val == 0) then
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
        integer val
   
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val /= 0) then
            print*, "FAILURE: failed LineLineIntersectionTestNoOverlap"
            stop
        endif

    END SUBROUTINE SphereSphereIntersectionTestNoOverlap

END MODULE SphereSphereIntersection