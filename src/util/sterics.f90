#include "../defines.inc"

! npagane | risca lab | feb 2020 | fortran implementation of gjk algorithm 
! adapted from MATLAB code https://github.com/mws262/MATLAB-GJK-Collision-Detection/blob/master/GJK.m

! define GJK intersection module that contains calculation and test
MODULE GJKAlgorithm
    use precision, only: pi, dp
    implicit none
    contains

    ! calculation
    FUNCTION GJK(s1, s2, nVerts)
        implicit none
        integer, intent(in) :: nVerts
        integer, parameter :: iteration = WLC_P__GJK_POLYGON ! keep high, this MF was part of the bug (thxUS)
        ! by upping iteration, we get better at detecting barely penetrating objects so 
        ! does not provide a huge advantage 
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        real(dp), dimension(3) :: a, b, c, d, aout, bout, cout, dout
        real(dp), dimension(3), parameter :: v = (/0.8, 0.5, 1.0/) ! random yet determined
        integer GJK

        ! default value (i.e. no intersection)
        GJK = 0

        ! check for shape overlap
        if ( all(sum(s2-s1,1) == 0) ) then 
           GJK = 1!iteration
           return
        endif

        ! line segment
        call pickLine(s1, s2, nVerts, v, a, b)

        ! triangle
        call pickTriangle(s1, s2, nVerts, a, b, &
                            aout, bout, cout, GJK, iteration)
        a = aout
        b = bout
        c = cout

        ! tetrahedron
        if (GJK > 0) then 
            call pickTetrahedron(s1, s2, nVerts, a, b, c, &
                              aout, bout, cout, dout, GJK, iteration)
        endif

    END FUNCTION GJK

    ! pickLine
    SUBROUTINE pickLine(s1, s2, nVerts, v, a, b)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        real(dp), dimension(3), intent(in) :: v 
        real(dp), dimension(3), intent(out) :: a, b
        
        ! make first line of simplex
        a = support(s2, s1, nVerts, v)
        b = support(s2, s1, nVerts, -1*v)

    END SUBROUTINE pickLine

    SUBROUTINE pickTriangle(s1, s2, nVerts, a, b, &
                            aout, bout, cout, flag, iteration)
        use vector_utils, only: cross
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        integer, intent(out) :: flag ! success
        real(dp), dimension(3), intent(out) :: aout, bout, cout
        real(dp), dimension(3) :: ab, ao, ac, abc, abp, acp, v
        integer i

        ! default value (i.e. no success)
        flag = 0
        aout = a
        bout = b

        ! first try
        ab = bout - aout
        ao = -1*aout
        v = cross(cross(ab, ao), ab) ! v is perp to ab and pointing towards origin
        cout = bout
        bout = aout
        aout = support(s2, s1, nVerts, v)

        ! iterate until convergence
        do i = 1,iteration
            ! check if found
            ab = bout - aout
            ao = -1*aout
            ac = cout - aout
            ! find normal to face
            abc = cross(ab, ac)
            ! perp to ab and ac away from triangle
            abp = cross(ab, abc)
            acp = cross(abc, ac)
            ! check if triangle contains origin
            if (dot_product(abp, ao) > 0 ) then 
                cout = bout 
                bout = aout
                v = abp
            else if (dot_product(acp, ao) > 0 ) then 
                bout = aout 
                v = acp
            else 
                flag = 1
                exit
            endif
            aout = support(s2, s1, nVerts, v)
        enddo

    END SUBROUTINE pickTriangle

    SUBROUTINE pickTetrahedron(s1, s2, nVerts, a, b, c, &
                              aout, bout, cout, dout, flag, iteration)
        use vector_utils, only: cross
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real(dp), dimension(3), intent(in) :: a, b, c
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        integer, intent(out) :: flag ! success
        real(dp), dimension(3), intent(out) :: aout, bout, cout, dout
        real(dp), dimension(3) :: ab, ao, ac, ad, abc, acd, adb, v
        real(dp) dot
        integer i

        ! default value (i.e. no success)
        flag = 0
        aout = a
        bout = b
        cout = c

        ! first try
        ab = bout - aout
        ac = cout - aout
        abc = cross(ab,ac)
        ao = -1*aout

        ! check if simplex is above or below origin
        dot = dot_product(abc, ao)
        if (dot > 0) then 
            dout = cout
            cout = bout
            bout = aout
            v = abc
            aout = support(s2, s1, nVerts, v)
        else
            dout = bout
            bout = aout
            v = -1*abc
            aout = support(s2, s1, nVerts, v)
        endif

        ! iterate until convergence
        do i = 1, iteration
            ab = bout - aout
            ao = -1*aout
            ac = cout - aout
            ad = dout - aout
            abc = cross(ab, ac)
            if (dot_product(abc, ao) <= 0) then
                acd = cross(ac, ad)
                if (dot_product(acd, ao) > 0) then 
                    bout = cout
                    cout = dout
                    ab = ac
                    ac = ad
                    abc = acd
                else if (dot_product(acd, ao) < 0) then 
                    adb = cross(ad, ab)
                    !print*, 'GJK', dot_product(adb, ao), dot_product(adb, ao) > 0
                    if (dot_product(adb, ao) > 0) then 
                        cout = bout
                        bout = dout
                        ac = ab
                        ab = ad
                        abc = adb
                    else
                        flag = 1!i
                        exit
                    endif
                endif
            endif
            ! try again
            if (dot_product(abc, ao) > 0) then
                dout = cout
                cout = bout
                bout = aout
                v = abc
                aout = support(s2, s1, nVerts, v)
            else
                dout = bout
                bout = aout
                v = -1*abc
                aout = support(s2, s1, nVerts, v)
            endif
        enddo

    END SUBROUTINE pickTetrahedron

    FUNCTION getExtremaPoint(s, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), dimension(nVerts, 3), intent(in) :: s
        real(dp), dimension(3), intent(in) :: v
        real(dp), dimension(3) :: getExtremaPoint
        real(dp) :: mag(nVerts)
        integer maxInd(1)

        ! find the furthest data point in v direction in shape s
        mag = s(:,1)*v(1) + s(:,2)*v(2) + s(:,3)*v(3)
        maxInd = maxloc(mag)
        getExtremaPoint = (/s(maxInd(1),1), s(maxInd(1),2), s(maxInd(1),3)/)

    END FUNCTION getExtremaPoint

    FUNCTION support(s2, s1, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        real(dp), dimension(3), intent(in) :: v
        real(dp), dimension(3) :: point1, point2, support

        ! support function for minkowski difference
        point1 = getExtremaPoint(s1, nVerts, v)
        point2 = getExtremaPoint(s2, nVerts, -1*v)
        support =  point1-point2

    END FUNCTION support

    FUNCTION constructPolygonPrism(pos1, pos2, wrap, u, v, s)
        use vector_utils, only: cross
        implicit none
        real(dp), dimension(3), intent(in) :: pos1 ! first bead position
        real(dp), dimension(3), intent(in) :: pos2 ! second bead position
        real(dp), intent(in) :: wrap ! num of basepairs wrapped
        real(dp), dimension(3), intent(in) :: u ! u angle of bead 1
        real(dp), dimension(3), intent(in) :: v ! v angle of bead 1
        integer, intent(in) :: s ! num sides of desired polygon
        real(dp) :: r, h ! radius and height of desired shape
        integer i
        real(dp), dimension(3,3) :: mtrx
        real(dp) space, spaceInit, incr, offset1, offset2
        real(dp), dimension(3) :: pos, center, vec
        real(dp), dimension(s,3) :: constructPolygonPrism

        incr = 2*pi/(s/2)

        ! construct material rotation matrix
        mtrx(:,1) = v
        mtrx(:,2) = cross(u,v)
        mtrx(:,3) = u

        ! determine if nucleosome or not
        if (wrap /= 0) then 
            center = [4.84550_DP, -2.44450_DP, 0.66940_DP]
            pos = pos1
            h = WLC_P__NUCLEOSOME_HEIGHT ! nm height
            r = WLC_P__NUCLEOSOME_RADIUS ! nm radius
            spaceInit = 0
            offset1 = -h/2
            offset2 = h/2
            ! create parametric t for first face
            space = spaceInit
            do i = 1, (s/2)
                ! rotate into material frame
                vec = [r*cos(space), offset1, r*sin(space)]
                ! construct polygon
                constructPolygonPrism(i,:) = pos + MATMUL(mtrx, center+vec)
                space = space + incr
            enddo
            ! create parametric t for second face
            space = spaceInit
            do i = (s/2)+1, s
                ! rotate into material frame
                vec = [r*cos(space), offset2, r*sin(space)]
                ! construct polygon
                constructPolygonPrism(i,:) = pos + MATMUL(mtrx, center+vec)
                space = space + incr
            enddo
        else ! dna 
            center = 0.0_dp
            pos = (pos2 + pos1) / 2.00_DP
            h = sqrt(dot_product(pos2 - pos1, pos2-pos1)) ! nm height
            r = WLC_P__DNA_RADIUS ! nm radius
            spaceInit = 2*pi/s
            offset1 = -h/2
            offset2 = h/2
            ! create parametric t for first face
            space = spaceInit
            do i = 1, (s/2)
                ! rotate into material frame
                vec = [r*sin(space), r*cos(space), offset1]
                ! construct polygon
                constructPolygonPrism(i,:) = pos + MATMUL(mtrx, center+vec)
                space = space + incr
            enddo
            ! create parametric t for second face
            space = spaceInit
            do i = (s/2)+1, s
                ! rotate into material frame
                vec = [r*sin(space), r*cos(space), offset2]
                ! construct polygon
                constructPolygonPrism(i,:) = pos + MATMUL(mtrx, center+vec)
                space = space + incr
            enddo
        endif 

    END FUNCTION constructPolygonPrism

    FUNCTION findCenterPolygonPrism(pos1, pos2, wrap, u, v)
        use vector_utils, only: cross
        implicit none
        real(dp), dimension(3), intent(in) :: pos1 ! first bead position
        real(dp), dimension(3), intent(in) :: pos2 ! second bead position
        real(dp), intent(in) :: wrap ! num of basepairs wrapped
        real(dp), dimension(3), intent(in) :: u ! u angle of bead 1
        real(dp), dimension(3), intent(in) :: v ! v angle of bead 1
        real(dp), dimension(3,3) :: mtrx
        real(dp), dimension(3) :: pos, center
        real(dp), dimension(3) :: findCenterPolygonPrism

        ! determine if nucleosome or not
        if (wrap /= 0) then 
            ! construct material rotation matrix
            mtrx(:,1) = v
            mtrx(:,2) = cross(u,v)
            mtrx(:,3) = u
            center = [4.84550_DP, -2.44450_DP, 0.66940_DP]
            pos = pos1
            ! find center of polygon
            findCenterPolygonPrism = pos + MATMUL(mtrx, center)
        else ! dna 
            pos = (pos2 + pos1) / 2.00_DP
            ! find center of polygon
            findCenterPolygonPrism = pos
        endif 

    END FUNCTION findCenterPolygonPrism


    SUBROUTINE sameShapeTest()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/ 4.76313972_DP,  0.0_DP , -4.76313972_DP, -4.76313972_DP,  0.0_DP ,4.76313972_DP, & 
                    4.76313972_DP,  0.0_DP , -4.76313972_DP, -4.76313972_DP,  0.0_DP ,4.76313972_DP/)
        s1(:,2) = (/2.750_DP,  5.5_DP,  2.750_DP, -2.750_DP, -5.5_DP,-2.750_DP, &
                    2.750_DP,  5.5_DP,  2.750_DP, -2.750_DP, -5.5_DP,-2.750_DP/)
        s1(:,3) = (/0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP/)
        s2(:,1) = (/ 4.76313972_DP,  0.0_DP , -4.76313972_DP, -4.76313972_DP,  0.0_DP ,4.76313972_DP, &
                    4.76313972_DP,  0.0_DP , -4.76313972_DP, -4.76313972_DP,  0.0_DP ,4.76313972_DP/)
        s2(:,2) = (/2.750_DP,  5.5_DP,  2.750_DP, -2.750_DP, -5.5_DP,-2.750_DP, &            
                    2.750_DP,  5.5_DP,  2.750_DP, -2.750_DP, -5.5_DP,-2.750_DP/)
        s2(:,3) = (/0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP,-3.0_DP/)
                            
        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: sameShapeTest"
            !stop
        endif

    END SUBROUTINE sameShapeTest

    SUBROUTINE noIntersectX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,1) = (/ 4.500_DP,    4.000_DP,    4.500_DP,    5.500_DP,    6.000_DP,    5.500_DP, &
                    4.500_DP,    4.000_DP,    4.500_DP,    5.500_DP,    6.000_DP,    5.500_DP/)
        s2(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectX"
            stop
        endif

    END SUBROUTINE noIntersectX

    SUBROUTINE intersectX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,1) = (/ 0.500_DP,    0.00_DP,    0.500_DP,    1.500_DP,    2.000_DP,    1.500_DP, &
                    0.500_DP,    0.00_DP,    0.500_DP,    1.500_DP,    2.000_DP,    1.500_DP/)
        s2(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectX"
            stop
        endif

    END SUBROUTINE intersectX

    SUBROUTINE tangentX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,1) = (/ 1.500_DP,    1.00_DP,    1.500_DP,    2.500_DP,    3.000_DP,    2.500_DP, &
                    1.500_DP,    1.00_DP,    1.500_DP,    2.500_DP,    3.000_DP,    2.500_DP/)
        s2(:,2) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentX"
            stop
        endif

    END SUBROUTINE tangentX

    SUBROUTINE noIntersectY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,2) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,2) = (/ 4.5000_DP,    4.000_DP,    4.500_DP,    5.500_DP,    6.000_DP,    5.500_DP, &
                    4.5000_DP,    4.000_DP,    4.500_DP,    5.500_DP,    6.000_DP,    5.500_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectY"
            stop
        endif

    END SUBROUTINE noIntersectY

    SUBROUTINE intersectY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,2) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,2) = (/ 0.500_DP,    0.00_DP,    0.500_DP,    1.500_DP,    2.000_DP,    1.500_DP, &
                    0.500_DP,    0.00_DP,    0.500_DP,    1.500_DP,    2.000_DP,    1.500_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
 
        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectY"
            stop
        endif

    END SUBROUTINE intersectY
 
    SUBROUTINE tangentY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,2) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,2) = (/ 1.500_DP,    1.00_DP,    1.500_DP,    2.500_DP,    3.000_DP,    2.500_DP, &
                    1.500_DP,    1.00_DP,    1.500_DP,    2.500_DP,    3.000_DP,    2.500_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,3) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentY"
            stop
        endif

    END SUBROUTINE tangentY

    SUBROUTINE noIntersectZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,3) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,3) = (/ 4.5000_DP,    4.0000_DP,    4.5000_DP,    5.5000_DP,    6.0000_DP,    5.5000_DP, &
                    4.5000_DP,    4.0000_DP,    4.5000_DP,    5.5000_DP,    6.0000_DP,    5.5000_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectZ"
            stop
        endif

    END SUBROUTINE noIntersectZ

    SUBROUTINE intersectZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res
        
        s1(:,3) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,3) = (/ 0.5000_DP,    0.000_DP,    0.5000_DP,    1.5000_DP,    2.0_DP,    1.5000_DP, &
                    0.5000_DP,    0.000_DP,    0.5000_DP,    1.5000_DP,    2.0_DP,    1.5000_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
 

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectZ"
            stop
        endif

    END SUBROUTINE intersectZ

    SUBROUTINE tangentZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,3) = (/-0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP, &
                    -0.5000_DP,   -1.0000_DP,   -0.5000_DP,    0.5000_DP,    1.0000_DP,    0.5000_DP/)
        s1(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s1(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)
        s2(:,3) = (/ 1.5000_DP,    1.0000_DP,    1.5000_DP,    2.5000_DP,    3.0000_DP,    2.5000_DP, &
                    1.5000_DP,    1.0000_DP,    1.5000_DP,    2.5000_DP,    3.0000_DP,    2.5000_DP/)
        s2(:,1) = (/-0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP, &
                    -0.8660_DP,         0.0_DP,    0.8660_DP,    0.8660_DP,         0.0_DP,   -0.8660_DP/)
        s2(:,2) = (/0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP, 5.5_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentZ"
            stop
        endif

    END SUBROUTINE tangentZ

SUBROUTINE runtimeTest1()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/4984.30367254335490_DP,4983.64113512541190_DP,4982.68133531636980_DP, 4982.38407292527060_DP, &
                 4983.04661034321360_DP,4984.00641015225570_DP, 4984.62377468239260_DP,  4983.96123726444970_DP, &
                  4983.00143745540750_DP,4982.70417506430840_DP, 4983.36671248225140_DP, 4984.32651229129350_DP /)
        s1(:,2) = (/4979.91428820798680_DP, 4980.27296032209510_DP,4980.37408775725450_DP,4980.11654307830490_DP, &
                 4979.75787096419660_DP, 4979.65674352903810_DP, 4981.51989683327430_DP, 4981.87856894738250_DP, &
                 4981.97969638254110_DP, 4981.72215170359230_DP,4981.36347958948410_DP,4981.26235215432460_DP/)
        s1(:,3) = (/5024.49835438601170_DP,5025.15592452248980_DP, 5024.89408952020490_DP,5023.97468438144280_DP, &
                5023.31711424496460_DP,5023.57894924724950_DP, 5023.94509417055450_DP, 5024.60266430703270_DP,&
                5024.34082930474780_DP,5023.42142416598560_DP,5022.76385402950750_DP,5023.02568903179240_DP/)
        s2(:,1) = (/ 4984.36825238218440_DP, 4984.08669571551580_DP, 4984.73699763229430_DP,4985.66885621574150_DP, &
                4985.95041288241100_DP,4985.30011096563250_DP, 4983.88349138438750_DP, 4983.60193471771800_DP, &
                 4984.25223663449650_DP,4985.18409521794460_DP, 4985.46565188461320_DP,4984.81534996783470_DP/)
        s2(:,2) = (/ 4980.83138307401800_DP,4981.76779965932660_DP, 4982.52040423158630_DP, 4982.33659221853760_DP, &
                 4981.40017563322910_DP,4980.64757106096930_DP, 4981.03535885760450_DP, 4981.97177544291300_DP, &
                 4982.72438001517280_DP, 4982.54056800212400_DP, 4981.60415141681550_DP,4980.85154684455570_DP/)
        s2(:,3) = (/5022.24983773471100_DP,5022.04043424330390_DP,5022.14384495284410_DP,5022.45665915379050_DP,&
                5022.66606264519760_DP,5022.56265193565740_DP,5023.81377535209410_DP, 5023.60437186068700_DP,&
                5023.70778257022720_DP,5024.02059677117360_DP, 5024.23000026258070_DP, 5024.12658955304050_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest1"
            stop
        endif

    END SUBROUTINE runtimeTest1

SUBROUTINE runtimeTest2()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s2(:,1) = (/4984.30367254335490_DP,4983.64113512541190_DP,4982.68133531636980_DP, 4982.38407292527060_DP, &
                 4983.04661034321360_DP,4984.00641015225570_DP, 4984.62377468239260_DP,  4983.96123726444970_DP, &
                  4983.00143745540750_DP,4982.70417506430840_DP, 4983.36671248225140_DP, 4984.32651229129350_DP /)
        s2(:,2) = (/4979.91428820798680_DP, 4980.27296032209510_DP,4980.37408775725450_DP,4980.11654307830490_DP, &
                 4979.75787096419660_DP, 4979.65674352903810_DP, 4981.51989683327430_DP, 4981.87856894738250_DP, &
                 4981.97969638254110_DP, 4981.72215170359230_DP, 4981.36347958948410_DP,4981.26235215432460_DP/)
        s2(:,3) = (/5024.49835438601170_DP,5025.15592452248980_DP, 5024.89408952020490_DP,5023.97468438144280_DP, &
                5023.31711424496460_DP,5023.57894924724950_DP, 5023.94509417055450_DP, 5024.60266430703270_DP,&
                5024.34082930474780_DP,5023.42142416598560_DP,5022.76385402950750_DP,5023.02568903179240_DP/)
        s1(:,1) = (/ 4984.36825238218440_DP, 4984.08669571551580_DP, 4984.73699763229430_DP,4985.66885621574150_DP , &
                4985.95041288241100_DP,4985.30011096563250_DP, 4983.88349138438750_DP, 4983.60193471771800_DP, &
                 4984.25223663449650_DP,4985.18409521794460_DP, 4985.46565188461320_DP,4984.81534996783470_DP /)
        s1(:,2) = (/ 4980.83138307401800_DP ,4981.76779965932660_DP, 4982.52040423158630_DP, 4982.33659221853760_DP, &
                 4981.40017563322910_DP,4980.64757106096930_DP, 4981.03535885760450_DP, 4981.97177544291300_DP, &
                 4982.72438001517280_DP, 4982.54056800212400_DP, 4981.60415141681550_DP,4980.85154684455570_DP/)
        s1(:,3) = (/5022.24983773471100_DP,5022.04043424330390_DP,5022.14384495284410_DP,5022.45665915379050_DP,&
                5022.66606264519760_DP,5022.56265193565740_DP,5023.81377535209410_DP, 5023.60437186068700_DP,&
                5023.70778257022720_DP,5024.02059677117360_DP, 5024.23000026258070_DP, 5024.12658955304050_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest2"
            stop
        endif

    END SUBROUTINE runtimeTest2

    SUBROUTINE runtimeTest3()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s2(:,1) = (/ -24.3671114456530410_DP,-22.9929409459261440_DP,-18.3468017955073110_DP,-15.0748331448153700_DP,&
                    -16.4490036445422660_DP,-21.0951427949610990_DP,-26.5483297118460020_DP, -25.1741592121191080_DP,&
                    -20.5280200617002750_DP,-17.2560514110083330_DP,-18.6302219107352260_DP,-23.2763610611540630_DP/)
        s2(:,2) = (/-9.98713843965114910_DP,-11.6006829331148700_DP,-13.8117633871436370_DP,-14.4092993477086820_DP,&
                    -12.7957548542449620_DP,-10.5846744002161940_DP, -14.9263499863614370_DP,-16.5398944798251580_DP, &
                    -18.7509749338539250_DP,-19.3485108944189720_DP,-17.7349664009552510_DP, -15.5238859469264840_DP/)
        s2(:,3) = (/29.0882289715605720_DP,33.8367160574957780_DP,34.5880576713820890_DP, 30.5909121993331960_DP,&
                    25.8424251133979900_DP,25.0910834995116760_DP,28.0411013899974580_DP, 32.7895884759326680_DP,&
                    33.5409300898189780_DP,29.5437846177700830_DP,24.7952975318348760_DP,24.0439559179485620_DP  /)
        s1(:,1) =  (/-22.7715144200787590_DP,-23.5427673846274350_DP, -23.7281526602769330_DP,-23.1422849713777590_DP ,&
                    -22.3710320068290790_DP,-22.1856467311795810_DP,-23.8026352241923470_DP,-24.5738881887410270_DP ,&
                    -24.7592734643905250_DP, -24.1734057754913500_DP,-23.4021528109426700_DP,-23.2167675352931720_DP /)
        s1(:,2) =  (/-19.5328639715821200_DP,-19.8404052511304800_DP,-19.4049578021917450_DP,-18.6619690737046500_DP ,&
                    -18.3544277941562940_DP,-18.7898752430950250_DP,-18.3762648073647310_DP,-18.6838060869130870_DP,&
                    -18.2483586379743520_DP,-17.5053699094872610_DP,-17.1978286299389000_DP,-17.6332760788776350_DP /)
        s1(:,3) = (/32.9398343482040870_DP, 32.3825308168585340_DP,31.5016120812269770_DP,31.1779968769409770_DP, &
                    31.7353004082865300_DP,32.6162191439180870_DP, 33.7285481556173960_DP, 33.1712446242718360_DP, &
                    32.2903258886402820_DP, 31.9667106843542830_DP, 32.5240142156998360_DP,33.4049329513313890_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest3"
            stop
        endif

    END SUBROUTINE runtimeTest3

    SUBROUTINE runtimeTest4()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/ -24.3671114456530410_DP,-22.9929409459261440_DP,-18.3468017955073110_DP,-15.0748331448153700_DP,&
                    -16.4490036445422660_DP,-21.0951427949610990_DP,-26.5483297118460020_DP, -25.1741592121191080_DP,&
                    -20.5280200617002750_DP,-17.2560514110083330_DP,-18.6302219107352260_DP,-23.2763610611540630_DP/)
        s1(:,2) = (/-9.98713843965114910_DP,-11.6006829331148700_DP,-13.8117633871436370_DP,-14.4092993477086820_DP,&
                    -12.7957548542449620_DP,-10.5846744002161940_DP, -14.9263499863614370_DP,-16.5398944798251580_DP, &
                    -18.7509749338539250_DP,-19.3485108944189720_DP,-17.7349664009552510_DP, -15.5238859469264840_DP/)
        s1(:,3) = (/29.0882289715605720_DP,33.8367160574957780_DP,34.5880576713820890_DP, 30.5909121993331960_DP,&
                    25.8424251133979900_DP,25.0910834995116760_DP,28.0411013899974580_DP, 32.7895884759326680_DP,&
                    33.5409300898189780_DP,29.5437846177700830_DP,24.7952975318348760_DP,24.0439559179485620_DP  /)
        s2(:,1) =  (/-22.7715144200787590_DP,-23.5427673846274350_DP, -23.7281526602769330_DP,-23.1422849713777590_DP ,&
                    -22.3710320068290790_DP,-22.1856467311795810_DP,-23.8026352241923470_DP,-24.5738881887410270_DP ,&
                    -24.7592734643905250_DP, -24.1734057754913500_DP,-23.4021528109426700_DP,-23.2167675352931720_DP /)
        s2(:,2) =  (/-19.5328639715821200_DP,-19.8404052511304800_DP,-19.4049578021917450_DP,-18.6619690737046500_DP ,&
                    -18.3544277941562940_DP,-18.7898752430950250_DP,-18.3762648073647310_DP,-18.6838060869130870_DP,&
                    -18.2483586379743520_DP,-17.5053699094872610_DP,-17.1978286299389000_DP,-17.6332760788776350_DP /)
        s2(:,3) = (/32.9398343482040870_DP, 32.3825308168585340_DP,31.5016120812269770_DP,31.1779968769409770_DP, &
                    31.7353004082865300_DP,32.6162191439180870_DP, 33.7285481556173960_DP, 33.1712446242718360_DP, &
                    32.2903258886402820_DP, 31.9667106843542830_DP, 32.5240142156998360_DP,33.4049329513313890_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest4"
            stop
        endif

    END SUBROUTINE runtimeTest4

    SUBROUTINE runtimeTest5()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/ -14.33210_DP,  -19.03230_DP,  -19.46260_DP,  -15.19270_DP, -10.49240_DP,  -10.06210_DP,&
                 -14.13830_DP, -18.83850_DP,  -19.26880_DP, -14.99890_DP, -10.29870_DP, -9.86830_DP/)
        s1(:,2) = (/ -20.35950_DP, -20.20520_DP  ,-19.47660_DP , -18.90230_DP , -19.05660_DP , -19.78520_DP,&
                  -25.79920_DP,-25.64490_DP  ,-24.91630_DP  ,-24.34200_DP , -24.49630_DP,  -25.22490_DP/)
        s1(:,3) = (/7.17300_DP,  9.39200_DP  , 14.52270_DP ,  17.43440_DP ,  15.21530_DP  , 10.08470_DP  , &
                  7.96170_DP, 10.18070_DP ,  15.31140_DP ,  18.22310_DP ,  16.00410_DP  , 10.87340_DP/)
        s2(:,1) =  (/-10.24010_DP, -15.12380_DP  ,-18.02480_DP , -16.04220_DP , -11.15850_DP ,  -8.25750_DP , &
                     -12.04400_DP,-16.92760_DP , -19.82870_DP , -17.84600_DP  ,-12.96230_DP , -10.06130_DP/)
        s2(:,2) =  (/-27.02840_DP,-25.32850_DP , -24.35540_DP , -25.08240_DP , -26.78230_DP,  -27.75540_DP,&
                    -32.22400_DP, -30.52410_DP  ,-29.55100_DP , -30.27800_DP , -31.97790_DP , -32.95100_DP/)
        s2(:,3) = (/14.26150_DP, 13.71400_DP,  17.91850_DP ,  22.67040_DP  , 23.21780_DP ,  19.01340_DP  , &
                     14.21930_DP, 13.67180_DP, 17.87630_DP ,  22.62820_DP ,  23.17560_DP ,  18.97120_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest5"
            stop
        endif

    END SUBROUTINE runtimeTest5
   
    SUBROUTINE runtimeTest6()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s2(:,1) = (/ -14.33210_DP,  -19.03230_DP,  -19.46260_DP,  -15.19270_DP, -10.49240_DP,  -10.06210_DP,&
                 -14.13830_DP, -18.83850_DP,  -19.26880_DP, -14.99890_DP, -10.29870_DP, -9.86830_DP/)
        s2(:,2) = (/ -20.35950_DP, -20.20520_DP  ,-19.47660_DP , -18.90230_DP , -19.05660_DP , -19.78520_DP,&
                  -25.79920_DP,-25.64490_DP  ,-24.91630_DP  ,-24.34200_DP , -24.49630_DP,  -25.22490_DP/)
        s2(:,3) = (/7.17300_DP,  9.39200_DP  , 14.52270_DP ,  17.43440_DP ,  15.21530_DP  , 10.08470_DP  , &
                  7.96170_DP, 10.18070_DP ,  15.31140_DP ,  18.22310_DP ,  16.00410_DP  , 10.87340_DP/)
        s1(:,1) =  (/-10.24010_DP, -15.12380_DP  ,-18.02480_DP , -16.04220_DP , -11.15850_DP ,  -8.25750_DP , &
                     -12.04400_DP,-16.92760_DP , -19.82870_DP , -17.84600_DP  ,-12.96230_DP , -10.06130_DP/)
        s1(:,2) =  (/-27.02840_DP,-25.32850_DP , -24.35540_DP , -25.08240_DP , -26.78230_DP,  -27.75540_DP,&
                    -32.22400_DP, -30.52410_DP  ,-29.55100_DP , -30.27800_DP , -31.97790_DP , -32.95100_DP/)
        s1(:,3) = (/14.26150_DP, 13.71400_DP,  17.91850_DP ,  22.67040_DP  , 23.21780_DP ,  19.01340_DP  , &
                     14.21930_DP, 13.67180_DP, 17.87630_DP ,  22.62820_DP ,  23.17560_DP ,  18.97120_DP/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest6"
            stop
        endif

    END SUBROUTINE runtimeTest6

END MODULE