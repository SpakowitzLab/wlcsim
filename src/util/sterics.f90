! npagane | risca lab | feb 2020 | fortran implementation of gjk algorithm 
! adapted from MATLAB code https://github.com/mws262/MATLAB-GJK-Collision-Detection/blob/master/GJK.m

! define sphere-sphere intersection module that contains calculation and test
MODULE GJKAlgorithm
    use precision, only: dp, pi
    implicit none
    contains

    ! calculation
    FUNCTION GJK(s1, s2, nVerts)
        implicit none
        integer, intent(in) :: nVerts
        integer, parameter :: iteration = 6 ! dont know optimal number
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        real(dp), dimension(3) :: a, b, c, d, aout, bout, cout, dout
        real(dp), dimension(3), parameter :: v = (/0.8, 0.5, 1.0/)
        integer GJK

        ! default value (i.e. no intersection)
        GJK = 0

        ! check for shape overlap
        if ( sum(s2(:,1)-s1(:,1)) + sum(s2(:,2)-s1(:,2)) + sum(s2(:,3)-s1(:,3)) == 0 ) then 
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
        integer, intent(in) :: wrap ! num of basepairs wrapped
        real(dp), dimension(3), intent(in) :: u ! u angle of bead 1
        real(dp), dimension(3), intent(in) :: v ! v angle of bead 1
        integer, intent(in) :: s ! num sides of desired polygon
        real(dp) :: r, h ! radius and height of desired shape
        integer i
        real(dp), dimension(3,3) :: mtrx, rot
        real(dp) space, spaceInit, incr, offset1, offset2
        real(dp), dimension(3) :: pos, center, vec
        real(dp), dimension(s,3) :: constructPolygonPrism
        real(dp) :: angle = 2.0*pi/10.5

        incr = 2*pi/(s/2)

        ! construct material rotation matrix
        mtrx(:,1) = v
        mtrx(:,2) = cross(u,v)
        mtrx(:,3) = u

        ! determine if nucleosome or not
        if (wrap /= 1) then 
            center = [4.8455, -2.4445, 0.6694]
            pos = pos1
            h = 5.5 ! nm height
            r = 5.2 ! nm radius
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
            pos = (pos2 + pos1) / 2.0
            h = sqrt(dot_product(pos2 - pos1, pos2-pos1)) ! nm height
            r = 1.0 ! nm radius 
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

    SUBROUTINE sameShapeTest()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/ 4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972, & 
                    4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972/)
        s1(:,2) = (/2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75, &
                    2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75/)
        s1(:,3) = (/0.,0.,0.,0.,0.,0.,-3.,-3.,-3.,-3.,-3.,-3./)
        s2(:,1) = (/ 4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972, &
                    4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972/)
        s2(:,2) = (/2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75, &            
                    2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75/)
        s2(:,3) = (/0.,0.,0.,0.,0.,0.,-3.,-3.,-3.,-3.,-3.,-3./)
                            
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

        s1(:,1) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,1) = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                    4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        s2(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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

        s1(:,1) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,1) = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                    .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        s2(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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

        s1(:,1) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,1) = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                    1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        s2(:,2) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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

        s1(:,2) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,2) = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                    4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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

        s1(:,2) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,2) = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                    .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
 
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

        s1(:,2) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,2) = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                    1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,3) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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

        s1(:,3) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,3) = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                    4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

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
        
        s1(:,3) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,3) = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                    .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
 

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

        s1(:,3) = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                    -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        s1(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s1(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        s2(:,3) = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                    1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        s2(:,1) = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                    -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        s2(:,2) = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentZ"
            stop
        endif

    END SUBROUTINE tangentZ


END MODULE

! test module
! PROGRAM GJKTest 
!     use GJKAlgorithm, only: GJK, sameShapeTest, noIntersectX, intersectX, tangentX, &
!                             noIntersectY, intersectY, tangentY, &
!                             noIntersectZ, intersectZ, tangentZ
!     implicit none

!     call sameShapeTest()
!     call noIntersectX()
!     call intersectX()
!     call tangentX()
!     call noIntersectY()
!     call intersectY()
!     call tangentY()
!     call noIntersectZ()
!     call intersectZ()
!     call tangentZ()

!     print*, "SUCCESS: successful completion of all GJK collision unit tests"

! END PROGRAM