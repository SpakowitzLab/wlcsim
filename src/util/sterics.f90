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
        integer, parameter :: iteration = 7 ! dont know optimal number (MATLAB used 6)
        ! by upping iteration, we get better at detecting barely penetrating objects so 
        ! does not provide a huge advantage 
        real(dp), dimension(nVerts, 3), intent(in) :: s1, s2
        real(dp), dimension(3) :: a, b, c, d, aout, bout, cout, dout
        real(dp), dimension(3), parameter :: v = (/0.8, 0.5, 1.0/) ! dont know optimum
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
        real(dp), parameter :: angle = 2.0*pi/10.5

        incr = 2*pi/(s/2)

        ! construct material rotation matrix
        mtrx(:,1) = v
        mtrx(:,2) = cross(u,v)
        mtrx(:,3) = u

        ! determine if nucleosome or not
        if (wrap /= 1) then 
            center = [4.8455, -2.4445, 0.6694]
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
            pos = (pos2 + pos1) / 2.0
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

SUBROUTINE runtimeTest1()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12,3) :: s1
        real(dp), dimension(12,3) :: s2
        integer res

        s1(:,1) = (/4984.3036725433549,4983.6411351254119,4982.6813353163698, 4982.3840729252706 , &
                 4983.0466103432136 ,4984.0064101522557 , 4984.6237746823926 ,  4983.9612372644497, &
                  4983.0014374554075 ,4982.7041750643084 , 4983.3667124822514 , 4984.3265122912935 /)
        s1(:,2) = (/4979.9142882079868  , 4980.2729603220951 ,4980.3740877572545 ,4980.1165430783049, &
                 4979.7578709641966 , 4979.6567435290381, 4981.5198968332743, 4981.8785689473825, &
                 4981.9796963825411, 4981.7221517035923,4981.3634795894841 ,4981.2623521543246/)
        s1(:,3) = (/5024.4983543860117 ,5025.1559245224898 , 5024.8940895202049 ,5023.9746843814428 , &
                5023.3171142449646 ,5023.5789492472495 , 5023.9450941705545, 5024.6026643070327 ,&
                5024.3408293047478 ,5023.4214241659856  ,5022.7638540295075 ,5023.0256890317924/)
        s2(:,1) = (/ 4984.3682523821844 , 4984.0866957155158 , 4984.7369976322943,4985.6688562157415 , &
                4985.9504128824110 ,4985.3001109656325, 4983.8834913843875, 4983.6019347177180, &
                 4984.2522366344965 ,4985.1840952179446 , 4985.4656518846132 ,4984.8153499678347 /)
        s2(:,2) = (/ 4980.8313830740180 ,4981.7677996593266 , 4982.5204042315863 , 4982.3365922185376, &
                 4981.4001756332291 ,4980.6475710609693 , 4981.0353588576045 , 4981.9717754429130 , &
                 4982.7243800151728 , 4982.5405680021240, 4981.6041514168155,4980.8515468445557 /)
        s2(:,3) = (/5022.2498377347110,5022.0404342433039 ,5022.1438449528441,5022.4566591537905 ,&
                5022.6660626451976,5022.5626519356574 ,5023.8137753520941, 5023.6043718606870 ,&
                5023.7077825702272 ,5024.0205967711736 , 5024.2300002625807, 5024.1265895530405/)

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

        s2(:,1) = (/4984.3036725433549,4983.6411351254119,4982.6813353163698, 4982.3840729252706 , &
                 4983.0466103432136 ,4984.0064101522557 , 4984.6237746823926 ,  4983.9612372644497, &
                  4983.0014374554075 ,4982.7041750643084 , 4983.3667124822514 , 4984.3265122912935 /)
        s2(:,2) = (/4979.9142882079868  , 4980.2729603220951 ,4980.3740877572545 ,4980.1165430783049, &
                 4979.7578709641966 , 4979.6567435290381, 4981.5198968332743, 4981.8785689473825, &
                 4981.9796963825411, 4981.7221517035923,4981.3634795894841 ,4981.2623521543246/)
        s2(:,3) = (/5024.4983543860117 ,5025.1559245224898 , 5024.8940895202049 ,5023.9746843814428 , &
                5023.3171142449646 ,5023.5789492472495 , 5023.9450941705545, 5024.6026643070327 ,&
                5024.3408293047478 ,5023.4214241659856  ,5022.7638540295075 ,5023.0256890317924/)
        s1(:,1) = (/ 4984.3682523821844 , 4984.0866957155158 , 4984.7369976322943,4985.6688562157415 , &
                4985.9504128824110 ,4985.3001109656325, 4983.8834913843875, 4983.6019347177180, &
                 4984.2522366344965 ,4985.1840952179446 , 4985.4656518846132 ,4984.8153499678347 /)
        s1(:,2) = (/ 4980.8313830740180 ,4981.7677996593266 , 4982.5204042315863 , 4982.3365922185376, &
                 4981.4001756332291 ,4980.6475710609693 , 4981.0353588576045 , 4981.9717754429130 , &
                 4982.7243800151728 , 4982.5405680021240, 4981.6041514168155,4980.8515468445557 /)
        s1(:,3) = (/5022.2498377347110,5022.0404342433039 ,5022.1438449528441,5022.4566591537905 ,&
                5022.6660626451976,5022.5626519356574 ,5023.8137753520941, 5023.6043718606870 ,&
                5023.7077825702272 ,5024.0205967711736 , 5024.2300002625807, 5024.1265895530405/)

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

        s2(:,1) = (/ -24.367111445653041 ,-22.992940945926144,-18.346801795507311,-15.074833144815370,&
                    -16.449003644542266,-21.095142794961099 ,-26.548329711846002 , -25.174159212119108 ,&
                    -20.528020061700275 ,-17.256051411008333 ,-18.630221910735226 ,-23.276361061154063/)
        s2(:,2) = (/-9.9871384396511491 ,-11.600682933114870,-13.811763387143637,-14.409299347708682 ,&
                    -12.795754854244962 ,-10.584674400216194, -14.926349986361437 ,-16.539894479825158, &
                    -18.750974933853925 ,-19.348510894418972 ,-17.734966400955251 , -15.523885946926484/)
        s2(:,3) = (/29.088228971560572,33.836716057495778 ,34.588057671382089 , 30.590912199333196,&
                    25.842425113397990 ,25.091083499511676,28.041101389997458 , 32.789588475932668 ,&
                    33.540930089818978 ,29.543784617770083,24.795297531834876 ,24.043955917948562  /)
        s1(:,1) =  (/-22.771514420078759,-23.542767384627435 , -23.728152660276933 ,-23.142284971377759 ,&
                    -22.371032006829079,-22.185646731179581,-23.802635224192347 ,-24.573888188741027 ,&
                    -24.759273464390525, -24.173405775491350 ,-23.402152810942670 ,-23.216767535293172 /)
        s1(:,2) =  (/-19.532863971582120 ,-19.840405251130480 ,-19.404957802191745 ,-18.661969073704650 ,&
                    -18.354427794156294 ,-18.789875243095025 ,-18.376264807364731 ,-18.683806086913087,&
                    -18.248358637974352 ,-17.505369909487261 ,-17.197828629938900 ,-17.633276078877635 /)
        s1(:,3) = (/32.939834348204087, 32.382530816858534,31.501612081226977 ,31.177996876940977, &
                    31.735300408286530 ,32.616219143918087 , 33.728548155617396, 33.171244624271836, &
                    32.290325888640282, 31.966710684354283 , 32.524014215699836 ,33.404932951331389/)

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

        s1(:,1) = (/ -24.367111445653041 ,-22.992940945926144,-18.346801795507311,-15.074833144815370,&
                    -16.449003644542266,-21.095142794961099 ,-26.548329711846002 , -25.174159212119108 ,&
                    -20.528020061700275 ,-17.256051411008333 ,-18.630221910735226 ,-23.276361061154063/)
        s1(:,2) = (/-9.9871384396511491 ,-11.600682933114870,-13.811763387143637,-14.409299347708682 ,&
                    -12.795754854244962 ,-10.584674400216194, -14.926349986361437 ,-16.539894479825158, &
                    -18.750974933853925 ,-19.348510894418972 ,-17.734966400955251 , -15.523885946926484/)
        s1(:,3) = (/29.088228971560572,33.836716057495778 ,34.588057671382089 , 30.590912199333196,&
                    25.842425113397990 ,25.091083499511676,28.041101389997458 , 32.789588475932668 ,&
                    33.540930089818978 ,29.543784617770083,24.795297531834876 ,24.043955917948562  /)
        s2(:,1) =  (/-22.771514420078759,-23.542767384627435 , -23.728152660276933 ,-23.142284971377759 ,&
                    -22.371032006829079,-22.185646731179581,-23.802635224192347 ,-24.573888188741027 ,&
                    -24.759273464390525, -24.173405775491350 ,-23.402152810942670 ,-23.216767535293172 /)
        s2(:,2) =  (/-19.532863971582120 ,-19.840405251130480 ,-19.404957802191745 ,-18.661969073704650 ,&
                    -18.354427794156294 ,-18.789875243095025 ,-18.376264807364731 ,-18.683806086913087,&
                    -18.248358637974352 ,-17.505369909487261 ,-17.197828629938900 ,-17.633276078877635 /)
        s2(:,3) = (/32.939834348204087, 32.382530816858534,31.501612081226977 ,31.177996876940977, &
                    31.735300408286530 ,32.616219143918087 , 33.728548155617396, 33.171244624271836, &
                    32.290325888640282, 31.966710684354283 , 32.524014215699836 ,33.404932951331389/)

        res = GJK(s1, s2, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: runtimeTest4"
            stop
        endif

    END SUBROUTINE runtimeTest4
   

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