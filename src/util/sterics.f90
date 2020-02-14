! npagane | risca lab | feb 2020 | fortran implementation of gjk algorithm 
! adapted from MATLAB code https://github.com/mws262/MATLAB-GJK-Collision-Detection/blob/master/GJK.m

! define sphere-sphere intersection module that contains calculation and test
MODULE GJKAlgorithm
    implicit none
    contains

    ! calculation
    FUNCTION GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        implicit none
        integer, intent(in) :: nVerts
        integer, parameter :: iteration = 6
        real(dp), intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real(dp), dimension(3) :: a, b, c, d, aout, bout, cout, dout
        real(dp), dimension(3), parameter :: v = (/0.8, 0.5, 1.0/)
        integer GJK

        ! default value (i.e. no intersection)
        GJK = 0

        ! check for shape overlap
        if ( sum(s2x-s1x) + sum(s2y-s1y) + sum(s2z-s1z) == 0 ) then 
           GJK = iteration
           return
        endif

        ! line segment
        call pickLine(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, v, a, b)

        ! triangle
        call pickTriangle(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, &
                            aout, bout, cout, GJK, iteration)
        a = aout
        b = bout
        c = cout

        ! tetrahedron
        if (GJK > 0) then 
            call pickTetrahedron(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, c, &
                              aout, bout, cout, dout, GJK, iteration)
        endif

    END FUNCTION GJK

    ! pickLine
    SUBROUTINE pickLine(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, v, a, b)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real(dp), dimension(3), intent(in) :: v 
        real(dp), dimension(3), intent(out) :: a, b
        
        ! make first line of simplex
        a = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        b = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, -1*v)

    END SUBROUTINE pickLine

    SUBROUTINE pickTriangle(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, &
                            aout, bout, cout, flag, iteration)
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
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
        aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)

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
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        enddo

    END SUBROUTINE pickTriangle

    SUBROUTINE pickTetrahedron(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, c, &
                              aout, bout, cout, dout, flag, iteration)
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real(dp), dimension(3), intent(in) :: a, b, c
        real(dp), intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
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
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        else
            dout = bout
            bout = aout
            v = -1*abc
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
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
                        flag = i
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
                aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
            else
                dout = bout
                bout = aout
                v = -1*abc
                aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
            endif
        enddo

    END SUBROUTINE pickTetrahedron

    FUNCTION getExtremaPoint(sx, sy, sz, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), intent(in) :: sx(nVerts), sy(nVerts), sz(nVerts)
        real(dp), dimension(3), intent(in) :: v
        real(dp), dimension(3) :: getExtremaPoint
        real(dp) :: mag(nVerts)
        integer maxInd(1)

        ! find the furthest data point in v direction in shape s
        mag = sx*v(1) + sy*v(2) + sz*v(3)
        maxInd = maxloc(mag)
        getExtremaPoint = (/sx(maxInd(1)), sy(maxInd(1)), sz(maxInd(1))/)

    END FUNCTION getExtremaPoint

    FUNCTION support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real(dp), intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real(dp), dimension(3), intent(in) :: v
        real(dp), dimension(3) :: point1, point2, support

        ! support function for minkowski difference
        point1 = getExtremaPoint(s1x, s1y, s1z, nVerts, v)
        point2 = getExtremaPoint(s2x, s2y, s2z, nVerts, -1*v)
        support =  point1-point2

    END FUNCTION support

    FUNCTION cross(a, b)
        real(dp), dimension(3) :: cross
        real(dp), dimension(3), intent(in) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    END FUNCTION cross

    FUNCTION constructPolygonPrism(cen, u, v, r, h, s, wrap)
        implicit none
        real(dp), dimension(3), intent(in) :: cen ! center position
        real(dp), dimension(3), intent(in) :: u ! u angle
        real(dp), dimension(3), intent(in) :: v ! v angle
        real(dp), intent(in) :: r, h ! radius and height of shape
        integer, intent(in) :: s ! num sides of polygon
        integer, intent(in) :: wrap ! num of basepairs wrapped
        real(dp), dimension(s) :: constructPolygonPrism
        real(dp), dimension(s) :: t, x, y, z
        integer i
        real(dp), dimension(3,3) :: mtrx
        real(dp) space, incr

        ! construct rotation matrix
        mtrx(:,1) = v
        mtrx(:,2) = cross(u,v)
        mtrx(:,3) = u

        ! create parametric t
        space = 1.0/(2*s)
        incr = 1.0/(s)
        do i = 1, s
            t(i) = space
            space = space + incr
        enddo
        x = r*sin(t) + MATMUL(mtrx, cen)(3)
        y = r*sin(t) + MATMUL(mtrx, cen)(2)
        z(1:nint(s/2)) = cen(3) - h/2*MATMUL(mtrx, cen)(1)
        z(nint(s/2):) = cen(3) + h/2*MATMUL(mtrx, cen)(1)
        ! if nucleosome
        if (wrap /= 1) then

        else
            print*, 'figure out'
        endif

    END FUNCTION

    SUBROUTINE sameShapeTest()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1x = (/ 4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972, & 
                                      4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972/)
        real(dp), dimension(12) :: s1y = (/2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75, &
                                      2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75/)
        real(dp), dimension(12) :: s1z = (/0.,0.,0.,0.,0.,0.,-3.,-3.,-3.,-3.,-3.,-3./)
        real(dp), dimension(12) :: s2x = (/ 4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972, &
                                      4.76313972,  0. , -4.76313972, -4.76313972,  0. ,4.76313972/)
        real(dp), dimension(12) :: s2y = (/2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75, &            
                                      2.75  ,  5.5,  2.75 , -2.75, -5.5,-2.75/)
        real(dp), dimension(12) :: s2z = (/0.,0.,0.,0.,0.,0.,-3.,-3.,-3.,-3.,-3.,-3./)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: sameShapeTest"
            !stop
        endif

    END SUBROUTINE sameShapeTest

    SUBROUTINE noIntersectX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1x = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2x = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                                      4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        real(dp), dimension(12) :: s2y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectX"
            stop
        endif

    END SUBROUTINE noIntersectX

    SUBROUTINE intersectX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1x = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2x = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                                      .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        real(dp), dimension(12) :: s2y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectX"
            stop
        endif

    END SUBROUTINE intersectX

    SUBROUTINE tangentX()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1x = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2x = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                                      1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        real(dp), dimension(12) :: s2y = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentX"
            stop
        endif

    END SUBROUTINE tangentX

    SUBROUTINE noIntersectY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1y = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2y = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                                      4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectY"
            stop
        endif

    END SUBROUTINE noIntersectY

    SUBROUTINE intersectY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1y = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2y = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                                      .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res
 
        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectY"
            stop
        endif

    END SUBROUTINE intersectY
 
    SUBROUTINE tangentY()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1y = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2y = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                                      1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2z = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentY"
            stop
        endif

    END SUBROUTINE tangentY

    SUBROUTINE noIntersectZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1z = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2z = (/ 4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000, &
                                      4.5000,    4.0000,    4.5000,    5.5000,    6.0000,    5.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res > 0 ) then
            print*, "failed test: noIntersectZ"
            stop
        endif

    END SUBROUTINE noIntersectZ

    SUBROUTINE intersectZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1z = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2z = (/ .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000, &
                                      .5000,    0.000,    .5000,    1.5000,    2.0000,    1.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: intersectZ"
            stop
        endif

    END SUBROUTINE intersectZ

    SUBROUTINE tangentZ()
        implicit none
        integer :: nVerts = 12
        real(dp), dimension(12) :: s1z = (/-0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000, &
                                      -0.5000,   -1.0000,   -0.5000,    0.5000,    1.0000,    0.5000/)
        real(dp), dimension(12) :: s1x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                     -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s1y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        real(dp), dimension(12) :: s2z = (/ 1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000, &
                                      1.5000,    1.000,    1.5000,    2.5000,    3.0000,    2.5000/)
        real(dp), dimension(12) :: s2x = (/-0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660, &
                                      -0.8660,         0.,    0.8660,    0.8660,         0.,   -0.8660/)
        real(dp), dimension(12) :: s2y = (/0., 0., 0., 0., 0., 0., 5.5, 5.5, 5.5, 5.5, 5.5, 5.5/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        !print*, res
        if (res == 0 ) then
            print*, "failed test: tangentZ"
            stop
        endif

    END SUBROUTINE tangentZ


END MODULE

! test module
PROGRAM GJKTest 
    use GJKAlgorithm, only: GJK, sameShapeTest, noIntersectX, intersectX, tangentX, &
                            noIntersectY, intersectY, tangentY, &
                            noIntersectZ, intersectZ, tangentZ
    implicit none

    call sameShapeTest()
    call noIntersectX()
    call intersectX()
    call tangentX()
    call noIntersectY()
    call intersectY()
    call tangentY()
    call noIntersectZ()
    call intersectZ()
    call tangentZ()

    print*, "SUCCESS: successful completion of all GJK collision unit tests"

END PROGRAM