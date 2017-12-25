! ---------------------------------------
!
!   This function calculates whether two cylinders with
!   hemispherical caps intersect.
!   Written by Quinn in Dec, 2017
!   General algorithm for closest approch of two lines by Dan Sunday
!
! ---------------------------------------
#include "../defines.inc"
function cylinders(x1,y1,x2,y2) result(collide)
    !use params, only : dp, eps
    use precision, only: dp

    implicit none

    real(dp), intent(in) :: x1(3) ! start of segment 1
    real(dp), intent(in) :: y1(3) ! end of segment 1
    real(dp), intent(in) :: x2(3) ! start of segment 2
    real(dp), intent(in) :: y2(3) ! end of segment 2
    logical  collide

    real(dp)  distance
    real(dp) z1(3) ! nearest point on segment 1
    real(dp) z2(3) ! nearest point on segment 2
    real(dp) v1(3)
    real(dp) v2(3)
    real(dp) ww(3)
    real(dp) t1 ! parameter for distance along line 1. 0-1 is segment
    real(dp) t2
    real(dp) aa,bb,cc,dd,ee
    real(dp) denom
    logical linePoint !function
    real(dp) eps
    real(dp) l1,l2,l3,l4
    real(dp) pp(3)
    real(dp) offset
    eps = 0.000000001_dp

    ! exclude interaction if it has no chance of happening
    if ( min(x1(1),y1(1)) > max(x2(1),y2(1)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    elseif ( min(x2(1),y2(1)) > max(x1(1),y1(1)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    elseif ( min(x1(2),y1(2)) > max(x2(2),y2(2)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    elseif ( min(x2(2),y2(2)) > max(x1(2),y1(2)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    elseif ( min(x1(3),y1(3)) > max(x2(2),y2(2)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    elseif ( min(x2(3),y2(3)) > max(x1(3),y1(3)) + WLC_P__CHAIN_D ) then
        collide = .FALSE.
        return
    endif

    v1 = y1 - x1
    v2 = y2 - x2
    ww = x1 - x2
    aa = dot_product(v1,v1)
    bb = dot_product(v1,v2)
    cc = dot_product(v2,v2)
    dd = dot_product(v1,ww)
    ee = dot_product(v2,ww)
    denom = aa*cc - bb*bb

    if (denom < eps) then
        !Nearly parallel
        pp = (y1-x1)/norm2(y1-x1)
        
        l1 =0.0_dp
        l2 = dot_product(y1-x1,pp)
        l3 = dot_product(x2-x1,pp)
        l4 = dot_product(y2-x1,pp)
        if (min(l1,l2) > max(l3,l4)) then
            offset = min(l1,l2) - max(l3,l4)
        elseif ( max(l1,l2) < min(l3,l4) ) then 
            offset = min(l3,l4) - max(l1,l2)
        else
            offset = 0.0_dp
        endif
        ww = x2 - x1
        ww = ww - pp*dot_product(ww,pp)
        if (dot_product(ww,ww) + offset**2 > WLC_P__CHAIN_D**2) then
            collide = .FALSE.
            return
        else
            !print*, "Nearly parallel collide. offset=",offset
            collide = .True.
            return
        endif
    endif

    t1 = (bb*ee - cc*dd)/denom
    t2 = (aa*ee - bb*dd)/denom
   
    !print*, "denom=",denom," t1",t1," t2",t2
    if (t1 < 0.0_dp .and. t2 < 0.0_dp) then
        distance = dot_product(x1 - x2,x1 - x2)
        if (distance > WLC_P__CHAIN_D**2) then 
            collide = .FALSE.
            return
        else
            !print*, "x ends collide"
            collide = .TRUE.
            return
        endif
    elseif (t1 < 0.0_dp .and. t2 > 1.0_dp) then
        distance = dot_product(x1 - y2,x1 - y2)
        if (distance > WLC_P__CHAIN_D**2) then 
            collide = .FALSE.
            return
        else
            !print*, "x1 collides with y2"
            collide = .TRUE.
            return
        endif
    elseif (t1 > 1.0_dp .and. t2 < 0.0_dp) then
        distance = dot_product(y1 - x2,y1 - x2)
        if (distance > WLC_P__CHAIN_D**2) then 
            collide = .FALSE.
            return
        else
            !print*, "y1 collides with x2"
            collide = .TRUE.
            return
        endif
    elseif (t1 > 1.0_dp .and. t2 > 1.0_dp) then
        distance = dot_product(y1 - y2,y1 - y2)
        if (distance > WLC_P__CHAIN_D**2) then 
            collide = .FALSE.
            return
        else
            !print*, "y1 collides with y2"
            collide = .TRUE.
            return
        endif
    elseif (t1 < 0.0_dp) then
        collide = linePoint(v2,x1-x2)
    elseif (t1 > 1.0_dp) then
        collide = linePoint(v2,y1-x2)
    elseif (t2 < 0.0_dp) then
        collide = linePoint(v1,x2-x1)
    elseif (t2 > 1.0_dp) then
        collide = linePoint(v1,y2-x1)
    else
        ! intermediate segment
        distance = norm2(x1 - x2 + ((bb*ee-cc*dd)*v1 - (aa*ee-bb*dd)*v2)/denom)
        if (distance > WLC_P__CHAIN_D) then
            collide = .FALSE.
            return
        else
            !print*, "two stick collision"
            collide = .TRUE.
            return
        endif
    endif


end function cylinders

#include "../defines.inc"
function linePoint(vv,xx)
    use precision, only: dp

    implicit none

    real(dp), intent(in) :: vv(3) ! vector
    real(dp), intent(in) :: xx(3) ! point
    logical  linePoint
    real(dp) dd(3)
    real(dp) pp(3)

    pp = vv/norm2(vv)
    dd = pp*dot_product(pp,xx) - xx
    
    if (dot_product(dd,dd) > WLC_P__CHAIN_D**2) then
        linePoint = .FALSE.
        return
    else
        print*, "Line point collision"
        linePoint = .TRUE.
        return
    endif
    return
end function linePoint


