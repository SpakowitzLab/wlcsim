subroutine calc_dist(r,nt,fpt_dist,in_rxn_rad,k1,k2)
    implicit none
    integer :: i
    double precision :: d, distance
    integer, intent(in) :: nt, in_rxn_rad(nt,nt), k1
    double precision, intent(in) :: r(nt,3), fpt_dist
    integer,intent(out) :: k2
    ! when more than 2 beads are within reaction distance, find the 2 that are closest together
    distance = fpt_dist 
    do i = 1, nt
        if ((in_rxn_rad(k1,i).eq.1) .and. (k1.ne.i)) then
            d = sqrt(((r(k1,1)-r(i,1))**2)+((r(k1,2)-r(i,2))**2)+((r(k1,3)-r(i,3))**2))
            if (d.lt.distance) then
                distance = d
                k2 = i
            end if
        end if
    end do
end
   
