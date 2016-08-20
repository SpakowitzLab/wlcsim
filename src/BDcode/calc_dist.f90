subroutine calc_dist(r,nt,c1,fpt_dist,in_rxn_rad)
    implicit none
    integer :: ind, c2
    double precision :: d, distance
    integer, intent(in) :: nt, c1, in_rxn_rad(nt,nt)
    double precision, intent(in) :: r(nt,3), fpt_dist
    ! when more than 2 beads are within reaction distance, find the 2 that are closest together
    distance = fpt_dist 
    do ind = 1, nt
        if ((in_rxn_rad(c1,ind).eq.1) .and. (c1.ne.ind)) then
            d = sqrt(((r(c1,1)-r(ind,1))**2)+((r(c1,2)-r(ind,2))**2)+((r(c1,3)-r(ind,3))**2))
            if (d.lt.distance) then
                distance = d
                c2 = ind
            end if
        end if
    end do
end
   
