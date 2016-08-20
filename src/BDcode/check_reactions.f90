subroutine check_reactions(r, nt, meth_status, in_rxn_rad, could_react, fpt_dist)
    implicit none
    integer, intent(in) :: in_rxn_rad(nt,nt), meth_status(nt)
    double precision, intent(inout) :: r(nt,3), fpt_dist
    double precision :: check_pair(nt,nt)
    integer, intent(inout) :: could_react, nt
    integer :: k1, k2, c1, c2, pairs(2,nt)
    ! initialize variables
    do k1 = 1, nt
        do k2 = 1, nt
            check_pair(k1,k2) = 0
        end do
     end do

    do k1 = 1, nt
        pairs(1,k1) = 0
        pairs(2,k1) = 0
    end do
     
    ! for pairs of beads that are close enough to react, check that one
    ! is methylated and one unmethylated
    do k1 = 1, nt
        if (sum(in_rxn_rad(k1,:)).eq.1) then
            k2 = maxloc(in_rxn_rad(k1,:),1)
            if (check_pair(k1,k2).eq.0) then
                check_pair(k1,k2) = 1
                check_pair(k2,k1) = 1
                if (meth_status(k1).eq.1 .and. meth_status(k2).eq.0) then
                    could_react = could_react + 1
                    pairs(1,could_react) = k1
                    pairs(2,could_react) = k2
                else if (meth_status(k1).eq.0 .and. meth_status(k2).eq.1) then
                    could_react = could_react + 1
                    pairs(1,could_react) = k2
                    pairs(2,could_react) = k1
                end if
            end if
        else if (sum(in_rxn_rad(k1,:)).gt.1) then
            c1 = k1
            call calc_dist(r,nt,c1,fpt_dist,in_rxn_rad)
            k2 = c2
            if (check_pair(k1,k2).eq.0) then
                check_pair(k1,k2) = 1
                check_pair(k2,k1) = 1
                if (meth_status(k1).eq.1 .and. meth_status(k2).eq.0) then
                    could_react = could_react + 1
                    pairs(1,could_react) = k1
                    pairs(2,could_react) = k2
                else if (meth_status(k1).eq.0 .and. meth_status(k2).eq.1) then
                    could_react = could_react + 1
                    pairs(1,could_react) = k2
                    pairs(2,could_react) = k1
                end if
            end if
        end if
    end do
end 
     

