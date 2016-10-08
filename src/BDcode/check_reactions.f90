subroutine check_reactions(r, nt, meth_status, in_rxn_rad, could_react, fpt_dist, pairs)
    implicit none
    integer, intent(in) :: meth_status(nt)
    double precision, intent(inout) :: r(nt,3), fpt_dist
    double precision :: check_pair(nt,nt)
    integer, intent(inout) :: could_react, nt, in_rxn_rad(nt,nt), pairs(2,nt)
    integer :: k1,k2

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
    do k1 = 1, nt - 1
        do k2 = k1 + 1, nt
            if ((in_rxn_rad(k1,k2).eq.1).and.(check_pair(k1,k2).eq.0)) then
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
        end do
    end do
end 
     

