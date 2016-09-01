subroutine methyl_profile(nt,meth_status,ktot,km,kd,num_methylated,time,rxn_happen,pairs,dt,nuc_site,num_spread,num_decay)
    use mt19937, only : grnd
    implicit none
    integer, intent(in) :: nt, pairs(2,nt), nuc_site
    integer, intent(inout) :: meth_status(nt), rxn_happen, num_spread, num_methylated, num_decay
    double precision, intent(in) :: km, kd, ktot, time, dt
    double precision :: time_rxn, rn1, rn2, rn3, dt_mod, prob_no_rxn, pr
    integer :: site_rxn, count, i

    ! for pairs of beads that could transfer a methyl mark,
    ! perform Gillespie algorithm to determine if reaction happens and then update methyl profile
    dt_mod = dt
    if (rxn_happen.eq.1) then
        ! does a reaction occur?
        rn1 = grnd()
        prob_no_rxn = exp(-ktot*dt_mod)
        if (rn1.gt.prob_no_rxn) then
            ! which reaction occurred?
            rn2 = grnd()
            if (rn2.lt.((kd/ktot)*num_methylated)) then ! one site is demethylated
                site_rxn = ceiling(rn2/(kd/ktot))
                if (site_rxn.eq.nuc_site) then
                    rxn_happen = 0
                    go to 16
                end if
                count = 0
                i = 1
                do while (count.lt.site_rxn)
                    count = count + meth_status(i)
                    i = i+1
                end do
                meth_status(i-1) = 0 
                num_decay = num_decay + 1
            else ! one site is methylated 
                pr = rn2 - ((kd/ktot)*num_methylated)
                site_rxn = ceiling(pr/(km/ktot))
                meth_status(pairs(2,site_rxn)) = 1
                num_spread = num_spread + 1
            end if
            ! at what time did it occur?
            rn3 = grnd()
            time_rxn = time - (1/ktot)*log(rn2*(prob_no_rxn-1)+1)
            dt_mod = time + dt - time_rxn
        else
            rxn_happen = 0
        end if
        16 continue
    end if

    num_methylated = sum(meth_status)

end 
   
     
    
 
