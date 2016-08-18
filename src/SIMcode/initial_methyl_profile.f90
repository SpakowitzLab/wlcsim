subroutine initial_methyl_profile(nt)
    implicit none
    integer, intent(in) :: nt
    integer :: ind, nuc_site, meth_status(nt)

    nuc_site = ceiling(real(nt/2.0))

    do ind = 1, nuc_site - 1
        meth_status(ind) = 0
    end do

    meth_status(nuc_site) = 1

    do ind = nuc_site + 1, nt
        meth_status(ind) = 0
    end do
    
end 
   
     
    
 
