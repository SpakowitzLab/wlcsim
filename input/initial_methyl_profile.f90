subroutine methyl_profile(nt)
  implicit none
    integer :: ind
    integer, intent(in) :: nt
    integer, intent(out) :: meth_status(nt)
    ! set initial methylation profile
    
    ! current version: place one nucleation site in the center of the chain
    do ind = 1,nt/2
        meth_status(ind) = 0
    end do
    meth_status((nt/2)+1) = 1
    do ind = (nt/2)+2,nt
        meth_status(ind) = 0
    end do

   
    



    

     
    
        
  
   
     
    
 
