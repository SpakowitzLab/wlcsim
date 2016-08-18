subroutine tot_rate_constant(nt,could_react,meth_status,km,kd)
    implicit none
    integer, intent(in) :: nt, could_react, meth_status(nt)
    double precision, intent(in) :: km, kd
    double precision :: ktot
    integer :: num_methylated
    ! determine total rate constant for all possible reactions
    num_methylated = sum(meth_status)
    ktot = num_methylated*kd + could_react*km
end
  
    
    
   
     
