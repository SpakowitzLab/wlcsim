subroutine tot_rate_constant(nt,could_react,meth_status,km,kd,ktot,num_methylated)
    implicit none
    integer, intent(in) :: nt, could_react, meth_status(nt)
    double precision, intent(in) :: km, kd
    double precision, intent(inout) :: ktot
    integer, intent(inout) :: num_methylated

    ! determine total rate constant for all possible reactions
    num_methylated = sum(meth_status)
    ktot = num_methylated*kd + could_react*km
end
  
    
    
   
     
