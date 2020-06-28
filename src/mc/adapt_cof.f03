Subroutine adapt_cof(downSuccess, nPTReplicas, S_val, N_average, &
                     lowerRepExe, upperRepExe, &
                     lowerRail, upperRail, &
                     repAnnealSpeed, replicaBounds)
   use params, only: dp
   implicit none

! inputs
   integer, intent(in) :: nPTReplicas
   integer, intent(in) :: downSuccess(nPTReplicas)
   integer, intent(in) :: N_average
   real(dp), intent(in) :: lowerRepExe
   real(dp), intent(in) :: upperRepExe
   real(dp), intent(in) :: lowerRail
   real(dp), intent(in) :: upperRail
   real(dp), intent(in) :: repAnnealSpeed
   logical, intent(in) :: replicaBounds ! output between 0 and 1
! input/output
   real(dp), intent(inout) :: S_val(nPTReplicas)

! internal variables
   real(dp), allocatable :: newS_val(:)
   integer rep
   real(dp) successRate

   allocate (newS_val(1:nPTReplicas))

   if (.false.) then
      if ((downSuccess(2) .lt. 0.99_dp) .and. (newS_val(1) .lt. 1.0_dp)) then
         newS_val(1) = S_val(1) + repAnnealSpeed
      else
         newS_val(1) = S_val(1)
      endif
   else
      newS_val(1) = S_val(1)
   endif
   do rep = 2, nPTReplicas
      successRate = dble(downSuccess(rep))/dble(N_average)
      if (S_val(rep) .lt. S_val(rep - 1)) then
         print *, "Error in adapt_cof!"
         print *, "S_val(", rep, ") = ", S_val(rep), " < S_val(", rep - 1, ") = ", S_val(rep - 1)
         print *, S_val
         stop 1
      endif
      ! addapt spacing
      if (successRate .lt. lowerRepExe) then
         newS_val(rep) = newS_val(rep - 1) + (S_val(rep) - S_val(rep - 1))*0.95_dp
      elseif (successRate .gt. upperRepExe) then
         newS_val(rep) = newS_val(rep - 1) + (S_val(rep) - S_val(rep - 1))*1.05_dp
      else
         newS_val(rep) = newS_val(rep - 1) + S_val(rep) - S_val(rep - 1)
      endif
      ! don't addapt too fast.  a.k.a. anneal
      if (newS_val(rep) .gt. repAnnealSpeed + S_val(rep)) then
         newS_val(rep) = S_val(rep) + repAnnealSpeed
      elseif (newS_val(rep) .lt. S_val(rep) - repAnnealSpeed) then
         newS_val(rep) = S_val(rep) - repAnnealSpeed
      endif

      ! inforce limits on addaption
      if ((newS_val(rep) - newS_val(rep - 1)) .lt. lowerRail) then
         newS_val(rep) = newS_val(rep - 1) + lowerRail
      endif
      if ((newS_val(rep) - newS_val(rep - 1)) .gt. upperRail) then
         newS_val(rep) = newS_val(rep - 1) + upperRail
      endif
      if (newS_val(rep) .lt. newS_val(rep - 1)) then
         print *, "Error in adaptS_val!!!"
         stop 1
      endif

      if (replicaBounds) then
         ! Warning, I think this causes problems with the lower rail
         if (newS_val(rep) .lt. 0.0_dp) newS_val(rep) = 0.0_dp
         if (newS_val(rep) .gt. 1.0_dp) newS_val(rep) = 1.0_dp
      endif
   enddo
   S_val = newS_val
   deallocate (newS_val)
end subroutine
