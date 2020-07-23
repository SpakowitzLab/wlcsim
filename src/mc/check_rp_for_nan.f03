subroutine check_rp_for_nan(success, MCTYPE)
   use params, only: wlc_RP, ERROR_UNIT, moveNames, wlc_nPointsMoved, &
                     wlc_pointsMoved, dp, wlc_spider_id, wlc_spider_dr
   implicit none
   logical, intent(out) :: success
   integer, intent(in) :: MCTYPE
   integer I, J
   success = .True.
   do I = 1, wlc_nPointsMoved
      do J = 1, 3
         if (isnan(wlc_RP(J, wlc_pointsMoved(I)))) then
            success = .False.
         endif
      enddo
   enddo
   if (.not. success) then
      write (ERROR_UNIT, *) "--------------------------------"
      write (ERROR_UNIT, *) "WARNING: wlc_RP is NAN at moved bead"
      write (ERROR_UNIT, *) "MCTYPE", MCTYPE, " ", moveNames(MCTYPE)
      do I = 1, wlc_nPointsMoved
         write (ERROR_UNIT, *) "Bead", I, " is ", wlc_pointsMoved(I)
         write (ERROR_UNIT, *) wlc_RP(:, wlc_pointsMoved(I))
      enddo
      if (MCTYPE == 12) then
         call mc_describe_spider(wlc_spider_id, wlc_spider_dr, ERROR_UNIT)
      endif
      write (ERROR_UNIT, *) "--------------------------------"
   endif

end subroutine
