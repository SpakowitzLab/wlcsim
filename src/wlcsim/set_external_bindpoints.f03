#include "../defines.inc"
subroutine set_external_bindpoints(rand_stat)
   use params, only: dp, wlc_external_bind_points
   use mersenne_twister
   implicit none
   type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
   integer ii
   integer test(1)
   wlc_external_bind_points = .False.
   if (WLC_P__N_EXTERNAL_BIND_POINTS > WLC_P__NT) then
      print *, "Can bind more thant WLC_P__NT"
      stop
   endif
   open (unit=1, file='data/bind_points', status='new')
   do ii = 1, WLC_P__N_EXTERNAL_BIND_POINTS
      do while (.True.)
         call random_index(WLC_P__NT, test, rand_stat)
         if (.not. wlc_external_bind_points(test(1))) then
            wlc_external_bind_points(test(1)) = .True.
            write (1, *) test(1)
            exit
         endif
      enddo
   enddo
   close (1)
end subroutine

