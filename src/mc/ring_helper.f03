#include "../defines.inc"
module ringHelper
   implicit none
contains

! Add a bend point which is to the left of the selected segment
   subroutine bend_points_left_ring(IT1)
      use polydispersity, only: is_left_end, rightmost_from
      use params, only: wlc_nBend, wlc_bendPoints
      implicit none
      integer, intent(in) :: IT1

      wlc_nBend = wlc_nBend + 1
      if (.not. is_left_end(IT1)) then
         wlc_bendPoints(wlc_nBend) = IT1 - 1
      else
         wlc_bendPoints(wlc_nBend) = rightmost_from(IT1)
      endif
   end subroutine

! Add a bend point which is to the right of the selected segment
   subroutine bend_points_right_ring(IT2)
      use params, only: wlc_nBend, wlc_bendPoints
      implicit none
      integer, intent(in) :: IT2

      wlc_nBend = wlc_nBend + 1
      wlc_bendPoints(wlc_nBend) = IT2
   end subroutine

! Add a point moved which is to the left  of the selected segment
   subroutine points_moved_left_ring(IT1)
      use polydispersity, only: is_left_end, rightmost_from
      use params, only: wlc_RP, wlc_R, wlc_U, wlc_VP, wlc_V, wlc_UP, &
                        wlc_pointsMoved, wlc_nPointsMoved
      implicit none
      integer, intent(in) :: IT1

      integer I

      wlc_nPointsMoved = wlc_nPointsMoved + 1
      if (.not. is_left_end(IT1)) then
         I = IT1 - 1
      else
         I = rightmost_from(IT1)
      endif
      wlc_RP(:, I) = wlc_R(:, I)
      wlc_UP(:, I) = wlc_U(:, I)
      if (WLC_P__LOCAL_TWIST) wlc_VP(:, I) = wlc_V(:, I)
      wlc_pointsMoved(wlc_nPointsMoved) = I
   end subroutine

! Add a point moved which is to the right of the selected segment
   subroutine points_moved_right_ring(IT2)
      use polydispersity, only: is_right_end, leftmost_from
      use params, only: wlc_RP, wlc_R, wlc_U, wlc_VP, wlc_V, wlc_UP, &
                        wlc_pointsMoved, wlc_nPointsMoved
      implicit none
      integer, intent(in) :: IT2

      integer I

      wlc_nPointsMoved = wlc_nPointsMoved + 1
      if (.not. is_right_end(IT2)) then
         I = IT2 + 1
      else
         I = leftmost_from(IT2)
      endif
      wlc_RP(:, I) = wlc_R(:, I)
      wlc_UP(:, I) = wlc_U(:, I)
      if (WLC_P__LOCAL_TWIST) wlc_VP(:, I) = wlc_V(:, I)
      wlc_pointsMoved(wlc_nPointsMoved) = I
   end subroutine

! Find next bead in the ring from the input index.
! This subroutine handles the case when the current bead is at the end
! of the chain array.
   function nextBead(IT) result(ITNext)
      use polydispersity, only: is_right_end, leftmost_from
      implicit none
      integer, intent(in) :: IT
      integer ITNext

      if (WLC_P__NP == 1) then
         if (IT < WLC_P__NB) then
            ITNext = IT + 1
         else
            ITNext = 1
         endif
      else
         if (.not. is_right_end(IT)) then
            ITNext = IT + 1
         else
            ITNext = leftmost_from(IT)
         endif
      endif
   end function

! Find previous bead in the ring from the input index.
   function prevBead(IT) result(ITPrev)
      use polydispersity, only: is_left_end, rightmost_from
      implicit none
      integer, intent(in) :: IT
      integer ITPrev

      if (WLC_P__NP == 1) then
         if (IT > 1) then
            ITPrev = IT - 1
         else
            ITPrev = WLC_P__NB
         endif
      else
         if (.not. is_left_end(IT)) then
            ITPrev = IT - 1
         else
            ITPrev = rightmost_from(IT)
         endif
      endif
   end function

end module ringHelper
