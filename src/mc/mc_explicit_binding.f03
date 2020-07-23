#include "../defines.inc"
!---------------------------------------------------------------!
!
!    This subroutine calculates the energy of "explicit_binding"
!    which is a spring potential between specified beads meant to
!    represent the effect of things like Loop Extrusion Factors.
!
!---------------------------------------------------------------!
subroutine mc_explicit_binding()
! values from wlcsim_data
   use params, only: wlc_R, wlc_RP, &
                     wlc_ExplicitBindingPair, wlc_nPointsMoved, wlc_pointsMoved, &
                     wlc_other_beads, wlc_network_start_index
   use params, only: dp
   use energies, only: energyOf, explicitBinding_
   implicit none

!   Internal variables

   integer ii, jj, indx
   integer otherEnd
   real(dp) r(3)

   energyOf(explicitBinding_)%dx = 0.0_dp

   do jj = 1, wlc_nPointsMoved
      ii = wlc_pointsMoved(jj)
      if (WLC_P__NETWORK) then
         do indx = wlc_network_start_index(ii), &
            wlc_network_start_index(ii + 1) - 1
            otherEnd = wlc_other_beads(indx)
            ! plus new
            if (isnan(wlc_RP(1, otherEnd))) then  ! assumes that RP is up to date or nan
               r = wlc_RP(:, ii) - wlc_R(:, otherEnd)
            else
               if (otherEnd .le. ii) cycle ! don't doule count
               r = wlc_RP(:, ii) - wlc_RP(:, otherEnd)
            endif
            energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx + &
                                            (r(1)**2 + r(2)**2 + r(3)**2)

            ! minus old
            r = wlc_R(:, ii) - wlc_R(:, otherEnd)
            energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx - &
                                            (r(1)**2 + r(2)**2 + r(3)**2)
         enddo
      else
         otherEnd = wlc_ExplicitBindingPair(ii)
         if (otherEnd .le. 0) cycle

         ! plus new
         if (isnan(wlc_RP(1, otherEnd))) then  ! assumes that RP is up to date or nan
            r = wlc_RP(:, ii) - wlc_R(:, otherEnd)
         else
            if (otherEnd .le. ii) cycle ! don't doule count
            r = wlc_RP(:, ii) - wlc_RP(:, otherEnd)
         endif
         energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx + &
                                         (r(1)**2 + r(2)**2 + r(3)**2)

         ! minus old
         r = wlc_R(:, ii) - wlc_R(:, otherEnd)
         energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx - &
                                         (r(1)**2 + r(2)**2 + r(3)**2)
      endif

   enddo
   RETURN
END

!---------------------------------------------------------------!
subroutine mc_explicit_binding_from_scratch()
! values from wlcsim_data
   use params, only: wlc_R, wlc_ExplicitBindingPair, wlc_network_start_index, &
                     wlc_other_beads
   use params, only: dp
   use energies, only: energyOf, explicitBinding_
   implicit none

!   iputs

!   Internal variables

   integer ii, indx
   integer otherEnd
   real(dp) r(3)

   energyOf(explicitBinding_)%dx = 0.0_dp

   do ii = 1, WLC_P__NT
      if (WLC_P__NETWORK) then
         do indx = wlc_network_start_index(ii), &
            wlc_network_start_index(ii + 1) - 1
            otherEnd = wlc_other_beads(indx)
            if (otherEnd .le. ii) cycle ! don't double count

            ! plus new
            r = wlc_R(:, ii) - wlc_R(:, otherEnd)
            energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx + &
                                            (r(1)**2 + r(2)**2 + r(3)**2)
         enddo
      else
         otherEnd = wlc_ExplicitBindingPair(ii)
         if (otherEnd .le. 0) cycle
         if (wlc_ExplicitBindingPair(otherEnd) .ne. ii) then
            print *, "Reflexive not found", ii, otherEnd
         endif
         if (otherEnd .le. ii) cycle ! Exclude <0 and don't double count

         ! plus new
         r = wlc_R(:, ii) - wlc_R(:, otherEnd)
         energyOf(explicitBinding_)%dx = energyOf(explicitBinding_)%dx + &
                                         (r(1)**2 + r(2)**2 + r(3)**2)
      endif
   enddo
   RETURN
END

!---------------------------------------------------------------!
