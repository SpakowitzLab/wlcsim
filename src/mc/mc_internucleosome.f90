#include "../defines.inc"
subroutine mc_internucleosome()
!Calculate the change in the polymer internucleosome attraction energy
! the actual energy function is in the nucleosome module file and is 
! based off of work from the de pablo group

! values from wlcsim_data
   use params, only: dp, NAN, wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V&
                     , wlc_R, wlc_UP, wlc_RP, wlc_nPointsMoved, wlc_pointsMoved
   use nucleosome, only: internucleosome_energy
   use energies, only: energyOf, internucleosome_

   implicit none
   real(dp) delInt     ! change in internucleosome attraction energy
   integer i, j, k
   integer ignore(WLC_P__NT) 

   ignore = NAN
   delInt = 0
   do k = 1, wlc_nPointsMoved
      i = wlc_pointsMoved(k)
      ignore(k) = i
      if (wlc_nucleosomeWrap(i) == 1) cycle
      do j = 1, WLC_P__NT
         if (wlc_nucleosomeWrap(j) == 1 .or. ANY(ignore == j)) cycle
         ! old config
         delInt = delInt - internucleosome_energy(wlc_R(:, i), wlc_R(:, j), &
                                                  wlc_U(:, i), wlc_U(:, j), &
                                                  wlc_V(:, i), wlc_V(:, j))
         ! new config
         if (isnan(wlc_RP(1, i)) .eqv. .false.) then 
            if (isnan(wlc_RP(1, j)) .eqv. .false.) then ! i in moved, j in moved
               delInt = delInt + internucleosome_energy(wlc_RP(:, i), wlc_RP(:, j), &
                                                        wlc_UP(:, i), wlc_UP(:, j), &
                                                        wlc_VP(:, i), wlc_VP(:, j))
            else ! i in moved, j not in moved
               delInt = delInt + internucleosome_energy(wlc_RP(:, i), wlc_R(:, j),&
                                                        wlc_UP(:, i), wlc_U(:, j),&
                                                        wlc_VP(:, i), wlc_V(:, j))
            endif
         else 
            if (isnan(wlc_RP(1, j)) .eqv. .false.) then ! i not in moved, j in moved
               delInt = delInt + internucleosome_energy(wlc_R(:, i), wlc_RP(:, j), &
                                                        wlc_U(:, i), wlc_UP(:, j), &
                                                        wlc_V(:, i), wlc_VP(:, j))
            else ! i not in moved, j not in moved
               delInt = delInt + internucleosome_energy(wlc_R(:, i), wlc_R(:, j), &
                                                        wlc_U(:, i), wlc_U(:, j), &
                                                        wlc_V(:, i), wlc_V(:, j))
            endif
         endif
      enddo
   enddo
   energyOf(internucleosome_)%dx = delInt

END subroutine 