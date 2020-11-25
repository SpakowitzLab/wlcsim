#include "../defines.inc"
!--------------------------------------------------------------*
!
!            Allows nucleosomes to breathe in MC move
!
!                  implemented by NP, oct 2020
!
!---------------------------------------------------------------

subroutine mc_nucleosome_breathe(IB1, IB2, IT1, IT2, rand_stat, success)
! this move will "breathe" a nucleosome to allow for under and overwrapping on
! the histone octamer. ill write more here soon.
   use params, only: wlc_V, wlc_R, wlc_RP, wlc_AB, wlc_U &
                     , wlc_UP, wlc_ABP, wlc_VP, wlc_pointsMoved, wlc_nPointsMoved, &
                     wlc_nucleosomeWrap, wlc_basepairs, wlc_nBend, wlc_bendPoints, &
                     wlc_basepairs_prop, wlc_nucleosomeWrap_prop

   use mersenne_twister
   use params, only: dp
   use windowTools, only: exponential_random_int
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain, get_IB, length_of_chain
   use nucleosome, only: nucleosome_prop

   implicit none
   integer, intent(out) :: IB1   ! Test bead position 1
   integer, intent(out) :: IT1   ! Index of test bead 1
   integer, intent(out) :: IB2   ! Test bead position 2
   integer, intent(out) :: IT2   ! Index of test bead 2
   logical, intent(out) :: success ! success of move to take place

! Things for random number generator
   type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
   real(dp) urand(1) ! single random number
   real(dp) DR    ! Displacement for breathe move (integer BPs)
   integer I, J, K ! test bead
   integer prevNuc, nextNuc, linkerSum
   integer nNucs
   integer nucArray(WLC_P__NT)
   real(dp) tempR(3), tempU(3), tempV(3)
   real(dp), parameter :: eps = 0.00001 ! rescale to avoid urand vals of 0
   integer, parameter :: max_wrap = 160
   integer, parameter :: min_base = 5
   integer, parameter :: max_base = 15

! initialize
   success = .false.

! find nucs
   K = 1
   do I = 1, WLC_P__NT
      if (wlc_nucleosomeWrap(I) == 0) cycle
      nucArray(K) = I
      K = K + 1
   enddo
   nNucs = K - 1

! select nuc to move
   call random_number(urand, rand_stat)
   K = ceiling(nNucs*(urand(1) + eps)/(1 + 1.1*eps))
   I = nucArray(K)

! select bp to unwrap 
   DR = 0.0_dp
   do while (DR == 0) 
      call random_gauss(urand, rand_stat)
      DR = urand(1)
   enddo

   wlc_basepairs_prop = wlc_basepairs
   wlc_nucleosomeWrap_prop = wlc_nucleosomeWrap


   ! change JUST entry-exit angle
   if (wlc_nucleosomeWrap(I) + DR <= max_wrap &
      .and. wlc_basepairs(I) - DR <= max_base &
      .and. wlc_basepairs(I) - DR >= min_base ) then 
      wlc_nucleosomeWrap_prop(I) = wlc_nucleosomeWrap(I) + DR
      wlc_basepairs_prop(I) = wlc_basepairs(I) - DR
      success = .true.
   endif

   if (success) then
      IT1 = I
      IT2 = I + 1
      IB1 = get_IB(IT1)
      IB2 = get_IB(IT2)
      if (IB1 > 1) then
         wlc_nBend = wlc_nBend + 1
         wlc_bendPoints(wlc_nBend) = IT1
         wlc_RP(:, IT1) = wlc_R(:, IT1)
         wlc_UP(:, IT1) = wlc_U(:, IT1)
         wlc_VP(:, IT1) = wlc_V(:, IT1)
         wlc_nPointsMoved = wlc_nPointsMoved + 1
         wlc_pointsMoved(wlc_nPointsMoved) = IT1
      endif
      if (IB2 < length_of_chain(get_IP(IT2))) then
         wlc_nBend = wlc_nBend + 1
         wlc_bendPoints(wlc_nBend) = IT2
         J = IT2 + 1
         wlc_RP(:, J) = wlc_R(:, J)
         wlc_UP(:, J) = wlc_U(:, J)
         wlc_VP(:, J) = wlc_V(:, J)
         wlc_nPointsMoved = wlc_nPointsMoved + 1
         wlc_pointsMoved(wlc_nPointsMoved) = J
      endif
      ! rotate
      call nucleosome_prop(wlc_U(:, IT1), wlc_V(:, IT1), wlc_R(:, IT1), &
                          wlc_basepairs_prop(IT1), wlc_nucleosomeWrap_prop(IT1), &
                          tempU, tempV, tempR)
      ! update nuc
      tempR = tempR + WLC_P__LENGTH_PER_BP*tempU*wlc_basepairs_prop(IT1) 
      wlc_RP(:, IT2) = tempR
      wlc_UP(:, IT2) = tempU
      wlc_VP(:, IT2) = tempV
      wlc_nPointsMoved = wlc_nPointsMoved + 1
      wlc_pointsMoved(wlc_nPointsMoved) = IT2
   else
      wlc_basepairs_prop = wlc_basepairs
      wlc_nucleosomeWrap_prop = wlc_nucleosomeWrap
   endif

end subroutine
