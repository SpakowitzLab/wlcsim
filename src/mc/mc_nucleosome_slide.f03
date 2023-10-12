#include "../defines.inc"
!--------------------------------------------------------------*
!
!              Slides nucleosomes in MC move
!
!               implemented by NP may 2020
!        TODO: THIS IS NOT TESTED NOR FULLY DEVELOPED!
!
!---------------------------------------------------------------

subroutine mc_nucleosome_slide(IB1, IB2, IT1, IT2, rand_stat, success)
! this move will "slide" a nucleosome bead some fractional basepair
! this is implemented by changed the discretization of the neighboring beads
! i.e. if you move a nucleosome +1bp along the chain, then the discretization to
! downstream the nucleosome will shorten to adjust for this and the discretization
! of the linker upstream will increase. likewise the nucleosome is then shifted
! along its U-vector in the direction of this change. this move was intended for
! simulations with linker discretization (i.e. WLC_P__INCLUDE_DISCRETIZE_LINKER==TRUE),
! but it should in *theory* work for simulations that discretize per nucleosome
   use params, only: wlc_V, wlc_R, wlc_RP, wlc_AB, wlc_U &
                     , wlc_UP, wlc_ABP, wlc_VP, wlc_pointsMoved, wlc_nPointsMoved, &
                     wlc_nucleosomeWrap, wlc_basepairs, wlc_nBend, wlc_bendPoints, &
                     wlc_basepairs_prop

   use mersenne_twister
   use params, only: dp
   use windowTools, only: exponential_random_int
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain, &
                             get_IB, length_of_chain
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
   real(dp) :: MCAMP ! Amplitude of random change
   real(dp) DR    ! Displacement for slide move (integer BPs)
   integer I ! test bead
   integer II, JJ, KK, J ! test indices
   integer prevNuc, nextNuc, linkerSum
   integer nNucs
   integer nucArray(WLC_P__NT)
   real(dp), parameter :: eps = 0.00001 ! rescale to avoid urand vals of 0
   integer, parameter :: min_base = 3
   integer, parameter :: max_base = 15

! initialize
   success = .false.

! find nucs
   KK = 1
   do II = 1, WLC_P__NT
      if (wlc_nucleosomeWrap(II) == 0) cycle
      nucArray(KK) = II
      KK = KK + 1
   enddo
   nNucs = KK - 1

! select nuc to move
   call random_number(urand, rand_stat)
   KK = ceiling(nNucs*(urand(1) + eps)/(1 + 1.1*eps))
   I = nucArray(KK)

! select distance to move (in bp)
   DR = 0.0_dp
   MCAMP = max(nint(0.25*WLC_P__LL/WLC_P__NBPL), 1)


do while (DR == 0.0_dp)
   call random_number(urand, rand_stat)
   DR = MCAMP*(urand(1) - 0.5_dp)
enddo

! find neighboring nuclesomes
   prevNuc = KK - 1
   if (prevNuc > 0) then
      prevNuc = nucArray(prevNuc)
      if (get_IP(prevNuc) /= get_IP(I)) then
         prevNuc = first_bead_of_chain(get_IP(I))
      endif
   else
      prevNuc = 1
   endif
   nextNuc = KK + 1
   if (nextNuc <= nNucs) then
      nextNuc = nucArray(nextNuc)
      if (get_IP(nextNuc) /= get_IP(I)) then
         nextNuc = last_bead_of_chain(get_IP(I)) - 1
      endif
   else
      nextNuc = WLC_P__NT - 1
   endif

   wlc_basepairs_prop = wlc_basepairs
! change distance between beads
   outer1: do II = 1, I - prevNuc ! explore the previous linker space
      inner1: do JJ = 0, (nextNuc - 1) - I ! explore the next linker space
         if ((wlc_basepairs(I - II) + DR > min_base) .AND. (wlc_basepairs(I - II) + DR < max_base) .AND. &
             (wlc_basepairs(I + JJ) - DR > min_base) .AND. (wlc_basepairs(I + JJ) - DR < max_base)) then
            wlc_basepairs_prop(I - II) = wlc_basepairs(I - II) + DR
            wlc_basepairs_prop(I + JJ) = wlc_basepairs(I + JJ) - DR
            success = .true.
            exit outer1
         endif
      enddo inner1
   enddo outer1

   if (success) then
      IT1 = I - II
      IT2 = I + JJ
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
      do KK = IT1 + 1, IT2
         wlc_RP(:, KK) = wlc_R(:, KK) + (wlc_U(:, KK - 1)+wlc_U(:, KK))*WLC_P__LENGTH_PER_BP*DR
         wlc_UP(:, KK) = wlc_U(:, KK)
         wlc_VP(:, KK) = wlc_V(:, KK)
         wlc_nPointsMoved = wlc_nPointsMoved + 1
         wlc_pointsMoved(wlc_nPointsMoved) = KK
      enddo
   else
      wlc_basepairs_prop = wlc_basepairs
   endif

end subroutine
