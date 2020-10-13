#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Slides nucleosomes in MC move
!
!    inspired by Quinn (chemMove/crank), implemented by NP 2020
!
!---------------------------------------------------------------

subroutine mc_nucleosome_wrap(IB1, IB2, IT1, IT2, rand_stat, success)
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
                     wlc_basepairs_prop, wlc_nucleosomeWrap_prop

   use mersenne_twister
   use params, only: dp
   use windowTools, only: exponential_random_int
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
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
   real(dp), parameter :: optRatio = 0.25 ! max ratio of average discretization to allow for bp slide
   integer I ! test bead
   integer II, JJ, KK, J ! test indices
   integer max_bp
   integer prevNuc, nextNuc, linkerSum
   integer nNucs
   integer nucArray(WLC_P__NT)
   real(dp) tempR(3), tempU(3), tempV(3)
   real(dp), parameter :: eps = 0.00001 ! rescale to avoid urand vals of 0

! initialize
   success = .false.

! find nucs
   KK = 1
   do II = 1, WLC_P__NT
      if (wlc_nucleosomeWrap(II) == 1) cycle
      nucArray(KK) = II
      KK = KK + 1
   enddo
   nNucs = KK - 1

! select nuc to move
   call random_number(urand, rand_stat)
   KK = ceiling(nNucs*(urand(1) + eps)/(1 + 1.1*eps))
   I = nucArray(KK)

! select bp to unwrap 
   DR = 0.0_dp
   call random_gauss(urand, rand_stat)

   DR = nint(urand(1))

! ! find neighboring nuclesomes
!    prevNuc = KK - 1
!    if (prevNuc > 0) then
!       prevNuc = nucArray(prevNuc)
!       if (get_IP(prevNuc) /= get_IP(I)) then
!          prevNuc = first_bead_of_chain(get_IP(I))
!       endif
!    else
!       prevNuc = 1
!    endif
!    nextNuc = KK + 1
!    if (nextNuc <= nNucs) then
!       nextNuc = nucArray(nextNuc)
!       if (get_IP(nextNuc) /= get_IP(I)) then
!          nextNuc = last_bead_of_chain(get_IP(I)) - 1
!       endif
!    else
!       nextNuc = WLC_P__NT - 1
!    endif

   wlc_basepairs_prop = wlc_basepairs
   wlc_nucleosomeWrap_prop = wlc_nucleosomeWrap
! ! change distance between beads
!    outer1: do II = 1, I - prevNuc ! explore the previous linker space
!       inner1: do JJ = 0, (nextNuc - 1) - I ! explore the next linker space
!          if ((wlc_basepairs(I - II) + DR > 3) .AND. (wlc_basepairs(I - II) + DR < max_bp) .AND. &
!              (wlc_basepairs(I + JJ) - DR > 3) .AND. (wlc_basepairs(I + JJ) - DR < max_bp)) then
!             wlc_basepairs_prop(I - II) = wlc_basepairs(I - II) + DR
!             wlc_basepairs_prop(I + JJ) = wlc_basepairs(I + JJ) - DR
!             wlc_nucleosomeWrap_prop(I) = wlc_nucleosomeWrap(I) + DR
!             success = .true.
!             exit outer1
!          endif
!       enddo inner1
!    enddo outer1

! change JUST entry-exit angle
if (wlc_nucleosomeWrap_prop(I) + DR <= 150 .and. wlc_nucleosomeWrap_prop(I) + DR >= 144 &
   .and. wlc_basepairs_prop(I) + DR <= 15 .and. wlc_basepairs_prop(I) + DR >= 3 ) then 
   wlc_nucleosomeWrap_prop(I) = wlc_nucleosomeWrap(I) + DR
   wlc_basepairs_prop(I) = wlc_basepairs(I) - DR
   success = .true.
endif

   if (success) then
      IB1 = I !- 2 !- II
      IB2 = I + 1!+ 2 !+ JJ
      IT1 = IB1
      IT2 = IB2 
      ! update nuc
      wlc_RP(:, I) = wlc_R(:, I)
      wlc_UP(:, I) = wlc_U(:, I)
      wlc_VP(:, I) = wlc_V(:, I)
      wlc_nPointsMoved = wlc_nPointsMoved + 1
      wlc_pointsMoved(wlc_nPointsMoved) = I
      ! rotate
      call nucleosome_prop(wlc_U(:, I), wlc_V(:, I), wlc_R(:, I), &
                          wlc_basepairs(I), wlc_nucleosomeWrap(I), &
                          tempU, tempV, wlc_RP(:, I + 1))
      call nucleosome_prop(wlc_U(:, I), wlc_V(:, I), wlc_R(:, I), &
                          wlc_basepairs_prop(I), wlc_nucleosomeWrap_prop(I), &
                          wlc_UP(:, I + 1), wlc_VP(:, I + 1), tempR)
      ! after nuc 
      wlc_RP(:, I + 1) = wlc_RP(:, I + 1) + WLC_P__LENGTH_PER_BP*tempU*wlc_basepairs_prop(I)
      !wlc_UP(:, I + 1) = wlc_U(:, I + 1)
      !wlc_VP(:, I + 1) = wlc_V(:, I + 1)
      wlc_nPointsMoved = wlc_nPointsMoved + 1
      wlc_pointsMoved(wlc_nPointsMoved) = I + 1
   else
      wlc_basepairs_prop = wlc_basepairs
      wlc_nucleosomeWrap_prop = wlc_nucleosomeWrap
   endif

end subroutine
