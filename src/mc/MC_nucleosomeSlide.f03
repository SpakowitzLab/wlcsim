#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Slides nucleosomes in MC move
!
!    inspired by Quinn (chemMove/crank), implemented by NP 2020
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_nucleosomeSlide(IT1,IT2,rand_stat,success)
! values from wlcsim_data
use params, only: wlc_V, wlc_R, wlc_RP, wlc_AB, wlc_U&
    , wlc_UP, wlc_ABP, wlc_VP, wlc_pointsMoved, wlc_nPointsMoved, &
    wlc_nucleosomeWrap, wlc_basepairs, wlc_nBend, wlc_bendPoints

use mersenne_twister
use params, only: dp
use windowTools, only: exponential_random_int
use polydispersity, only: get_I
use nucleosome, only: nucleosomeProp

implicit none
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) ::success

! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urand(1) ! single random number
real(dp) :: window = 1 ! move 1 bps 
integer DR    ! Displacement for slide move (integer BPs)
integer I ! test bead
real(dp), dimension(3) :: tempLoc ! temporary cartesian location
integer II,JJ,KK ! test indices
integer prevNuc(1), nextNuc(1)
integer, parameter :: nNucs = nint((WLC_P__L/WLC_P__LENGTH_PER_BP-WLC_P__LL)/(147+WLC_P__LL)) ! assuming all octasomes
integer nucArray(nNucs)
integer loc(1), newloc(1)

! find nucs
loc = 1
do KK = 1, nNucs
    newloc = findloc(wlc_nucleosomeWrap(loc(1):),147)
    nucArray(KK) = loc(1)+newloc(1)-1
    loc=loc+newloc
enddo

! select nuc to move
call random_number(urand,rand_stat)
I = nucArray(ceiling(nNucs*(urand(1))))

! select distance to move (in bp)
call random_number(urand,rand_stat)
DR = exponential_random_int(window, rand_stat)*((-1)**nint(urand(1)))
if (DR == 0) then 
    return
endif

! find neighboring nuclesomes
prevNuc = findloc(nucArray, I) - 1
if (prevNuc(1) > 0) then
    prevNuc = nucArray(prevNuc(1))
endif
nextNuc = findloc(nucArray, I) + 1
if (nextNuc(1) < nNucs+1) then 
    nextNuc = nucArray(nextNuc(1))
else
    nextNuc = WLC_P__NB
endif

! shorten distance between beads
success = .false.
outer: do II = 1, I-(prevNuc(1)+1) ! explore the previous linker space
    inner: do JJ = 1, (nextNuc(1)-1)-I ! explore the next linker space
        if ((wlc_basepairs(I-II) + DR > 1) .AND. (wlc_basepairs(I-II) + DR < 2*(sum(wlc_basepairs)/WLC_P__NB)) .AND. &
        (wlc_basepairs(I+JJ) - DR > 1) .AND. (wlc_basepairs(I+JJ) - DR < 2*(sum(wlc_basepairs)/WLC_P__NB)) ) then
            wlc_basepairs(I-II) = wlc_basepairs(I-II) + DR
            wlc_basepairs(I+JJ) = wlc_basepairs(I+JJ) - DR
            success = .true.
            exit outer
        endif
    enddo inner
enddo outer

if (success) then
    !print*, 'move', wlc_basepairs
    ! adjust positions of beads between adjusted length beads
    if (DR > 0) then 
        IT1 = I-II+1
        IT2 = I+JJ
    else
        IT1 = I-II
        IT2 = I+JJ+1
    endif
    do KK = IT1, IT2
       !wlc_RP(:,KK) = wlc_R(:,KK) + wlc_U(:,KK)*WLC_P__LENGTH_PER_BP*DR
       wlc_RP(:,KK) = wlc_R(:,KK)
       wlc_UP(:,KK) = wlc_U(:,KK)
       if (WLC_P__LOCAL_TWIST) wlc_VP(:,I) = wlc_V(:,KK)
    enddo
    do KK = IT1,IT2
       wlc_nPointsMoved=wlc_nPointsMoved+1
       wlc_pointsMoved(wlc_nPointsMoved)=KK
    enddo
endif

end subroutine
