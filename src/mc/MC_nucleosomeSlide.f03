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
subroutine MC_nucleosomeSlide(IB1,IB2,IT1,IT2,rand_stat,success)
! values from wlcsim_data
use params, only: wlc_V, wlc_R, wlc_RP, wlc_AB, wlc_U&
    , wlc_UP, wlc_ABP, wlc_VP, wlc_pointsMoved, wlc_nPointsMoved, &
    wlc_nucleosomeWrap, wlc_basepairs, wlc_nBend, wlc_bendPoints, &
    wlc_basepairs_prop

use mersenne_twister
use params, only: dp
use windowTools, only: exponential_random_int
use polydispersity, only: get_I
use nucleosome, only: nucleosomeProp

implicit none
integer, intent(out) :: IB1   ! Test bead position 1
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Test bead position 2
integer, intent(out) :: IT2   ! Index of test bead 2
logical, intent(out) :: success

! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real(dp) urand(1) ! single random number
real(dp) :: window = 10 ! window of exponential
real(dp) :: MCAMP ! Amplitude of random change
integer DR    ! Displacement for slide move (integer BPs)
integer I ! test bead
integer II,JJ,KK,J ! test indices
real(dp) tempR(3), tempU(3), tempV(3)
integer, parameter :: max_bp = 11 ! 11 bp is maximum allowed bp size
integer prevNuc, nextNuc
integer nNucs
integer nucArray(WLC_P__NT)

! initialize
success = .false.

! find nucs
KK = 1
do II = 1, WLC_P__NT
    if (wlc_nucleosomeWrap(II)==1) cycle
    nucArray(KK) = II
    KK = KK+1
enddo
nNucs = KK-1

! select nuc to move
call random_number(urand,rand_stat)
KK = ceiling(nNucs*(urand(1)))
I = nucArray(KK)

! select distance to move (in bp)
DR = 0.0_dp
MCAMP = exponential_random_int(window,rand_stat) + 3
if (MCAMP > window) then ! 10 bp slide (~36% of moves)
    call random_number(urand,rand_stat)
    DR = 10.0_dp*(-1)**(urand(1)-0.5_dp)
else ! few bp slide
    do while (DR == 0)
        call random_number(urand,rand_stat)
        DR = MCAMP*(urand(1)-0.5_dp)
    enddo
endif
!print*, DR

! find neighboring nuclesomes
prevNuc = KK-1
if (prevNuc > 0) then
    prevNuc = nucArray(prevNuc)
else
    prevNuc = 0
endif
nextNuc = KK+1
if (nextNuc <= nNucs) then 
    nextNuc = nucArray(nextNuc)
else
    nextNuc = WLC_P__NB
endif

wlc_basepairs_prop = wlc_basepairs
! change distance between beads
outer: do II = 1, I-prevNuc ! explore the previous linker space
    inner: do JJ = 0, (nextNuc-1)-I ! explore the next linker space
        if ((wlc_basepairs(I-II) + DR > 1) .AND. (wlc_basepairs(I-II) + DR < max_bp) .AND. &
        (wlc_basepairs(I+JJ) - DR > 1) .AND. (wlc_basepairs(I+JJ) - DR < max_bp) ) then
            wlc_basepairs_prop(I-II) = wlc_basepairs(I-II) + DR
            wlc_basepairs_prop(I+JJ) = wlc_basepairs(I+JJ) - DR
            success = .true.
            exit outer
        endif
    enddo inner
enddo outer

if (success) then
    IB1 = I-II
    IB2 = I+JJ
    IT1 = IB1+1
    IT2 = IB2
    if (IB1>=1) then 
        wlc_nBend = wlc_nBend + 1
        wlc_bendPoints(wlc_nBend)=IB1
        J=IB1
        wlc_RP(:,J)=wlc_R(:,J)
        wlc_UP(:,J)=wlc_U(:,J)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,J) = wlc_V(:,J)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=J
    endif
    if (IB2<WLC_P__NT) then 
        wlc_nBend = wlc_nBend + 1
        wlc_bendPoints(wlc_nBend)=IB2
        J=IB2+1
        wlc_RP(:,J)=wlc_R(:,J)
        wlc_UP(:,J)=wlc_U(:,J)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,J) = wlc_V(:,J)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=J
    endif
    do KK = IT1, IT2
        !call nucleosomeProp(wlc_U(:,KK-1), wlc_V(:,KK-1), wlc_R(:,KK-1), &
        !                    wlc_basepairs_prop(KK-1),wlc_nucleosomeWrap(KK-1), &
        !                    tempU, tempV, tempR)
        wlc_RP(:,KK) = wlc_RP(:,KK-1) + wlc_UP(:,KK-1)*WLC_P__LENGTH_PER_BP*wlc_basepairs_prop(KK-1)
        wlc_UP(:,KK) = wlc_U(:,KK)
        wlc_VP(:,KK) = wlc_V(:,KK)
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=KK
    enddo
endif

end subroutine
