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
real(dp) :: MCAMP ! Amplitude of random change
real(dp) DR    ! Displacement for slide move (integer BPs)
real(dp), parameter :: optRatio = 0.25 ! max ratio of average discretization to allow for bp slide 
integer I ! test bead
integer II,JJ,KK,J ! test indices
real(dp) tempR(3), tempU(3), tempV(3)
integer max_bp 
integer prevNuc, nextNuc, linkerSum
integer nNucs
integer nucArray(WLC_P__NT)

! initialize
success = .false.
max_bp = 2*sum(wlc_basepairs)/WLC_P__NT

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
MCAMP = max(nint(optRatio*sum(wlc_basepairs)/WLC_P__NT),3)

!call random_number(urand,rand_stat)
!if (urand(1)>=0.5) then ! 10 bp slide (half of moves)
!    call random_number(urand,rand_stat)
!    DR = 10*(-1)**nint(urand(1))
!else ! few bp slide
    call random_number(urand,rand_stat)
    DR = MCAMP*(urand(1)-0.5_dp)
!endif

! find neighboring nuclesomes
prevNuc = KK-1
if (prevNuc > 0) then
    prevNuc = nucArray(prevNuc)
else
    prevNuc = 1
endif
nextNuc = KK+1
if (nextNuc <= nNucs) then 
    nextNuc = nucArray(nextNuc)
else
    nextNuc = WLC_P__NB
endif

wlc_basepairs_prop = wlc_basepairs
! change distance between beads
outer1: do II = 1, I-prevNuc ! explore the previous linker space
    inner1: do JJ = 0, (nextNuc-1)-I ! explore the next linker space
        if ((wlc_basepairs(I-II) + DR > 3) .AND. (wlc_basepairs(I-II) + DR < max_bp) .AND. &
        (wlc_basepairs(I+JJ) - DR > 3) .AND. (wlc_basepairs(I+JJ) - DR < max_bp) ) then
            wlc_basepairs_prop(I-II) = wlc_basepairs(I-II) + DR
            wlc_basepairs_prop(I+JJ) = wlc_basepairs(I+JJ) - DR
            success = .true.
            exit outer1
        endif
    enddo inner1
enddo outer1

! piece-wise movement of linker rather than just on one bead
linkerSum = 0
if (success .eqv. .false.) then 
    outer2: do II = 1, I-prevNuc ! explore the previous linker space
        inner2: do JJ = 0, (nextNuc-1)-I ! explore the next linker space
            if ((wlc_basepairs(I-II)+DR/(I-prevNuc) > 3) .AND. (wlc_basepairs(I-II)+DR/(I-prevNuc) < max_bp) .AND. &
            (wlc_basepairs(I+JJ)-DR/(I-prevNuc) > 3) .AND. (wlc_basepairs(I+JJ)-DR/(I-prevNuc) < max_bp) ) then
                wlc_basepairs_prop(I-II) = wlc_basepairs(I-II) + DR/(I-prevNuc)
                wlc_basepairs_prop(I+JJ) = wlc_basepairs(I+JJ) - DR/(I-prevNuc)
                linkerSum = linkerSum + abs(DR/(I-prevNuc))
                if (linkerSum>=10) then 
                    success = .true.
                    !print*, 'WOO'
                    exit outer2
                endif
            endif
        enddo inner2
    enddo outer2
endif

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
        call nucleosomeProp(wlc_UP(:,KK-1), wlc_VP(:,KK-1), wlc_RP(:,KK-1), &
                            wlc_basepairs_prop(KK-1),wlc_nucleosomeWrap(KK-1), &
                            tempU, tempV, tempR)
        wlc_RP(:,KK) = tempR + tempU*WLC_P__LENGTH_PER_BP*wlc_basepairs_prop(KK)
        wlc_UP(:,KK) = tempU
        !if (wlc_nucleosomeWrap(KK) == 1) then 
        !    wlc_VP(:,KK) = tempV
        !else
            wlc_VP(:,KK) = wlc_V(:,KK) ! dont change twist since incompatible
        !endif 
        wlc_nPointsMoved=wlc_nPointsMoved+1
        wlc_pointsMoved(wlc_nPointsMoved)=KK
    enddo
else
    wlc_basepairs_prop = wlc_basepairs
endif

end subroutine
