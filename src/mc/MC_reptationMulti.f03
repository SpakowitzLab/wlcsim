#include "../defines.inc"
!--------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! variables that need to be allocated only on certain branches moved into MD to prevent segfaults
! please move other variables in as you see fit
subroutine MC_reptationMulti(IT1,IT2,IB1,IB2,rand_stat,forward,super,nsteps)
! values from wlcsim_data
use params, only: wlc_U, wlc_RP, wlc_R, wlc_V, wlc_UP&
    , wlc_ABP, wlc_AB, wlc_VP
use mersenne_twister
use params, only: dp
use vector_utils, only: random_perp, cross

implicit none
integer, intent(in) :: nsteps
logical, intent(in) :: super
integer, intent(out) :: IT1   ! Index of test bead 1
integer, intent(out) :: IT2   ! Index of test bead 2
integer, intent(out) :: IB1   ! Index of test bead 1
integer, intent(out) :: IB2   ! Index of test bead 2

integer IP    ! Test polymer
! Things for random number generator
type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
real urnd(1) ! single random number
integer irnd(1)
real(dp) DR(3)    ! Displacement for slide move
real(dp) Uvec(3) ! parallel component of triad
real(dp) pDir(3) ! perp component of triad
real(dp) tDir(3) ! twist component of triad
real(dp) r_relative(3) ! r in new coordinate system
real(dp) u_relative(3) ! u in new coordinate system
real(dp) v_relative(3) ! v in new coordinate system
logical, intent(out) :: forward
integer beadi
real(dp), allocatable :: Rtemp(:,:)
real(dp), allocatable :: Utemp(:,:)
real(dp), allocatable :: Vtemp(:,:)

print*, "Proceed with caution.  The multirep move has never been tested"
stop

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
    wlc_RP = wlc_R
    wlc_UP = wlc_U
endif

allocate( Rtemp(3,nsteps) )
allocate( Utemp(3,nsteps) )
if (WLC_P__LOCAL_TWIST) allocate( Utemp(3,nsteps) )
! single bead reptation
call random_index(WLC_P__NP,irnd,rand_stat)
IP=irnd(1)
IT1 = WLC_P__NB*(IP-1) + 1
IT2 = WLC_P__NB*(IP-1) + WLC_P__NB
IB1 = 1
IB2 = WLC_P__NB
Rtemp=wlc_R(:,IT1:IT2)
Utemp=wlc_U(:,IT1:IT2)
Vtemp=wlc_V(:,IT1:IT2)
! move forward or backward
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    forward = .true.
    ! note, this is not the most efficient way of doing things
    ! it would be more efficient to move all nbpm at a time.
    do beadi =1,nsteps
        dR = Rtemp(:,2)-Rtemp(:,1)

        Uvec = Utemp(:,1)
        ! chose coordinate system
        if (WLC_P__LOCAL_TWIST) then
            pDir = Vtemp(:,1)
            tDir = cross(Uvec,pDir)
        else
            call random_perp(Uvec,pDir,tDir,rand_stat)
        endif
        ! find next r and u in new coordinate system
        u_relative(1) = dot_product(Uvec,Utemp(:,2))
        u_relative(2) = dot_product(pDir,Utemp(:,2))
        u_relative(3) = dot_product(tDir,Utemp(:,2))
        if (WLC_P__LOCAL_TWIST) then
            v_relative(1) = dot_product(Uvec,Vtemp(:,2))
            v_relative(2) = dot_product(pDir,Vtemp(:,2))
            v_relative(3) = dot_product(tDir,Vtemp(:,2))
        endif
        r_relative(1) = dot_product(Uvec,dR)
        r_relative(2) = dot_product(pDir,dR)
        r_relative(3) = dot_product(tDir,dR)

        ! orient coordinate system with end of chain
        Uvec = Utemp(:,nsteps)
        if (WLC_P__LOCAL_TWIST) then
            pDir = Vtemp(:,nsteps)
            tDir = cross(Uvec,pDir)
        else
            call random_perp(Uvec,pDir,tDir,rand_stat)
        endif
        ! update UP and RP
        wlc_UP(:,IT2) = Uvec*u_relative(1) + pDir*u_relative(2) + tDir*u_relative(3)
        wlc_UP(:,IT2) = wlc_UP(:,IT2)/norm2(wlc_UP(:,IT2))
        if (WLC_P__LOCAL_TWIST) then
            wlc_VP(:,IT2) = Uvec*v_relative(1) + pDir*v_relative(2) + tDir*v_relative(3)
            wlc_VP(:,IT2) = wlc_VP(:,IT2)/norm2(wlc_VP(:,IT2))
        endif
        wlc_RP(:,IT2) = Rtemp(:,nsteps) + Uvec*r_relative(1) + pDir*r_relative(2) + tDir*r_relative(3)

        wlc_RP(:,IT1:IT2-1)=Rtemp(:,2:nsteps)
        wlc_UP(:,IT1:IT2-1)=Utemp(:,2:nsteps)
        if (WLC_P__LOCAL_TWIST) wlc_VP(:,IT1:IT2-1)=Vtemp(:,2:nsteps)

        Rtemp=wlc_RP(:,IT1:IT2)
        Utemp=wlc_UP(:,IT1:IT2)
        if (WLC_P__LOCAL_TWIST) Vtemp=wlc_VP(:,IT1:IT2)
    enddo
    if (super) then
        wlc_ABP(IT1:IT2-nsteps)=wlc_ABP(IT1+nsteps:IT2)
        wlc_ABP(IT2-nsteps+1:IT2)=wlc_AB(IT1:IT1+nsteps-1)
    endif
else
    forward = .false.
    do beadi =1,nsteps
        dR = Rtemp(:,nsteps)-Rtemp(:,nsteps-1)

        Uvec = Utemp(:,nsteps)
        ! chose coordinate system
        if (WLC_P__LOCAL_TWIST) then
            pDir = Vtemp(:,IT2)
            tDir = cross(Uvec,pDir)
        else
            call random_perp(Uvec,pDir,tDir,rand_stat)
        endif
        ! find next r and u in new coordinate system
        u_relative(1) = dot_product(Uvec,Utemp(:,nsteps - 1))
        u_relative(2) = dot_product(pDir,Utemp(:,nsteps - 1))
        u_relative(3) = dot_product(tDir,Utemp(:,nsteps - 1))
        if (WLC_P__LOCAL_TWIST) then
            v_relative(1) = dot_product(Uvec,Vtemp(:,nsteps - 1))
            v_relative(2) = dot_product(pDir,Vtemp(:,nsteps - 1))
            v_relative(3) = dot_product(tDir,Vtemp(:,nsteps - 1))
        endif
        r_relative(1) = dot_product(Uvec,dR)
        r_relative(2) = dot_product(pDir,dR)
        r_relative(3) = dot_product(tDir,dR)

        ! orient coordinate system with end of chain
        Uvec = Utemp(:,1)
        if (WLC_P__LOCAL_TWIST) then
            pDir = Vtemp(:,nsteps)
            tDir = cross(Uvec,pDir)
        else
            call random_perp(Uvec,pDir,tDir,rand_stat)
        endif
        ! update UP and RP
        wlc_UP(:,IT1) = Uvec*u_relative(1) + pDir*u_relative(2) + tDir*u_relative(3)
        wlc_UP(:,IT1) = wlc_UP(:,IT1)/norm2(wlc_UP(:,IT1))
        if (WLC_P__LOCAL_TWIST) then
            wlc_VP(:,IT1) = Uvec*v_relative(1) + pDir*v_relative(2) + tDir*v_relative(3)
            wlc_VP(:,IT1) = wlc_VP(:,IT1)/norm2(wlc_VP(:,IT1))
        endif
        wlc_RP(:,IT1) = Rtemp(:,1)-Uvec(:)*r_relative(1)-pDir(:)*r_relative(2)-tDir(:)*r_relative(3)

        wlc_RP(:,IT1 + 1:IT2) = Rtemp(:,1:nsteps - 1)
        wlc_UP(:,IT1 + 1:IT2) = Utemp(:,1:nsteps - 1)
        if (WLC_P__LOCAL_TWIST) then
            wlc_VP(:,IT1 + 1:IT2) = Vtemp(:,1:nsteps-1)
        endif

        Rtemp=wlc_RP(:,IT1:IT2)
        Utemp=wlc_UP(:,IT1:IT2)
        if (WLC_P__LOCAL_TWIST) Vtemp=wlc_VP(:,IT1:IT2)
    enddo
    if (super) then
        wlc_ABP(IT1+nsteps:IT2)=wlc_ABP(IT1:IT2-nsteps)
        wlc_ABP(IT1:IT1+nsteps-1)=wlc_AB(IT2-nsteps+1:IT2)
    endif
endif
deallocate(Rtemp)
deallocate(Utemp)
if (WLC_P__LOCAL_TWIST) deallocate(Vtemp)
end subroutine
