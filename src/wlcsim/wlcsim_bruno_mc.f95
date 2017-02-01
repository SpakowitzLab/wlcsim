!---------------------------------------------------------------*

subroutine wlcsim_bruno_mc(save_ind, wlc_d, wlc_p)

!     WLC Simulation Package:
!     Simulation Package for Brownian dynamics and
!     Monte Carlo Simulation

use params, only: dp, wlcsim_params, wlcsim_data, pack_as_para, &
    printEnergies, printSimInfo

implicit none

integer, intent(in) :: save_ind
type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d

real(dp), save, allocatable, dimension(:,:):: R0 ! Conformation of polymer chains
real(dp), save, allocatable, dimension(:,:):: U0 ! Conformation of polymer chains

!     Energy variables

real(dp) EELAS(3) ! Elastic energy
real(dp) EPONP    ! Poly-poly energy

!     Structure analysis

real(dp) SIG(3,3)
real(dp) COR

integer NUM_POSSIBLE_COLLISIONS

! Exit early if all first passage times have been recorded and the relevant flag is set
if (wlc_p%collisionDetectionType /= 0 .AND. wlc_p%exitWhenCollided) then
    NUM_POSSIBLE_COLLISIONS = wlc_p%NT*wlc_p%NT - wlc_p%NT
    if (COUNT(wlc_d%coltimes.NE.-1.0d+0) == NUM_POSSIBLE_COLLISIONS) then
        ! we've already exited this function previously, giving us the
        ! opportunity to save, and are now reentering it, so we can just quit
        stop
    endif
endif

if (save_ind == 1) then
    ! perform initialization mc if applicable
    !brown always true
    call InitializeEnergiesForVerifier(wlc_p, wlc_d)
    allocate(R0(wlc_p%NT,3))
    allocate(U0(wlc_p%NT,3))
endif

call MCsim(wlc_p, wlc_d, wlc_p%stepsPerSave)

call VerifyEnergiesFromScratch(wlc_p, wlc_d)

call stress(SIG, wlc_d%R, wlc_d%U, wlc_p%NT, wlc_p%NB, wlc_p%NP, &
            pack_as_para(wlc_p), wlc_p%INTERP_BEAD_LENNARD_JONES, wlc_p%SIMtype)
call stressp(COR, wlc_d%R, wlc_d%U, R0, U0, wlc_p%NT, wlc_p%NB, &
             wlc_p%NP, pack_as_para(wlc_p), wlc_p%INTERP_BEAD_LENNARD_JONES, wlc_p%SIMtype)

call energy_elas(EELAS, wlc_d%R, wlc_d%U, wlc_p%NT, wlc_p%NB, &
                 wlc_p%NP, pack_as_para(wlc_p), wlc_p%ring, wlc_p%twist, &
                 wlc_p%lk, wlc_p%lt, wlc_p%l)
EPONP=0.
if (wlc_p%INTERP_BEAD_LENNARD_JONES) then
    ! ring is always false for me
    call energy_self_chain(EPONP, wlc_d%R, wlc_p%NT, wlc_p%NB, &
                     pack_as_para(wlc_p), .FALSE.)
endif

print*, '________________________________________'
call printSimInfo(save_ind, wlc_p, wlc_d)
call printEnergies(wlc_d)

end


!---------------------------------------------------------------*

