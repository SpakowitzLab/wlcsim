!---------------------------------------------------------------*

subroutine wlcsim_bruno_looping_events(save_ind, wlc_d, wlc_p, outfile_base)

!     WLC Simulation Package:
!     Simulation Package for Brownian dynamics and
!     Monte Carlo Simulation

use params, only: dp, wlcsim_params, wlcsim_data, pack_as_para, &
    printEnergies, printSimInfo, save_simulation_state, MAXFILENAMELEN, &
    outFileUnit

implicit none

integer, intent(in) :: save_ind
type(wlcsim_params), intent(in) :: wlc_p
type(wlcsim_data), intent(inout) :: wlc_d
character(MAXFILENAMELEN) :: outfile_base
character(MAXFILENAMELEN) :: loop_file_name

real(dp), save, allocatable, dimension(:,:):: R0 ! Conformation of polymer chains
real(dp), save, allocatable, dimension(:,:):: U0 ! Conformation of polymer chains

integer, save, allocatable, dimension(:,:) :: events ! which events happened this time step
integer, save, allocatable, dimension(:) :: num_events ! to safely index into events
integer, save, allocatable, dimension(:,:) :: col_state ! 1 or 0 for looped/unlooped

real(dp) TSAVE    ! Time of save point

!     Energy variables

real(dp) EELAS(3) ! Elastic energy
real(dp) EPONP    ! Poly-poly energy

!     Structure analysis

real(dp) SIG(3,3)
real(dp) COR

!     Variables for the random number generators

integer IDUM              ! Seed for the generator

integer :: k1, k2

if (save_ind == 1) then
    ! perform initialization mc if applicable, save before and after
    !brown always true
    call save_simulation_state(-1, wlc_d, wlc_p, outfile_base, 'NEW')
    call InitializeEnergiesForVerifier(wlc_p, wlc_d)
    call MCsim(wlc_p, wlc_d, wlc_p%nInitMCSteps)
    call VerifyEnergiesFromScratch(wlc_p, wlc_d)
    call save_simulation_state(0, wlc_d, wlc_p, outfile_base, 'REPLACE')
    ! for BDsim
    allocate(R0(wlc_p%NT,3))
    allocate(U0(wlc_p%NT,3))
    ! for get_looping events
    loop_file_name = trim(adjustL(outfile_base)) // 'loop_times'
    allocate(events(wlc_p%NT,wlc_p%NT))
    allocate(num_events(wlc_p%NT))
    allocate(col_state(wlc_p%NT, wlc_p%NT))
    do k2 = 1, wlc_p%NT
        do k1 = 1, k2 - 1
            if (abs(wlc_d%r(k1,1) - wlc_d%r(k2,1)) < wlc_p%collisionRadius &
                    .and. abs(wlc_d%r(k1,2) - wlc_d%r(k2,2)) < wlc_p%collisionRadius &
                    .and. abs(wlc_d%r(k1,3) - wlc_d%r(k2,3)) < wlc_p%collisionRadius) then
                col_state(k1, k2) = 1
            else
                col_state(k1, k2) = 0
            end if
        end do
    end do
endif

! run brownian dynamics until next save point
! when called for save point save_ind, we should run from times
! (save_ind-1)*stepsPerSave*dt to save_ind*stepsPerSave*dt
! so that for save_ind 1, we run from time = 0 to time = timePerSave
TSAVE = save_ind*wlc_p%stepsPerSave*wlc_p%dt
do while (wlc_d%time < TSAVE)
    !brown always true
    call BDsim(wlc_d%R, wlc_d%U, wlc_p%NT, wlc_p%NB, wlc_p%NP, wlc_d%TIME, wlc_d%time + wlc_p%dt, &
            wlc_p%DT, .true., wlc_p%INTERP_BEAD_LENNARD_JONES, IDUM, pack_as_para(wlc_p), wlc_p%SIMtype, &
            wlc_d%coltimes, wlc_p%collisionRadius, wlc_p%collisionDetectionType)
    call get_looping_events(wlc_d%R, wlc_p%NT, wlc_p%collisionRadius, &
            col_state, num_events, events)
    call print_loop_events(loop_file_name, wlc_d%time, wlc_p%nt, events, &
            num_events)
enddo


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

end subroutine wlcsim_bruno_looping_events


subroutine print_loop_events(file_name, time, nt, events, num_events)
    use params, only : dp, outFileUnit, MAXFILENAMELEN
    implicit none
    integer, intent(in) :: nt
    integer, intent(in) :: events(nt, nt), num_events(nt)
    real(dp), intent(in) :: time
    character(MAXFILENAMELEN) :: file_name

    integer :: k1, k2

    open(unit=outFileUnit, file=file_name, action='write', position='append')
    do k2 = 1, nt
        do k1 = 1, num_events(k2)
            write(outFileUnit, *) events(k1,k2), k2, time
        end do
    end do
    close(unit=outFileUnit)
end subroutine print_loop_events

