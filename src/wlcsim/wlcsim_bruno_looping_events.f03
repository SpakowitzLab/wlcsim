#include "../defines.inc"
!---------------------------------------------------------------*

subroutine wlcsim_bruno_looping_events(save_ind, wlc_p, outfile_base)
! values from wlcsim_data
use params, only: wlc_time, wlc_U, wlc_coltimes, wlc_r, wlc_TIME&
    , wlc_R

!     WLC Simulation Package:
!     Simulation Package for Brownian dynamics and
!     Monte Carlo Simulation

use params, only: dp, wlcsim_params,  pack_as_para, &
    printEnergies, printSimInfo, save_simulation_state, MAXFILENAMELEN, &
    outFileUnit

implicit none

integer, intent(in) :: save_ind
type(wlcsim_params), intent(in) :: wlc_p
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
    call save_simulation_state(-1, wlc_p, outfile_base, 'NEW')
    call InitializeEnergiesForVerifier(wlc_p)
    call MCsim(wlc_p, .false., WLC_P__NINITMCSTEPS)
    call VerifyEnergiesFromScratch(wlc_p)
    call save_simulation_state(0, wlc_p, outfile_base, 'REPLACE')
    ! for BDsim
    allocate(R0(3,WLC_P__NT))
    allocate(U0(3,WLC_P__NT))
    ! for get_looping events
    loop_file_name = trim(adjustL(outfile_base)) // 'loop_times'
    allocate(events(WLC_P__NT,WLC_P__NT))
    allocate(num_events(WLC_P__NT))
    allocate(col_state(WLC_P__NT, WLC_P__NT))
    do k2 = 1, WLC_P__NT
        do k1 = 1, k2 - 1
            if (abs(wlc_r(1,k1) - wlc_r(1,k2)) < WLC_P__COLLISIONRADIUS &
                    .and. abs(wlc_r(2,k1) - wlc_r(2,k2)) < WLC_P__COLLISIONRADIUS &
                    .and. abs(wlc_r(3,k1) - wlc_r(3,k2)) < WLC_P__COLLISIONRADIUS) then
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
TSAVE = save_ind*WLC_P__STEPSPERSAVE*wlc_p%DT
do while (wlc_time < TSAVE)
    !brown always true
    call BDsim(wlc_R, wlc_U, wlc_TIME, wlc_time + wlc_p%DT, &
            wlc_p%DT, .true., WLC_P__INTERP_BEAD_LENNARD_JONES, IDUM, pack_as_para(wlc_p), wlc_p%SIMTYPE, &
            wlc_coltimes, WLC_P__COLLISIONRADIUS, WLC_P__COLLISIONDETECTIONTYPE)
    call get_looping_events(wlc_R, WLC_P__NT, WLC_P__COLLISIONRADIUS, &
            col_state, num_events, events)
    call print_loop_events(loop_file_name, wlc_time, WLC_P__NT, events, &
            num_events)
enddo


call stress(SIG, wlc_R, wlc_U, &
            pack_as_para(wlc_p), WLC_P__INTERP_BEAD_LENNARD_JONES, wlc_p%SIMTYPE)
call stressp(COR, wlc_R, wlc_U, R0, U0, &
             pack_as_para(wlc_p), WLC_P__INTERP_BEAD_LENNARD_JONES, wlc_p%SIMTYPE)

call energy_elas(EELAS, wlc_p)
EPONP = 0.
if (WLC_P__INTERP_BEAD_LENNARD_JONES) then
    ! ring is always false for me
    call energy_self_chain(EPONP, wlc_R, WLC_P__NT, WLC_P__NB, &
                     pack_as_para(wlc_p), .FALSE.)
endif

print*, '________________________________________'
call printSimInfo(save_ind)
call printEnergies()

end subroutine wlcsim_bruno_looping_events


subroutine print_loop_events(file_name, time, nt, events, num_events)
    use params, only : dp, outFileUnit, MAXFILENAMELEN
    implicit none
    integer, intent(in) :: nt
    integer, intent(in) :: events(nt, nt), num_events(nt)
    real(dp), intent(in) :: time
    character(MAXFILENAMELEN) :: file_name

    integer :: k1, k2

    open(unit = outFileUnit, file = file_name, action = 'write', position = 'append')
    do k2 = 1, nt
        do k1 = 1, num_events(k2)
            write(outFileUnit, *) events(k1,k2), k2, time
        end do
    end do
    close(unit = outFileUnit)
end subroutine print_loop_events

