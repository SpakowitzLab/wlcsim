#include "../defines.inc"
!---------------------------------------------------------------*

subroutine wlcsim_bruno(save_ind, wlc_p)
! values from wlcsim_data
   use params, only: wlc_R, wlc_TIME, wlc_coltimes, wlc_U

!     WLC Simulation Package:
!     Simulation Package for Brownian dynamics and
!     Monte Carlo Simulation

   use params, only: dp, wlcsim_params, pack_as_para, &
                     printEnergies, printSimInfo

   implicit none

   integer, intent(in) :: save_ind
   type(wlcsim_params), intent(in) :: wlc_p

   real(dp), save, allocatable, dimension(:, :):: R0 ! Conformation of polymer chains
   real(dp), save, allocatable, dimension(:, :):: U0 ! Conformation of polymer chains

   real(dp) TSAVE    ! Time of save point

!     Energy variables

   real(dp) EELAS(3) ! Elastic energy
   real(dp) EPONP    ! Poly-poly energy

!     Structure analysis

   real(dp) SIG(3, 3)
   real(dp) COR

!     Variables for the random number generators

   integer IDUM              ! Seed for the generator

   integer NUM_POSSIBLE_COLLISIONS
   integer NT_temp ! for avoiding overflow on compile

! Exit early if all first passage times have been recorded and the relevant flag is set
   if (WLC_P__COLLISIONDETECTIONTYPE /= 0 .AND. WLC_P__EXITWHENCOLLIDED) then
      NT_temp = WLC_P__NT
      NUM_POSSIBLE_COLLISIONS = WLC_P__NT*NT_temp - WLC_P__NT
      if (COUNT(wlc_coltimes /= -1.0d+0) == NUM_POSSIBLE_COLLISIONS) then
         ! we've already exited this function previously, giving us the
         ! opportunity to save, and are now reentering it, so we can just quit
         stop
      endif
   endif

   if (save_ind == 1) then
      ! perform initialization mc if applicable
      !brown always true
      call initialize_energies_from_scratch(wlc_p)
      call mcsim(wlc_p, WLC_P__NINITMCSTEPS)
      call verify_energies_from_scratch(wlc_p)
      allocate (R0(3, WLC_P__NT))
      allocate (U0(3, WLC_P__NT))
   endif

! run brownian dynamics until next save point
! when called for save point save_ind, we should run from times
! (save_ind-1)*stepsPerSave*dt to save_ind*stepsPerSave*dt
! so that for save_ind 1, we run from time = 0 to time = timePerSave
   TSAVE = save_ind*WLC_P__STEPSPERSAVE*wlc_p%DT
!brown always true
   call BDsim(wlc_R, wlc_U, wlc_TIME, TSAVE, &
              wlc_p%DT, .true., WLC_P__INTERP_BEAD_LENNARD_JONES, IDUM, pack_as_para(wlc_p), wlc_p%SIMTYPE, &
              wlc_coltimes, WLC_P__COLLISIONRADIUS, WLC_P__COLLISIONDETECTIONTYPE)

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

   print *, '________________________________________'
   call printSimInfo(save_ind)
   call printEnergies()

end

!---------------------------------------------------------------*

