#include "../defines.inc"
!-----------------------------------------------------------!
!
!         Calculate Energy dou to exteranl potential
!
!            Started by Quinn, Dec 2017
!     Sets: dEField and dx_Field based on R, RP, and, HA
!-----------------------------------------------------------

function NSTCS_plus_binding(x, ii) result(energy)
   ! Non-Specific to cuge side plus binding
   use params, only: dp
   implicit none
   real(dp), intent(in) :: x
   integer, intent(in) :: ii
   real(dp) energy
   real(dp) offset
   energy = 0.0_dp
   if (x < WLC_P__BINDING_R) then
      energy = energy + 1.0_dp
   endif
   if (WLC_P__BOUNDARY_TYPE == 'ExtnedBinsPast') then
      offset = WLC_P__DBIN
   else
      offset = 0.0_dp
   endif
   if (WLC_P__EXTERNAL_FIELD_TYPE == 'NSTCS_plus_binding') then
      if (MOD(ii, 3000) == 0 .and. x > WLC_P__BINDING_R + offset) then
         energy = energy - 100_dp
         energy = energy - 10.0_dp*x
      endif
   endif

end function
function random_to_cube(x, ii) result(energy)
   ! Non-Specific to cuge side plus binding
   use params, only: dp, wlc_external_bind_points
   implicit none
   real(dp), intent(in) :: x
   integer, intent(in) :: ii
   real(dp) energy
   real(dp) offset
   energy = 0.0_dp
   if (WLC_P__BOUNDARY_TYPE == 'ExtnedBinsPast') then
      offset = WLC_P__DBIN
   else
      offset = 0.0_dp
   endif
   if (WLC_P__EXTERNAL_FIELD_TYPE == 'Random_to_cube_side') then
      if (wlc_external_bind_points(ii) .and. x > WLC_P__BINDING_R + offset) then
         energy = energy - 100_dp
         energy = energy - 10.0_dp*x
      endif
   endif

end function

subroutine mc_external_field()
! values from wlcsim_data
   use params, only: wlc_RP, wlc_R, wlc_nPointsMoved, wlc_pointsMoved, nan
   use params, only: dp
   use energies, only: energyOf, external_
   implicit none
   integer ii, jj
   real(dp) vv(3)
   real(dp), parameter :: center(3) = [WLC_P__LBOX_X/2.0_dp, &
                                       WLC_P__LBOX_Y/2.0_dp, &
                                       WLC_P__LBOX_Z/2.0_dp]
   real(dp) centers(3)
   integer ix, iy, iz
   real(dp) NSTCS_plus_binding, random_to_cube

   do jj = 1, wlc_nPointsMoved
      ii = wlc_pointsMoved(jj)
      if (WLC_P__EXTERNAL_FIELD_TYPE == 'nonSpecificToCubeSide') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  NSTCS_plus_binding(wlc_RP(1, ii), ii)

         energyOf(external_)%dx = energyOf(external_)%dx - &
                                  NSTCS_plus_binding(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'NSTCS_plus_binding') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  NSTCS_plus_binding(wlc_RP(1, ii), ii)

         energyOf(external_)%dx = energyOf(external_)%dx - &
                                  NSTCS_plus_binding(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'Random_to_cube_side') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  random_to_cube(wlc_RP(1, ii), ii)

         energyOf(external_)%dx = energyOf(external_)%dx - &
                                  random_to_cube(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'toExcludedSphereInPeriodic') then
         do ix = 1, WLC_P__N_SPHERES_TO_SIDE
         do iy = 1, WLC_P__N_SPHERES_TO_SIDE
         do iz = 1, WLC_P__N_SPHERES_TO_SIDE
            centers(1) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_X/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            centers(2) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Y/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            centers(3) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Z/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            vv(1) = modulo(wlc_RP(1, ii), WLC_P__LBOX_X) - centers(1)
            vv(2) = modulo(wlc_RP(2, ii), WLC_P__LBOX_Y) - centers(2)
            vv(3) = modulo(wlc_RP(3, ii), WLC_P__LBOX_Z) - centers(3)
            if (dot_product(vv, vv) < &
                (WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
               energyOf(external_)%dx = energyOf(external_)%dx + 1.0_dp
            endif
            vv(1) = modulo(wlc_R(1, ii), WLC_P__LBOX_X) - centers(1)
            vv(2) = modulo(wlc_R(2, ii), WLC_P__LBOX_Y) - centers(2)
            vv(3) = modulo(wlc_R(3, ii), WLC_P__LBOX_Z) - centers(3)
            if (dot_product(vv, vv) < &
                (WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
               energyOf(external_)%dx = energyOf(external_)%dx - 1.0_dp
            endif
         enddo
         enddo
         enddo
      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'nonSpecificToSphere') then
         vv(1) = modulo(wlc_RP(1, ii), WLC_P__LBOX_X) - center(1)
         vv(2) = modulo(wlc_RP(2, ii), WLC_P__LBOX_Y) - center(2)
         vv(3) = modulo(wlc_RP(3, ii), WLC_P__LBOX_Z) - center(3)
         if (dot_product(vv, vv) > &
             (-WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
            energyOf(external_)%dx = energyOf(external_)%dx + 1.0_dp
         endif
         vv(1) = modulo(wlc_R(1, ii), WLC_P__LBOX_X) - center(1)
         vv(2) = modulo(wlc_R(2, ii), WLC_P__LBOX_Y) - center(2)
         vv(3) = modulo(wlc_R(3, ii), WLC_P__LBOX_Z) - center(3)
         if (dot_product(vv, vv) > &
             (-WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
            energyOf(external_)%dx = energyOf(external_)%dx - 1.0_dp
         endif

      endif
   enddo

END

subroutine mc_external_field_from_scratch()
! values from wlcsim_data
   use params, only: wlc_R
   use params, only: dp, nan
   use energies, only: energyOf, external_
   implicit none
   integer ii
   real(dp) vv(3)
   real(dp), parameter :: center(3) = [WLC_P__LBOX_X/2.0_dp, &
                                       WLC_P__LBOX_Y/2.0_dp, &
                                       WLC_P__LBOX_Z/2.0_dp]
   real(dp) centers(3)
   integer ix, iy, iz
   real(dp) NSTCS_plus_binding, random_to_cube

   do ii = 1, WLC_P__NT
      if (WLC_P__EXTERNAL_FIELD_TYPE == 'nonSpecificToCubeSide') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  NSTCS_plus_binding(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'NSTCS_plus_binding') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  NSTCS_plus_binding(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'Random_to_cube_side') then
         energyOf(external_)%dx = energyOf(external_)%dx + &
                                  random_to_cube(wlc_R(1, ii), ii)

      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'toExcludedSphereInPeriodic') then
         do ix = 1, WLC_P__N_SPHERES_TO_SIDE
         do iy = 1, WLC_P__N_SPHERES_TO_SIDE
         do iz = 1, WLC_P__N_SPHERES_TO_SIDE
            centers(1) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_X/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            centers(2) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Y/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            centers(3) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Z/real(WLC_P__N_SPHERES_TO_SIDE, dp)
            vv(1) = modulo(wlc_R(1, ii), WLC_P__LBOX_X) - centers(1)
            vv(2) = modulo(wlc_R(2, ii), WLC_P__LBOX_Y) - centers(2)
            vv(3) = modulo(wlc_R(3, ii), WLC_P__LBOX_Z) - centers(3)
            if (dot_product(vv, vv) < &
                (WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
               energyOf(external_)%dx = energyOf(external_)%dx + 1.0_dp
            endif
         enddo
         enddo
         enddo
      elseif (WLC_P__EXTERNAL_FIELD_TYPE == 'nonSpecificToSphere') then
         vv(1) = modulo(wlc_R(1, ii), WLC_P__LBOX_X) - center(1)
         vv(2) = modulo(wlc_R(2, ii), WLC_P__LBOX_Y) - center(2)
         vv(3) = modulo(wlc_R(3, ii), WLC_P__LBOX_Z) - center(3)
         if (dot_product(vv, vv) > &
             (-WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2) then
            energyOf(external_)%dx = energyOf(external_)%dx + 1.0_dp
         endif

      endif
   enddo

END
