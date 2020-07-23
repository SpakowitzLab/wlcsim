#include "../defines.inc"
pure function in_confinement(RP, NT, IT1, IT2)
   use params, only: dp, wlcsim_params
   use inputparams, only: MAXPARAMLEN

   implicit none
   Interface
      pure function inside_confinement(RR)
         use params, only: dp
         real(dp), intent(in) :: RR(3)
         logical inside_confinement
      end function
   end Interface
   integer, intent(in) :: IT1, IT2, NT
   real(dp), intent(in) :: RP(3, NT)
   logical in_confinement
   integer i

   in_confinement = .TRUE.
   do I = IT1, IT2
      in_confinement = inside_confinement(RP(:, I))
      if (in_confinement .eqv. .False.) return
   enddo
end function
pure function list_confinement()
   use params, only: dp, wlcsim_params, nan, wlc_RP, wlc_pointsMoved, wlc_nPointsMoved
   use inputparams, only: MAXPARAMLEN

   implicit none
   Interface
      pure function inside_confinement(RR)
         use params, only: dp
         real(dp), intent(in) :: RR(3)
         logical inside_confinement
      end function
   end Interface

   logical list_confinement
   integer i, j

   list_confinement = .TRUE.
   do I = 1, wlc_nPointsMoved
      J = wlc_pointsMoved(I)
      list_confinement = inside_confinement(wlc_RP(:, J))
      if (list_confinement .eqv. .False.) return
   enddo

end function
pure function inside_confinement(RR)
   use params, only: dp, wlcsim_params, nan
   use inputparams, only: MAXPARAMLEN

   implicit none

   real(dp), intent(in) :: RR(3)
   logical inside_confinement
   real(dp) rad, length, r2
   real(dp), parameter :: center(3) = [WLC_P__LBOX_X/2.0_dp, &
                                       WLC_P__LBOX_Y/2.0_dp, &
                                       WLC_P__LBOX_Z/2.0_dp]

   real(dp) centers(3)
   integer ix, iy, iz
   inside_confinement = .True.
   if (WLC_P__CONFINETYPE == 'platesInZperiodicXY') then
      ! Confinement only in the z-direction
      ! limits: 0 and LBox(3)
      if (WLC_P__BOUNDARY_TYPE == "ExtendBinsPast") then
         if ((RR(3) < 0.0_dp + WLC_P__DBIN) .or. &
             (RR(3) > WLC_P__CONFINEMENT_SLIT_WIDTH - WLC_P__DBIN)) then
            inside_confinement = .False.
         endif
      elseif (WLC_P__BOUNDARY_TYPE == "SolidEdgeBin") then
         if ((RR(3) < 0.0_dp) .or. &
             (RR(3) > WLC_P__CONFINEMENT_SLIT_WIDTH)) then
            inside_confinement = .False.
         endif
      endif
   elseif (WLC_P__CONFINETYPE == 'cube') then
      if (WLC_P__BOUNDARY_TYPE == "ExtendBinsPast") then
         if ((RR(1) < 0.0_dp + WLC_P__DBIN) &
             .or. (RR(1) > WLC_P__CONFINEMENT_CUBE_LENGTH - WLC_P__DBIN) &
             .or. (RR(2) < 0.0_dp + WLC_P__DBIN) &
             .or. (RR(2) > WLC_P__CONFINEMENT_CUBE_LENGTH - WLC_P__DBIN) &
             .or. (RR(3) < 0.0_dp + WLC_P__DBIN) &
             .or. (RR(3) > WLC_P__CONFINEMENT_CUBE_LENGTH - WLC_P__DBIN)) then
            inside_confinement = .False.
         endif
      elseif (WLC_P__BOUNDARY_TYPE == "SolidEdgeBin") then
         if ((RR(1) < 0.0_dp) &
             .or. (RR(1) > WLC_P__CONFINEMENT_CUBE_LENGTH) &
             .or. (RR(2) < 0.0_dp) &
             .or. (RR(2) > WLC_P__CONFINEMENT_CUBE_LENGTH) &
             .or. (RR(3) < 0.0_dp) &
             .or. (RR(3) > WLC_P__CONFINEMENT_CUBE_LENGTH)) then
            inside_confinement = .False.
         endif
      endif
   elseif (WLC_P__CONFINETYPE == 'sphere') then
      ! sphere with given diameter
      rad = (WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2
      if ((RR(1) - center(1))**2 + (RR(2) - center(2))**2 + &
          (RR(3) - center(3))**2 > rad) then
         inside_confinement = .False.
      endif
   elseif (WLC_P__CONFINETYPE == 'excludedShpereInPeriodic') then
      ! Periodic boundary conditions with an excluded sphere
      do ix = 1, WLC_P__N_SPHERES_TO_SIDE
      do iy = 1, WLC_P__N_SPHERES_TO_SIDE
      do iz = 1, WLC_P__N_SPHERES_TO_SIDE
         centers(1) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_X/real(WLC_P__N_SPHERES_TO_SIDE, dp)
         centers(2) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Y/real(WLC_P__N_SPHERES_TO_SIDE, dp)
         centers(3) = (real(ix, dp) - 0.5_DP)*WLC_P__LBOX_Z/real(WLC_P__N_SPHERES_TO_SIDE, dp)
         rad = (WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp)**2
         if ((modulo(RR(1), WLC_P__LBOX_X) - centers(1))**2 + &
             (modulo(RR(2), WLC_P__LBOX_Y) - centers(2))**2 + &
             (modulo(RR(3), WLC_P__LBOX_Z) - centers(3))**2 < rad) then
            inside_confinement = .False.
         endif
      enddo
      enddo
      enddo

   elseif (WLC_P__CONFINETYPE == 'ecoli') then
      ! cylinder with hemispherical caps, one tip at origin
      ! full length - lbox(1)/CONFINEMENT_ECOLI_LENGTH
      ! diameter - lbox(2:3)/CONFINEMENT_ECOLI_DIAMETER
      length = WLC_P__CONFINEMENT_ECOLI_LENGTH
      rad = WLC_P__CONFINEMENT_ECOLI_DIAMETER/2
      r2 = RR(2)**2 + RR(3)**2
      if (r2 > rad &
          .or. RR(1) > length &
          .or. RR(1) < 0.0_dp) then
         inside_confinement = .False.
         return
      elseif (RR(1) >= rad .and. RR(1) <= length - rad) then
         ! we're inside the main cylinder section, no need to check
         ! intersection with caps
         return
         ! if we are inside cap touching origin
      elseif (RR(1) < rad) then
         r2 = r2 + (RR(1) - rad)**2
         ! we are inside cap far from origin
      else
         r2 = r2 + (RR(1) - length + rad)**2
      endif
      if (r2 > rad*rad) then
         inside_confinement = .False.
      endif
      ! always inside confinement otherwise
      ! print*, "Undefined comfone Type"
      ! stop 1
   endif
   return

end function

