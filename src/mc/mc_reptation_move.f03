#include "../defines.inc"

subroutine mc_reptation(IT1, IT2, IB1, IB2, rand_stat, forward, super)
!--------------------------------------------------------------*
!
!           Makes Monti Carlo reptation move
!
!    Quinn separated out this file on 8/9/17
!
!---------------------------------------------------------------

! values from wlcsim_data
   use params, only: wlc_VP, wlc_ABP, wlc_V, wlc_AB, wlc_R &
                     , wlc_RP, wlc_U, wlc_UP, wlc_nPointsMoved, wlc_pointsMoved
   use mersenne_twister
   use params, only: dp
   use vector_utils, only: random_perp, cross
   use polydispersity, only: first_bead_of_chain, last_bead_of_chain, get_IB, length_of_chain

   implicit none
   logical, intent(in) :: super
   integer, intent(out) :: IT1   ! Index of test bead 1
   integer, intent(out) :: IT2   ! Index of test bead 2
   integer, intent(out) :: IB1   ! Index of test bead 1
   integer, intent(out) :: IB2   ! Index of test bead 2

   integer IP    ! Test polymer
   integer I  ! Test indices
! Things for random number generator
   type(random_stat), intent(inout) :: rand_stat  ! status of random number generator
   real(dp) urnd(1) ! single random number
   integer irnd(1)
   real(dp) DR(3)    ! Displacement for slide move
   real(dp) Uvec(3) ! parallel component of triad
   real(dp) pDir(3) ! perp component of triad
   real(dp) tDir(3) ! twist component of triad
   real(dp) r_relative(3) ! r in new coordinate system
   real(dp) u_relative(3) ! u in new coordinate system
   real(dp) v_relative(3) ! v in new coordinate system
   logical, intent(out) :: forward

!TOdo saving RP is not actually needed, even in these cases, but Brad's code assumes that we have RP.
   if (WLC_P__RING .OR. WLC_P__INTERP_BEAD_LENNARD_JONES) then
      wlc_RP = wlc_R
      wlc_UP = wlc_U
   endif

! single bead reptation
   call random_index(WLC_P__NP, irnd, rand_stat)
   IP = irnd(1)
   IT1 = first_bead_of_chain(IP)
   IT2 = last_bead_of_chain(IP)
   IB1 = 1
   IB2 = length_of_chain(IP)
! move forward or backward
   call random_number(urnd, rand_stat)
   if (urnd(1) .lt. 0.5_dp) then
      forward = .true.
      dR = wlc_R(:, IT1 + 1) - wlc_R(:, IT1)

      Uvec = wlc_U(:, IT1)
      ! chose coordinate system
      if (WLC_P__LOCAL_TWIST) then
         pDir = wlc_V(:, IT1)
         tDir = cross(Uvec, pDir)
      else
         call random_perp(Uvec, pDir, tDir, rand_stat)
      endif
      ! find next r and u in new coordinate system
      u_relative(1) = dot_product(Uvec, wlc_U(:, IT1 + 1))
      u_relative(2) = dot_product(pDir, wlc_U(:, IT1 + 1))
      u_relative(3) = dot_product(tDir, wlc_U(:, IT1 + 1))
      if (WLC_P__LOCAL_TWIST) then
         v_relative(1) = dot_product(Uvec, wlc_V(:, IT1 + 1))
         v_relative(2) = dot_product(pDir, wlc_V(:, IT1 + 1))
         v_relative(3) = dot_product(tDir, wlc_V(:, IT1 + 1))
      endif
      r_relative(1) = dot_product(Uvec, dR)
      r_relative(2) = dot_product(pDir, dR)
      r_relative(3) = dot_product(tDir, dR)

      ! orient coordinate system with end of chain
      Uvec = wlc_U(:, IT2)
      if (WLC_P__LOCAL_TWIST) then
         pDir = wlc_V(:, IT2)
         tDir = cross(Uvec, pDir)
      else
         call random_perp(Uvec, pDir, tDir, rand_stat)
      endif
      ! update UP and RP
      wlc_UP(:, IT2) = Uvec*u_relative(1) + pDir*u_relative(2) + tDir*u_relative(3)
      wlc_UP(:, IT2) = wlc_UP(:, IT2)/norm2(wlc_UP(:, IT2))
      if (WLC_P__LOCAL_TWIST) then
         wlc_VP(:, IT2) = Uvec*v_relative(1) + pDir*v_relative(2) + tDir*v_relative(3)
         wlc_VP(:, IT2) = wlc_VP(:, IT2)/norm2(wlc_VP(:, IT2))
      endif
      wlc_RP(:, IT2) = wlc_R(:, IT2) + Uvec*r_relative(1) + pDir*r_relative(2) + tDir*r_relative(3)

      do I = IT1, IT2 - 1
         wlc_RP(:, I) = wlc_R(:, I + 1)
         wlc_UP(:, I) = wlc_U(:, I + 1)
         if (WLC_P__LOCAL_TWIST) wlc_VP(:, I) = wlc_V(:, I + 1)
      enddo
      if (super) then
         do I = IT1, IT2 - 1
            wlc_ABP(I) = wlc_AB(I + 1)
         enddo
         wlc_ABP(IT2) = wlc_AB(IT1)
      endif

      ! RperpMag = sqrt(r_relative(2)**2 + r_relative(3)**2)
      ! RparaMag = r_relative(1)
      ! call test_equiv_forward(U,R,UP,RP,WLC_P__NT,IT1,IT2,RparaMag,RperpMag)

   else
      forward = .false.
      dR = wlc_R(:, IT2) - wlc_R(:, IT2 - 1)

      Uvec = wlc_U(:, IT2)
      ! chose coordinate system
      if (WLC_P__LOCAL_TWIST) then
         pDir = wlc_V(:, IT2)
         tDir = cross(Uvec, pDir)
      else
         call random_perp(Uvec, pDir, tDir, rand_stat)
      endif
      ! find next r and u in new coordinate system
      u_relative(1) = dot_product(Uvec, wlc_U(:, IT2 - 1))
      u_relative(2) = dot_product(pDir, wlc_U(:, IT2 - 1))
      u_relative(3) = dot_product(tDir, wlc_U(:, IT2 - 1))
      if (WLC_P__LOCAL_TWIST) then
         v_relative(1) = dot_product(Uvec, wlc_V(:, IT2 - 1))
         v_relative(2) = dot_product(pDir, wlc_V(:, IT2 - 1))
         v_relative(3) = dot_product(tDir, wlc_V(:, IT2 - 1))
      endif
      r_relative(1) = dot_product(Uvec, dR)
      r_relative(2) = dot_product(pDir, dR)
      r_relative(3) = dot_product(tDir, dR)

      ! orient coordinate system with end of chain
      Uvec = wlc_U(:, IT1)
      if (WLC_P__LOCAL_TWIST) then
         pDir = wlc_V(:, IT2)
         tDir = cross(Uvec, pDir)
      else
         call random_perp(Uvec, pDir, tDir, rand_stat)
      endif
      ! update UP and RP
      wlc_UP(:, IT1) = Uvec*u_relative(1) + pDir*u_relative(2) + tDir*u_relative(3)
      wlc_UP(:, IT1) = wlc_UP(:, IT1)/norm2(wlc_UP(:, IT1))
      if (WLC_P__LOCAL_TWIST) then
         wlc_VP(:, IT1) = Uvec*v_relative(1) + pDir*v_relative(2) + tDir*v_relative(3)
         wlc_VP(:, IT1) = wlc_VP(:, IT1)/norm2(wlc_VP(:, IT1))
      endif
      wlc_RP(:, IT1) = wlc_R(:, IT1) - Uvec(:)*r_relative(1) - pDir(:)*r_relative(2) - tDir(:)*r_relative(3)

      do I = IT1 + 1, IT2
         wlc_RP(:, I) = wlc_R(:, I - 1)
         wlc_UP(:, I) = wlc_U(:, I - 1)
         if (WLC_P__LOCAL_TWIST) wlc_VP(:, I) = wlc_V(:, I - 1)
      enddo
      if (super) then
         do I = IT1 + 1, IT2
            wlc_ABP(I) = wlc_AB(I - 1)
         enddo
         wlc_ABP(IT1) = wlc_AB(IT2)
      endif
   endif
   do I = IT1, IT2
      wlc_nPointsMoved = wlc_nPointsMoved + 1
      wlc_pointsMoved(wlc_nPointsMoved) = I
   enddo
end subroutine
