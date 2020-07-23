#include "../defines.inc"
!-----------------------------------------------------------!
!
!         Calculate HP1 Binding Energy
!
!            Started by Quinn 12/16/15
!
!
!  sign convention: WLC_P__EM and WLC_P__EU are more positive for favorable binding
!  Typical Values: WLC_P__EU = -1.52 and WLC_P__EM = 0.01

subroutine mc_bind(IT1, IT2, AB, ABP, METH)
   use params, only: dp, wlcsim_params
   use energies, only: energyOf, mu_, bind_
   implicit none
   integer, intent(in) :: IT1    ! Start test bead
   integer, intent(in) :: IT2    ! Final test bead
   integer, intent(in) :: AB(WLC_P__NT)   ! Chemical identity (a.k.a. binding state)
   integer, intent(in) :: ABP(WLC_P__NT)  ! Test Chemical identity
   integer, intent(in) :: METH(WLC_P__NT) ! Methalation state (unerlyin chamical type)
   integer I      ! Index of bead being compared

   real(dp), parameter :: selfInt = WLC_P__HP1_BIND* &
                          ((WLC_P__BEADVOLUME/(WLC_P__DBIN**3))**2)
   real(dp), parameter :: EM_cor = -1.0_dp*(WLC_P__EM) - selfInt
   real(dp), parameter :: EU_cor = -1.0_dp*(WLC_P__EU) - selfInt
   real(dp), parameter :: intr = WLC_P__COUPLING_ENERGY - 2.0_dp*selfInt ! Interaction stregth between two HP1 on same bead

! index(meth,bind)
! bind 0=0,0  1=0,1  2=1,0 3=1,1
! meth 0=0,0  1=0,1  2=1,1

!       0,0          0,1          1,1
! 0,0
! 0,1
! 1,0
! 1,1

   real(dp), parameter, dimension(0:2, 0:3) :: dEBind_table = &
                                               reshape( &
                                               [0.0_dp, 0.0_dp, 0.0_dp, &
                                                EU_cor, EM_cor, EM_cor, &
                                                EU_cor, EU_cor, EM_cor, &
                                                2.0_dp*EU_cor + intr, EU_cor + EM_cor + intr, 2.0_dp*EM_cor + intr] &
                                               , [3, 4])

   real(dp), parameter, dimension(0:3) :: dxMu_table = &
                                          [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]

   if (WLC_P__TWO_TAIL) then
      do I = IT1, IT2, WLC_P__NBPM
         energyOf(bind_)%dx = energyOf(bind_)%dx + dEBind_table(METH(I), ABP(I)) &
                              - dEBind_table(METH(I), AB(I))
         energyOf(mu_)%dx = energyOf(mu_)%dx - (dxMu_table(ABP(I)) - dxMu_table(AB(I)))
      ENDdo
   else
      do I = IT1, IT2, WLC_P__NBPM
         if (METH(I) == 1) then
            energyOf(bind_)%dx = energyOf(bind_)%dx + EM_cor*real(ABP(I) - AB(I), dp)
         else
            energyOf(bind_)%dx = energyOf(bind_)%dx + EU_cor*real(ABP(I) - AB(I), dp)
         endif
         energyOf(mu_)%dx = energyOf(mu_)%dx - real(ABP(I) - AB(I), dp)
      ENDdo
   endif

END
