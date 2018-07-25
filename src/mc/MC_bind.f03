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

subroutine MC_bind(wlc_p,IT1,IT2,AB,ABP,METH,DEBind,dx_mu,demu)
use params, only: dp,wlcsim_params
implicit none
type(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: IT1    ! Start test bead
integer, intent(in) :: IT2    ! Final test bead
integer, intent(in) :: AB(WLC_P__NT)   ! Chemical identity (a.k.a. binding state)
integer, intent(in) :: ABP(WLC_P__NT)  ! Test Chemical identity
integer, intent(in) :: METH(WLC_P__NT) ! Methalation state (unerlyin chamical type)
real(dp), intent(out) :: DEBind    ! Change in binding energy
real(dp), intent(out) :: DEMu    ! Change in chemcial potential energy
real(dp), intent(out) :: dx_mu ! -n_bound
integer I      ! Index of bead being compared

real(dp), parameter :: selfInt= WLC_P__HP1_BIND* &
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

real(dp), parameter, dimension(0:2,0:3) :: dEBind_table = &
    reshape(&
    [0.0_dp,            0.0_dp,                 0.0_dp, &
     EU_cor,            EM_cor,                 EM_cor, &
     EU_cor,            EU_cor,                 EM_cor, &
     2.0_dp*EU_cor+intr, EU_cor+EM_cor+intr, 2.0_dp*EM_cor+intr] &
     ,[3,4])

real(dp), parameter, dimension(0:3) :: dxMu_table = &
    [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]

DEBind = 0.0_dp
Dx_mu = 0.0_dp
if (WLC_P__TWO_TAIL) then
    do I = IT1,IT2,WLC_P__NBPM
        DEBind = DEBind + dEBind_table(METH(I),ABP(I)) &
                        - dEBind_table(METH(I), AB(I))
        Dx_mu = Dx_mu - (dxMu_table(ABP(I)) - dxMu_table(AB(I)))
    ENDdo
else
    do I = IT1,IT2,WLC_P__NBPM
        if(METH(I) == 1) then
            DEBind = DEBind + EM_cor*real(ABP(I)-AB(I),dp)
        else
            DEBind = DEBind + EU_cor*real(ABP(I)-AB(I),dp)
        endif
        Dx_mu = Dx_mu - real(ABP(I)-AB(I),dp)
    ENDdo
endif
DEMu=Dx_mu*wlc_p%MU
END
