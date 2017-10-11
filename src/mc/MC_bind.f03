!-----------------------------------------------------------!
!
!         Calculate HP1 Binding Energy
!
!            Started by Quinn 12/16/15
!
!
!  sign convention: wlc_p%EM and wlc_p%EU are more positive for favorable binding
!  Typical Values: wlc_p%EU = -1.52 and wlc_p%EM = 0.01

subroutine MC_bind(wlc_p,IT1,IT2,AB,ABP,METH,DEBind,dx_mu)
use params, only: dp,wlcsim_params
implicit none
type(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: IT1    ! Start test bead
integer, intent(in) :: IT2    ! Final test bead
integer, intent(in) :: AB(wlc_p%NT)   ! Chemical identity (a.k.a. binding state)
integer, intent(in) :: ABP(wlc_p%NT)  ! Test Chemical identity
integer, intent(in) :: METH(wlc_p%NT) ! Methalation state (unerlyin chamical type)
real(dp), intent(out) :: DEBind    ! Change in binding energy
real(dp), intent(out) :: dx_mu ! -n_bound
integer I      ! Index of bead being compared
DEBind = 0.0_dp
Dx_mu = 0.0_dp
do I = IT1,IT2,wlc_p%NBPM
    if(METH(I) == 1) then
        DEBind = DEBind + (-wlc_p%mu-wlc_p%EM)*real(ABP(I)-AB(I))
        Dx_mu = Dx_mu-real(ABP(I)-AB(I))
    else
        DEBind = DEBind + (-wlc_p%mu-wlc_p%EU)*real(ABP(I)-AB(I))
        Dx_mu = Dx_mu-real(ABP(I)-AB(I))
        !print*, 'In MC_bind wlc_p%EU:',EU,' wlc_p%EM:',EM
    endif
ENDdo

END
