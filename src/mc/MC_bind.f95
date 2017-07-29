!-----------------------------------------------------------!
!
!         Calculate HP1 Binding Energy
!
!            Started by Quinn 12/16/15
!
!
!  sign convention: EM and EU are more positive for favorable binding
!  Typical Values: EU = -1.52 and EM = 0.01

subroutine MC_bind(NT,BPM,IT1,IT2,AB,ABP,METH,EU,EM,DEBind,mu,dx_mu)
use params, only: dp
implicit none
integer, intent(in) :: NT     ! Total number of beads in simulation
integer, intent(in) :: BPM    ! Number of beads per monomer
integer, intent(in) :: IT1    ! Start test bead
integer, intent(in) :: IT2    ! Final test bead
integer, intent(in) :: AB(NT)   ! Chemical identity (a.k.a. binding state)
integer, intent(in) :: ABP(NT)  ! Test Chemical identity
integer, intent(in) :: METH(NT) ! Methalation state (unerlyin chamical type)
real(dp), intent(in) :: EU        ! Binding energy of Unemethalted state
real(dp), intent(in) :: EM        ! Binding energy of methalated state
real(dp), intent(in) :: mu   ! Chemical potential of HP1
real(dp), intent(out) :: DEBind    ! Change in binding energy
real(dp), intent(out) :: dx_mu ! -n_bound
integer I      ! Index of bead being compared
DEBind = 0.0_dp
Dx_mu = 0.0_dp
do I = IT1,IT2,BPM
    if(METH(I) == 1) then
        DEBind = DEBind + (-mu-EM)*real(ABP(I)-AB(I))
        Dx_mu = Dx_mu-real(ABP(I)-AB(I))
    else
        DEBind = DEBind + (-mu-EU)*real(ABP(I)-AB(I))
        Dx_mu = Dx_mu-real(ABP(I)-AB(I))
        !print*, 'In MC_bind EU:',EU,' EM:',EM
    endif
ENDdo

if (abs(DEBind).gt.100000) then

    do I = IT1,IT2,BPM
        print*, ABP(I), AB(I)
    ENDdo
    print*, "range:",IT1,IT2
    print*, "error in MC_bind"
    print*, "DEBind",DEBind
    print*, "mu", mu
    print*, "EM",EM,"   EU",EU
    stop 1
endif


END
