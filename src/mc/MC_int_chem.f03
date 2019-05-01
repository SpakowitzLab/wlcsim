#include "../defines.inc"
!---------------------------------------------------------------!
!
!     Written by Quinn 8/4/2017
!--------------------------------------------------------------!
subroutine MC_int_chem(wlc_p,I1,I2)
! values from wlcsim_data
use params, only: wlc_inDPHI, wlc_DPHIA, wlc_DPHIB, wlc_U, wlc_AB&
    , wlc_ABP, wlc_R, wlc_DPHI_l2, wlc_NPHI
use params, only: dp, wlcsim_params
use energies, only: energyOf, maierSaupe_
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
LOGICAL initialize   ! if true, calculate absolute energy
integer, intent(in) :: I1           ! Test bead position 1
integer, intent(in) :: I2           ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer IB                ! Bead index
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ

! Copy so I don't have to type wlc_p% everywhere
real(dp) contribution
integer m_index ! m from spherical harmonics
real(dp), dimension(-2:2) :: phi2
real(dp) AminusB ! +1 if A and -1 if B

real(dp), parameter, dimension(0:3) :: number_bound_table = [0.0_dp, 1.0_dp, &
                                                             1.0_dp, 2.0_dp]
if (WLC_P__FIELDINTERACTIONTYPE == 'chromatin2') then
    print*, "chemical move not set up for chromatin2."
    stop
endif

wlc_NPHI = 0
do IB = I1,I2
   RBin(1) = wlc_R(1,IB)
   RBin(2) = wlc_R(2,IB)
   RBin(3) = wlc_R(3,IB)
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (energyOf(maierSaupe_)%isOn) then
       call Y2calc(wlc_U(:,IB),phi2)
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif


   AminusB = number_bound_table(wlc_ABP(IB))-&
             number_bound_table(wlc_AB(IB))

   do ISX = 1,2
      do ISY = 1,2
         do ISZ = 1,2
            WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
            inDBin = IX(ISX) + (IY(ISY)-1)*WLC_P__NBIN_X + (IZ(ISZ)-1)*WLC_P__NBIN_X*WLC_P__NBIN_Y
            contribution = AminusB*WTOT*WLC_P__BEADVOLUME/&
                              (WLC_P__DBIN**3)

            ! Generate list of which phi's change and by how much
            I = wlc_NPHI
            do
               if (I.eq.0) then
                  wlc_NPHI = wlc_NPHI + 1
                  wlc_inDPHI(wlc_NPHI) = inDBin
                  wlc_DPHIA(wlc_NPHI) = contribution
                  wlc_DPHIB(wlc_NPHI) = -1.0*contribution
                  if(energyOf(maierSaupe_)%isOn) then
                      do m_index = -2,2
                          wlc_DPHI_l2(m_index,wlc_NPHI) = &
                                     phi2(m_index)*contribution
                      enddo
                  endif
                  exit
               elseif (inDBin == wlc_inDPHI(I)) then
                  wlc_DPHIA(I) = wlc_DPHIA(I) +  contribution
                  wlc_DPHIB(I) = wlc_DPHIB(I) -  contribution
                  if(energyOf(maierSaupe_)%isOn) then
                      do m_index = -2,2
                          wlc_DPHI_l2(m_index,I) = wlc_DPHI_l2(m_index,I) + &
                                                     phi2(m_index)*contribution
                      enddo
                  endif
                  exit
               else
                  I = I-1
               endif
            enddo
         enddo
      enddo
   enddo
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
initialize = .False.
call hamiltonian(wlc_p,initialize)

RETURN
END

!---------------------------------------------------------------!
