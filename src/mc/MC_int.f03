#include "../defines.inc"
!---------------------------------------------------------------!
!
!     This subroutine calculates the field energy from scratch.
!     Also sets PhiA and PhiB to 0.
!     Set x_chi, x
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!     Edited by Quinn in 2016

subroutine MC_int_initialize(wlc_p)
! values from wlcsim_data
use params, only: wlc_PHIB, wlc_NPHI, wlc_PHI_l2, wlc_U, wlc_AB&
    , wlc_PHIA, wlc_R, wlc_phi_l2, wlc_dphi_l2, wlc_inDPHI, wlc_DPHIA&
    , wlc_DPHIB
use params, only: dp, wlcsim_params
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
LOGICAL,parameter :: initialize = .TRUE.   ! if true, calculate absolute energy

!   Internal variables
integer IB                ! Bead index
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) contribution
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
real(dp), dimension(-2:2) :: phi2
integer m_index ! m for spherical harmonics

! Copy so I don't have to type wlc_p% everywhere
integer NBinX(3)
NBinX = wlc_p%NBINX

! -------------------------------------------------------------
!
!  Calculate change (or value if initialize) of phi for A and B
!
!--------------------------------------------------------------
wlc_PHIA = 0.0_dp
wlc_DPHIA = 0.0_dp
wlc_PHIB = 0.0_dp
wlc_DPHIB = 0.0_dp
wlc_inDPHI = 0
if (wlc_p%CHI_L2_ON) then
    wlc_phi_l2 = 0.0_dp
    wlc_dphi_l2 = 0.0_dp
endif

wlc_NPHI = 0
do IB = 1,WLC_P__NT
   RBin(1) = wlc_R(1,IB)
   RBin(2) = wlc_R(2,IB)
   RBin(3) = wlc_R(3,IB)

   if (wlc_p%CHI_L2_ON .and. wlc_AB(IB).eq.1 ) then
       call Y2calc(wlc_U(:,IB),phi2)
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif


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
   do ISX = 1,2
      do ISY = 1,2
         do ISZ = 1,2
            if (((IX(ISX).le.0).OR.(IX(ISX).ge.(NBinX(1) + 1))) .or.&
                ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBinX(2) + 1))) .or.&
                ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBinX(3) + 1)))) then
                print*, "Out of Bounds!"
                stop
            endif

            WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
            inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
            contribution =  WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)

            if (wlc_AB(IB) == 1 .or. wlc_AB(IB) == 2) then! A, chrystal, singally bound
                ! Set all phi values on initialize
                wlc_PHIA(inDBin) = wlc_PHIA(inDBin) + contribution
                if(wlc_p%CHI_L2_ON) then
                    do m_index = -2,2
                        wlc_PHI_l2(m_index,indBin) = wlc_PHI_l2(m_index,indBin) + phi2(m_index)*contribution
                    enddo
                endif
            else if (wlc_AB(IB) == 0) then
                ! Set all phi values on initialize
                wlc_PHIB(inDBin) = wlc_PHIB(inDBin) + contribution
            else if (wlc_AB(IB) == 3) then
                ! phiA + phiB = phi_poly
                ! phiA = n_hp1*v/Vol**3
                wlc_PHIA(inDBin) = wlc_PHIA(inDBin) + 2.0_dp*contribution
                wlc_PHIB(inDBin) = wlc_PHIB(inDBin) - contribution
            endif
         enddo
      enddo
   enddo
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte abosolute energy
!
!---------------------------------------------------------------------
call hamiltonian(wlc_p,initialize)

RETURN
END

!---------------------------------------------------------------!
!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!     Edited by Quinn in 2016
!--------------------------------------------------------------!
subroutine MC_int_update(wlc_p,I1,I2)
! values from wlcsim_data
use params, only: wlc_NPHI, wlc_inDPHI,wlc_ind_in_list
use params, only: wlcsim_params
implicit none
TYPE(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: I1           ! Test bead position 1
integer, intent(in) :: I2           ! Test bead position 2
LOGICAL, parameter :: initialize = .False.  ! if true, calculate absolute energy
integer I,J

wlc_NPHI = 0

call CalcDphi(wlc_p,I1,I2)
do I = 1,wlc_NPHI
   J = wlc_inDPHI(I)
   wlc_ind_in_list(J) = -1
enddo
call hamiltonian(wlc_p,initialize) ! calculate change in energy based on density change
end subroutine MC_int_update

subroutine MC_int_update_spider(wlc_p,spider_id)
! values from wlcsim_data
use params, only: wlc_NPHI, wlc_spiders, wlc_inDPHI, wlc_ind_in_list
use params, only: wlcsim_params
implicit none
TYPE(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: spider_id
LOGICAL, parameter :: initialize = .False.  ! if true, calculate absolute energy
integer section_n,I1,I2,I,J

wlc_NPHI = 0
do section_n = 1, wlc_spiders(spider_id)%nSections
    I1 = wlc_spiders(spider_id)%moved_sections(1,section_n)
    I2 = wlc_spiders(spider_id)%moved_sections(2,section_n)
    call CalcDphi(wlc_p,I1,I2)
enddo

do I = 1,wlc_NPHI
   J = wlc_inDPHI(I)
   wlc_ind_in_list(J) = -1
enddo


call hamiltonian(wlc_p,initialize) ! calculate change in energy based on density change
end subroutine MC_int_update_spider

subroutine CalcDphi(wlc_p,I1,I2)
!  Note: This subroutine assumes you have set wlc_bin_in_list=FALSE
!  some time before the start of the move.

! values from wlcsim_data
use params, only: wlc_NPHI, wlc_RP, wlc_U, wlc_AB, wlc_R&
    , wlc_UP, wlc_DPHI_l2, wlc_inDPHI, wlc_DPHIA, wlc_DPHIB&
    , wlc_ind_in_list
use params, only: dp, wlcsim_params
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p
integer, intent(in) :: I1           ! Test bead position 1
integer, intent(in) :: I2           ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer IB                ! Bead index
integer rrdr ! -1 if r, 1 if r + dr
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ

! Copy so I don't have to type wlc_p% everywhere
integer NBinX(3)
integer m_index ! m for spherical harmonics
real(dp), dimension(-2:2) :: phi2
real(dp) contribution
real(dp) change

NBinX = wlc_p%NBINX

do IB = I1,I2
  do rrdr = -1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_R(1,IB)
       RBin(2) = wlc_R(2,IB)
       RBin(3) = wlc_R(3,IB)
   else
       RBin(1) = wlc_RP(1,IB)
       RBin(2) = wlc_RP(2,IB)
       RBin(3) = wlc_RP(3,IB)
   endif
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
   if (wlc_p%CHI_L2_ON .and. wlc_AB(IB).eq.1) then
       if (rrdr == -1) then
           call Y2calc(wlc_U(:,IB),phi2)
       else
           call Y2calc(wlc_UP(:,IB),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif
   change = real(rrdr,dp)*WLC_P__BEADVOLUME*(WLC_P__DBIN**3)

   if (wlc_AB(IB) == 1 .or. wlc_AB(IB) == 2) then ! A, chrystal, singally bound
       do ISZ = 1,2
          do ISY = 1,2
             do ISX = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = WTOT*change
                I = wlc_ind_in_list(indBin)
                if (I == -1) then
                    wlc_NPHI = wlc_NPHI + 1
                    wlc_ind_in_list(indBin) = wlc_NPHI
                    wlc_inDPHI(wlc_NPHI) = inDBin
                    wlc_DPHIA(wlc_NPHI) = contribution
                    wlc_DPHIB(wlc_NPHI) = 0.0_dp
                    if(wlc_p%CHI_L2_ON) then
                        do m_index = -2,2
                            wlc_DPHI_l2(m_index,wlc_NPHI) = &
                                + phi2(m_index)*contribution
                        enddo
                    endif
                else
                    wlc_DPHIA(I) = wlc_DPHIA(I) + contribution
                    if(wlc_p%CHI_L2_ON) then
                        do m_index = -2,2
                            wlc_DPHI_l2(m_index,I) = wlc_DPHI_l2(m_index,I) &
                                + phi2(m_index)*contribution
                        enddo
                    endif
                endif
             enddo
          enddo
       enddo
   else if (wlc_AB(IB) == 0) then
       do ISZ = 1,2
          do ISY = 1,2
             do ISX = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = WTOT*change
                I = wlc_ind_in_list(indBin)
                if (I == -1) then
                    wlc_NPHI = wlc_NPHI + 1
                    wlc_ind_in_list(indBin) = wlc_NPHI
                    wlc_inDPHI(wlc_NPHI) = inDBin
                    wlc_DPHIA(wlc_NPHI) = 0.0_dp
                    if(wlc_p%CHI_L2_ON) then
                        do m_index = -2,2
                            wlc_DPHI_l2(m_index,wlc_NPHI) = 0.0_dp
                        enddo
                    endif
                    wlc_DPHIB(wlc_NPHI) = contribution
                else
                    wlc_DPHIB(I) = wlc_DPHIB(I) + contribution
                endif
             enddo
          enddo
       enddo
   else if (wlc_AB(IB) == 3) then
       do ISZ = 1,2
          do ISY = 1,2
             do ISX = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = WTOT*change
                I = wlc_ind_in_list(indBin)
                if (I == -1) then
                    wlc_NPHI = wlc_NPHI + 1
                    wlc_ind_in_list(indBin) = wlc_NPHI
                    wlc_inDPHI(wlc_NPHI) = inDBin
                    wlc_DPHIA(wlc_NPHI) = 2.0_dp*contribution
                    if(wlc_p%CHI_L2_ON) then
                        do m_index = -2,2
                            wlc_DPHI_l2(m_index,wlc_NPHI) = 0.0_dp
                        enddo
                    endif
                    wlc_DPHIB(wlc_NPHI) = -1.0_dp*contribution
                else
                    wlc_DPHIA(I) = wlc_DPHIA(I) + 2.0_dp*contribution
                    wlc_DPHIB(I) = wlc_DPHIB(I) - contribution
                endif
             enddo
          enddo
       enddo
   else
       print*, "AB state",wlc_AB(IB)," not defined"
       stop
   endif
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads

RETURN
END

!---------------------------------------------------------------!
