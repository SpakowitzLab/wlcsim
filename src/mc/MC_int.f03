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

subroutine MC_int_initialize(wlc_p,wlc_d)
use params, only: dp, wlcsim_params, wlcsim_data
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
LOGICAL initialize   ! if true, calculate absolute energy

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
wlc_d%PHIA = 0.0_dp
wlc_d%DPHIA = 0.0_dp
wlc_d%PHIB = 0.0_dp
wlc_d%DPHIB = 0.0_dp
wlc_d%inDPHI = 0
if (wlc_p%CHI_L2_ON) then
    wlc_d%phi_l2 = 0.0_dp
    wlc_d%dphi_l2 = 0.0_dp
endif

wlc_d%NPHI = 0
do IB = 1,WLC_P__NT
   RBin(1) = wlc_d%R(1,IB)
   RBin(2) = wlc_d%R(2,IB)
   RBin(3) = wlc_d%R(3,IB)

   if (wlc_p%CHI_L2_ON .and. wlc_d%AB(IB).eq.1 ) then
       call Y2calc(wlc_d%U(:,IB),phi2)
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

            if (wlc_d%AB(IB) == 1 .or. wlc_d%AB(IB) == 2) then! A, chrystal, singally bound
                ! Set all phi values on initialize
                wlc_d%PHIA(inDBin) = wlc_d%PHIA(inDBin) + contribution
                if(wlc_p%CHI_L2_ON) then
                    do m_index = -2,2
                        wlc_d%PHI_l2(m_index,indBin) = wlc_d%PHI_l2(m_index,indBin) + phi2(m_index)*contribution
                    enddo
                endif
            else if (wlc_d%AB(IB) == 0) then
                ! Set all phi values on initialize
                wlc_d%PHIB(inDBin) = wlc_d%PHIB(inDBin) + contribution
            else if (wlc_d%AB(IB) == 3) then 
                ! phiA + phiB = phi_poly
                ! phiA = n_hp1*v/Vol**3
                wlc_d%PHIA(inDBin) = wlc_d%PHIA(inDBin) + 2.0_dp*contribution
                wlc_d%PHIB(inDBin) = wlc_d%PHIB(inDBin) - contribution
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
initialize = .True.
call hamiltonian(wlc_p,wlc_d,initialize)

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
subroutine MC_int_update(wlc_p,wlc_d,I1,I2)
use params, only: dp, wlcsim_params, wlcsim_data
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p
TYPE(wlcsim_data), intent(inout) :: wlc_d
LOGICAL initialize   ! if true, calculate absolute energy
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

NBinX = wlc_p%NBINX

wlc_d%NPHI = 0
do IB = I1,I2
  do rrdr = -1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
   else
       RBin(1) = wlc_d%RP(1,IB)
       RBin(2) = wlc_d%RP(2,IB)
       RBin(3) = wlc_d%RP(3,IB)
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
   if (wlc_p%CHI_L2_ON .and. wlc_d%AB(IB).eq.1) then
       if (rrdr == -1) then
           call Y2calc(wlc_d%U(:,IB),phi2)
       else
           call Y2calc(wlc_d%UP(:,IB),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0_dp
   endif

   if (wlc_d%AB(IB) == 1 .or. wlc_d%AB(IB) == 2) then ! A, chrystal, singally bound
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = rrdr*WTOT*WLC_P__BEADVOLUME/&
                                  (WLC_P__DBIN**3)
                I = wlc_d%NPHI
                ! Generate list of which phi's change and by how much
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = contribution
                      wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = &
                                  + phi2(m_index)*contribution
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + contribution
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,I) = wlc_d%DPHI_l2(m_index,I) &
                                  + phi2(m_index)*contribution
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
   else if (wlc_d%AB(IB) == 0) then
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = rrdr*WTOT*WLC_P__BEADVOLUME/&
                                  (WLC_P__DBIN**3)
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              ! This is somewhat wastefull, could eliminate for speedup by having another NPHI for L=2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = 0.0_dp
                          enddo
                      endif
                      wlc_d%DPHIB(wlc_d%NPHI) = contribution
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + contribution
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else if (wlc_d%AB(IB) == 3) then
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                contribution = rrdr*WTOT*WLC_P__BEADVOLUME/&
                                  (WLC_P__DBIN**3)
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = 2.0_dp*contribution
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              ! This is somewhat wastefull, could eliminate for speedup by having another NPHI for L=2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = 0.0_dp
                          enddo
                      endif
                      wlc_d%DPHIB(wlc_d%NPHI) = -1.0_dp*contribution
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + 2.0*contribution
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) - contribution
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else
       print*, "AB state",wlc_d%AB(IB)," not defined"
       stop
   endif
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
initialize = .False.
call hamiltonian(wlc_p,wlc_d,initialize)

RETURN
END

!---------------------------------------------------------------!
