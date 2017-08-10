!---------------------------------------------------------------!
!
!     Written by Quinn 8/4/2017
!--------------------------------------------------------------!
subroutine MC_int_chem(wlc_p,wlc_d,I1,I2)
use params, only: dp, wlcsim_params, wlcsim_data
implicit none

TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
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
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type wlc_p% everywhere
integer NBinX(3)
real(dp) temp
real(dp) contribution 
integer m_plus3 ! m from spherical harmonics
real(dp) phi2(5)
real(dp) AminusB ! +1 if A and -1 if B

NBinX = wlc_p%NBinX

wlc_d%NPHI = 0
do IB = I1,I2
   RBin(1) = wlc_d%R(1,IB)
   RBin(2) = wlc_d%R(2,IB)
   RBin(3) = wlc_d%R(3,IB)
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p%confineType,RBin,wlc_p%LBOX,wlc_p%NBinX,wlc_p%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (wlc_p%chi_l2_on) then
       call Y2calc(wlc_d%U(:,IB),phi2)
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0
   endif
   AminusB = real(wlc_d%ABP(IB)-wlc_d%AB(IB))

   if (wlc_p%confineType == 'none' .or. wlc_p%confineType == 'periodicUnequal') then
       ! If periodic than you can assume that all bins are included and have a volume
       ! of dbin**3
       temp = wlc_p%beadVolume/(wlc_p%dbin**3)
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                contribution=temp*WTOT
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = contribution*AminusB
                      wlc_d%DPHIB(wlc_d%NPHI) = -1.0*contribution*AminusB
                      if(wlc_p%chi_l2_on) then
                          do m_plus3 =1,5
                              wlc_d%DPHI_l2(m_plus3,wlc_d%NPHI) = &
                                     AminusB*phi2(m_plus3)*contribution
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) +  contribution*AminusB
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) -  contribution*AminusB
                      if(wlc_p%chi_l2_on) then
                          do m_plus3 =1,5
                              wlc_d%DPHI_l2(m_plus3,I) = wlc_d%DPHI_l2(m_plus3,I) + &
                                  AminusB*phi2(m_plus3)*contribution
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
    else ! not periodic
       ! if not periodic you need to check wheter bin is outside and
       ! also need to divide by volume of bin
       temp = wlc_p%beadVolume
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)

                contribution=temp*WTOT/wlc_d%Vol(indBin)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = contribution*AminusB
                      wlc_d%DPHIB(wlc_d%NPHI) = -1.0*contribution*AminusB
                      if(wlc_p%chi_l2_on) then
                          do m_plus3 =1,5
                              wlc_d%DPHI_l2(m_plus3,wlc_d%NPHI) =  +&
                                  AminusB*phi2(m_plus3)*contribution
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIA(wlc_d%NPHI) = wlc_d%DPHIA(wlc_d%NPHI) +  contribution*AminusB
                      wlc_d%DPHIB(wlc_d%NPHI) = wlc_d%DPHIB(wlc_d%NPHI) -  contribution*AminusB
                      if(wlc_p%chi_l2_on) then
                          do m_plus3 =1,5
                              wlc_d%DPHI_l2(m_plus3,I) = wlc_d%DPHI_l2(m_plus3,I) + &
                                  AminusB*phi2(m_plus3)*contribution
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
    endif
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
