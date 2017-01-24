!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016

SUBROUTINE MC_int_rep(mc,md,I1,I2,forward)
use params
IMPLICIT NONE

!   iputs
TYPE(wlcsim_params), intent(in) :: mc   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: md
INTEGER, intent(in) :: I1  ! Test bead position 1
INTEGER, intent(in) :: I2  ! Test bead position 2

!   Internal variables
INTEGER I                 ! For looping over bins
INTEGER II                ! For looping over IB
INTEGER IB                ! Bead index
INTEGER rrdr ! -1 if r, 1 if r+dr
INTEGER IX(2),IY(2),IZ(2)
DOUBLE PRECISION WX(2),WY(2),WZ(2)
DOUBLE PRECISION WTOT       ! total weight ascribed to bin
DOUBLE PRECISION RBIN(3)    ! bead position
INTEGER INDBIN              ! index of bin
INTEGER ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type mc% everywhere
INTEGER NBINX(3)
DOUBLE PRECISION temp    !for speeding up code
LOGICAL, intent(in) :: forward ! move forward

NBINX=mc%NBINX

md%NPHI=0
! -------------------------------------------------------------
!
!  Calculate end beads
!
!--------------------------------------------------------------

do II=1,2
  if (II.eq.1) then
      IB=I1
      if (forward) then
          rrdr=-1
      else
          rrdr=1
      endif
  elseif (II.eq.2) then
      IB=I2
      if (forward) then
          rrdr=1
      else
          rrdr=-1
      endif
  else
      print*, "Error in MC_int_rep, II={1,2}"
      stop 1
  endif
   ! subract current and add new
   if (rrdr.eq.-1) then
       RBIN(1)=md%R(IB,1)
       RBIN(2)=md%R(IB,2)
       RBIN(3)=md%R(IB,3)
   else
       RBIN(1)=md%RP(IB,1)
       RBIN(2)=md%RP(IB,2)
       RBIN(3)=md%RP(IB,3)
   endif
   isA=md%AB(IB).eq.1
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(mc%confineType,RBIN,mc%LBOX,mc%NBINX,mc%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=md%NPHI
                do
                   if (I.eq.0) then
                      md%NPHI=md%NPHI+1
                      md%INDPHI(md%NPHI)=INDBIN
                      md%DPHIA(md%NPHI)=rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                      md%DPHIB(md%NPHI)=0.0_dp
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      md%DPHIA(I)=md%DPHIA(I)+rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                      exit
                   else
                      I=I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=md%NPHI
                do
                   if (I.eq.0) then
                      md%NPHI=md%NPHI+1
                      md%INDPHI(md%NPHI)=INDBIN
                      md%DPHIA(md%NPHI)=0.0_dp
                      md%DPHIB(md%NPHI)=rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      md%DPHIB(I)=md%DPHIB(I)+rrdr*WTOT*mc%beadVolume/md%Vol(INDBIN)
                      exit
                   else
                      I=I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!



!----------------------------------------------------------
!
!  Now do intermediate Beads
!
!-----------------------------------------------------------
do IB=I1,I2-1
   if (md%AB(IB).eq.md%AB(IB+1)) CYCLE
   if (forward) then
       RBIN(1)=md%RP(IB,1)
       RBIN(2)=md%RP(IB,2)
       RBIN(3)=md%RP(IB,3)
       isA=md%AB(IB).eq.1
   else
       RBIN(1)=md%R(IB,1)
       RBIN(2)=md%R(IB,2)
       RBIN(3)=md%R(IB,3)
       isA=md%AB(IB+1).eq.1
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(mc%confineType,RBIN,mc%LBOX,mc%NBINX,mc%dbin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=md%NPHI
                do
                   if (I.eq.0) then
                      md%NPHI=md%NPHI+1
                      md%INDPHI(md%NPHI)=INDBIN
                      temp=WTOT*mc%beadVolume/md%Vol(INDBIN)
                      md%DPHIA(md%NPHI)=temp
                      md%DPHIB(md%NPHI)=-temp
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      temp=WTOT*mc%beadVolume/md%Vol(INDBIN)
                      md%DPHIA(I)=md%DPHIA(I)+temp
                      md%DPHIB(I)=md%DPHIB(I)-temp
                      exit
                   else
                      I=I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=md%NPHI
                do
                   if (I.eq.0) then
                      md%NPHI=md%NPHI+1
                      md%INDPHI(md%NPHI)=INDBIN
                      temp=WTOT*mc%beadVolume/md%Vol(INDBIN)
                      md%DPHIA(md%NPHI)=-temp
                      md%DPHIB(md%NPHI)=temp
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      temp=WTOT*mc%beadVolume/md%Vol(INDBIN)
                      md%DPHIA(I)=md%DPHIA(I)-temp
                      md%DPHIB(I)=md%DPHIB(I)+temp
                      exit
                   else
                      I=I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
call hamiltonian(mc,md,.false.)

RETURN
END

!---------------------------------------------------------------!
