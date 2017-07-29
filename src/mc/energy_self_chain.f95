!---------------------------------------------------------------*

!
!     This subroutine calcualtes the self-interaction of a single
!     chain with itself. Chain may be linear or a ring.
!
!     The interaction energy is determined used the distance of closest
!     approach between line segments. Energy of segments adjacent along the
!     chain is not calculated.

!     Interaction is calculated for a repulsive lennard jones potential.
!     Arbitrary pair-wise potentials may be substituted

subroutine ENERGY_SELF_CHAin(EPONP,R,NT,NB,PARA,RinG)
    !TOdo change so works with multiple chains

    use params, only: dp

    implicit none

  integer NB,NT            ! Current number of beads
  real(dp) R(NT,3)   ! Bead positions
  real(dp) EPONP ! Self-interaction force
  real(dp) FMAG     ! Mag of force

  !     Variables for the calculation

  real(dp) U1(3),U2(3),U1U2
  real(dp) D1,D2
  real(dp) R12(3),D12,R12T(3),R12C1(3),R12C2(3)
  real(dp) S1,S2
  real(dp) GI(3)
  integer IT1,IT2,IT1P1,IT2P1

  !     Parameters in the simulation

  real(dp) PARA(10)
  real(dp) LHC      ! HC length
  real(dp) VHC      ! Potential strengths
  real(dp) GAM
  real(dp) LBOX     ! Box edge length
  real(dp) XIR
  logical RinG              ! Is polymer a ring?
  integer NMAX


  ! EB = PARA(1)
  ! EPAR = PARA(2)
  ! EPERP = PARA(3)
  GAM = PARA(4)
  ! ETA = PARA(5)
  XIR = PARA(6)
  ! XIU = PARA(7)
  LBOX = PARA(8)
  LHC = PARA(9)
  VHC = PARA(10)


  !     Calculate the self-interaction forces


  if (RinG) then
     NMAX = NB
  else
     NMAX = NB-1
  ENDif


  EPONP = 0.
  do IT1 = 3,NMAX
     if (IT1 == NB.AND.RinG) then
        IT1P1 = 1
     else
        IT1P1 = IT1 + 1
     ENDif
     do IT2 = 1,IT1-2
        if (IT2 == NB.AND.RinG) then
           IT2P1 = 1
        else
           IT2P1 = IT2 + 1
        ENDif

        if (IT1P1 == IT2.OR.IT2P1 == IT1) then
           GOTO 70
        ENDif
        R12(1) = R(IT2,1)-R(IT1,1)
        R12(2) = R(IT2,2)-R(IT1,2)
        R12(3) = R(IT2,3)-R(IT1,3)
        ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

        D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
        U1(1) = R(IT1P1,1)-R(IT1,1)
        U1(2) = R(IT1P1,2)-R(IT1,2)
        U1(3) = R(IT1P1,3)-R(IT1,3)
        D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
        U1(1) = U1(1)/D1
        U1(2) = U1(2)/D1
        U1(3) = U1(3)/D1

        U2(1) = R(IT2P1,1)-R(IT2,1)
        U2(2) = R(IT2P1,2)-R(IT2,2)
        U2(3) = R(IT2P1,3)-R(IT2,3)
        D2 = sqrt(U2(1)**2. + U2(2)**2. + U2(3)**2.)
        U2(1) = U2(1)/D2
        U2(2) = U2(2)/D2
        U2(3) = U2(3)/D2

        U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
        if (U1U2 == 1..OR.U1U2 == -1.) then
           D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
           GOTO 60
        endif

        GI(1) = U1(1)-U1U2*U2(1)
        GI(2) = U1(2)-U1U2*U2(2)
        GI(3) = U1(3)-U1U2*U2(3)

        S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

        if (S1 > D1.OR.S1 < 0.) then
           R12T = R(IT2P1,:)-R(IT1P1,:)
           R12C1 = R(IT2P1,:)-R(IT1,:)
           R12C2 = R(IT1P1,:)-R(IT2,:)
           D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2,R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2,R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
           GOTO 60
        endif

        GI(1) = U2(1)-U1U2*U1(1)
        GI(2) = U2(2)-U1U2*U1(2)
        GI(3) = U2(3)-U1U2*U1(3)

        S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2 > D2.OR.S2 < 0.) then
           R12T = R(IT2P1,:)-R(IT1P1,:)
           R12C1 = R(IT2P1,:)-R(IT1,:)
           R12C2 = R(IT1P1,:)-R(IT2,:)
           D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2,R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2,R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
           GOTO 60
        endif

        R12(1) = R12(1) + S2*U2(1)-S1*U1(1)
        R12(2) = R12(2) + S2*U2(2)-S1*U1(2)
        R12(3) = R12(3) + S2*U2(3)-S1*U1(3)

        D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)

60      if (D12 > LHC) then
           goto 70
        endif

        FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6. + 1.)/12.


        EPONP = EPONP + FMAG
70      continue
     ENDdo
  ENDdo



  RETURN
END subroutine ENERGY_SELF_CHAin

!---------------------------------------------------------------*

