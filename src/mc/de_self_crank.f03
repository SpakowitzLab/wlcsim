!---------------------------------------------------------------*

!
!     This subroutine calculates the change in self-interaction of a DNA-like
!     molecule as a chain that interacts through a repulsive lennard-jones potential
!
!     The change in self-interaction energy associated with a crankshaft move
!     is calculated.
!
!     The program takes the terminal bead indices of the segment that is
!     rotated by the crankshaft move.
!
!
!     The interaction energy is determined used the distance of closest
!     approach between line segments. Energy of segments adjacent along the
!     chain is not calculated.


subroutine de_self_crank(DE,R,RP,NT,N,NP,PARA,RinG,IB1,IB2)
  use params, only : dp

  implicit none
  integer N,NT,NP            ! Current number of beads
  real(dp) R(3,NT)   ! Bead positions
  real(dp) RP(3,NT)  ! Test bead positions
  real(dp) DE       ! Change in self-energy
  real(dp) E        ! Self-energy before move
  real(dp) EP       ! Self-energy of test polymer
  real(dp) FMAG     ! Mag of force
  real(dp) RIJ      ! Interbead dist
  real(dp) EIJ(3)   ! Interbead unit vector
  integer I, J              ! Index holders
  integer SKIP              ! Bead skip index

  !     Variables for the calculation

  real(dp) U1(3),U2(3),U1U2
  real(dp) D1,D2
  real(dp) R12(3),D12,E12(3),R12T(3),R12C1(3),R12C2(3)
  real(dp) S1,S2
  real(dp) GI(3)
  integer I1,J1,I2,J2
  integer IB1,IB2
  integer IO,II,IOP1,IIP1
  integer DIO,DII

  !     Parameters in the simulation

  real(dp) PARA(10)
  real(dp) LHC      ! HC length
  real(dp) VHC 	! Potential strengths
  real(dp) GAM
  real(dp) LBOX     ! Box edge length
  real(dp) SUM
  real(dp) DT
  logical RinG              ! Is polymer a ring?
  integer NMAX

  real(dp) D12Min,FMAGMin


  GAM = PARA(4)
  LBOX = PARA(8)
  LHC = PARA(9)
  VHC = PARA(10)


  ! Determine number of segments inside segment moved and outside

  if (RinG) then
     if (IB2 >= IB1) then
        DII = IB2-IB1
     else
        DII = (N-IB1) + IB2
     ENDif
  else
     DII = IB2-IB1
  ENDif
  DIO = N-DII



  !///////////////////////////////////////////////////////////////////////
  !//////////////////////////////////////////////////////////////////////
  !Calculate the interaction energy for original polymer (before move)
  ! Calculate only the portion of the energy that changes during the move
  ! i.e. the interaction between the segment moved and the segment not moved
  !////////////////////////////////////////////////////////////////////////
  !///////////////////////////////////////////////////////////////////////

  E = 0.
  II = IB1
  do I = 1,DII
     if (II == N.AND.RinG) then
        IIP1 = 1
     elseif (II == N + 1.AND.RinG) then
        II = 1
        IIP1 = II + 1
     else
        IIP1 = II + 1
     ENDif
     IO = IB2
     do J = 1,DIO
        if (IO == N.AND.RinG) then
           IOP1 = 1
        elseif (IO == N + 1.AND.RinG) then
           IO = 1
           IOP1 = IO + 1
        elseif (IO == N.AND.(.not.RinG)) then
           IO = 0
           GOTO 70
        else
           IOP1 = IO + 1
        ENDif

        if (IIP1 == IO.OR.IOP1 == II) then
           GOTO 70
        ENDif
        R12(1) = R(1,IO)-R(1,II)
        R12(2) = R(2,IO)-R(2,II)
        R12(3) = R(3,IO)-R(3,II)
        ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

        D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
        ! if (D12 > (3.*GAM)) then
        !    goto 70
        ! endif

        U1(1) = R(1,IIP1)-R(1,II)
        U1(2) = R(2,IIP1)-R(2,II)
        U1(3) = R(3,IIP1)-R(3,II)
        D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
        U1(1) = U1(1)/D1
        U1(2) = U1(2)/D1
        U1(3) = U1(3)/D1

        U2(1) = R(1,IOP1)-R(1,IO)
        U2(2) = R(2,IOP1)-R(2,IO)
        U2(3) = R(3,IOP1)-R(3,IO)
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
           R12T = R(:,IOP1)-R(:,IIP1)
           R12C1 = R(:,IOP1)-R(:,II)
           R12C2 = R(:,IIP1)-R(:,IO)
           D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2,R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2,R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
           GOTO 60
        endif

        GI(1) = U2(1)-U1U2*U1(1)
        GI(2) = U2(2)-U1U2*U1(2)
        GI(3) = U2(3)-U1U2*U1(3)

        S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2 > D2.OR.S2 < 0.) then
           R12T = R(:,IOP1)-R(:,IIP1)
           R12C1 = R(:,IOP1)-R(:,II)
           R12C2 = R(:,IIP1)-R(:,IO)
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


        E = E + FMAG
70      continue
        IO = IO + 1
     ENDdo

     II = II + 1
  ENDdo

  !/////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  !Calculation the energy for the test polymer (polymer after the move)
  !Calculate only the portion of the energy that changes during the move.
  !////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  EP = 0.
  II = IB1
  do I = 1,DII
     if (II == N.AND.RinG) then
        IIP1 = 1
     elseif (II == N + 1.AND.RinG) then
        II = 1
        IIP1 = II + 1
     else
        IIP1 = II + 1
     ENDif
     IO = IB2
     do J = 1,DIO
        if (IO == N.AND.RinG) then
           IOP1 = 1
        elseif (IO == N + 1.AND.RinG) then
           IO = 1
           IOP1 = IO + 1
        elseif (IO == N.AND.(.not.RinG)) then
           IO = 0
           GOTO 90
        else
           IOP1 = IO + 1
        ENDif

        if (IIP1 == IO.OR.IOP1 == II) then
           GOTO 90
        ENDif
        R12(1) = RP(1,IO)-RP(1,II)
        R12(2) = RP(2,IO)-RP(2,II)
        R12(3) = RP(3,IO)-RP(3,II)
        ! R12(1) = R12(1)-nint(R12(1)/LBOX)*LBOX
        ! R12(2) = R12(2)-nint(R12(2)/LBOX)*LBOX
        ! R12(3) = R12(3)-nint(R12(3)/LBOX)*LBOX

        D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)
        ! if (D12 > (3.*GAM)) then
        !    goto 90
        ! endif

        U1(1) = RP(1,IIP1)-RP(1,II)
        U1(2) = RP(2,IIP1)-RP(2,II)
        U1(3) = RP(3,IIP1)-RP(3,II)
        D1 = sqrt(U1(1)**2. + U1(2)**2. + U1(3)**2.)
        U1(1) = U1(1)/D1
        U1(2) = U1(2)/D1
        U1(3) = U1(3)/D1

        U2(1) = RP(1,IOP1)-RP(1,IO)
        U2(2) = RP(2,IOP1)-RP(2,IO)
        U2(3) = RP(3,IOP1)-RP(3,IO)
        D2 = sqrt(U2(1)**2. + U2(2)**2. + U2(3)**2.)
        U2(1) = U2(1)/D2
        U2(2) = U2(2)/D2
        U2(3) = U2(3)/D2

        U1U2 = U1(1)*U2(1) + U1(2)*U2(2) + U1(3)*U2(3)
        if (U1U2 == 1..OR.U1U2 == -1.) then
           D12 = SQRT(R12(1)**2 + R12(2)**2 + R12(3)**1)
           GOTO 80
        endif

        GI(1) = U1(1)-U1U2*U2(1)
        GI(2) = U1(2)-U1U2*U2(2)
        GI(3) = U1(3)-U1U2*U2(3)

        S1 = (R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

        if (S1 > D1.OR.S1 < 0.) then
           R12T = RP(:,IOP1)-RP(:,IIP1)
           R12C1 = RP(:,IOP1)-RP(:,II)
           R12C2 = RP(:,IIP1)-RP(:,IO)
           D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2,R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2,R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
           GOTO 80
        endif

        GI(1) = U2(1)-U1U2*U1(1)
        GI(2) = U2(2)-U1U2*U1(2)
        GI(3) = U2(3)-U1U2*U1(3)

        S2 = -(R12(1)*GI(1) + R12(2)*GI(2) + R12(3)*GI(3))/(1.-U1U2**2.)

        if (S2 > D2.OR.S2 < 0.) then
           R12T = RP(:,IOP1)-RP(:,IIP1)
           R12C1 = RP(:,IOP1)-RP(:,II)
           R12C2 = RP(:,IIP1)-RP(:,IO)
           D12 = SQRT(Min(R12(1)**2 + R12(2)**2 + R12(3)**2,R12T(1)**2 + R12T(2)**2 + R12T(3)**2,&
                & R12C1(1)**2 + R12C1(2)**2 + R12C1(3)**2,R12C2(1)**2 + R12C2(2)**2 + R12C2(3)**2))
           GOTO 80
        endif

        R12(1) = R12(1) + S2*U2(1)-S1*U1(1)
        R12(2) = R12(2) + S2*U2(2)-S1*U1(2)
        R12(3) = R12(3) + S2*U2(3)-S1*U1(3)

        D12 = sqrt(R12(1)**2. + R12(2)**2. + R12(3)**2.)

80      if (D12 > LHC) then
           goto 90
        endif

        FMAG = VHC*((LHC/D12)**12.-2.*(LHC/D12)**6. + 1.)/12.


        EP = EP + FMAG
90      continue
        IO = IO + 1
     ENDdo

     II = II + 1
  ENDdo


  !Get the energy difference

  DE = EP-E


  RETURN
ENDsubroutine de_self_crank

!---------------------------------------------------------------*

