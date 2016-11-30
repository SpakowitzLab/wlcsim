!   Generates the initial distribution of "A"'s and "B"'s for simulations of
!   copolymers.
!
!   Quinn updated on 5/22/16 to use thread safe randum number generator
!

subroutine initchem(AB,NT,N,G,NP,FA,LAM,rand_stat)

!  use mt19937, only : grnd, init_genrand, rnorm, mt, mti
    use mersenne_twister
    use params, only: dp

  integer, intent(out) :: AB(NT)     ! Chemical identity of beads
  integer, intent(in) :: N           ! Number of monomers per polymer
  integer, intent(in) :: G           ! Number of beads per monomer
  integer, intent(in) :: NP          ! Number of polymer chains
  integer, intent(in) :: NT          ! Total number of beads

  integer I,J,K,IB
  real TEST(1)   ! changed to real by Quinn
  type(random_stat), intent(inout) ::rand_stat    ! status of random number generator
  !integer ABVAL

  real(dp), intent(in) :: FA   ! Fraction of A beads
  real(dp), intent(in) :: LAM  ! Chemical correlation parameter
  real(dp) PAA,PBB,PAB,PBA ! Chemical identity statistics


  !		Translate LAM and FA to probabilities

  PAA=FA*(1.0_dp-LAM)+LAM
  PBB=FA*(LAM-1.0_dp)+1.0_dp
  PBA=1.0_dp-PAA
  PAB=1.0_dp-PBB

  !		Determine the bead identities

  IB=1
  do 10 I=1,NP
     !TEST=grnd()
     call random_number(TEST,rand_stat)
     if (dble(TEST(1)).lt.FA) then
        AB(IB)=1
     else
        AB(IB)=0
     endif
     IB=IB+1
     do 15 K=2,G
        AB(IB)=AB(IB-1)
        IB=IB+1
15   continue

     do 20 J=2,N
        !TEST=grnd()
        call random_number(TEST,rand_stat)
        if (AB(IB-1).EQ.1) then
           if (TEST(1).LE.PAA) then
              AB(IB)=1
           else
              AB(IB)=0
           endif
        else
           if (TEST(1).LE.PAB) then
              AB(IB)=1
           else
              AB(IB)=0
           endif
        endif
        IB=IB+1

     do 30 K=2,G
        AB(IB)=AB(IB-1)
        IB=IB+1
30      continue
20      continue
10      continue

RETURN
end

!---------------------------------------------------------------*
