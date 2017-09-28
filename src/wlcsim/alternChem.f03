!   Generates the initial distribution of "A"'s and "B"'s for simulations of
!   copolymers.
!
!   Quinn updated on 5/22/16 to use thread safe randum number generator
!

subroutine alternChem(AB,NT,N,G,NP,FA,LAM,rand_stat)

!  use mt19937, only : grnd, init_genrand, rnorm, mt, mti
    use mersenne_twister
    use params, only: dp
  implicit none
  integer, intent(in) :: NT          ! Total number of beads
  integer, intent(out) :: AB(NT)     ! Chemical identity of beads
  integer, intent(in) :: N           ! Number of monomers per polymer
  integer, intent(in) :: G           ! Number of beads per monomer
  integer, intent(in) :: NP          ! Number of polymer chains
  real(dp), intent(in) :: LAM !Currently unused

  integer I,J,K,IB
  real TEST(1)   ! changed to real by Quinn
  type(random_stat), intent(inout) ::rand_stat    ! status of random number generator
  !integer ABVAL

  real(dp), intent(in) :: FA   ! Fraction of A beads
  integer nA ! number of A beads per block

  !     idiot check
  if (NT .ne. NP*N*G) then
      call stop_if_err(1, "initchem: nt \= nP*nBpM*nMpP")
  endif

  !		Translate LAM and FA to probabilities


  nA = nint(float(G)*FA)

!		Determine the bead identities

IB = 1
do I = 1,NP
    !TEST = grnd()
    call random_number(TEST,rand_stat)
    if (dble(TEST(1)) < 0.5) then
        do J = 2,N
            do K=1,nA
                AB(IB) = 1
                IB = IB + 1
            enddo
            do K=1,G-nA
                AB(IB) = 0
                IB = IB + 1
            enddo
        enddo
    else
        do J = 2,N
            do K=1,G-nA
                AB(IB) = 0
                IB = IB + 1
            enddo
            do K=1,nA
                AB(IB) = 1
                IB = IB + 1
            enddo
        enddo
    endif
enddo

RETURN
end

!---------------------------------------------------------------*
