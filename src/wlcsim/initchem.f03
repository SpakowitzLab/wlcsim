!   Generates the initial distribution of "A"'s and "B"'s for simulations of
!   copolymers.
!
!   Quinn updated on 5/22/16 to use thread safe randum number generator
!

subroutine initchem(AB,NT,N,G,NP,FA,LAM,rand_stat)

!  use mt19937, only : grnd, init_genrand, rnorm, mt, mti
    use mersenne_twister
    use params, only: dp
  implicit none
  integer, intent(in) :: NT          ! Total number of beads
  integer, intent(out) :: AB(NT)     ! Chemical identity of beads
  integer, intent(in) :: N           ! Number of monomers per polymer
  integer, intent(in) :: G           ! Number of beads per monomer
  integer, intent(in) :: NP          ! Number of polymer chains

  integer I,J,K,IB
  real(dp) TEST(1)   ! changed to real by Quinn
  type(random_stat), intent(inout) ::rand_stat    ! status of random number generator
  !integer ABVAL

  real(dp), intent(in) :: FA   ! Fraction of A beads
  real(dp), intent(in) :: LAM  ! Chemical correlation parameter
  real(dp) PAA,PBB,PAB,PBA ! Chemical identity statistics

  !     idiot check
  if (NT .ne. NP*N*G) then
      call stop_if_err(1, "initchem: nt \= nP*nBpM*nMpP")
  endif

  !		Translate LAM and FA to probabilities

  PAA = FA*(1.0_dp-LAM) + LAM
  PBB = FA*(LAM-1.0_dp) + 1.0_dp
  PBA = 1.0_dp-PAA
  PAB = 1.0_dp-PBB

  !		Determine the bead identities

IB = 1
do I = 1,NP
    !TEST = grnd()
    call random_number(TEST,rand_stat)
    if (dble(TEST(1)) < FA) then
        AB(IB) = 1
    else
        AB(IB) = 0
    endif
    IB = IB + 1
    do K = 2,G
        AB(IB) = AB(IB-1)
        IB = IB + 1
    enddo

    do J = 2,N
        !TEST = grnd()
        call random_number(TEST,rand_stat)
        if (AB(IB-1) == 1) then
           if (TEST(1) <= PAA) then
              AB(IB) = 1
           else
              AB(IB) = 0
           endif
        else
           if (TEST(1) <= PAB) then
              AB(IB) = 1
           else
              AB(IB) = 0
           endif
        endif
        IB = IB + 1

        do K = 2,G
            AB(IB) = AB(IB-1)
            IB = IB + 1
        enddo
    enddo
enddo

RETURN
end

!---------------------------------------------------------------*
