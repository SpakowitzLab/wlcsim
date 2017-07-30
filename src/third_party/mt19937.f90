! A Fortran-program for MT19937: Real number version
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-11-26  Time: 17:09:23
! Latest revision - 5 February 2002
! A new seed initialization routine has been added based upon the new
! C version dated 26 January 2002.
! This version assumes that integer overflows do NOT cause crashes.
! This version is compatible with Lahey's ELF90 compiler,
! and should be compatible with most full Fortran 90 or 95 compilers.
! Notice the strange way in which umask is specified for ELF90.
 
!   genrand() generates one pseudorandom real number (double) which is
! uniformly distributed on [0,1]-interval, for each call.
! sgenrand(seed) set initial values to the working area of 624 words.
! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
! integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.

! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Library General Public License as published by
! the Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.   This library is distributed in the hope that
! it will be useful, but WITHout ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> real(dp) function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
!               If n < 0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

module mt19937
implicit none
integer, PARAMETER :: dp = SELECTED_REAL_KinD(12, 60)

! Period parameters
integer, PARAMETER :: n = 624, n1 = n + 1, m = 397, mata = -1727483681
!                                    constant vector a
integer, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
integer, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
integer, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
integer, SAVE      :: mt(0:n-1), mti = n1
!                     mti = =N + 1 means mt[N] is not initialized
LOGICAL, SAVE :: RNORMRESTART = .FALSE.
REAL, SAVE :: RNORMVSAVE, RNORMSLNSAVE

private
PUBLIC :: dp, sgrnd, grnd, init_genrand, rnorm, MT,mti
! efk: global variables to allow restarting of rnorm() calculation from savefile
PUBLIC :: RNORMRESTART, RNORMSLNSAVE, RNORMVSAVE

contains


subroutine sgrnd(seed)
! This is the original version of the seeding routine.
! It was replaced in the Japanese version in C on 26 January 2002
! It is recommended that routine init_genrand is used instead.

integer, intent(in)   :: seed

!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

mt(0)= IAND(seed, -1)
do  mti = 1,n-1
  mt(mti) = IAND(69069 * mt(mti-1), -1)
END do

RETURN
END subroutine sgrnd
!***********************************************************************

subroutine init_genrand(seed)
! This initialization is based upon the multiplier given on p.106 of the
! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.

! This version assumes that integer overflow does NOT cause a crash.

integer, intent(in)  :: seed

integer  :: latest

mt(0) = seed
latest = seed
do mti = 1, n-1
  latest = IEOR( latest, ISHFT( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
END do

RETURN
END subroutine init_genrand
!***********************************************************************

FUNCTION grndORIG() RESULT(fn_val)
REAL (dp) :: fn_val

integer, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x = 0,1
integer       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)

if(mti >= n) then
!                       generate N words at one time
  if(mti == n + 1) then
!                            if sgrnd() has not been called,
    CALL sgrnd(4357)
!                              a default initial seed is used
  END if
  
  do  kk = 0, n-m-1
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk + 1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk + m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END do
  do  kk = n-m, n-2
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk + 1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk + (m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END do
  y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
  mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti = 0
END if

y = mt(mti)
mti = mti + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

if(y < 0) then
  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
else
  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
END if

RETURN
END FUNCTION grndORIG


FUNCTION tshftu(y) RESULT(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl

! EFK: 2009/08/01 converted to use MT19937 random number generator
FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

implicit none
REAL  :: fn_val

! Local variables

REAL            :: u, sum
REAL, SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TinY( one )

!print*, 'testx0:', second, rnormrestart, v, sln, rnormvsave, rnormslnsave
 
if (second) then
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

else if (RNORMRESTART) then
   ! efk: restart from global variables in save file
   SECOND = .FALSE.
   v = rnormvsave; sln = rnormslnsave
   FN_Val = v*sln
   RNORMRESTART = .FALSE.
else
! First call; generate a pair of random normals

  second = .true.
  do
    U = GRND( )
    V = GRND()
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    if(sum < one) EXIT
  END do
  sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
  fn_val = u*sln
END if

RNORMVSAVE = V; RNORMSLNSAVE = SLN; RNORMRESTART = SECOND

RETURN
END FUNCTION rnorm

! EFK: 2009/09/03 wrap around GRND so that it does not ever give exactly 0 or 1
FUNCTION grnd() RESULT(fn_val)
  REAL (dp) :: fn_val
  
  FN_VAL = 0D0
  do while (FN_VAL == 0D0.OR.FN_VAL == 1D0)
     FN_VAL = GRNdoRIG()     
  ENDdo
END FUNCTION grnd

END module mt19937



