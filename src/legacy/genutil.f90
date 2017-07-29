module GENUTIL
  ! generally useful utilities
  use MT19937 ! mersenne random number generator
  use QUATUTIL ! utilities for dealing with quaternions

  implicit none

contains

  integer FUNCTION STRinG2NUM(STRinGin,APPENDNUM)
    ! convert the string to a unique number based on ascii characters
    ! the characters SPACE, {,},(,),[,],",`,<,> and all nonprintable characters are ignored
    ! at most the last five characters (ignoring the unacceptable characters above) at the end of the string are used
    ! any leading "!" do not affect the final number (these map to 0)
    ! if APPENDNUM is specificied, only use the last 4 characters of the string as well as the additional number modulo 84

    implicit none
    CHARACTER(LEN = *) :: STRinGin
    CHARACTER(LEN = 5) :: STRinG
    integer, OPTIONAL :: APPENDNUM
    integer :: DIGARRAY(5)
    integer :: ALLOWED(84)
    integer :: N, I, D, COUNT
    CHARACTER*84 :: ALLOWEDSTR

    ! set the allowed characters
    ALLOWED(1:6) = (/33,35,36,37,38,39/)
    ALLOWED(7:24) = (/(I,I = 42,59)/)
    ALLOWED(25:27) = (/61,63,64/)
    ALLOWED(28:53) = (/(I, I = 65,90)/)
    ALLOWED(54:56) = (/92,94,95/)
    ALLOWED(57:82) = (/(I, I = 97,122)/)
    ALLOWED(83:84) = (/124,126/)

    N = LEN(STRinGin)
    if (PRESENT(APPENDNUM)) then
       STRinG(1:4) = STRinGin(N-3:N)
       STRinG(5:5) = ACHAR(ALLOWED(MOD(APPENDNUM,84) + 1))
    else
       STRinG = STRinGin(N-4:N)
    ENDif
    N =  5


    do I = 1,84
       ALLOWEDSTR(I:I) = ACHAR(ALLOWED(I))
    ENDdo

    DIGARRAY = 0
    COUNT = 0
    do I = 0,N-1
       D = inDEX(ALLOWEDSTR,STRinG(N-I:N-I),.FALSE.)
       if (D == 0) then
          print*, 'Ignoring character:', D
          CYCLE
       ENDif

       DIGARRAY(5-COUNT) = D-1
       COUNT = COUNT + 1
       if (COUNT >= 5) EXIT
    ENDdo

    STRinG2NUM = BASE2DEC(DIGARRAY,5,84)
  END FUNCTION STRinG2NUM

  integer FUNCTION BASE2DEC(DIGARRAY,N,BASE)
  ! given a number in some integer base (specified as a list of digits)
  ! convert that number to a decimal integer
  ! N is the size of the list
  ! if resulting number is too large, wrap around to negative numbers
  ! starting from the right, only use as many of the digits as
  ! will fit into the resulting integer between -HUGE and HUGE
  ! if any digit is greater than base-1, print error and stop

  implicit none
  integer, DIMENSION(N) :: DIGARRAY
  integer, intent(in) :: N, BASE
  integer :: MAXDIG, I, D

  MAXDIG = inT(LOG(2*DBLE(HUGE(BASE)) + 2)/LOG(DBLE(BASE)))

  BASE2DEC = 0
  do I = 0, Min(N-1,MAXDIG-1)
     D = DIGARRAY(N-I)
     if (D == 0) CYCLE
     if (D > BASE-1) then
        PRinT*, 'ERROR in BASE2DEC: digit is bigger than base.', I, D, BASE
        STOP 1
     ENDif

     BASE2DEC = BASE2DEC + D*BASE**I
  ENDdo

  END FUNCTION BASE2DEC

  subroutine inTERPARRAY(ARRAY,NA,COL,VAL,inD,inTERP)
    ! for an 2D array with dimensions NA
    ! use the values in column COL to interpolate for the value VAL
    ! return the index inD such that ARRAY(inD,COL)<VAL<ARRAY(inD + 1,COL)
    ! and the interpolation of all other columns in inTERP
    ! If val is out of bounds returns inD = 0 or inD = NL
    implicit none
    integer, intent(in) :: NA(2), COL
    real(dp), intent(in) :: ARRAY(NA(1),NA(2))
    real(dp), intent(in) :: VAL
    integer, intent(out) :: inD
    real(dp), intent(out) :: inTERP(NA(2))
    real(dp) :: FRAC

    CALL inTERP1(ARRAY(:,COL),NA(1),VAL,inD)

    if (inD == 0.OR.inD >= NA(1)) RETURN

    FRAC = (VAL-ARRAY(inD,COL))/(ARRAY(inD + 1,COL)-ARRAY(inD,COL))
    inTERP = (1-FRAC)*ARRAY(inD,:) + FRAC*ARRAY(inD + 1,:)

  END subroutine inTERPARRAY

  subroutine inTERP1(LIST,NL,VAL,inD)
    ! for a monotonically increasing, real(dp) list
    ! find the index where LIST(inD) <VAL<LIST(inD + 1)
    integer, intent(in) :: NL
    real(dp), intent(in) :: LIST(NL)
    real(dp), intent(in) :: VAL
    integer, intent(out) :: inD
    real(dp) :: MinL, MAXL,PREVMAXL
    integer :: MinI, MAXI,PREVMAXI
    LOGICAL :: VERBOSE = .FALSE.

    MinI = 1; MAXI = NL; MinL = LIST(1); MAXL = LIST(NL)
    if (VAL < MinL) then
       inD = 0; RETURN
    elseif (VAL == MinL) then
       inD = 1; RETURN
    elseif (VAL > MAXL) then
       inD = NL; RETURN
    elseif (VAL == MAXL) then
       inD = NL-1; RETURN
    ENDif

    do while (MAXI-MinI > 1.OR.MAXL < VAL)
       if (MAXL > VAL) then
          PREVMAXL = MAXL; PREVMAXI = MAXI
          MAXI = MinI + (MAXI-MinI)/2
          MAXL = LIST(MAXI)
       else
          MinI = MAXI; MAXI = PREVMAXI
          MinL = MAXL; MAXL = PREVMAXL
       ENDif
       if (VERBOSE) PRinT*, 'MinI, MAXI, MinL, MAXL', MinI, MAXI, MinL, MAXL,VAL
       if (maxi.eq.mini) then
          print*, 'something weird in interp1:', list(1), list(nl), val
          stop 1
       endif
    ENDdo

    if (.NOT.(MAXI == MinI + 1.AND.LIST(MinI) <= VAL.AND.LIST(MAXI) >= VAL)) then
       PRinT*, 'SOMETHinG IS WEIRD in inTERP1', val, mini, maxi, list(mini), list(maxi)
       STOP 1
    ENDif

    inD = MinI
  END subroutine inTERP1

  subroutine REPLACESUBSTR(inSTRinG,C,REPL)
    ! replace the last instance of the substring C in inSTRinG with REPL
    implicit none
    CHARACTER*100 :: inSTRinG
    CHARACTER(LEN = *) :: C
    CHARACTER(LEN = *) :: REPL
    integer :: LENC, inD

    inSTRinG = ADJUSTL(inSTRinG)

    LENC = LEN_TRIM(C)

    inD = inDEX(inSTRinG,C,.TRUE.)
    if (inD > 0) then! if * was found in the string

       inSTRinG = inSTRinG(1:inD-1) // TRIM(ADJUSTL(REPL)) // inSTRinG(inD + LENC:100)
    END if
  END subroutine REPLACESUBSTR

  subroutine NORMALIZE(X)
    ! normalize a 3 dimensional vector

    implicit none

    real(dp), intent(inout) :: X(3)
    real(dp) :: DX

    DX = SQRT(doT_PRODUCT(X,X))
    X(:) = X(:)/DX

    RETURN
  END subroutine NORMALIZE

  real(dp) FUNCTION NORM(X)
    ! norm of 3D vector X

    implicit none
    real(dp), intent(in) :: X(3)

    NORM = sqrt(doT_PRODUCT(X,X))

  END FUNCTION NORM

  subroutine CROSS_PRODUCT(A, B, C)
    ! take the cross product of 3D vectors A and B; return result in C

    implicit none

    real(dp), intent(in) :: A(3), B(3)
    real(dp), intent(out) :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1)-A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END subroutine CROSS_PRODUCT

  subroutine RANdoMAXIS(REFAX,CTRANGE,RANDAX)
    ! generate a random axis, within a certain range in cos(theta)
    ! relative to the reference axis
    real(dp), intent(in) :: REFAX(3), CTRANGE
    real(dp), intent(out) :: RANDAX(3)
    real(dp) :: THETA, pHI, E1(3), E2(3), E3(3), X, Y, Z

    THETA = 1.0D0
    do while (THETA == 1.0D0)
       ! get random number; MT19937 algorithm uses closed interval [0,1],
       ! so ignore when R is exactly 1
       THETA = GRND() !get a random number
    ENDdo
    !CALL RANdoM_NUMBER(THETA)
    THETA = acos(1D0 - THETA*MAX(CTRANGE,2D0))

    PHI = 1.0D0
    do while (PHI == 1.0D0)
       PHI = GRND() !get a random number
    ENDdo
    !CALL RANdoM_NUMBER(PHI)
    PHI = PHI*2*PI

    ! axis system relative to which angles are defined
    E3 = REFAX
    if (E3(2) == 0 .AND. E3(3) == 0) then
       E2 = (/0D0,1D0,0D0/)
    else
       CALL CROSS_PRODUCT(E3, (/1D0,0D0,0D0/),E2)
    END if
    CALL CROSS_PRODUCT(E2, E3, E1)

    CALL NORMALIZE(E1); CALL NORMALIZE(E2); CALL NORMALIZE(E3)

    ! generate the axis around which to rotate
    X = sin(THETA)*cos(PHI)
    Y = sin(THETA)*sin(PHI)
    Z = cos(THETA)

    RANDAX = X*E1 + Y*E2 + Z*E3

  END subroutine RANdoMAXIS

  subroutine GETPERP(V1,V2)
    ! get some unit vector perpendicular to V1 and store it in V2
    implicit none
    real(dp), intent(in) :: V1(3)
    real(dp), intent(out) :: V2(3)

    if (V1(2) == 0.AND.V1(3) == 0) then
       V2 = (/0D0,1D0,0D0/)
    else
       CALL CROSS_PRODUCT(V1,(/0D0,1D0,0D0/),V2)
       CALL NORMALIZE(V2)
    ENDif
  END subroutine GETPERP
END module GENUTIL
