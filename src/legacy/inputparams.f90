module inPUTPARAMS
  ! various subroutines for reading input from a parameter file

  implicit none
  integer, PARAMETER :: MAXLinELEN = 500 ! maximum length of an input line
  ! current line read in from the input file
  CHARACTER(LEN = MAXLinELEN) :: CURLinE
  integer :: CURPOS ! current position in that line

  ! special characters
  CHARACTER, PARAMETER :: SPACE = ' ', COMMENT = '#'

  ! special strings
  CHARACTER(LEN = 26), PARAMETER :: UPLETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(LEN = 26), PARAMETER :: LOWLETTERS = 'abcdefghijklmnopqrstuvwxyz'

contains

  subroutine READLinE(FILEUNIT,FILEEND,NITEMS)
    ! read a line from the given file unit
    ! store it in the global variable CURLinE for access by all subroutines in this module
    ! Returns END = .true. if end-of-file is hit
    ! Also returns the number of items (delimited by spaces) on the line

    implicit none
    integer, intent(in) :: FILEUNIT
    LOGICAL, intent(out) :: FILEEND
    integer, intent(out) :: NITEMS
    integer :: IOS, SPACEinD, ADDLEN, LinELEN
    CHARACTER (LEN = MAXLinELEN) :: ADDLinE

    FILEEND = .FALSE.

    CURLinE = ''
    LinELEN = 0
    do while (LinELEN < MAXLinELEN) ! read until no longer hit line continuation
       READ(FILEUNIT,'(A)',IOSTAT = IOS) ADDLinE
       if (IOS > 0) then
          PRinT*, 'ERROR in READLinE: some problem with reading in line from input file. IOS = ', IOS
          STOP 1
       elseif (IOS < 0) then
          FILEEND= .TRUE.
          EXIT
       ENDif

       ! Check for line continuation, then add on the line
       ADDLEN = LEN(TRIM(ADDLinE))

       if (ADDLEN > 2) then
          if (ADDLinE(ADDLEN-2:ADDLEN) == ' + + + ') then
             CURLinE = CURLinE(1:LinELEN) // ADDLinE(1:ADDLEN-3)
             LinELEN = LinELEN + ADDLEN-3
          else
             CURLinE = CURLinE(1:LinELEN) // ADDLinE(1:ADDLEN)
             LinELEN = LinELEN + ADDLEN
             EXIT
          ENDif
       else
          CURLinE = CURLinE(1:LinELEN) // ADDLinE(1:ADDLEN)
          LinELEN = LinELEN + ADDLEN
          EXIT
       ENDif

    ENDdo

    NITEMS = 0
    CURPOS = 1
    CALL ADVANCECURPOS
    do while (CURPOS <= MAXLinELEN)
       NITEMS = NITEMS + 1 ! found an item

       ! index of subsequent space
       SPACEinD = inDEX(CURLinE(CURPOS:MAXLinELEN),SPACE,.FALSE.)
       if (SPACEinD == 0) then
          ! no more spaces; this was last item
          EXIT
       else
          CURPOS = SPACEinD + CURPOS-1
          CALL ADVANCECURPOS
       ENDif
    ENDdo

    ! reset current position, ignore all initial spaces
    CURPOS = 1
    CALL ADVANCECURPOS

  END subroutine READLinE

  subroutine ADVANCECURPOS
    ! advance the current position in the current line to the next character that is not a space
    ! if everything from curpos onwards is spaces, then curpos becomes maxlinelen + 1
    implicit none
    CHARACTER :: CHAR

    CHAR = CURLinE(CURPOS:CURPOS)
    do while (CHAR == SPACE)
       CURPOS = CURPOS + 1
       if (CURPOS > MAXLinELEN) RETURN
       CHAR = CURLinE(CURPOS:CURPOS)
    ENDdo
  END subroutine ADVANCECURPOS

  subroutine READA(STR, CASESET)
    ! read in a string from the current position on the current line, up until the next space
    ! if CASESET is supplied then CASESET = -1 converts to lowercase
    ! and CASESET = 1 converts to uppercase

    implicit none
    CHARACTER(LEN = *), intent(out) :: STR
    integer, intent(in) , OPTIONAL :: CASESET
    integer :: SPACEinD

    ! index of next space
    SPACEinD = inDEX(CURLinE(CURPOS:MAXLinELEN),' ',.FALSE.)
    SPACEinD = SPACEinD + CURPOS-1

    STR = CURLinE(CURPOS:SPACEinD-1)

    STR = ADJUSTL(STR)

    ! move current position to next non-space character
    CURPOS = SPACEinD
    CALL ADVANCECURPOS

    if (PRESENT(CASESET)) then
       if (CASESET < 0) then
          ! convert to lowercase
          CALL LOWERCASE(STR)
       elseif (CASESET > 0) then
          ! convert to uppercase
          CALL UPPERCASE(STR)
       ENDif
    ENDif
  END subroutine READA

  subroutine READU(STR)
    ! uppercase read-in routine for backwards compatibility to old input code
    implicit none
    CHARACTER(LEN = *), intent(out) :: STR

    CALL READA(STR,CASESET = 1)
  END subroutine READU

 subroutine READL(STR)
    ! uppercase read-in routine for backwards compatibility to old input code
    implicit none
    CHARACTER(LEN = *), intent(out) :: STR

    CALL READA(STR,CASESET = -1)
  END subroutine READL

  subroutine LOWERCASE(STR)
    ! convert all uppercase characters in a string to lowercase
    implicit none

    CHARACTER(LEN = *), intent(inout) :: STR
    integer :: LSTR, I, LETinD

    LSTR = LEN(TRIM(STR))
    do I = 1,LSTR
       ! is this character present in the uppercase letters? If so where?
       LETinD = inDEX(UPLETTERS, STR(I:I),.FALSE.)
       if (LETinD > 0) then
          STR(I:I) = LOWLETTERS(LETinD:LETinD)
       ENDif
    ENDdo
  END subroutine LOWERCASE

  subroutine UPPERCASE(STR)
    ! convert all lowercase characters in a string to uppercase
    implicit none

    CHARACTER(LEN = *), intent(inout) :: STR
    integer :: LSTR, I, LETinD

    LSTR = LEN(TRIM(STR))
    do I = 1,LSTR
       ! is this character present in the lowercase letters? If so where?
       LETinD = inDEX(LOWLETTERS, STR(I:I),.FALSE.)
       if (LETinD > 0) then
          STR(I:I) = UPLETTERS(LETinD:LETinD)
       ENDif
    ENDdo
  END subroutine UPPERCASE

  subroutine READF(ANS)
    ! read in a double-precision floating point number
    ! exponential notation (1e0, 2d-1,3D + 08, 4.5E-01) is allowed

    implicit none
    real(dp), intent(out) :: ANS
    CHARACTER*100 :: BASESTR, EXPSTR, STR, FMTSTR, FMTSTR1, FMTSTR2
    integer :: EXPinD, PTinD, STRLEN, EXPNUM, err
    real(dp) :: BASENUM

    ! read in as a string, making exponential designators uppercase
    CALL READA(STR,CASESET = 1)

    ! find index of exponential designator
    EXPinD = inDEX(STR,'D',.FALSE.)
    if (EXPinD == 0) then
       EXPinD = inDEX(STR,'E',.FALSE.)
    ENDif

    ! get the base number (without the exponential part)
    ! and the exponential part separately
    if (EXPinD == 0) then
       BASESTR= STR
       EXPSTR = '0'
    else
       BASESTR = STR(1:EXPinD-1)
       EXPSTR = STR(EXPinD + 1:LEN(STR))
    ENDif

    ! find decimal point in base number
    PTinD = inDEX(BASESTR,'.',.FALSE.)

    ! format specifier for reading in decimal base number
    STRLEN = LEN(TRIM(BASESTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    if (PTinD > 0) then
       WRITE(FMTSTR2,'(I10)')  STRLEN-PTinD
    else
       FMTSTR2 = '0'
    ENDif

    FMTSTR = '(F' // TRIM(ADJUSTL(FMTSTR1))// '.' // TRIM(ADJUSTL(FMTSTR2)) // ')'

    ! read in the base number
    READ(BASESTR,FMTSTR,IOSTAT = ERR) BASENUM

    if (ERR /= 0) then
       PRinT*, 'ERROR in READF: something wrong with reading base number  as float:.'
       PRinT*, 'attempted to read: ', BASESTR
       PRinT*, 'ERROR CODE:', ERR
       STOP 1
    ENDif

    ! read in the exponential number
    STRLEN = LEN(TRIM(EXPSTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'

    READ(EXPSTR,FMTSTR,IOSTAT = ERR) EXPNUM

    if (ERR /= 0) then
       PRinT*, 'ERROR in READF: something wrong with reading exponent as integer:.'
       PRinT*, 'attempted to read: ', EXPSTR
       PRinT*, 'ERROR CODE:', ERR
       STOP 1
    ENDif

    ANS = BASENUM*10D0**(EXPNUM)

  END subroutine READF

  subroutine READI(ANS)
    ! read in an integer
    ! exponential notation (2D5, -1D3) is allowed
    ! returns error if number is too large to fit in integer kind

    implicit none
    integer, intent(out) :: ANS
    CHARACTER*100 :: BASESTR, EXPSTR, STR, FMTSTR, FMTSTR1, FMTSTR2
    integer :: EXPinD, STRLEN, BASENUM,EXPNUM, err
    real(dp) :: EXPLIM

    ! read in as a string, making exponential designators uppercase
    CALL READA(STR,CASESET = 1)

    ! find index of exponential designator
    EXPinD = inDEX(STR,'D',.FALSE.)
    if (EXPinD == 0) then
       EXPinD = inDEX(STR,'E',.FALSE.)
    ENDif

    ! get the base number (without the exponential part)
    ! and the exponential part separately
    if (EXPinD == 0) then
       BASESTR= STR
       EXPSTR = '0'
    else
       BASESTR = STR(1:EXPinD-1)
       EXPSTR = STR(EXPinD + 1:LEN(STR))
    ENDif

     ! read in the base number
    STRLEN = LEN(TRIM(BASESTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'
    READ(BASESTR,FMTSTR,IOSTAT = ERR) BASENUM

    if (ERR /= 0) then
       PRinT*, 'ERROR in READI: something wrong with reading base number as integer:.'
       PRinT*, 'attempted to read: ', BASESTR
       PRinT*, 'ERROR CODE:', ERR
       STOP 1
    ENDif

    ! read in the exponent
    STRLEN = LEN(TRIM(EXPSTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'
    READ(EXPSTR,FMTSTR,IOSTAT = ERR) EXPNUM

    if (ERR /= 0) then
       PRinT*, 'ERROR in READI: something wrong with reading exponent as integer:.'
       PRinT*, 'attempted to read: ', EXPSTR
       PRinT*, 'ERROR CODE:', ERR
       STOP 1
    ENDif

    ! cannot have negative exponent
    if(EXPNUM < 0) then
       PRinT*, 'ERROR in READI: negative exponents not allowed for integer input'
       PRinT*, 'Attempted to read as integer: ', TRIM(STR)
       STOP 1
    ENDif

    if (ABS(BASENUM) > 0) then
       ! double check that the number is not too big to fit as integer
       ! and that the exponent is not negative
       EXPLIM = DBLE(HUGE(ANS))/DBLE(ABS(BASENUM))
       EXPLIM = LOG(EXPLIM)/LOG(10D0)
       if (EXPNUM > FLOOR(EXPLIM)) then
          PRinT*, 'ERROR in READI: number is too big to fit as integer.'
          PRinT*, 'Attempted to read as integer: ', TRIM(STR)
          STOP 1
       ENDif
    ENDif

    ANS = BASENUM*10**(EXPNUM)

  END subroutine READI

  subroutine REAdo(ANS)
    ! read in a logical argument
    ! must be either T, F, TRUE, FALSE, '1', or '0' (1 for true, 0 for false)
    ! not case sensitive
    ! returns error if some other input does not match one of these strings

    implicit none
    LOGICAL, intent(out) :: ANS
    CHARACTER*100 :: STR

    ! read in as a string, making uppercase
    CALL READA(STR,CASESET = 1)

    if (STR == 'T'.OR.STR == '1'.OR.STR == 'TRUE') then
       ANS = .TRUE.
    elseif (STR == 'F'.OR.STR == '0'.OR.STR == 'FALSE')  then
       ANS = .FALSE.
    else
       PRinT*, 'ERROR in REAdo: failed to read string as logical. Logical string must be T, F, TRUE, FALSE, 0, or 1'
       print*, 'Attempted to read in: ', str
       stop 1
    ENDif

  END subroutine REAdo

END module inPUTPARAMS
