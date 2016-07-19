MODULE INPUTPARAMS
  ! various subroutines for reading input from a parameter file

  IMPLICIT NONE
  INTEGER, PARAMETER :: MAXLINELEN = 500 ! maximum length of an input line
  ! current line read in from the input file
  CHARACTER(LEN=MAXLINELEN) :: CURLINE
  INTEGER :: CURPOS ! current position in that line

  ! special characters
  CHARACTER, PARAMETER :: SPACE = ' ', COMMENT = '#'

  ! special strings
  CHARACTER(LEN=26), PARAMETER :: UPLETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(LEN=26), PARAMETER :: LOWLETTERS = 'abcdefghijklmnopqrstuvwxyz'

CONTAINS

  SUBROUTINE READLINE(FILEUNIT,FILEEND,NITEMS)
    ! read a line from the given file unit
    ! store it in the global variable CURLINE for access by all subroutines in this module
    ! Returns END=.true. if end-of-file is hit
    ! Also returns the number of items (delimited by spaces) on the line

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FILEUNIT
    LOGICAL, INTENT(OUT) :: FILEEND
    INTEGER, INTENT(OUT) :: NITEMS
    INTEGER :: IOS, SPACEIND, ADDLEN, LINELEN
    CHARACTER (LEN=MAXLINELEN) :: ADDLINE

    FILEEND = .FALSE.

    CURLINE = ''
    LINELEN = 0
    DO WHILE (LINELEN.LT.MAXLINELEN) ! read until no longer hit line continuation
       READ(FILEUNIT,'(A)',IOSTAT=IOS) ADDLINE       
       IF (IOS.GT.0) THEN
          PRINT*, 'ERROR IN READLINE: some problem with reading in line from input file. IOS = ', IOS
          STOP 1
       ELSEIF (IOS.LT.0) THEN
          FILEEND= .TRUE.
          EXIT
       ENDIF

       ! Check for line continuation, then add on the line
       ADDLEN = LEN(TRIM(ADDLINE))
       
       IF (ADDLEN.GT.2) THEN
          IF (ADDLINE(ADDLEN-2:ADDLEN).EQ.'+++') THEN             
             CURLINE = CURLINE(1:LINELEN) // ADDLINE(1:ADDLEN-3)
             LINELEN = LINELEN + ADDLEN-3
          ELSE
             CURLINE = CURLINE(1:LINELEN) // ADDLINE(1:ADDLEN)
             LINELEN = LINELEN + ADDLEN
             EXIT
          ENDIF
       ELSE
          CURLINE = CURLINE(1:LINELEN) // ADDLINE(1:ADDLEN)
          LINELEN = LINELEN + ADDLEN
          EXIT
       ENDIF
       
    ENDDO   

    NITEMS = 0
    CURPOS = 1
    CALL ADVANCECURPOS    
    DO WHILE (CURPOS.LE.MAXLINELEN)
       NITEMS = NITEMS + 1 ! found an item

       ! index of subsequent space
       SPACEIND = INDEX(CURLINE(CURPOS:MAXLINELEN),SPACE,.FALSE.)
       IF (SPACEIND.EQ.0) THEN
          ! no more spaces; this was last item
          EXIT
       ELSE
          CURPOS = SPACEIND+CURPOS-1
          CALL ADVANCECURPOS
       ENDIF
    ENDDO

    ! reset current position, ignore all initial spaces
    CURPOS = 1
    CALL ADVANCECURPOS

  END SUBROUTINE READLINE

  SUBROUTINE ADVANCECURPOS
    ! advance the current position in the current line to the next character that is not a space
    ! if everything from curpos onwards is spaces, then curpos becomes maxlinelen+1
    IMPLICIT NONE
    CHARACTER :: CHAR

    CHAR = CURLINE(CURPOS:CURPOS)
    DO WHILE (CHAR.EQ.SPACE)
       CURPOS = CURPOS + 1
       IF (CURPOS.GT.MAXLINELEN) RETURN
       CHAR = CURLINE(CURPOS:CURPOS)
    ENDDO
  END SUBROUTINE ADVANCECURPOS

  SUBROUTINE READA(STR, CASESET)
    ! read in a string from the current position on the current line, up until the next space
    ! if CASESET is supplied then CASESET=-1 converts to lowercase
    ! and CASESET=1 converts to uppercase

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(OUT) :: STR
    INTEGER, INTENT(IN) , OPTIONAL :: CASESET
    INTEGER :: SPACEIND

    ! index of next space
    SPACEIND = INDEX(CURLINE(CURPOS:MAXLINELEN),' ',.FALSE.)
    SPACEIND = SPACEIND+CURPOS-1

    STR = CURLINE(CURPOS:SPACEIND-1)
    
    STR = ADJUSTL(STR)

    ! move current position to next non-space character
    CURPOS = SPACEIND
    CALL ADVANCECURPOS
   
    IF (PRESENT(CASESET)) THEN
       IF (CASESET.LT.0) THEN
          ! convert to lowercase
          CALL LOWERCASE(STR)
       ELSEIF (CASESET.GT.0) THEN
          ! convert to uppercase
          CALL UPPERCASE(STR)
       ENDIF
    ENDIF
  END SUBROUTINE READA
  
  SUBROUTINE READU(STR)
    ! uppercase read-in routine for backwards compatibility to old input code
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(OUT) :: STR

    CALL READA(STR,CASESET=1)
  END SUBROUTINE READU

 SUBROUTINE READL(STR)
    ! uppercase read-in routine for backwards compatibility to old input code
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(OUT) :: STR

    CALL READA(STR,CASESET=-1)
  END SUBROUTINE READL

  SUBROUTINE LOWERCASE(STR)
    ! convert all uppercase characters in a string to lowercase
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(INOUT) :: STR
    INTEGER :: LSTR, I, LETIND

    LSTR = LEN(TRIM(STR))
    DO I = 1,LSTR
       ! is this character present in the uppercase letters? If so where?
       LETIND = INDEX(UPLETTERS, STR(I:I),.FALSE.)
       IF (LETIND.GT.0) THEN
          STR(I:I) = LOWLETTERS(LETIND:LETIND)
       ENDIF
    ENDDO
  END SUBROUTINE LOWERCASE

  SUBROUTINE UPPERCASE(STR)
    ! convert all lowercase characters in a string to uppercase
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(INOUT) :: STR
    INTEGER :: LSTR, I, LETIND

    LSTR = LEN(TRIM(STR))
    DO I = 1,LSTR
       ! is this character present in the lowercase letters? If so where?
       LETIND = INDEX(LOWLETTERS, STR(I:I),.FALSE.)
       IF (LETIND.GT.0) THEN
          STR(I:I) = UPLETTERS(LETIND:LETIND)
       ENDIF
    ENDDO
  END SUBROUTINE UPPERCASE

  SUBROUTINE READF(ANS)
    ! read in a double-precision floating point number
    ! exponential notation (1e0, 2d-1,3D+08, 4.5E-01) is allowed

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: ANS
    CHARACTER*100 :: BASESTR, EXPSTR, STR, FMTSTR, FMTSTR1, FMTSTR2
    INTEGER :: EXPIND, PTIND, STRLEN, EXPNUM, err
    DOUBLE PRECISION :: BASENUM

    ! read in as a string, making exponential designators uppercase
    CALL READA(STR,CASESET=1)
        
    ! find index of exponential designator
    EXPIND = INDEX(STR,'D',.FALSE.)
    IF (EXPIND.EQ.0) THEN
       EXPIND = INDEX(STR,'E',.FALSE.)
    ENDIF

    ! get the base number (without the exponential part)
    ! and the exponential part separately
    IF (EXPIND.EQ.0) THEN
       BASESTR= STR
       EXPSTR = '0'
    ELSE
       BASESTR = STR(1:EXPIND-1)
       EXPSTR = STR(EXPIND+1:LEN(STR))
    ENDIF    
    
    ! find decimal point in base number
    PTIND = INDEX(BASESTR,'.',.FALSE.)    

    ! format specifier for reading in decimal base number
    STRLEN = LEN(TRIM(BASESTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    IF (PTIND.GT.0) THEN
       WRITE(FMTSTR2,'(I10)')  STRLEN-PTIND
    ELSE
       FMTSTR2 = '0'
    ENDIF

    FMTSTR = '(F' // TRIM(ADJUSTL(FMTSTR1))// '.' // TRIM(ADJUSTL(FMTSTR2)) // ')'
  
    ! read in the base number
    READ(BASESTR,FMTSTR,IOSTAT=ERR) BASENUM
    
    IF (ERR.NE.0) THEN
       PRINT*, 'ERROR IN READF: something wrong with reading base number  as float:.'
       PRINT*, 'attempted to read: ', BASESTR
       PRINT*, 'ERROR CODE:', ERR
       STOP 1
    ENDIF
    
    ! read in the exponential number
    STRLEN = LEN(TRIM(EXPSTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'
    
    READ(EXPSTR,FMTSTR,IOSTAT=ERR) EXPNUM

    IF (ERR.NE.0) THEN
       PRINT*, 'ERROR IN READF: something wrong with reading exponent as integer:.'
       PRINT*, 'attempted to read: ', EXPSTR
       PRINT*, 'ERROR CODE:', ERR
       STOP 1
    ENDIF    

    ANS = BASENUM*10D0**(EXPNUM)
        
  END SUBROUTINE READF

  SUBROUTINE READI(ANS)
    ! read in an integer
    ! exponential notation (2D5, -1D3) is allowed
    ! returns error if number is too large to fit in integer kind

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ANS
    CHARACTER*100 :: BASESTR, EXPSTR, STR, FMTSTR, FMTSTR1, FMTSTR2
    INTEGER :: EXPIND, STRLEN, BASENUM,EXPNUM, err
    DOUBLE PRECISION :: EXPLIM
    
    ! read in as a string, making exponential designators uppercase
    CALL READA(STR,CASESET=1)
        
    ! find index of exponential designator
    EXPIND = INDEX(STR,'D',.FALSE.)
    IF (EXPIND.EQ.0) THEN
       EXPIND = INDEX(STR,'E',.FALSE.)
    ENDIF

    ! get the base number (without the exponential part)
    ! and the exponential part separately
    IF (EXPIND.EQ.0) THEN
       BASESTR= STR
       EXPSTR = '0'
    ELSE
       BASESTR = STR(1:EXPIND-1)
       EXPSTR = STR(EXPIND+1:LEN(STR))
    ENDIF    

     ! read in the base number
    STRLEN = LEN(TRIM(BASESTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'
    READ(BASESTR,FMTSTR,IOSTAT=ERR) BASENUM

    IF (ERR.NE.0) THEN
       PRINT*, 'ERROR IN READI: something wrong with reading base number as integer:.'
       PRINT*, 'attempted to read: ', BASESTR
       PRINT*, 'ERROR CODE:', ERR
       STOP 1
    ENDIF

    ! read in the exponent
    STRLEN = LEN(TRIM(EXPSTR))
    WRITE(FMTSTR1,'(I10)') STRLEN
    FMTSTR = '(I'//TRIM(ADJUSTL(FMTSTR1))//')'    
    READ(EXPSTR,FMTSTR,IOSTAT=ERR) EXPNUM

    IF (ERR.NE.0) THEN
       PRINT*, 'ERROR IN READI: something wrong with reading exponent as integer:.'
       PRINT*, 'attempted to read: ', EXPSTR
       PRINT*, 'ERROR CODE:', ERR
       STOP 1
    ENDIF

    ! cannot have negative exponent
    IF(EXPNUM.LT.0) THEN
       PRINT*, 'ERROR IN READI: negative exponents not allowed for integer input'
       PRINT*, 'Attempted to read as integer: ', TRIM(STR)
       STOP 1
    ENDIF

    IF (ABS(BASENUM).GT.0) THEN
       ! double check that the number is not too big to fit as integer
       ! and that the exponent is not negative
       EXPLIM = DBLE(HUGE(ANS))/DBLE(ABS(BASENUM))
       EXPLIM = LOG(EXPLIM)/LOG(10D0)
       IF (EXPNUM.GT.FLOOR(EXPLIM)) THEN
          PRINT*, 'ERROR IN READI: number is too big to fit as integer.'
          PRINT*, 'Attempted to read as integer: ', TRIM(STR)
          STOP 1
       ENDIF
    ENDIF

    ANS = BASENUM*10**(EXPNUM)
        
  END SUBROUTINE READI

  SUBROUTINE READO(ANS)
    ! read in a logical argument
    ! must be either T, F, TRUE, FALSE, '1', or '0' (1 for true, 0 for false)
    ! not case sensitive
    ! returns error if some other input does not match one of these strings

    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: ANS
    CHARACTER*100 :: STR

    ! read in as a string, making uppercase
    CALL READA(STR,CASESET=1)    

    IF (STR.EQ.'T'.OR.STR.EQ.'1'.OR.STR.EQ.'TRUE') THEN
       ANS = .TRUE.
    ELSEIF (STR.EQ.'F'.OR.STR.EQ.'0'.OR.STR.EQ.'FALSE')  THEN
       ANS = .FALSE.
    ELSE
       PRINT*, 'ERROR IN READO: failed to read string as logical. Logical string must be T, F, TRUE, FALSE, 0, or 1'
       print*, 'Attempted to read in: ', str
       stop 1
    ENDIF      

  END SUBROUTINE READO

END MODULE INPUTPARAMS
