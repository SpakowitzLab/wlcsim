MODULE INPUTUTIL
  ! utilities for reading in structures from a file
  IMPLICIT NONE

CONTAINS

   SUBROUTINE READSNAPSHOTS(CGRP,INFILE,GETPARAM,NUM,SUCCESS)
    ! read in a snapshot file
    ! set the chain configuration to the last complete snapshot
    ! NUM is the additional saved number for the last snapshot (ie: MC step)
    ! chainp must be initialized and have all params set already
    ! COMPLETE indicates whether a complete structure was successfully read
    ! if not successful, coords remain unchanged
    USE MANYCHAINS, ONLY : CHAINGROUP, SETUPCHAINGROUP, COPYCHAINGROUP, CLEANUPCHAINGROUP
    USE CHAINUTIL, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP, CGRP2
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*), INTENT(IN) :: INFILE
    LOGICAL, INTENT(IN) :: GETPARAM
    INTEGER, INTENT(OUT) :: NUM
    LOGICAL, INTENT(OUT) :: SUCCESS
    TYPE(CHAINGROUP), TARGET :: GROUP2
    INTEGER :: NUMSAVE, NCHAIN, MAXNPT
    LOGICAL :: COMPLETE
    INTEGER :: B, C
    CHARACTER :: CDUM
    DOUBLE PRECISION :: DUM1, DUM2, TMPVEC(6)
    INTEGER :: NBEAD, ICHECK

    ! make a copy group
    CGRP2=>GROUP2
    MAXNPT = CGRP%CHAINS(1)%MAXNPT
    CALL SETUPCHAINGROUP(CGRP2,CGRP%NCHAIN,CGRP%NCONNECT,CGRP%NFORCE,MAXNPT)
    CALL COPYCHAINGROUP(CGRP,CGRP2)

    
    NUMSAVE=0

    OPEN(UNIT=11,FILE=INFILE,STATUS='OLD')
    COMPLETE = .FALSE.; SUCCESS = .FALSE.
    DO
       READ(11,'(A,1X,2I12)',IOSTAT=ICHECK) CDUM,NCHAIN,NUM
      
       IF (ICHECK.LT.0) exit ! end of file

       IF (NCHAIN.NE.CGRP%NCHAIN) THEN
          PRINT*, 'ERROR IN READSNAPSHOT: wrong number of chains', NCHAIN, CGRP%NCHAIN
          STOP 1
       ENDIF

       DO C = 1,NCHAIN
          CHAINP=>CGRP%CHAINS(C)
          READ(11,'(A,1X,I12)',IOSTAT=ICHECK) CDUM,NBEAD
          IF (NBEAD.GT.CHAINP%MAXNPT) THEN
             PRINT*, 'ERROR IN READSNAPSHOT: too many beads', NBEAD, CHAINP%MAXNPT
             STOP 1
          ELSE
             CHAINP%NPT = NBEAD
          ENDIF
       
          IF (ICHECK.LT.0) exit ! end of file


          COMPLETE = .TRUE.
          DO B = 1,CHAINP%NPT
             READ(11, '(A,1X,12G20.10)',IOSTAT=ICHECK) CDUM, CHAINP%POS(:,B),CHAINP%UVEC(:,B), TMPVEC          

             IF (ICHECK.LT.0) THEN
                print*, 'bad structure: stopping read.', b
                COMPLETE = .FALSE.
                EXIT ! end of file
             ENDIF
             IF (CDUM.NE.'A') THEN
                PRINT*, 'ERROR IN READSNAPSHOT: bad beadline', CDUM
                stop 2
             ENDIF

             IF (B.LT.CHAINP%NPT.AND.GETPARAM) THEN
                CHAINP%LS(B) = TMPVEC(1)
                CHAINP%LP(B) = TMPVEC(2)
                CHAINP%GAM(B) = TMPVEC(3)
                CHAINP%EPAR(B) = TMPVEC(4)
                CHAINP%EPERP(B) = TMPVEC(5)
                CHAINP%EC(B) = TMPVEC(6)
             ENDIF
          ENDDO
          IF (.NOT.COMPLETE) EXIT
       ENDDO

       IF (COMPLETE) THEN
          SUCCESS = .TRUE.
          CALL COPYCHAINGROUP(CGRP,CGRP2)
       ELSE
          ! incomplete structure; revert to previous one
          CALL COPYCHAINGROUP(CGRP2,CGRP)
          NUM = NUMSAVE
          EXIT
       ENDIF
    ENDDO

    IF (.NOT.SUCCESS) NUM=0
    CLOSE(11)
    
    CALL CLEANUPCHAINGROUP(CGRP2)
  END SUBROUTINE READSNAPSHOTS

  SUBROUTINE READSNAPSHOTSold(CHAINP,INFILE,NUM,SUCCESS)
    ! read in a snapshot file
    ! set the chain configuration to the last complete snapshot
    ! NUM is the additional saved number for the last snapshot (ie: MC step)
    ! chainp must be initialized and have all params set already
    ! COMPLETE indicates whether a complete structure was successfully read
    ! if not successful, coords remain unchanged
    USE CHAINUTIL, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*), INTENT(IN) :: INFILE
    INTEGER, INTENT(OUT) :: NUM
    LOGICAL, INTENT(OUT) :: SUCCESS
    DOUBLE PRECISION :: POSSAVE(3,CHAINP%MAXNPT), UVECSAVE(3,CHAINP%MAXNPT)
    INTEGER :: NUMSAVE
    LOGICAL :: COMPLETE
    INTEGER :: B
    CHARACTER :: CDUM
    DOUBLE PRECISION :: DUM1, DUM2, TMPVEC(6)
    INTEGER :: NBEAD, ICHECK

    POSSAVE = CHAINP%POS; UVECSAVE = CHAINP%UVEC
    NUMSAVE=0

    OPEN(UNIT=11,FILE=INFILE,STATUS='OLD')
    COMPLETE = .FALSE.; SUCCESS = .FALSE.
    DO
       READ(11,'(A,1X,2I12,2G20.10)',IOSTAT=ICHECK) CDUM,NBEAD,NUM,DUM1,DUM2
      
       IF (ICHECK.LT.0) exit ! end of file

       IF (NBEAD.GT.CHAINP%MAXNPT) THEN
          PRINT*, 'ERROR IN READSNAPSHOT: too many beads', NBEAD, CHAINP%MAXNPT
          STOP 1
       ELSE
          CHAINP%NPT = NBEAD
       ENDIF

       COMPLETE = .TRUE.
       DO B = 1,CHAINP%NPT
          READ(11, '(A,1X,12G20.10)',IOSTAT=ICHECK) CDUM, CHAINP%POS(:,B),CHAINP%UVEC(:,B), TMPVEC          
!          print*, 'testx1:', icheck, b
          IF (ICHECK.LT.0) THEN
             print*, 'bad structure: stopping read.', b
             COMPLETE = .FALSE.
             EXIT ! end of file
          ENDIF
          IF (CDUM.NE.'A') THEN
             PRINT*, 'ERROR IN READSNAPSHOT: bad beadline', CDUM
             stop 2
          ENDIF

          IF (B.LT.CHAINP%NPT) THEN
             CHAINP%LS(B) = TMPVEC(1)
             CHAINP%LP(B) = TMPVEC(2)
             CHAINP%GAM(B) = TMPVEC(3)
             CHAINP%EPAR(B) = TMPVEC(4)
             CHAINP%EPERP(B) = TMPVEC(5)
             CHAINP%EC(B) = TMPVEC(6)
          ENDIF
       ENDDO

       IF (COMPLETE) THEN
          SUCCESS = .TRUE.
          POSSAVE = CHAINP%POS; UVECSAVE = CHAINP%UVEC; NUMSAVE = NUM
       ELSE
          ! incomplete structure; revert to previous one
          CHAINP%POS = POSSAVE
          CHAINP%UVEC = UVECSAVE
          NUM = NUMSAVE
          EXIT
       ENDIF
    ENDDO

    IF (.NOT.SUCCESS) NUM=0
    CLOSE(11)

  END SUBROUTINE READSNAPSHOTSOLD
END MODULE INPUTUTIL
