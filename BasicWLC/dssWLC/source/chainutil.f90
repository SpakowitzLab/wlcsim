MODULE CHAINUTIL
  USE KEYS, ONLY : VERBOSE
  IMPLICIT NONE

  ! utilities for defining and dealing with chain object

  TYPE CHAIN 
     INTEGER :: NPT ! current number of beads
     INTEGER :: MAXNPT ! maximal possible number of beads

     ! positions and quaternion orientations for each bead
     ! 1st index is coordinate, 2nd index is bead number
     DOUBLE PRECISION, POINTER :: POS(:,:)
     DOUBLE PRECISION, POINTER :: UVEC(:,:)
     ! SHEARABLE, STRETCHABLE: can the chain stretch and shear?
     ! COUPLED: is there bend-shear coupling?
     ! FINITEXT: prevent segments from stretching beyond contour length
     LOGICAL :: SHEARABLE, STRETCHABLE, COUPLED, FINITEXT
     DOUBLE PRECISION,POINTER :: LS(:) ! segment length
     ! -------- used for brownian dynamics only
     ! DELS is not currently used????
     DOUBLE PRECISION, POINTER :: DELS(:) !different lengths for each segment
     ! single vector of coordinates for the chain
     ! for each bead, lists 3 coords for position and 3 for the U vector
!     DOUBLE PRECISION, POINTER :: COORDS(:) 
!     INTEGER :: NCRD
     ! ---------

     ! energy parameters (FOR EACH SEGMENT)
     ! LP = persistence length
     ! GAM = natural speed
     ! EPERP,EPAR = shear and stretch moduli
     ! EC = bend-shear coupling     
     DOUBLE PRECISION, POINTER :: LP(:),GAM(:),EPERP(:),EPAR(:),EC(:)

     ! FORCE = external force in kT/nm
     ! FINITSHEAR = how much of shear energy component is corrected for finite extension
     ! anything greater than 0 will limit overall segment extension
     DOUBLE PRECISION :: FINITSHEAR, FORCE
     ! spring modulus for U vector normalization constraint
     DOUBLE PRECISION :: CONSTMOD
     DOUBLE PRECISION :: MU
     ! does the chain object have a force associated with it?
     LOGICAL :: HASFORCE

     ! sterics
     LOGICAL :: STERICS ! use sterics?
     DOUBLE PRECISION :: STERRAD,STERRAD2 ! steric radius squared
     DOUBLE PRECISION :: STERMOD, STERICENERGY
     ! when calculating sterics skip this many atoms on either side
     INTEGER :: STERSKIP 
     DOUBLE PRECISION :: MINSEGLEN, MAXSEGLEN
    
     ! store the junction energy associated with each bead
     DOUBLE PRECISION, POINTER :: BEADENERGY(:)
     ! and the energy associated with the external force
     DOUBLE PRECISION :: FORCEENERGY

     ! friction coefficients for bead position and orientation vector
     DOUBLE PRECISION, POINTER :: FRICTR(:), FRICTU(:)

     ! fixed beads (both positions and orientations fixed for now)
     INTEGER :: NFIXBEAD
     INTEGER, POINTER :: FIXBEAD(:)
     LOGICAL, POINTER :: ISFIXED(:)

     ! have arrays been allocated yet?
     LOGICAL :: ARRAYSET = .FALSE.
  END type CHAIN

  TYPE OBSTACLE
     ! object corresponding to a spherical obstacle to the chain
     ! potential is quadratic within a cutoff radius
     DOUBLE PRECISION :: COORDS(3)

     ! cutoff radius and modulus for the potential
     DOUBLE PRECISION :: RAD, MOD
     ! friction coefficient for the obstacle
     DOUBLE PRECISION :: FRICT
     
  END type OBSTACLE

CONTAINS    
  SUBROUTINE RESETLINKERS(CHAINP,CUTOFF,RESET)
    ! reset all linker lengths in the chain to the ground state value
    ! while keeping relative linker orientations unchanged
    ! only reset if any linker length is off by more than CUTOFF Fraction
    USE GENUTIL, ONLY : NORM
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: CUTOFF
    LOGICAL, INTENT(OUT) :: RESET
    INTEGER :: B
    DOUBLE PRECISION :: TANV(3,CHAINP%NPT-1), ntv, COM(3), NEWCOM(3)

    RESET = .FALSE.
    ! get normalized linker vectors
    DO B = 1,CHAINP%NPT-1
       TANV(:,B) = CHAINP%POS(:,B+1)-CHAINP%POS(:,B)
       NTV = NORM(TANV(:,B));
       IF ((NTV-CHAINP%LS(B))/CHAINP%LS(B).GT.CUTOFF) THEN
          RESET = .TRUE.
       ENDIF
       TANV(:,B) = TANV(:,B)/NTV
    END DO

    ! reset all beads; keep center of mass constant
    IF (RESET) THEN
       COM = SUM(CHAINP%POS,2)/CHAINP%NPT
       DO B = 2,CHAINP%NPT
          CHAINP%POS(:,B) = CHAINP%POS(:,B-1)+TANV(:,B-1)*CHAINP%LS(B-1)*CHAINP%GAM(B-1)
       ENDDO
       NEWCOM = SUM(CHAINP%POS,2)/CHAINP%NPT
      
       DO B = 1,CHAINP%NPT
          CHAINP%POS(:,B) = CHAINP%POS(:,B) + COM-NEWCOM
       ENDDO
    ENDIF

  END SUBROUTINE RESETLINKERS
  
  ! -------------- statistics functions ------------
  SUBROUTINE GETCHAINRG(CHAINP,RG)
    ! get radius of gyration of the chain, using the POS coordinates
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: RG
    DOUBLE PRECISION :: RC(3),DIFF(3)
    INTEGER :: B,I

    DO I = 1,3
       RC(I) = SUM(CHAINP%POS(I,1:CHAINP%NPT))/CHAINP%NPT
    ENDDO
    
    RG = 0D0
    DO B = 1,CHAINP%NPT
       DIFF = CHAINP%POS(:,B)-RC
       RG = RG + DOT_PRODUCT(DIFF,DIFF)
    ENDDO
    RG = SQRT(RG/CHAINP%NPT)

  END SUBROUTINE GETCHAINRG


! --------- INPUT/OUTPUT
  SUBROUTINE INPUTSNAPSHOT(CHAINLIST,NCHAIN,FILENAME,NSKIP,NREAD)
    ! read in NCHAIN snapshots from a file
    ! and put coordinates into a list of chains
    ! NREAD returns the number of chains successfully read 
    !(in case there are less than NCHAIN chains specified in the input file)
    ! NSKIP: skip the first few configurations
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NCHAIN,NSKIP
    TYPE(CHAIN), TARGET :: CHAINLIST(NCHAIN)
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(OUT) :: NREAD

    CHARACTER(LEN=*) :: FILENAME
    INTEGER :: TMPI1, TMPI2, NPTREAD, C, B, CHECK
    CHARACTER :: TMPC
    DOUBLE PRECISION :: LSREAD, LPREAD, GAMREAD, EPARREAD, EPERPREAD, ECREAD    
    DOUBLE PRECISION :: POSREAD(3), UVECREAD(3), CHECKPARAMS

    NREAD = 0
    OPEN(UNIT=55,FILE=FILENAME,STATUS='OLD')    
    DO C = 1,NSKIP
       READ(UNIT=55,FMT='(A,1X,2I12)',IOSTAT=CHECK) TMPC,TMPI1,TMPI2
       IF (CHECK.LT.0) THEN
          EXIT
       ELSEIF (CHECK.GT.0) THEN
          PRINT*, 'ERROR in reading in snapshot data. I line', C, TMPC, TMPI1, TMPI2
          print*, CHECK
          STOP 1
       ENDIF

       READ(55,'(A,1X,I12)',IOSTAT=CHECK) TMPC,NPTREAD
       !print*, 'testx3:', c, nptread
       IF (CHECK.LT.0) THEN
          EXIT
       ELSEIF (CHECK.GT.0) THEN
          PRINT*, 'ERROR in reading in snapshot data. C line', C, TMPC, NPTREAD
          STOP 1
       ENDIF

       DO B = 1,NPTREAD
          READ(55,'(A,1X,12G20.10)') TMPC, POSREAD,UVECREAD,&
               & LSREAD, LPREAD, GAMREAD, EPARREAD, &
               & EPERPREAD, ECREAD           
       ENDDO
       PRINT*, 'SKIPPING CONFIG:', C
    ENDDO

    DO C = 1,NCHAIN
       CHAINP=>CHAINLIST(C)
       READ(UNIT=55,FMT='(A,1X,2I12)',IOSTAT=CHECK) TMPC,TMPI1,TMPI2
       IF (CHECK.LT.0) THEN
          EXIT
       ELSEIF (CHECK.GT.0) THEN
          PRINT*, 'ERROR in reading in snapshot data. I line', C, TMPC, TMPI1, TMPI2
          STOP 1
       ENDIF

       READ(55,'(A,1X,I12)',IOSTAT=CHECK) TMPC,NPTREAD
       IF (CHECK.LT.0) THEN
          EXIT
       ELSEIF (CHECK.GT.0) THEN
          PRINT*, 'ERROR in reading in snapshot data. C line', C, TMPC, NPTREAD
          STOP 1
       ENDIF
       IF (NPTREAD.NE.CHAINP%NPT) THEN
          PRINT*, 'ERROR in reading snapshot data. Wrong number of beads expected.', C, NPTREAD, CHAINP%NPT
          STOP 1
       ENDIF

       DO B = 1,CHAINP%NPT
          READ(55,'(A,1X,12G20.10)') TMPC, CHAINP%POS(:,B),CHAINP%UVEC(:,B),&
               & LSREAD, LPREAD, GAMREAD, EPARREAD, &
               & EPERPREAD, ECREAD          
          if (B.LT.CHAINP%NPT) THEN
            ! print*, 'testx1:', LSREAD,chainp%ls(b)
             CHECKPARAMS = MaxVAL([ABS(LSREAD-CHAINP%LS(B)),ABS(LPREAD-CHAINP%LP(B)),&
                  & ABS(GAMREAD-CHAINP%GAM(B)),ABS(EPARREAD-CHAINP%EPAR(B)),&
                  & ABS(EPERPREAD-CHAINP%EPERP(B)),ABS(ECREAD-CHAINP%EC(B))])
             IF (CHECKPARAMS.GT.1D-9) THEN
                PRINT*, 'ERROR in reading snapshot. Wrong segment parameters', B
                print*, LSREAD, LPREAD, GAMREAD, EPARREAD, EPERPREAD, ECREAD
                PRINT*, CHAINP%LS(B), CHAINP%LP(B), CHAINP%EPAR(B), CHAINP%EPERP(B), CHAINP%EC(B)
                STOP 1
             ENDIF
          ENDIF

       ENDDO

       NREAD = NREAD + 1
       !print*, 'SUCCESSFULLY READ CONFIG:', C
    ENDDO
    CLOSE(55)
  END SUBROUTINE INPUTSNAPSHOT

  SUBROUTINE OUTPUTSNAPSHOT(CHAINP,FILENAME,NUM,APPEND,NEXTRA,EXTRADATA)
    ! output a snapshot to file, appending if file already exists
    ! NUM is an extra number to add to the info line
    ! optionally, attach an extra array of data to each bead
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*) :: FILENAME
    INTEGER, INTENT(IN) :: NUM
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER, INTENT(IN), OPTIONAL :: NEXTRA
    DOUBLE PRECISION, POINTER, OPTIONAL :: EXTRADATA(:,:)
    INTEGER :: B 
    CHARACTER*10 :: INTSTR
    CHARACTER*50 :: FMT
    LOGICAL :: WRITEEXTRA

    WRITEEXTRA = PRESENT(NEXTRA).AND.PRESENT(EXTRADATA)

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,POSITION='APPEND')
    ELSE
        OPEN(UNIT=99,FILE=FILENAME,POSITION='REWIND')
     ENDIF

    ! write information line
    WRITE(99,'(A,1X,2I12)') 'I',1,NUM
    WRITE(99,'(A,1X,I12)') 'C',CHAINP%NPT

    ! write bead lines
    IF (WRITEEXTRA) THEN
       WRITE(INTSTR,'(I10)') NEXTRA+12
       FMT = '(A,1X,' // TRIM(ADJUSTL(INTSTR)) // 'G20.10)'
    ELSE
       FMT = '(A,1X,12G20.10)'
    ENDIF

   ! PRINT*, 'TESTX1:', FMT

    DO B = 1,CHAINP%NPT-1  
       IF (WRITEEXTRA) THEN
          WRITE(99,FMT) 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B),&
               & CHAINP%LS(B), CHAINP%LP(B), CHAINP%GAM(B), CHAINP%EPAR(B), &
               & CHAINP%EPERP(B), CHAINP%EC(B),EXTRADATA(1:NEXTRA,B)
       ELSE
           WRITE(99,FMT) 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B),&
               & CHAINP%LS(B), CHAINP%LP(B), CHAINP%GAM(B), CHAINP%EPAR(B), &
               & CHAINP%EPERP(B), CHAINP%EC(B)
       ENDIF
    ENDDO
    B = CHAINP%NPT
    IF (WRITEEXTRA) THEN
       WRITE(99,FMT) 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B), 0,0,0,0,0,0, EXTRADATA(1:NEXTRA,B)
    ELSE
       WRITE(99,'(A,1X,12G20.10)') 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B), 0,0,0,0,0,0  
    ENDIF

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE OUTPUTCHAINOBST(CHAINP,OBP,FILENAMe)
    ! output chain and obstacle coords into a file named FILENAME
    ! based on COORDS array of chain (as used in BD sims)

    IMPLICIT NONE

    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(OBSTACLE), POINTER :: OBP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    INTEGER :: B

    OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')

    ! for each bead dump position coords followed by orientation coords

    DO B = 1,CHAINP%NPT
       WRITE(99,'(A1,1x,6G20.10)') 'B',CHAINP%POS(:,B),CHAINP%UVEC(:,B)
    ENDDO


    ! then dump obstacle coords
    WRITE(99,'(A1,1X,3G20.10)') 'O',OBP%COORDS

    CLOSE(99)
  END SUBROUTINE OUTPUTCHAINOBST


  ! ---------- energy functions for MC -----------
  SUBROUTINE GETENERGY(CHAINP,ENERGY)
    ! get the overall energy associated with the current chain configuration
    ! fill in the full chainp%beadenergy array with bead-specific energies
    ! also include the energy due to external force    
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    INTEGER :: B

    ENERGY = 0D0

    ! get bead-specific energies
    DO B = 2,CHAINP%NPT
       CALL GETBEADENERGY(CHAINP,B,CHAINP%BEADENERGY(B))
       !IF (VERBOSE) PRINT*, 'TESTX3:', B, CHAINP%BEADENERGY(B), &
       !     & sqrt(dot_product(chainp%pos(:,b)-chainp%pos(:,b-1),chainp%pos(:,b)-chainp%pos(:,b-1)))
    ENDDO

    ! get energy from external force
    IF (CHAINP%HASFORCE) THEN
       CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY )
    ELSE
          CHAINP%FORCEENERGY = 0D0
    END IF
    IF (CHAINP%STERICS) CALL GETSTERICENERGY(CHAINP,CHAINP%STERICENERGY)

    ! total energy
    ENERGY = CHAINP%STERICENERGY + SUM(CHAINP%BEADENERGY(1:CHAINP%NPT)) -CHAINP%NPT*CHAINP%MU +CHAINP%FORCEENERGY

  END SUBROUTINE GETENERGY

  SUBROUTINE GETSTERICENERGY(CHAINP,ENERGY)
    ! get steric energy from all pairs of beads
    IMPLICIT NONE
    
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    INTEGER :: B1, B2
    DOUBLE PRECISION :: DIFF(3), ND2

    ENERGY = 0D0
    DO B1 = 1,CHAINP%NPT
       DO B2 = 1,CHAINP%NPT
          IF (ABS(B1-B2).LE.CHAINP%STERSKIP) CYCLE
          DIFF = CHAINP%POS(:,B2)-CHAINP%POS(:,B1)
          ND2 = DOT_PRODUCT(DIFF,DIFF)
          IF (ND2.LT.CHAINP%STERRAD2) THEN
             ENERGY = ENERGY + CHAINP%STERMOD*(CHAINP%STERRAD2-ND2)**2
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE GETSTERICENERGY

  SUBROUTINE CHECKSTERICCLASH(CHAINP,NB1,SET1,NB2,SET2, CLASH)
    ! check for steric clashes between any bead in set1 and any bead in set2
    ! ignore pairs of identical beads

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: NB1, NB2, SET1(NB1), SET2(NB2)
    LOGICAL, INTENT(OUT) :: CLASH
    INTEGER :: C1, C2, B1, B2
    DOUBLE PRECISION :: DIFF(3), D2

    CLASH = .FALSE.

    DO C1 = 1,NB1
       B1 = SET1(C1)
       DO C2 = 1,NB2
          B2 = SET2(C2)
          IF (ABS(B1-B2).LE.CHAINP%STERSKIP) CYCLE
          
          DIFF = CHAINP%POS(:,B1)-CHAINP%POS(:,B2)
          D2 = DOT_PRODUCT(DIFF,DIFF)

          IF (D2.LT.CHAINP%STERRAD2) THEN
             CLASH = .TRUE.
             RETURN
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE CHECKSTERICCLASH

  SUBROUTINE GETFORCEENERGY(CHAINP,ENERGY)
    ! THIS IS OUT OF DATE AND INCOMPATIBLE WITH MANYCHAIN SETUP
    ! get the energy associated with external force
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY

    ENERGY = 0D0
    IF (CHAINP%HASFORCE) THEN
       ENERGY = -(CHAINP%POS(3,CHAINP%NPT)-CHAINP%POS(3,1)-SUM(CHAINP%LS(1:CHAINP%NPT-1)))*CHAINP%FORCE
    ENDIF
    !PRINT*, 'TESTX1:', CHAINP%FORCE, ENERGY
  END SUBROUTINE GETFORCEENERGY

  SUBROUTINE GETBEADQUINT(CHAINP,B,NRHO,NPHI,QINT)
    ! get the partition function for a particular bead, integrated over U vector
    ! assuming the bead is located between two beads with fixed R and U
    USE GENUTIL, ONLY : NORM, PI, CROSS_PRODUCT
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B, NRHO, NPHI
    DOUBLE PRECISION, INTENT(OUT) :: QINT
    DOUBLE PRECISION :: UVECSAVE(3), BEADESAVE
    DOUBLE PRECISION :: XAX(3), YAX(3), ZAX(3), DELRHO, DELPHI, ADD, SR
    DOUBLE PRECISION :: RHO, PHI, ENERGY1, ENERGY2, ENERGY0
    INTEGER :: RC, PC
    
    IF (B.LE.1.OR.B.GE.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN GETBEADQUINT: bad bead', B
       stop 1
    ENDIF
    IF (.NOT.(CHAINP%SHEARABLE.AND.CHAINP%STRETCHABLE)) THEN
       PRINT*, 'ERROR IN GETBEADQUINT: chain must be stretchable and shearable'
       STOP 1
    ENDIF

    ! save current uvector value
    UVECSAVE = CHAINP%UVEC(:,B)
    BEADESAVE = CHAINP%BEADENERGY(B)

    ! set up previous bead coordinate system
    ZAX = CHAINP%UVEC(:,B-1)
    IF (ZAX(2)*ZAX(2)+ZAX(3)*ZAX(3).EQ.0) THEN
       XAX = (/0D0,0D0,1D0/)
       YAX = (/0D0,-1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(ZAX,(/1D0,0D0,0D0/),XAX)
       XAX = XAX/NORM(XAX)
       CALL CROSS_PRODUCT(ZAX,XAX,YAX)       
    ENDIF
    
    ! get current energy to rescale everything by that
    CALL GETBEADENERGY(CHAINP,B,ENERGY1)
    CALL GETBEADENERGY(CHAINP,B+1,ENERGY2)
    ENERGY0 = ENERGY1+ENERGY2
    !print*, 'testx1:', energy0

    !UVEC2 = CHAINP%UVEC(:,B+1)
    !RHO12 = DOT_PRODUCT(UVEC1,UVEC2)

    DELRHO = 2d0/(NRHO-1);
    DELPHI = 2*PI/(NPHI-1);

    QINT = 0D0
    OPEN(UNIT=44,FILE='qintmat.out')
    DO RC = 1,NRHO
       RHO = 1-2*DBLE(RC-1)/(NRHO-1)
       SR = SQRT(1-RHO*RHO)
       DO PC = 1,NPHI
          PHI = 2*PI*DBLE(PC-1)/(NPHI-1)
          
          CHAINP%UVEC(:,B) = SR*COS(PHI)*XAX + SR*SIN(PHI)*YAX + RHO*ZAX
          CALL GETBEADENERGY(CHAINP,B,ENERGY1)
          CALL GETBEADENERGY(CHAINP,B+1,ENERGY2)
          
          ADD = EXP(-(ENERGY1+ENERGY2-ENERGY0))
         !  IF (B.EQ.8) THEN
         !     ! print*, 'testx1:', rc, pc, energy1+energy2-energy0, add*delrho*delphi, qint*delrho*delphi
         !     WRITE(44,*) RHO, PHI, energy1+energy2-energy0, add*delrho*delphi, qint*delrho*delphi
         ! ENDIF
          IF (RC.EQ.1.OR.RC.EQ.NRHO) THEN
             ADD = ADD*0.5
          ENDIF
          IF (PC.EQ.1.OR.PC.EQ.NRHO) THEN
             ADD = ADD*0.5
          ENDIF
          QINT = QINT + ADD
       ENDDO
    ENDDO
    CLOSE(44)

  ! print*, 'testxb:', energy0, delrho, delphi
   ! print*, 'testx2:', qint*delrho*delphi
    QINT = QINT*DELRHO*DELPHI*EXP(-ENERGY0)

    ! restore to original values
    CHAINP%UVEC(:,B) = UVECSAVE
    CHAINP%BEADENERGY(B) = BEADESAVE
  END SUBROUTINE GETBEADQUINT

  SUBROUTINE GETBEADQRINT(CHAINP,B,Npt,QINT)
    ! get the partition function for a particular bead, integrated over R position
    ! assuming the bead is located between two beads with fixed R and U
    USE GENUTIL, ONLY : NORM, PI, CROSS_PRODUCT
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B, NPT
    DOUBLE PRECISION, INTENT(OUT) :: QINT
    DOUBLE PRECISION :: POSSAVE(3), BEADESAVE
    DOUBLE PRECISION :: XAX(3), YAX(3), ZAX(3), ADD, ZRANGE,XYRANGE,DELZ,DELXY
    DOUBLE PRECISION :: X,Y,Z, ENERGY1, ENERGY2, ENERGY0, X0, Y0, Z0
    INTEGER :: XC,ZC,YC

    IF (B.LE.1.OR.B.GE.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN GETBEADQUINT: bad bead', B
       stop 1
    ENDIF
    IF (.NOT.(CHAINP%SHEARABLE.AND.CHAINP%STRETCHABLE)) THEN
       PRINT*, 'ERROR IN GETBEADQUINT: chain must be stretchable and shearable'
       STOP 1
    ENDIF

    ! save current position value
    POSSAVE = CHAINP%POS(:,B)
    BEADESAVE = CHAINP%BEADENERGY(B)

    ! set up previous bead coordinate system
    ZAX = CHAINP%UVEC(:,B-1)
    IF (ZAX(2)*ZAX(2)+ZAX(3)*ZAX(3).EQ.0) THEN
       XAX = (/0D0,0D0,1D0/)
       YAX = (/0D0,-1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(ZAX,(/1D0,0D0,0D0/),XAX)
       XAX = XAX/NORM(XAX)
       CALL CROSS_PRODUCT(ZAX,XAX,YAX)       
    ENDIF

    ! get current energy to rescale everything by that
    CALL GETBEADENERGY(CHAINP,B,ENERGY1)
    CALL GETBEADENERGY(CHAINP,B+1,ENERGY2)
    ENERGY0 = ENERGY1+ENERGY2
    !print*, 'testx1:', b, energy0

    X0 = DOT_PRODUCT(CHAINP%POS(:,B)-CHAINP%POS(:,B-1),XAX)
    Y0 = DOT_PRODUCT(CHAINP%POS(:,B)-CHAINP%POS(:,B-1),YAX)
    Z0 =  DOT_PRODUCT(CHAINP%POS(:,B)-CHAINP%POS(:,B-1),ZAX)

    !how far out to integrate
    ZRANGE = 3/SQRT(CHAINP%EPAR(B))
    XYRANGE = 3/SQRT(CHAINP%EPERP(B))

    DELZ = 2*ZRANGE/(NPT-1);
    DELXY = 2*XYRANGE/(NPT-1);

    QINT = 0D0
    !OPEN(UNIT=44,FILE='qintRmat.out')
    DO XC = 1,NPT
!       X = (-1+2*DBLE(XC-1)/(NPT-1))*XYRANGE
       X = (-1+2*DBLE(XC-1)/(NPT-1))*XYRANGE + Y0
       DO YC = 1,NPT
          !Y = (-1+2*DBLE(YC-1)/(NPT-1))*XYRANGE
          Y = (-1+2*DBLE(YC-1)/(NPT-1))*XYRANGE + Y0
          DO ZC= 1,NPT
             !Z = (-1+2*DBLE(ZC-1)/(NPT-1))*ZRANGE + CHAINP%GAM(B)*CHAINP%LS(B)
             Z = (-1+2*DBLE(ZC-1)/(NPT-1))*ZRANGE + Z0

             CHAINP%POS(:,B) = CHAINP%POS(:,B-1) + X*XAX + Y*YAX+Z*ZAX
             CALL GETBEADENERGY(CHAINP,B,ENERGY1)
             CALL GETBEADENERGY(CHAINP,B+1,ENERGY2)

             ADD = EXP(-(ENERGY1+ENERGY2-ENERGY0))
             !IF (B.EQ.8) THEN
             !print*, 'testx1:', rc, pc, energy1+energy2-energy0, add*delrho*delphi, qint*delrho*delphi
            ! WRITE(44,*) X,Y,Z, energy1+energy2-energy0, add*delXY*DELXY*DELZ, qint*delXY*DELXY*DELZ
             !ENDIF
             IF (ZC.EQ.1.OR.ZC.EQ.NPT) THEN
                ADD = ADD*0.5
             ENDIF
             IF (XC.EQ.1.OR.XC.EQ.NPT) THEN
                ADD = ADD*0.5
             ENDIF
             IF (YC.EQ.1.OR.YC.EQ.NPT) THEN
                ADD = ADD*0.5
             ENDIF
             QINT = QINT + ADD
          ENDDO
       ENDDO
    ENDDO
    !CLOSE(44)

    ! print*, 'testxb:', energy0, delrho, delphi
    QINT = QINT*DELZ*DELXY*DELXY*EXP(-ENERGY0)

    ! restore to original values
    CHAINP%POS(:,B) = POSSAVE
    CHAINP%BEADENERGY(B) = BEADESAVE
  END SUBROUTINE GETBEADQRINT

  SUBROUTINE GETBEADENERGY(CHAINP,B,ENERGY)
    ! get the energy associated with a particular bead
    ! (for the segment preceeding the bead)

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY    
    DOUBLE PRECISION :: DPOS(3),DU(3),UU,RPAR,RPERP(3),RPERP2
    DOUBLE PRECISION :: ROTMAT(3,3), NDR,NDRSHIFT,RMAX,RMAX2
    DOUBLE PRECISION :: DEL, GAM, LP,EPAR,EPERP,EC
    DOUBLE PRECISION :: D1D2,ND1,ND2,DPOS2(3)
    ! DOUBLE PRECISION :: FINSHIFT
    !FINSHIFT = 1D-3

    IF(B.EQ.1) THEN
       ENERGY=0D0
       RETURN
    END IF

    DEL = CHAINP%LS(B-1); LP = CHAINP%LP(B-1); EPAR = CHAINP%EPAR(B-1)
    EPERP = CHAINP%EPERP(B-1); GAM = CHAINP%GAM(B-1); EC = CHAINP%EC(B-1)

    ! bend angle
    DU = CHAINP%UVEC(:,B)-CHAINP%UVEC(:,B-1)
    UU = DOT_PRODUCT(CHAINP%UVEC(:,B),CHAINP%UVEC(:,B-1))    

    ! displacement relative to previous coord system

    DPOS = CHAINP%POS(:,B)-CHAINP%POS(:,B-1)

    IF (CHAINP%SHEARABLE) THEN
       RPAR = DOT_PRODUCT(DPOS,CHAINP%UVEC(:,B-1))
       RPERP = DPOS - RPAR*CHAINP%UVEC(:,B-1)
       RPERP2 = DOT_PRODUCT(RPERP,RPERP)
       ENERGY = LP/DEL*(1-UU)
    ELSEIF (CHAINP%STRETCHABLE) THEN
       IF (B.LT.CHAINP%NPT) THEN
          DPOS2 = CHAINP%POS(:,B+1)-CHAINP%POS(:,B)
          ND1 = SQRT(DOT_PRODUCT(DPOS,DPOS))
          ND2 = SQRT(DOT_PRODUCT(DPOS2,DPOS2))
          D1D2 = DOT_PRODUCT(DPOS,DPOS2)
          ENERGY = LP/DEL*(1-D1D2/ND1/ND2)
       ENDIF
       RPAR = ND1
    ELSE
       IF (B.LT.CHAINP%NPT) THEN
          DPOS2 = CHAINP%POS(:,B+1)-CHAINP%POS(:,B)
          D1D2 = DOT_PRODUCT(DPOS,DPOS2)
          ENERGY = LP/DEL*(1-D1D2/(DEL*GAM)**2)
       ELSE
          ENERGY = 0D0
       ENDIF
    ENDIF


    IF (CHAINP%STRETCHABLE) THEN
       IF (CHAINP%FINITEXT.AND.RPAR.GT.DEL*GAM) THEN
          IF ((RPAR-DEL*GAM)**2.GE.(DEL-DEL*GAM)**2) THEN
             ENERGY = ENERGY + HUGE(1D0)
          ELSE
             ENERGY = ENERGY - EPAR*CHAINP%FINITSHEAR/2*DEL*(1-GAM)**2* & 
                  & LOG(1-(RPAR-DEL*GAM)**2/(DEL-DEL*GAM)**2)             
             ENERGY = ENERGY + EPAR*(1-CHAINP%FINITSHEAR)/2/DEL* &
                  & (RPAR-DEL*GAM)**2        
          ENDIF
       ELSE
          ENERGY = ENERGY + EPAR/2/DEL*(RPAR - DEL*GAM)**2
       ENDIF
    ENDIF

    IF (CHAINP%SHEARABLE) THEN
       IF (CHAINP%FINITEXT) THEN
          IF (RPAR**2+RPERP2.GT.DEL**2) THEN
             ENERGY = ENERGY + HUGE(1D0)
             RETURN
          ENDIF
          RMAX2 = DEL**2-RPAR**2;
          ENERGY = ENERGY - EPERP*CHAINP%FINITSHEAR/2/DEL*RMAX2*LOG(1-(RPERP2)/RMAX2) 
          ENERGY = ENERGY + EPERP*(1-CHAINP%FINITSHEAR)/2/DEL*(RPERP2)          
       ELSE
          ENERGY = ENERGY + EPERP/2/DEL*(RPERP2) 
       ENDIF

       IF (CHAINP%COUPLED) THEN
          ENERGY = ENERGY + EC/DEL*DOT_PRODUCT(RPERP,DU)
       ENDIF
    ENDIF

  END SUBROUTINE GETBEADENERGY

! ------------ SETUP FUNCTIONS ---------------
  SUBROUTINE INITIALIZECHAIN(CHAINP,RANDOMIZE,RANDMAG)
    ! initialize chain configuration, straight in the z direction or random
    USE mt19937, ONLY : GRND
    USE QUATUTIL, ONLY : ROTANGAX
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, INTENT(IN) :: RANDOMIZE
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: RANDMAG(2)
    INTEGER :: N    , B
    DOUBLE PRECISION :: ENERGY
    DOUBLE PRECISION :: RVEC(3), AX(3), THETA,ROTMAT(3,3)

    IF (RANDOMIZE.AND..NOT.PRESENT(RANDMAG)) THEN
       PRINT*, 'ERROR IN INITIALIZECHAIN: no randomize magnitude supplied'
       STOP 1
    ENDIF

    IF (RANDOMIZE) THEN
       ! Randomly perturb configuration
       DO B = 1,CHAINP%NPT
          RVEC = (/GRND(),GRND(),GRND()/)*2*RANDMAG(1) - (/1,1,1/)*RANDMAG(1)
          IF (B.EQ.1) THEN
             CHAINP%POS(:,B) = RVEC
          ELSE
             IF (.NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
                ! make sure segment length is correct
                RVEC = (/0D0,0D0,CHAINP%LS(B-1)*CHAINP%GAM(B-1)/)+RVEC
                RVEC = RVEC/SQRT(SUM(RVEC**2))*CHAINP%LS(B-1)*CHAINP%GAM(B-1)
                CHAINP%POS(:,B) = CHAINP%POS(:,B-1)+RVEC
             ELSE
                CHAINP%POS(:,B) = CHAINP%POS(:,B-1)+&
                     &(/0D0,0D0,CHAINP%LS(B-1)*CHAINP%GAM(B-1)/)+RVEC
             ENDIF
          ENDIF

          ! random axis
          AX = (/GRND(),GRND(),GRND()/)
          AX = AX/SQRT(DOT_PRODUCT(AX,AX))

          ! random theta
          THETA = GRND()*RANDMAG(2)
          CALL ROTANGAX(THETA,AX,(/0D0,0D0,1D0/),CHAINP%UVEC(:,B),.TRUE.,ROTMAT)
       END DO
    ELSE
       ! set up a straight chain
       CHAINP%POS(1:2,:) = 0D0;
       CHAINP%POS(3,1) = 0D0;
       DO B = 2,CHAINP%NPT
          CHAINP%POS(3,B) = CHAINP%POS(3,B-1)+CHAINP%LS(B-1)*CHAINP%GAM(B-1)
       ENDDO
       CHAINP%UVEC = 0D0
       CHAINP%UVEC(3,:) = 1D0
    END IF

    ! Coords for brownian dynamics
    ! DO B = 1,CHAINP%NPT
    !    CHAINP%COORDS(6*(B-1)+1:6*(B-1)+3) = CHAINP%POS(:,B)
    !    CHAINP%COORDS(6*(B-1)+4:6*(B-1)+6) = CHAINP%UVEC(:,B)
    ! END DO


    CALL GETENERGY(CHAINP,ENERGY)

  END SUBROUTINE INITIALIZECHAIN

  SUBROUTINE COPYCHAIN(CHAINP1,CHAINP2)
    ! copy all parameters over from chainp1 to chainp2
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP1, CHAINP2
    
    IF (.NOT.CHAINP1%ARRAYSET.OR..NOT.CHAINP2%ARRAYSET) THEN
       PRINT*, 'ERROR: CANNOT COPY CHAINS UNLESS BOTH HAVE INITIALIZED ARRAYS'
       PRINT*, CHAINP1%ARRAYSET, CHAINP2%ARRAYSET
       STOP 1
    ENDIF

    IF (CHAINP1%MAXNPT.NE.CHAINP2%MAXNPT) THEN
       PRINT*, 'ERROR: CANNOT COPY CHAINS OF DIFFERENT SIZE'
       stop 1 
    ENDIF

    CHAINP2%LS = CHAINP1%LS
    CHAINP2%STRETCHABLE = CHAINP1%STRETCHABLE
    CHAINP2%SHEARABLE = CHAINP1%SHEARABLE
    CHAINP2%COUPLED = CHAINP1%COUPLED
    CHAINP2%LP = CHAINP1%LP
    CHAINP2%GAM = CHAINP1%GAM
    CHAINP2%EC = CHAINP1%EC
    CHAINP2%EPERP = CHAINP1%EPERP
    CHAINP2%EPAR = CHAINP1%EPAR
    CHAINP2%FORCE = CHAINP1%FORCE
    CHAINP2%HASFORCE = CHAINP1%HASFORCE

    CHAINP2%POS = CHAINP1%POS
    CHAINP2%UVEC = CHAINP1%UVEC
    CHAINP2%BEADENERGY = CHAINP1%BEADENERGY
    CHAINP2%FORCEENERGY = CHAINP1%FORCEENERGY
    CHAINP2%FINITEXT = CHAINP1%FINITEXT
    CHAINP2%FINITSHEAR = CHAINP1%FINITSHEAR

    CHAINP2%CONSTMOD = CHAINP1%CONSTMOD
    CHAINP2%DELS = CHAINP1%DELS
!    CHAINP2%COORDS = CHAINP1%COORDS
    CHAINP2%FRICTR = CHAINP1%FRICTR
    CHAINP2%FRICTU = CHAINP1%FRICTU

    CHAINP2%STERICS = CHAINP1%STERICS
    CHAINP2%STERRAD2 = CHAINP1%STERRAD2
    CHAINP2%STERRAD = CHAINP1%STERRAD
    CHAINP2%STERSKIP = CHAINP2%STERSKIP

    CHAINP2%NPT = CHAINP1%NPT
    CHAINP2%MAXNPT = CHAINP1%MAXNPT
    CHAINP2%MINSEGLEN = CHAINP1%MINSEGLEN
    CHAINP2%MAXSEGLEN = CHAINP1%MAXSEGLEN

    CHAINP2%STERMOD = CHAINP1%STERMOD
    CHAINP2%STERICENERGY = CHAINP1%STERICENERGY
    CHAINP2%MU = CHAINP1%MU
  END SUBROUTINE COPYCHAIN

  SUBROUTINE SETOBSTACLEPARAMS(OBP)
    ! set parameters for obstacle object
    USE KEYS, ONLY : FRICTOB,RADOB,MODOB
    IMPLICIT NONE
    TYPE(OBSTACLE), POINTER :: OBP

    OBP%RAD = RADOB
    OBP%MOD = MODOB
    OBP%FRICT= FRICTOB
  END SUBROUTINE SETOBSTACLEPARAMS

  SUBROUTINE SETCHAINPARAMS(CHAINP,FORCE)
    ! set parameters for the chain using keywords from KEYS module
    USE KEYS, ONLY : LS, SHEARABLE,STRETCHABLE,COUPLED,LP,EC,EPERP,&
         & EPAR,GAM,FINITEXT, FINITSHEAR, CONSTMOD,FRICTR,FRICTU,&
         & USESTERICS, STERRAD, STERSKIP, STARTNPT,MINSEGLEN,MAXSEGLEN,STERMOD,&
         & MU,NFIXBEAD,FIXBEAD, NEDGESEG, EDGELS, EDGELP, EDGEGAM, EDGEEPERP,&
         & EDGEEPAR,EDGEEC, FRICTPERLEN

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: FORCE
    INTEGER :: I, b

    CHAINP%NPT = STARTNPT
    CHAINP%NFIXBEAD = NFIXBEAD
    IF (.NOT.CHAINP%ARRAYSET.AND.NFIXBEAD.GT.0) THEN
       PRINT*, 'ERROR IN SETCHAINPARAMS. ARRAYS HAVE NOT BEEN SET.'
       STOP 1
    ENDIF
    IF (NFIXBEAD.GT.0) THEN
       CHAINP%FIXBEAD(1:NFIXBEAD) = FIXBEAD(1:NFIXBEAD,1)
       DO I = 1,NFIXBEAD
          IF (CHAINP%FIXBEAD(I).LT.1.OR.CHAINP%FIXBEAD(I).GT.CHAINP%NPT) THEN
             PRINT*, 'ERROR: BAD FIXBEAD', CHAINP%FIXBEAD(I)
             STOP 1
          ENDIF
          CHAINP%ISFIXED(CHAINP%FIXBEAD(I)) = .TRUE.
       ENDDO
    ENDIF
    CHAINP%LS = LS
    CHAINP%DELS = LS
    CHAINP%STRETCHABLE = STRETCHABLE
    CHAINP%SHEARABLE = SHEARABLE
    CHAINP%COUPLED = COUPLED
    CHAINP%LP = LP
    CHAINP%GAM = GAM
    CHAINP%EC = EC
    CHAINP%EPERP = EPERP
    CHAINP%EPAR = EPAR
    IF (NEDGESEG.GT.0) THEN
       DO B = 1,NEDGESEG
          CHAINP%LS(B) = EDGELS
          CHAINP%DELS(B) = EDGELS
          CHAINP%LP(B) = EDGELP
          CHAINP%GAM(B) = EDGEGAM
          CHAINP%EC(B) = EDGEEC
          CHAINP%EPERP(B) = EDGEEPERP
          CHAINP%EPAR(B) = EDGEEPAR
       END DO
       DO B = CHAINP%NPT-NEDGESEG,CHAINP%NPT-1
          CHAINP%LS(B) = EDGELS
          CHAINP%DELS(B) = EDGELS
          CHAINP%LP(B) = EDGELP
          CHAINP%GAM(B) = EDGEGAM
          CHAINP%EC(B) = EDGEEC
          CHAINP%EPERP(B) = EDGEEPERP
          CHAINP%EPAR(B) = EDGEEPAR
       END DO
    ENDIF

    IF (PRESENT(FORCE)) THEN
       CHAINP%FORCE = FORCE
       CHAINP%HASFORCE = .TRUE.
    ELSE
       CHAINP%FORCE = 0D0
    ENDIF
    CHAINP%FINITEXT = FINITEXT
    CHAINP%FINITSHEAR = FINITSHEAR
    CHAINP%CONSTMOD = CONSTMOD
    IF (FRICTPERLEN) THEN
       ! keyword parameters are friction per unit length
       ! different friction for edge beads; different friction if different segment lengths
       CHAINP%FRICTR(1) = FRICTR*CHAINP%LS(1)/2
       CHAINP%FRICTU(1) = FRICTU*CHAINP%LS(1)
       CHAINP%FRICTR(CHAINP%NPT) = FRICTR*CHAINP%LS(CHAINP%NPT-1)/2
       CHAINP%FRICTU(CHAINP%NPT) = FRICTU*CHAINP%LS(CHAINP%NPT-1)
       DO B = 2,CHAINP%NPT-1
          CHAINP%FRICTR(B) = FRICTR*(CHAINP%LS(B-1)/2+CHAINP%LS(B)/2)
          CHAINP%FRICTU(B) = FRICTU*CHAINP%LS(B)           
       ENDDO
    ELSE
       ! keyword parameters are just the friction for each bead
       CHAINP%FRICTR = FRICTR
       CHAINP%FRICTU = FRICTU
    ENDIF
    CHAINP%STERICS = USESTERICS
    CHAINP%STERRAD = STERRAD
    CHAINP%STERRAD2 = STERRAD*STERRAD
    CHAINP%STERSKIP = STERSKIP
    CHAINP%MINSEGLEN = MINSEGLEN
    CHAINP%MAXSEGLEN = MAXSEGLEN
    CHAINP%STERMOD = STERMOD
    CHAINP%MU = MU

  END SUBROUTINE SETCHAINPARAMS

  SUBROUTINE SETUPCHAIN(CHAINP,MAXNPT)
    USE KEYS, ONLY : MAXFIXBEAD
    ! set up the arrays for a chain object
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: maxNPT

    CHAINP%MAXNPT = MAXNPT
    CHAINP%HASFORCE = .FALSE.

    ALLOCATE(CHAINP%POS(3,MAXNPT),CHAINP%UVEC(3,MAXNPT),CHAINP%BEADENERGY(MAXNPT))

    ALLOCATE(CHAINP%LP(MAXNPT-1),CHAINP%LS(MAXNPT-1), CHAINP%EPAR(MAXNPT-1), &
         & CHAINP%EPERP(MAXNPT-1), CHAINP%EC(MAXNPT-1), CHAINP%GAM(MAXNPT-1))

!    CHAINP%NCRD = MAXNPT*6
!    ALLOCATE(CHAINP%COORDS(CHAINP%NCRD))
    ALLOCATE(CHAINP%DELS(MAXNPT-1),CHAINP%FIXBEAD(MAXFIXBEAD),CHAINP%ISFIXED(MAXNPT))

    ALLOCATE(CHAINP%FRICTR(MAXNPT),CHAINP%FRICTU(MAXNPT))

    CHAINP%BEADENERGY = 0D0
    CHAINP%STERICENERGY=0D0
    CHAINP%ISFIXED = .FALSE.
    CHAINP%NFIXBEAD = 0
    CHAINP%ARRAYSET = .TRUE.    
  END SUBROUTINE SETUPCHAIN

  SUBROUTINE CLEANUPCHAIN(CHAINP)
    ! clean up allocatable arrays for the object 
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP

    DEALLOCATE(CHAINP%POS,CHAINP%UVEC,CHAINP%DELS,&
         & CHAINP%BEADENERGY, CHAINP%LP, CHAINP%LS, CHAINP%gam, CHAINP%EPAR, &
         & CHAINP%EPERP, CHAINP%EC,CHAINP%FIXBEAD,CHAINP%ISFIXED, &
         & CHAINP%FRICTR, CHAINP%FRICTU)
    
  END SUBROUTINE CLEANUPCHAIN
END MODULE CHAINUTIL
