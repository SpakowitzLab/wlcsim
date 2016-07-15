MODULE MANYCHAINS
  USE CHAINUTIL
  IMPLICIT NONE

  TYPE CHAINGROUP
     ! object corresponding to an interconnected group (mesh) of chains
     ! number of chains in group
     INTEGER :: NCHAIN=0
     ! array of chain objects
     TYPE(CHAIN), POINTER :: CHAINS(:)

     ! list of interconnections between chains
     ! CONNECT(n,:) = bead # 1, chain # 1, bead # 2, chain # 2 (for nth connection)
     INTEGER :: NCONNECT = 0
     INTEGER, POINTER :: CONNECT(:,:)
     ! strength of position connections (harmonic springs)
     ! and strength of orientation connections (harmonic springs in angle)
     DOUBLE PRECISION :: CONPOSMOD,CONUVECMOD 
     ! for each connection are the positions and/or the orientations connected?
     LOGICAL :: CONNECTPOS, CONNECTUVEC
     ! are connections rigidly fixed?
     LOGICAL :: FIXCONNECT
     ! for each bead and chain, list ONE other bead and chain it's connected to
     ! -1 if not connected    
     INTEGER, POINTER :: CONBEAD(:,:,:)

     ! list of fixed beads
     ! last index 1 = fixed position; last index 2 = fixed orientation
     LOGICAL, POINTER :: FIXBEAD(:,:,:)

     ! list of forces on the group
     INTEGER :: NFORCE
     INTEGER, POINTER :: FORCEBEAD(:,:) ! what bead on what chain
     DOUBLE PRECISION, POINTER :: FORCE(:,:) ! force vector
  END TYPE CHAINGROUP

CONTAINS
  SUBROUTINE GROUPSNAPSHOT(CGRP,FILENAME,NUM,APPEND)
    ! output snapshot containing all the different chains
    ! NUM is an extra number given in the file (used to keep track of MC step)
    ! APPEND is whether to append to file or rewrite it

    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    CHARACTER(LEN=*) :: FILENAME
    INTEGER, INTENT(IN) :: NUM
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER :: B,c
    TYPE(CHAIN), POINTER :: CHAINP

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,POSITION='APPEND')
    ELSE
       OPEN(UNIT=99,FILE=FILENAME,POSITION='REWIND')
    ENDIF

    ! write information line
    WRITE(99,'(A,1X,2I12)') 'I',CGRP%NCHAIN,NUM

    DO C = 1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)
       ! chain information line
       WRITE(99,'(A,1X,I12)') 'C',CHAINP%NPT

       ! write bead lines
       DO B = 1,CHAINP%NPT-1       
          WRITE(99,'(A,1X,12G20.10)') 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B),&
               & CHAINP%LS(B), CHAINP%LP(B), CHAINP%GAM(B), CHAINP%EPAR(B), &
               & CHAINP%EPERP(B), CHAINP%EC(B)
       ENDDO
       B = CHAINP%NPT
       WRITE(99,'(A,1X,12G20.10)') 'A', CHAINP%POS(:,B),CHAINP%UVEC(:,B), 0,0,0,0,0,0  
    ENDDO
    CLOSE(99)

  END SUBROUTINE GROUPSNAPSHOT

  SUBROUTINE BEADALLENERGY(CGRP,C,B,ENERGY,RECALC)
    ! get all energy components involving bead B on chain C
    ! if recalc is true recalculate everything; otherwise use saved BEADENERGY arrays
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: B, C
    LOGICAL, INTENT(IN) :: RECALC
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: CONENERGY
    INTEGER :: F, CON

    CHAINP=>CGRP%CHAINS(C)

    ENERGY = 0
    
    ! Elastic energies
    IF (B.LT.CHAINP%NPT) THEN
       IF (RECALC) THEN
          CALL GETBEADENERGY(CHAINP,B+1,CHAINP%BEADENERGY(B+1))
       ENDIF
       ENERGY = ENERGY+CHAINP%BEADENERGY(B+1)
    ENDIF    
    IF (B.GT.1) THEN
       IF (RECALC) THEN
          CALL GETBEADENERGY(CHAINP,B,CHAINP%BEADENERGY(B))
       ENDIF
       ENERGY = ENERGY + CHAINP%BEADENERGY(B)
    ENDIF
   
    ! connection energies
    IF (.NOT.CGRP%FIXCONNECT) THEN
       DO CON = 1,CGRP%NCONNECT
          IF (CGRP%CONNECT(CON,1).EQ.B.AND.CGRP%CONNECT(CON,2).EQ.C &
               & .OR.CGRP%CONNECT(CON,3).EQ.B.AND.CGRP%CONNECT(CON,4).EQ.C) THEN
             CALL GETCONNECTENERGY(CGRP,CON,CONENERGY)
             ENERGY = ENERGY + CONENERGY
          ENDIF
       ENDDO
    ENDIF

    ! force energy
    DO F = 1,CGRP%NFORCE
       IF (CGRP%FORCEBEAD(F,2).EQ.C.AND.CGRP%FORCEBEAD(F,1).EQ.B) THEN          
          ENERGY = ENERGY - DOT_PRODUCT(CHAINP%POS(:,B),CGRP%FORCE(F,:))
       ENDIF
    ENDDO

    IF (CHAINP%STERICS) THEN
       PRINT*, 'ERROR IN BEADDALLENERGY: not yet set up for sterics'
       stop 1
    endif
  END SUBROUTINE BEADALLENERGY

  SUBROUTINE GROUPENERGY(CGRP,TOTENERGY)
    ! get total energy for a group of chains
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    DOUBLE PRECISION, INTENT(OUT) :: TOTENERGY
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: C,B,F,C1,C2,B1,B2,CON
    DOUBLE PRECISION :: DIFF(3), DIST2,RHO, POS1(3), POS2(3), UVEC1(3), UVEC2(3),CONENERGY, FORCENERGY

    TOTENERGY = 0

    ! internal chain energies
    DO C = 1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)
       ! get bead-specific energies       
       DO B = 2,CHAINP%NPT
          CALL GETBEADENERGY(CHAINP,B,CHAINP%BEADENERGY(B))
       ENDDO

       IF (CHAINP%STERICS) THEN
          PRINT*, 'ERROR: GROUPENERGY WITH STERICS NOT YET SET UP'
          STOP 1
       ENDIF

       TOTENERGY = TOTENERGY + SUM(CHAINP%BEADENERGY(1:CHAINP%NPT))-CHAINP%NPT*CHAINP%MU
    ENDDO
    
    ! connection energies
    IF (.NOT.CGRP%FIXCONNECT) THEN
       DO CON = 1,CGRP%NCONNECT
          CALL GETCONNECTENERGY(CGRP,CON,CONENERGY)
          !       PRINT*, 'TESTX1:',CON, CGRP%CONNECT(CON,:),CONENERGY
          TOTENERGY = TOTENERGY + CONENERGY
       ENDDO
    ENDIF

    ! force-based energies
    FORCENERGY = 0D0
    DO F = 1,CGRP%NFORCE       
       C = CGRP%FORCEBEAD(F,2); B = CGRP%FORCEBEAD(F,1)
       TOTENERGY = TOTENERGY - DOT_PRODUCT(CGRP%CHAINS(C)%POS(:,B),CGRP%FORCE(F,:))
       FORCENERGY = FORCENERGY - DOT_PRODUCT(CGRP%CHAINS(C)%POS(:,B),CGRP%FORCE(F,:))
    ENDDO
    !PRINT*, 'TESTX1:', FORCENERGY
  END SUBROUTINE GROUPENERGY

  SUBROUTINE GETCONNECTENERGY(CGRP,CON,energy)
    ! get the connection energy for the given connection number
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: CON
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    INTEGER :: C1,B1,C2,B2
    DOUBLE PRECISION :: POS1(3), POS2(3), DIFF(3), DIST2, RHO, UVEC1(3), UVEC2(3)
    B1 = CGRP%CONNECT(CON,1); B2 = CGRP%CONNECT(CON,3)
    C1 = CGRP%CONNECT(CON,2); C2 = CGRP%CONNECT(CON,4)

    ENERGY=0D0
    IF (CGRP%CONNECTPOS) THEN
       POS1 = CGRP%CHAINS(C1)%POS(:,B1)
       POS2 = CGRP%CHAINS(C2)%POS(:,B2)
       DIFF = POS2-POS1
       DIST2 =DOT_PRODUCT(DIFF,DIFF)       
       ENERGY = ENERGY + CGRP%CONPOSMOD/2*DIST2
    ENDIF
    IF (CGRP%CONNECTUVEC) THEN
       UVEC1 = CGRP%CHAINS(C1)%UVEC(:,B1)
       UVEC2 = CGRP%CHAINS(C2)%UVEC(:,B2)
       RHO = DOT_PRODUCT(UVEC1,UVEC2)
       ENERGY = ENERGY + CGRP%CONUVECMOD*(1-RHO)          
    ENDIF

  END SUBROUTINE GETCONNECTENERGY

  SUBROUTINE APPLYDEFORM(CGRP,SHEARMAT)   
    ! Apply deformation according to given matrix
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    DOUBLE PRECISION, INTENT(IN)  :: SHEARMAT(3,3)
    INTEGER :: C,B,I
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: TMP(3)

    DO C = 1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)
       DO B = 1,CHAINP%NPT
          DO I = 1,3
             TMP(I) = DOT_PRODUCT(SHEARMAT(I,:),CHAINP%POS(:,B))
          ENDDO
          CHAINP%POS(:,B) = TMP
       ENDDO
    ENDDO
  END SUBROUTINE APPLYDEFORM

  SUBROUTINE INITIALIZESQUARELATTICE(CGRP,DIST)
    ! Initialize a group of chains on a square lattice
    ! with spacing DIST
    ! subjected to shear deformation defined by SHEARMAT
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    DOUBLE PRECISION, INTENT(IN) :: DIST
    DOUBLE PRECISION :: TMP(3)
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: C, B, NC, I

    IF (CGRP%NCHAIN.EQ.0) THEN
       PRINT*, 'ERROR: No chains to initialize.'
       STOP 1
    ENDIF
       
    IF (MOD(CGRP%NCHAIN,2).NE.0) THEN
       PRINT*, 'ERROR: cannot have square lattice with odd number of chains.', CGRP%NCHAIN
       STOP 1
    ENDIF
    NC = CGRP%NCHAIN/2

    ! vertical chains
    DO C = 1,NC
       CHAINP=>CGRP%CHAINS(C)
       CHAINP%POS(1,:) = 0D0;
       CHAINP%POS(2,:) = (C-1)*DIST
       CHAINP%POS(3,1) = 0D0;
       DO B = 2,CHAINP%NPT
          CHAINP%POS(3,B) = CHAINP%POS(3,B-1) + CHAINP%LS(B-1)*CHAINP%GAM(B-1)
       ENDDO
       CHAINP%UVEC = 0D0
       CHAINP%UVEC(3,:) = 1D0
    ENDDO

    ! horizontal chains
    DO C = NC+1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)
       CHAINP%POS(1,:) = 0D0;
       CHAINP%POS(3,:) = (C-NC-1)*DIST
       CHAINP%POS(2,:) = 0D0;
       DO B = 2,CHAINP%NPT
          CHAINP%POS(2,B) = CHAINP%POS(2,B-1) + CHAINP%LS(B-1)*CHAINP%GAM(B-1)
       ENDDO
       CHAINP%UVEC = 0D0
       IF (CGRP%CONNECTUVEC) THEN
          CHAINP%UVEC(3,:) = 1D0
       ELSE
          CHAINP%UVEC(2,:) = 1D0
       ENDIF
    ENDDO
   
  END SUBROUTINE INITIALIZESQUARELATTICE

  SUBROUTINE INITIALIZEDIAMONDLATTICE(CGRP,SEGLEN,NDIAMOND,LENDIAMOND,WIDTHDIAMOND)
    ! initialize positions/orientations of beads in a diamond lattice
    ! lendiamond is the length of a diamond side in number of segments
    ! widthdiamond is the total width of a diamond in space
    ! seglen is the width of a segment
    ! NDIAMOND is the width and height of the mesh in terms of number of diamonds
    USE KEYS, ONLY : STARTCOLLAPSE
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: NDIAMOND(2),LENDIAMOND
    DOUBLE PRECISION :: SEGLEN,WIDTHDIAMOND
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: NL, C, B
    DOUBLE PRECISION :: HEIGHT, YSTP, ZSTP, STARTY, STARTZ, STARTPT(2)

    ! number of left leaning chains
    NL = NDIAMOND(1)+NDIAMOND(2)-1

    IF (CGRP%NCHAIN.NE.2*NL) THEN
       PRINT*, 'ERROR IN INITIALIZEDIAMONDLATTICE: wrong number of chains!'
       STOP 1
    ENDIF

    ! height of half a diamond
    HEIGHT = SQRT((SEGLEN*LENDIAMOND)**2 - (WIDTHDIAMOND/2)**2)
    YSTP = WIDTHDIAMOND/2/LENDIAMOND
    ZSTP = HEIGHT/LENDIAMOND

    ! place left leaning chains
    DO C = 1,NL
       CHAINP=>CGRP%CHAINS(C)

       IF (C.LE.NDIAMOND(1)) THEN
          STARTY = WIDTHDIAMOND*C
          STARTZ = 0D0
       ELSE
          STARTY = WIDTHDIAMOND*NDIAMOND(1)
          STARTZ = 2*HEIGHT*(C-NDIAMOND(1))
       ENDIF
       STARTPT = (/STARTY,STARTZ/)
       
       CHAINP%POS = 0D0
       DO B = 1,CHAINP%NPT          
          CHAINP%POS(2:3,B) = STARTPT + (B-1)*(/-YSTP,ZSTP/)
          CHAINP%UVEC(:,B) = (/0D0,0D0,1D0/)
       ENDDO
       IF (STARTCOLLAPSE) CHAINP%POS(2,:) = 0D0
    ENDDO

    ! place right leaning chains
    DO C = NL+1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)
       
       IF (C-NL.LE.NDIAMOND(2)) THEN
          STARTY = 0
          STARTZ = 2*HEIGHT*(NDIAMOND(2)-C+NL)
       ELSE
          STARTY = WIDTHDIAMOND*(C-NL-NDIAMOND(2))
          STARTZ = 0
       ENDIF
       STARTPT = (/STARTY,STARTZ/)

       CHAINP%POS = 0D0
       DO B = 1,CHAINP%NPT          
          CHAINP%POS(2:3,B) = STARTPT + (B-1)*(/YSTP,ZSTP/)
          CHAINP%UVEC(:,B) = (/0D0,0D0,1D0/)
       ENDDO
       IF (STARTCOLLAPSE) CHAINP%POS(2,:) = 0D0
    ENDDO
    
  END SUBROUTINE INITIALIZEDIAMONDLATTICE

  SUBROUTINE SETUPDIAMONDLATTICE(CGRP,STARTCONNECT,NDIAMOND,LENDIAMOND,WIDTHDIAMOND,FIXBOUNDARY)
    ! set up a diamond lattice of connections
    ! update chain lengths appropriately
    ! start with the specified connection index, to allow additional connections
    ! NDIAMOND is the width and height in # of diamonds
    ! LENDIAMOND is the length of a diamond side in chain segments
    ! WIDTHDIAMOND is the width of each diamond in space
    ! also initialize the chain
    ! if FIXBOUNDARY, also fix the beads on the boundary (top/bottom and/or sides) of the lattice    

    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: STARTCONNECT, NDIAMOND(2), LENDIAMOND
    DOUBLE PRECISION, INTENT(IN) :: WIDTHDIAMOND
    LOGICAL, INTENT(IN) :: FIXBOUNDARY(4)
    INTEGER :: NL, MI, MA, CC, C, C2, B, B2, D, NINTER, CON, CURB2
    LOGICAL :: TOP1, SIDE1, TOPEND, SIDEEND

    ! number of left leaning chains
    NL = NDIAMOND(1)+NDIAMOND(2)-1

    IF (CGRP%NCHAIN.NE.2*NL) THEN
       PRINT*, 'ERROR IN SETUPDIAMONDLATTICE: wrong number of chains!'
       STOP 1
    ENDIF
    
    MI = MINVAL(NDIAMOND); MA = MAXVAL(NDIAMOND)
    CON = STARTCONNECT-1
    DO CC = 1,CGRP%NCHAIN
       IF (CC.LE.NL) THEN
          C = CC
       ELSE
          C = CC - NL
       ENDIF

       ! length of chain
       IF (C.LT.MI) THEN
          CGRP%CHAINS(CC)%NPT = 2*C*LENDIAMOND+1          
       ELSE IF (C.LT.MA) THEN
          CGRP%CHAINS(CC)%NPT = 2*MI*LENDIAMOND+1
       ELSE
          CGRP%CHAINS(CC)%NPT = 2*(NDIAMOND(1)+NDIAMOND(2)-C)*LENDIAMOND+1          
       ENDIF     
    ENDDO
    
    ! set up connections
    DO C = 1,NL
       IF (C.LT.MI) THEN          
          NINTER = 2*C + 1 ! number of intersections for this chain
       ELSE IF (C.LT.MA) THEN
          NINTER = 2*MI+1
       ELSE
          NINTER = 2*(NDIAMOND(1)+NDIAMOND(2)-C)+1
       ENDIF    
       
       IF (C.LE.NDIAMOND(1)) THEN
          C2 = NDIAMOND(1)+2*NDIAMOND(2)+C          
       ELSE
          C2 = 2*NDIAMOND(2)+3*NDIAMOND(1)-C
       ENDIF

       ! keep track of bead on connected chain
       IF (C.LE.NDIAMOND(1)) THEN
          CURB2 = 1
       ELSE
          curB2 = CGRP%CHAINS(C2-1)%NPT
       ENDIF
                 
       DO D = 1,NINTER
          C2 = C2-1

          B = LENDIAMOND*(D-1)+1     
          B2 = CURB2
          IF (C.LE.NDIAMOND(1).AND.D.LT.C+1) THEN
             CURB2 = CURB2 + LENDIAMOND
          ELSE IF (C.LE.NDIAMOND(1)) then
             CURB2 = CURB2 - LENDIAMOND
          ELSE IF (C2.GT.NL+NDIAMOND(2)) THEN
             CURB2 = CURB2 + LENDIAMOND
          ELSE             
             CURB2 = CURB2 - LENDIAMOND
          ENDIF

          IF ((D.EQ.1.AND.C.EQ.NDIAMOND(1))&
               &.OR.(D.EQ.NINTER.AND.C.EQ.NDIAMOND(2))) CYCLE
          CON = CON+1
          
          IF (CON.GT.CGRP%NCONNECT) THEN
             PRINT*, 'ERROR IN SETUPDIAMONDLATTICE: too many connections'
             print*, CON, CGRP%NCONNECT, NDIAMOND
             STOP 1
          ENDIF

          CGRP%CONNECT(CON,1) = B; CGRP%CONNECT(CON,2) = C
          CGRP%CONNECT(CON,3) = B2; CGRP%CONNECT(CON,4) = C2             
       ENDDO
       
    ENDDO

    ! fix beads on the boundary
    IF (ANY(FIXBOUNDARY(1:2)).AND.ANY(FIXBOUNDARY(3:4))) THEN
       DO CC = 1,CGRP%NCHAIN
          ! start of chain is on top/bottom
          TOP1 = CC.LE.NDIAMOND(1).OR.CC.GE.NL+NDIAMOND(2)
          ! end of chain is on top/bottom
          TOPEND = CC.GE.NDIAMOND(2).AND.CC.LE.NL+NDIAMOND(1)
          ! start of chain is on sides
          SIDE1 = CC.GE.NDIAMOND(1).AND.CC.LE.NL+NDIAMOND(2)
          ! end of chain is on sides
          SIDEEND = CC.GE.NDIAMOND(2).AND.CC.LE.NL+NDIAMOND(1)

          IF ((FIXBOUNDARY(1).AND.TOP1).OR.(FIXBOUNDARY(2).AND.SIDE1)) THEN
             CGRP%FIXBEAD(1,CC,1) = FIXBOUNDARY(3)
             CGRP%FIXBEAD(1,CC,2) = FIXBOUNDARY(4)
          ENDIF
          IF ((FIXBOUNDARY(2).AND.SIDEEND).OR.(FIXBOUNDARY(1).AND.TOPEND))  THEN
             CGRP%FIXBEAD(CGRP%CHAINS(CC)%NPT,CC,1) = FIXBOUNDARY(3)
             CGRP%FIXBEAD(CGRP%CHAINS(CC)%NPT,CC,2) = FIXBOUNDARY(4)
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE SETUPDIAMONDLATTICE

  SUBROUTINE SETCHAINGROUPPARAMS(CGRP)
    ! set parameters for a group of chains based on KEYS
    ! including parameters for the chains themselves
    USE KEYS, ONLY : NCONNECT, CONNECTIONS, SQUARELATTICE,FORCEBEAD,FORCE,&
         & CONNECTPOS,CONNECTUVEC,CONPOSMOD,CONUVECMOD,FIXCONNECT,&
         & FIXBEAD,NFIXBEAD, FIXBOUNDARY, WIDTHDIAMOND,LENDIAMOND,NDIAMOND,&
         & DIAMONDLATTICE
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER :: C
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: B1,C1,B2,C2

    CGRP%CONNECTPOS = CONNECTPOS
    CGRP%CONNECTUVEC = CONNECTUVEC
    CGRP%CONPOSMOD = CONPOSMOD
    CGRP%CONUVECMOD = CONUVECMOD    
    CGRP%FIXCONNECT = FIXCONNECT
    CGRP%FIXBEAD = .FALSE.

    ! initialize chain parameters
    DO C = 1,CGRP%NCHAIN       
       CHAINP=>CGRP%CHAINS(C)
       IF (CGRP%NCHAIN.EQ.1.AND.CGRP%NFORCE.GE.1) THEN
          CALL SETCHAINPARAMS(CHAINP,SQRT(FORCE(1,1)**2 + FORCE(1,2)**2 + FORCE(1,3)**2))
       ELSE
          CALL SETCHAINPARAMS(CHAINP)
       ENDIF
    ENDDO


    ! explicitly specified connections
    DO C = 1,NCONNECT
       CGRP%CONNECT(C,:) = CONNECTIONS(C,:)
    ENDDO
    
    IF (DIAMONDLATTICE) THEN
       CALL SETUPDIAMONDLATTICE(CGRP,NCONNECT+1,NDIAMOND,LENDIAMOND,WIDTHDIAMOND,FIXBOUNDARY)
    ELSEIF (SQUARELATTICE) THEN
       ! set up square lattice of connections
       CALL CONNECTSQUARELATTICE(CGRP,NCONNECT+1,FIXBOUNDARY)
    ENDIF

    ! list of one connection for each bead
    CGRP%CONBEAD = -1
    DO C = 1,CGRP%NCONNECT
       C1 = CGRP%CONNECT(C,2); B1 = CGRP%CONNECT(C,1)
       C2 = CGRP%CONNECT(C,4); B2 = CGRP%CONNECT(C,3)
       CGRP%CONBEAD(B1,C1,:) = (/B2,C2/)
       CGRP%CONBEAD(B2,C2,:) = (/B1,C1/)
    ENDDO

    ! forces
    DO C = 1,CGRP%NFORCE
       CGRP%FORCEBEAD(C,:) = FORCEBEAD(C,:)
       CGRP%FORCE(C,:) = FORCE(C,:)
    ENDDO    
    
    ! fixed beads    
    DO C = 1,NFIXBEAD
       B1 = FIXBEAD(C,1)
       C1 = FIXBEAD(C,2)
       CGRP%FIXBEAD(B1,C1,1) = (FIXBEAD(C,3).GT.0) ! position
       CGRP%FIXBEAD(B1,C1,2) = (FIXBEAD(C,4).GT.0) ! orientation
    ENDDO
        
    ! PRINT*, 'TESTX FIXBEADS:', CGRP%FIXBEAD(1,1,:)
    ! DO C1 = 1,CGRP%NCHAIN
    !    CHAINP=>CGRP%CHAINS(C1)
    !    DO B1 = 1,CHAINP%NPT
    !       PRINT*, B1, C1, CGRP%FIXBEAD(B1,C1,:)
    !    ENDDO
    ! ENDDO
  END SUBROUTINE SETCHAINGROUPPARAMS

  SUBROUTINE CONNECTSQUARELATTICE(CGRP,STARTCON,FIXBOUNDARY)
    ! connect up a group of chains in a square lattice
    ! must be a perfect square number of chains, 
    ! with lattice connections falling on beads
    ! start saving connections list at given index
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: STARTCON
    LOGICAL, INTENT(IN) :: FIXBOUNDARY(4)
    INTEGER :: NC, NPT, C,C2,CON,DIST
    TYPE(CHAIN), POINTER :: CHAINP
    
    IF (MOD(CGRP%NCHAIN,2).NE.0) THEN
       PRINT*, 'ERROR IN SETCHAINGROUPPARAMS: cannot have square lattice with odd number of chains', CGRP%NCHAIN 
       STOP 1
    ENDIF

    NC = CGRP%NCHAIN/2

    CHAINP=>CGRP%CHAINS(1)
    NPT = CHAINP%NPT
    DO C = 1,CGRP%NCHAIN
       IF (CGRP%CHAINS(C)%NPT.NE.NPT) THEN
          PRINT*, 'ERROR IN SETCHAINGROUPPARAMS: cannot have square lattice with chains of different length'
          STOP 1
       ENDIF
    ENDDO

    IF (MOD(NPT-1,NC-1).NE.0) THEN
       PRINT*, 'ERROR IN SETCHAINGROUPPARAMS: lattice size is not integer number of beads',NPT,NC
       STOP 1
    ENDIF

    DIST = (NPT-1)/(NC-1)

    CON = STARTCON-1
    DO C = 1,NC
       DO C2 = 1,NC
          CON = CON+1
          
          IF (CON.GT.CGRP%NCONNECT) THEN
             PRINT*, 'ERROR IN CONNECTSQUARELATTICE:, too many connections!', CON
             STOP 1
          ENDIF

          ! first bead and chain
          CGRP%CONNECT(CON,2) = C
          CGRP%CONNECT(CON,1) = DIST*(C2-1)+1
          ! second bead and chain
          CGRP%CONNECT(CON,4) = C2+NC
          CGRP%CONNECT(CON,3) = DIST*(C-1)+1

          IF (ANY(FIXBOUNDARY)) THEN
             IF (.NOT.ALL(FIXBOUNDARY)) THEN
                PRINT*, 'ERROR IN CONNECTSQUARELATTICE: not yet set up with partial fixed boundary'
                STOP 1
             ENDIF

             ! fix positions for boundary intersections
             IF (C.EQ.1.OR.C.EQ.NC.OR.C2.EQ.1.OR.C2.EQ.NC) THEN               
                CGRP%FIXBEAD(CGRP%CONNECT(CON,1),C,1) = .TRUE.
                CGRP%FIXBEAD(CGRP%CONNECT(CON,3),C2+NC,1) = .TRUE.                
             ENDIF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE CONNECTSQUARELATTICE

  SUBROUTINE SETUPCHAINGROUP(CGRP,NCHAIN,NCON,NFORCE,MAXNPT)
    ! set up a group of chains, allocating all arrays
    ! do not initialize positions, parameters, or connections
    ! NCON is maximal allowed number of connections
    ! MAXNPT is the maximum number of beads in each chain
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    INTEGER, INTENT(IN) :: NCHAIN, NCON,NFORCE,MAXNPT
    INTEGER :: C
    TYPE(CHAIN), POINTER :: CHAINP

    CGRP%NCHAIN = NCHAIN
    CGRP%NCONNECT = NCON
    CGRP%NFORCE = NFORCE

    ALLOCATE(CGRP%CHAINS(NCHAIN),CGRP%CONNECT(NCON,4),&
         & CGRP%CONBEAD(MAXNPT,NCHAIN,2), &
         & CGRP%FORCEBEAD(NFORCE,2),CGRP%FORCE(NFORCE,3),&
         & CGRP%FIXBEAD(MAXNPT,NCHAIN,2))
    
    ! allocate arrays for the chains
    DO C = 1,CGRP%NCHAIN
       CHAINP=>CGRP%CHAINS(C)    
       CALL SETUPCHAIN(CHAINP,MAXNPT)
    ENDDO
    
  END SUBROUTINE SETUPCHAINGROUP

  SUBROUTINE CLEANUPCHAINGROUP(CGRP)
    ! deallocate everything for a chain group, including the arrays within the chains

    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: C

    DO C = 1,CGRP%NCHAIN
       IF (CGRP%CHAINS(C)%ARRAYSET) THEN
          CHAINP=>CGRP%CHAINS(C)
          CALL CLEANUPCHAIN(CHAINP)
       ENDIF
    ENDDO

    DEALLOCATE(CGRP%CHAINS,CGRP%CONNECT,CGRP%FORCEBEAD, CGRP%FORCE,CGRP%CONBEAD,CGRP%FIXBEAD)
  END SUBROUTINE CLEANUPCHAINGROUP

  SUBROUTINE COPYCHAINGROUP(CGRP1,CGRP2)
    ! copy everything from one chain group to another
    IMPLICIT NONE
    TYPE(CHAINGROUP), POINTER :: CGRP1, CGRP2
    INTEGER :: C
    TYPE(CHAIN), POINTER :: CHAINp1, CHAINp2

    IF (CGRP2%NCHAIN.NE.CGRP1%NCHAIN) THEN
       PRINT*,'ERROR IN COPYCHAINGROUP: cannot copy if different number of chains', CGRP1%NCHAIN, CGRP2%NCHAIN
       STOP 1
    ENDIF
    IF (CGRP2%NCONNECT.NE.CGRP1%NCONNECT) THEN
       PRINT*, 'ERROR IN COPYCHAINGROUP: different number of connections', CGRP1%NCONNECT, CGRP2%NCONNECT
       STOP 1
    ENDIF
     IF (CGRP2%NFORCE.NE.CGRP1%NFORCE) THEN
       PRINT*, 'ERROR IN COPYCHAINGROUP: different number of forces', CGRP1%NFORCE, CGRP2%NFORCE
       STOP 1
    ENDIF

    DO C = 1,CGRP1%NCHAIN
       CHAINP1=>CGRP1%CHAINS(C)
       CHAINP2=>CGRP2%CHAINS(C)

       IF (CHAINP1%MAXNPT.NE.CHAINP2%MAXNPT) THEN
          PRINT*, 'ERROR IN COPYCHAINGROUP: cannot copy if chains have different max size', C, CHAINP1%MAXNPT, CHAINP2%MAXNPT
          STOP 1
       ENDIF
       
       CALL COPYCHAIN(CHAINP1,CHAINP2)
    ENDDO
    
    CGRP2%CONNECT = CGRP1%CONNECT
    CGRP2%CONPOSMOD = CGRP1%CONPOSMOD
    CGRP2%CONUVECMOD = CGRP1%CONUVECMOD
    CGRP2%CONNECTPOS = CGRP1%CONNECTPOS
    CGRP2%CONNECTUVEC = CGRP1%CONNECTUVEC
    CGRP2%FIXCONNECT = CGRP1%FIXCONNECT
    CGRP2%CONBEAD = CGRP1%CONBEAD
    CGRP2%FIXBEAD = CGRP1%FIXBEAD
    CGRP2%FORCEBEAD = CGRP1%FORCEBEAD
    CGRP2%FORCE = CGRP1%FORCE

  END SUBROUTINE COPYCHAINGROUP
END MODULE MANYCHAINS
