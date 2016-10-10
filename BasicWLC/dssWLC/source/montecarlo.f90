MODULE MONTECARLO
  USE CHAINUTIL, ONLY : CHAIN, GETBEADENERGY,GETFORCEENERGY, GETSTERICENERGY
  USE MT19937, ONLY : GRND
  USE QUATUTIL, ONLY : ROTANGAX, PI
  USE BROWNDYN, ONLY : GETCHAINFORCEINT
  IMPLICIT NONE
  
CONTAINS
  SUBROUTINE RUNMONTECARLO(CGRP,TOTSTEPS,STATSTEPS,INITSTEPS,STARTSTEP)
    ! run a monte carlo simulation for a group of chains
    ! local bead moves only
    ! all step counts are actually number of sweeps (ie: individual steps multiplied by total number of beads)
    ! TOTSTEPS is total number of steps to do
    ! first INITSTEPS steps do not contribute to statistics
    ! Take average statistics every STATSTEPS steps thereafter
    ! start with step given by STARTSTEP+1
    ! if STARTSTEP>0, append to output files rather than rewriting

    ! WARNING: currently not using BEADENERGY arrays; should implement this later for increased efficiency
    USE KEYS, ONLY : TRACKDIST,MCPRINTFREQ,OUTFILE,MCOUTPUTFREQ,ADJUSTEVERY,&
         & FACCTARGET,FACCTOL,ADJUSTSCL,INITRANGE,SNAPSHOTFILE,DUMPSNAPSHOTS,&
         & SNAPSHOTEVERY
    USE MANYCHAINS, ONLY : CHAINGROUP, GROUPSNAPSHOT, GROUPENERGY,BEADALLENERGY
    USE CHAINUTIL, ONLY : CHAIN
    IMPLICIT NONE

    TYPE(CHAINGROUP), POINTER :: CGRP
    TYPE(CHAIN), POINTER :: CHAINP, CHAINP2
    INTEGER, INTENT(IN) :: TOTSTEPS,INITSTEPS,STATSTEPS,STARTSTEP
    INTEGER :: STEP, NCT,B,C,NPT,SUBSTEP
    DOUBLE PRECISION :: PREVPOS(3), PREVUVEC(3), PREVENERGY, AVGENERGY
    DOUBLE PRECISION :: ARANGE, RRANGE, R2AVG, R2CUR, DELE, TMP, DR(3), FACC
    INTEGER :: SNAPCT, C2, B2
    LOGICAL :: ACCEPT, ADJUSTED
    DOUBLE PRECISION :: ENERGY1, ENERGY2, TESTENERGY1, TESTENERGY2
    DOUBLE PRECISION :: PREVPOS2(3), PREVUVEC2(3)
    LOGICAL :: MVOTHERBEAD
    DOUBLE PRECISION :: SAVEBEADENERGY(2), SAVEBEADENERGY2(2)
    LOGICAL :: DOROT1, DOROT2, DOSHIFT1, DOSHIFT2
    INTEGER, ALLOCATABLE :: ALLBEADS(:,:)
    INTEGER :: TOTBEADS, BEADCT
    DOUBLE PRECISION :: FCT, TOTACCEPT,RU

    ! Full list of all beads to use in sampling
    TOTBEADS = 0
    DO C = 1,CGRP%NCHAIN
       TOTBEADS = TOTBEADS + CGRP%CHAINS(C)%NPT
    ENDDO

    PRINT*, 'TOTAL NUMBER OF BEADS:', TOTBEADS

    ALLOCATE(ALLBEADS(TOTBEADS,2))

    BEADCT = 0
    DO C = 1,CGRP%NCHAIN
       DO B = 1,CGRP%CHAINS(C)%NPT
          BEADCT = BEADCT + 1
          ALLBEADS(BEADCT,:) = (/B,C/)
       ENDDO
    ENDDO

    ! NPT = CGRP%CHAINS(1)%NPT
    ! DO C = 2,CGRP%NCHAIN
    !    IF (CGRP%CHAINS(C)%NPT.NE.NPT) THEN
    !       PRINT*, 'ERROR IN RUNMONTECARLO: group MC not set up yet for chains of different lengths'
    !       STOP 1
    !    ENDIF
    ! ENDDO


    IF (STARTSTEP.GT.0) THEN
       PRINT*, 'STARTING MC FROM STEP:', STARTSTEP
    ENDIF

    IF (STARTSTEP.GT.0) THEN
       OPEN(UNIT=55,FILE=OUTFILE,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
       OPEN(UNIT=55,FILE=OUTFILE,STATUS='UNKNOWN')
    ENDIF

    R2AVG = 0D0; AVGENERGY = 0D0
    NCT = 0
    TOTACCEPT = 0;
    SNAPCT = 0
    FCT = 0 ! count of total steps for fraction accepted

    ARANGE = INITRANGE(1)
    RRANGE = INITRANGE(2)

    ! get initial energy
    CALL GROUPENERGY(CGRP,PREVENERGY)
    PRINT*, 'STARTENERGY:', PREVENERGY

    ! Dump starting snapshot
    IF (DUMPSNAPSHOTS) THEN
       SNAPCT = SNAPCT + 1
       CALL GROUPSNAPSHOT(CGRP,SNAPSHOTFILE,STARTSTEP,STARTSTEP.GT.0)
    ENDIF

    DO STEP = STARTSTEP+1,TOTSTEPS
       DO SUBSTEP = 1,TOTBEADS
          FCT = FCT + 1

          ! take a Monte Carlo step


          ! pick a bead at random
          BEADCT = FLOOR(GRND()*TOTBEADS)+1
          ! find the corresponding chain and bead
          C = ALLBEADS(BEADCT,2)
          B = ALLBEADS(BEADCT,1)
          CHAINP=>CGRP%CHAINS(C)

          PREVPOS = CHAINP%POS(:,B)
          PREVUVEC = CHAINP%UVEC(:,B)
          IF (B.LT.CHAINP%NPT) THEN
             SAVEBEADENERGY = (/CHAINP%BEADENERGY(B),CHAINP%BEADENERGY(B+1)/)
          ELSE
             SAVEBEADENERGY = (/CHAINP%BEADENERGY(B),0D0/)
          ENDIF

          ! move another bead along with this one
          MVOTHERBEAD =CGRP%FIXCONNECT.AND.CGRP%CONBEAD(B,C,1).GT.0
          IF (MVOTHERBEAD) THEN
             C2 = CGRP%CONBEAD(B,C,2); B2 = CGRP%CONBEAD(B,C,1)
             CHAINP2=>CGRP%CHAINS(C2)
             IF (CGRP%CONNECTPOS) PREVPOS2 = CHAINP2%POS(:,B2)
             IF (CGRP%CONNECTUVEC) PREVUVEC2 = CHAINP2%UVEC(:,B2)

             IF (B2.LT.CHAINP2%NPT) THEN
                SAVEBEADENERGY2 = (/CHAINP2%BEADENERGY(B2),CHAINP2%BEADENERGY(B2+1)/)
             ELSE
                SAVEBEADENERGY2 = (/CHAINP2%BEADENERGY(B2),0D0/)
             ENDIF
          ENDIF

          DOROT1 = .NOT.CGRP%FIXBEAD(B,C,2)
          DOSHIFT1 = .NOT.CGRP%FIXBEAD(B,C,1)

          !CALL GROUPENERGY(CGRP,ENERGY1)       
          CALL BEADALLENERGY(CGRP,C,B,ENERGY1,.FALSE.)       
          IF (MVOTHERBEAD) THEN
             ! simultaneously move 2 joined beads
             CALL BEADALLENERGY(CGRP,C2,B2,TMP,.FALSE.)
             ENERGY1 = ENERGY1+TMP      
             DOROT1 = .NOT.(CGRP%FIXBEAD(B,C,2).OR.CGRP%FIXBEAD(B2,C2,2))
             DOSHIFT1 = .NOT.(CGRP%FIXBEAD(B,C,1).OR.CGRP%FIXBEAD(B2,C2,1))
             DOROT2 = DOROT1.AND.CGRP%CONNECTUVEC
             DOSHIFT2 = DOSHIFT1.AND.CGRP%CONNECTPOS
             CALL LOCALMOVE2BEAD(CHAINP,CHAINP2,B,B2,ARANGE,RRANGE,&
                  & DOROT1,DOSHIFT1,DOROT2,DOSHIFT2)
             CALL BEADALLENERGY(CGRP,C2,B2,TMP,.TRUE.)
             ENERGY2 = TMP
          ELSE          
             DOROT1 = .NOT.CGRP%FIXBEAD(B,C,2)
             DOSHIFT1 = .NOT.CGRP%FIXBEAD(B,C,1)
             CALL LOCALMOVE(CHAINP,ARANGE,RRANGE,B,DOROT1,DOSHIFT1)
             ENERGY2 = 0
          ENDIF
          CALL BEADALLENERGY(CGRP,C,B,TMP,.TRUE.)
          ENERGY2 = ENERGY2+TMP
          !CALL GROUPENERGY(CGRP,ENERGY2)
          !PRINT*, 'TESTX0:', ENERGY1, ENERGY2, PREVENERGY
          DELE = ENERGY2-ENERGY1

          ! CALL GROUPENERGY(CGRP,TESTENERGY2)
          ! IF (ABS(PREVENERGY+DELE - TESTENERGY2).GT.1D-12) THEN
          !    PRINT*, 'ERROR IN ENERGY:',PREVENERGY+DELE, ENERGY2, PREVENERGY, DELE, ENERGY1
          !    print*, dele+energy1, energy2, energy1, prevenergy
          !    PRINT*, C,B, CGRP%CONBEAD(B,C,:)
          !    STOP 1
          ! ENDIF


          ! decide whether to accept
          IF (DELE.LT.0) THEN
             ACCEPT = .TRUE.
          ELSE
             TMP = GRND()
             ACCEPT = (TMP.LT.EXP(-DELE))
          ENDIF

          IF (ACCEPT) THEN                    
             TOTACCEPT = TOTACCEPT + 1
             PREVENERGY = PREVENERGY + DELE               
          ELSE ! restore old coordinates 
             CHAINP%POS(:,B) = PREVPOS
             CHAINP%UVEC(:,B) = PREVUVEC
             IF (B.GT.1) THEN
                CHAINP%BEADENERGY(B) = SAVEBEADENERGY(1)
             ENDIF
             IF (B.LT.CHAINP%NPT) THEN
                CHAINP%BEADENERGY(B+1) = SAVEBEADENERGY(2)
             ENDIF
             IF (CGRP%FIXCONNECT.AND.CGRP%CONBEAD(B,C,1).GT.0) THEN   
                IF (CGRP%CONNECTPOS) CHAINP2%POS(:,B2) = PREVPOS2
                IF (CGRP%CONNECTUVEC) CHAINP2%UVEC(:,B2) = PREVUVEC2
                IF (B2.GT.1) THEN
                   CHAINP2%BEADENERGY(B2) = SAVEBEADENERGY2(1)
                ENDIF
                IF (B2.LT.CHAINP2%NPT) THEN
                   CHAINP2%BEADENERGY(B2+1) = SAVEBEADENERGY2(2)
                ENDIF
             ENDIF
          ENDIF


          IF (SUBSTEP.GT.1) CYCLE
          ! for the first substep, check whether to do various statistics and output

          IF (STEP.GT.INITSTEPS.AND.MOD(STEP,STATSTEPS).EQ.0) THEN
             ! update statistics
             NCT = NCT + 1      
             DR = CGRP%CHAINS(TRACKDIST(2))%POS(:,TRACKDIST(1)) &
                  - CGRP%CHAINS(TRACKDIST(4))%POS(:,TRACKDIST(3))       
             R2CUR = DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3)
             R2AVG = R2AVG + R2CUR
             AVGENERGY = AVGENERGY + PREVENERGY
          ENDIF

          FACC =  DBLE(TOTACCEPT)/FCT; 

          IF (MOD(STEP,MCPRINTFREQ).EQ.0) THEN
             PRINT '(A,G20.5,4G15.7,2I5,1x,L1,1x,2I4)', 'STEP, FACC, R2AVG, R2CUR:', DBLE(STEP), &
                  & FACC, R2AVG/NCT,  AVGENERGY/NCT, PREVENERGY, C, B,accept,CGRP%CONBEAD(B,C,:) 
          ENDIF
          IF (STEP.GT.INITSTEPS.AND.MOD(STEP,MCOUTPUTFREQ).EQ.0) THEN
             DR = CGRP%CHAINS(TRACKDIST(4))%POS(:,TRACKDIST(3)) &
                  - CGRP%CHAINS(TRACKDIST(2))%POS(:,TRACKDIST(1))    
             NPT = CGRP%CHAINS(1)%NPT
             RU = DOT_PRODUCT(CGRP%CHAINS(1)%POS(:,NPT)-CGRP%CHAINS(1)%POS(:,1)&
                  & ,CGRP%CHAINS(1)%UVEC(:,1))
             WRITE(55,*) STEP, FACC,R2AVG/NCT, DR, PREVENERGY,RU
             FLUSH(55)
          ENDIF

          IF (DUMPSNAPSHOTS.AND.STEP.GT.INITSTEPS.AND.MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
             SNAPCT = SNAPCT + 1
             CALL GROUPSNAPSHOT(CGRP,SNAPSHOTFILE,STEP,SNAPCT.GT.1.OR.STARTSTEP.GT.0)
          ENDIF

          IF (ADJUSTEVERY.GT.0.AND.MOD(STEP,ADJUSTEVERY).EQ.0) THEN
             ! adjust step range if needed


             ADJUSTED =.FALSE.
             IF (FACC.GT.FACCTARGET+FACCTOL) THEN
                IF (.NOT.(ARANGE.GE.2*PI.AND..NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE)) THEN
                   ARANGE = ARANGE*ADJUSTSCL; RRANGE = RRANGE*ADJUSTSCL
                   ADJUSTED = .TRUE.             
                ENDIF
             ELSE IF (FACC.LT.FACCTARGET-FACCTOL) THEN
                ARANGE = ARANGE/ADJUSTSCL; RRANGE = RRANGE/ADJUSTSCL
                ADJUSTED = .TRUE.        
             ENDIF
             ARANGE = MIN(ARANGE,2*PI)

             IF (ADJUSTED) PRINT*, 'ADJUSTED RANGES:', ARANGE, RRANGE
          ENDIF

       ENDDO
    ENDDO

    !    verbose = .true.
    CALL GROUPENERGY(CGRP,PREVENERGY)
    PRINT*, 'FINAL ENERGY: ', PREVENERGY

    CLOSE(55)

    DEALLOCATE(ALLBEADS)
  END SUBROUTINE RUNMONTECARLO

  SUBROUTINE RUNMONTECARLO1CHAIN(CHAINP,TOTSTEPS,STATSTEPS,INITSTEPS,STARTSTEP)
    ! run a monte carlo simulation for a single chain
    ! alternating the two types of crank moves
    ! TOTSTEPS is total number of steps to do
    ! first INITSTEPS steps do not contribute to statistics
    ! Take average statistics every STATSTEPS steps thereafter
    ! start with step given by STARTSTEP+1
    ! if STARTSTEP>0, append to output files rather than rewriting
    ! WARNING: this will not currently work with external force

    USE KEYS, ONLY : MCPRINTFREQ,VERBOSE,OUTFILE,MCOUTPUTFREQ, ADJUSTEVERY, &
         & FACCTARGET, FACCTOL, ADJUSTSCL, INITRANGE, SNAPSHOTFILE, &
         & DUMPSNAPSHOTS,SNAPSHOTEVERY, DOREDISC, APPENDSNAPSHOTS, &
         & USEBDENERGY, dolocalmoves, OUTPUTBEADWEIGHT, INTUWEIGHTNPT, INTRWEIGHTNPT
    USE CHAINUTIL, ONLY : SETUPCHAIN,COPYCHAIN,GETENERGY,CLEANUPCHAIN,&
         & GETCHAINRG,OUTPUTSNAPSHOT, GETBEADQUINT,GETBEADQRINT
    !USE REDISC, ONLY : REDISCREMOVE, REDISCADD
    IMPLICIT NONE

    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: TOTSTEPS,INITSTEPS,STATSTEPS,STARTSTEP
    TYPE(CHAIN), TARGET :: WLC
    TYPE(CHAIN), POINTER :: PREVCHAINP
    INTEGER :: STEP, NCT, B
    DOUBLE PRECISION :: ARANGE1, RRANGE1, ARANGE2, RRANGE2
    DOUBLE PRECISION :: DELE, PREVENERGY, TMP
    LOGICAL :: ACCEPT
    DOUBLE PRECISION :: DR(3), R2AVG, TESTENERGY,R2CUR, FACC1, FACC2,TANCORAVG
    INTEGER :: BPIVOT(2),  TOTACCEPT1, TOTACCEPT2, I
    INTEGER :: MVTYPE, SNAPCT, NFLEXEDGE, NMV1, NMV2
    LOGICAL :: ADJUSTED, ENERGYFROMBD
    DOUBLE PRECISION :: UU, RGAVG, RG, DELEADD, DELERM,ENERGY,RU
    DOUBLE PRECISION :: RFORCE(3,CHAINP%NPT),UFORCE(3,CHAINP%NPT), ENERGY1, ENERGY2
    DOUBLE PRECISION :: D1(3),D2(3),RHO1,RHO2,D3(3),D4(3),RHO3,RHO4
    INTEGER :: NFLEXBEAD, CHOOSEB, CT
    DOUBLE PRECISION :: ENERGYOLD, QINT
    DOUBLE PRECISION, TARGET :: BEADWEIGHT(2,CHAINP%NPT)
    DOUBLE PRECISION, POINTER :: BEADWEIGHTPT(:,:)

    BEADWEIGHTPT=>BEADWEIGHT

    ! use energy from BD calculations?
    ENERGYFROMBD = USEBDENERGY.OR.(CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE)

    IF (ENERGYFROMBD) THEN
       PRINT*, 'Using energy calculations from brownian dynamics code'
    ENDIF

    RHO1 = 0; RHO2 = 0; RHO3 = 0; RHO4 = 0;
    FACC1 = FACCTARGET; FACC2 = FACCTARGET
    NMV1 = 0; NMV2 = 0

    print*, 'Starting single-chain monte carlo, with force:', CHAINP%FORCE, CHAINP%HASFORCE

    IF (STARTSTEP.GT.0) THEN
       PRINT*, 'STARTING MC FROM STEP:', STARTSTEP
    ENDIF

    IF (STARTSTEP.GT.0) THEN
       OPEN(UNIT=55,FILE=OUTFILE,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
       OPEN(UNIT=55,FILE=OUTFILE,STATUS='UNKNOWN')
    ENDIF

    R2AVG = 0D0; RGAVG = 0D0
    NCT = 0
    TOTACCEPT1 = 0; TOTACCEPT2 = 0
    TANCORAVG = 0D0
    SNAPCT = 0

    ARANGE1 = INITRANGE(1)
    RRANGE1 = INITRANGE(2)
    ARANGE2 = INITRANGE(3)
    RRANGE2 = INITRANGE(4)

    ! get initial energy
    IF (ENERGYFROMBD) THEN
       CALL GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,prevenergy,.FALSE.)
    ELSE
       CALL GETENERGY(CHAINP,PREVENERGY)
    ENDIF
    PRINT*, 'STARTENERGY:', PREVENERGY, chainp%npt
    
    ! print*, 'testx1:', chainp%beadenergy
    ! do i = 1,chainp%npt
    !    print*, chainp%pos(:,i),chainp%uvec(:,i)
    ! enddo
    ! stop 1

    ! set up chain copy
    PREVCHAINP=>WLC

    CALL SETUPCHAIN(PREVCHAINP,CHAINP%MAXNPT)    
    CALL COPYCHAIN(CHAINP,PREVCHAINP)

    STEP = 0
    ! Dump starting snapshot
    IF (DUMPSNAPSHOTS) THEN
       SNAPCT = SNAPCT + 1
       IF (OUTPUTBEADWEIGHT) THEN
          ! output weight of each mobile bead by integrating over U and R vectors
          BEADWEIGHT = 0D0
          DO B = 1,CHAINP%NPT
             IF (.NOT.CHAINP%ISFIXED(B)) THEN
                CALL GETBEADQUINT(CHAINP,B,INTUWEIGHTNPT,INTUWEIGHTNPT,QINT)
                BEADWEIGHT(1,B) = QINT
                CALL GETBEADQRINT(CHAINP,B,INTrWEIGHTNPT,QINT)
                BEADWEIGHT(2,B) = QINT
             ENDIF
          ENDDO

          CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,STEP,STARTSTEP.GT.0.OR.APPENDSNAPSHOTS,2,BEADWEIGHTPT)
       ELSE
          CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,STEP,STARTSTEP.GT.0.OR.APPENDSNAPSHOTS)

       ENDIF

    ENDIF

    DO STEP = STARTSTEP+1,TOTSTEPS
       ! take a Monte Carlo step

       IF (DOLOCALMOVES) THEN
          ! number of moveable beads
          NFLEXBEAD = CHAINP%NPT - CHAINP%NFIXBEAD
          ! choose which moveable bead to move
          ! WARNING: this is really inefficient...
          CHOOSEB = FLOOR(GRND()*NFLEXBEAD)+1
          CT=0
          DO B = 1,CHAINP%NPT
             IF (.NOT.CHAINP%ISFIXED(B)) THEN
                CT = CT + 1
             ENDIF
             IF (CT.EQ.CHOOSEB) THEN
                BPIVOT(1) = B
                EXIT
             ENDIF
          ENDDO
          IF (B.GT.CHAINP%NPT) THEN
             PRINT*, 'ERROR IN LOCALMOVE MONTECARLO: BAD B', B, CT, CHOOSEB
             STOP 1
          ENDIF

          ! do the local bead move
          CALL LOCALMOVE(CHAINP,ARANGE1,RRANGE1,BPIVOT(1),.TRUE.,.TRUE.)       
          ! get the change in energy
          IF (BPIVOT(1).LT.CHAINP%NPT) THEN
             ENERGYOLD = CHAINP%BEADENERGY(BPIVOT(1))+CHAINP%BEADENERGY(BPIVOT(1)+1)
             CALL GETBEADENERGY(CHAINP,BPIVOT(1),CHAINP%BEADENERGY(BPIVOT(1)))
             CALL GETBEADENERGY(CHAINP,BPIVOT(1)+1,CHAINP%BEADENERGY(BPIVOT(1)+1))
             DELE = CHAINP%BEADENERGY(BPIVOT(1))+CHAINP%BEADENERGY(BPIVOT(1)+1) - ENERGYOLD
          ELSE
             ENERGYOLD =  CHAINP%BEADENERGY(BPIVOT(1))
             CALL GETBEADENERGY(CHAINP,BPIVOT(1),CHAINP%BEADENERGY(BPIVOT(1)))
             DELE = CHAINP%BEADENERGY(BPIVOT(1)) - ENERGYOLD
          ENDIF


          NMV1 = NMV1 + 1
          MVTYPE = 1
       ELSE
          ! decide what type of crank step to take
          MVTYPE = MOD(STEP,2)+1
          IF (CHAINP%NFIXBEAD.GT.0.AND.MVTYPE.EQ.1) THEN
             ! find how many flexible edge beads there are
             NFLEXEDGE = 0
             IF (CHAINP%FIXBEAD(1).GT.1) &
                  & NFLEXEDGE = NFLEXEDGE + CHAINP%FIXBEAD(1)+1
             IF (CHAINP%FIXBEAD(CHAINP%NFIXBEAD).LT.CHAINP%NPT) &
                  & NFLEXEDGE = NFLEXEDGE + CHAINP%NPT-CHAINP%FIXBEAD(CHAINP%NFIXBEAD)+1

             ! use that to get fraction of type 1 crank moves (tail pivots)       
             TMP = GRND()
             IF (TMP.GT.DBLE(NFLEXEDGE)/CHAINP%NPT) THEN
                MVTYPE = 2
             ENDIF
          ENDIF

          !IF (ENERGYFROMBD) THEN
          !   CALL GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,ENERGY1,.FALSE.)
          !ENDIF

          IF (MVTYPE.EQ.1) THEN
             CALL CRANKMOVE1(CHAINP,ARANGE1,RRANGE1,BPIVOT(1),DELE)     
             NMV1 = NMV1 + 1
          ELSEIF (MVTYPE.EQ.2) THEN
             CALL CRANKMOVE2(CHAINP,ARANGE2,RRANGE2,BPIVOT,DELE)
             NMV2 = NMV2 + 1
          ELSEIF (MVTYPE.EQ.3) THEN
             PRINT*, 'MOVE 3 IS NOT SET UP FOR 1CHAIN MONTECARLO'
             STOP 1
             !CALL LOCALMOVE(CHAINP,ARANGE1,RRANGE1,BPIVOT(1),DELE)
          ELSE
             PRINT*, 'ERROR IN MONTECARLO: UNKNOWN MOVETYPE'
             STOP 1
          ENDIF
       ENDIF

       IF (ENERGYFROMBD) THEN
          CALL GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,ENERGY2,.FALSE.)
          !          DELE = ENERGY2-ENERGY1
          DELE = ENERGY2-PREVENERGY
       ENDIF

       ! Rediscretize chain
       ! IF (DOREDISC.AND.DELE.LT.HUGE(1D0)) THEN      
       !    CALL REDISCREMOVE(CHAINP,DELERM)
       !    CALL REDISCADD(CHAINP,DELEADD)
       !    CALL GETENERGY(CHAINP,ENERGY)
       !    DELE = ENERGY-PREVENERGY
       ! ENDIF

       ! decide whether to accept
       IF (DELE.LT.0) THEN
          ACCEPT = .TRUE.
       ELSE
          TMP = GRND()
          ACCEPT = (TMP.LT.EXP(-DELE))
       ENDIF

       IF (ACCEPT) THEN          
          !D1 = CHAINP%POS(:,2)-CHAINP%POS(:,1)
          !D2 = CHAINP%POS(:,CHAINP%NPT) - CHAINP%POS(:,CHAINP%NPT-1)
          !UU = DOT_PRODUCT(D1,D2)/SQRT(DOT_PRODUCT(D1,D1)*DOT_PRODUCT(D2,D2))
          !print*, 'testx1:', uu, energy2
          IF (MVTYPE.EQ.2) THEN
             !PREVCHAINP%BEADENERGY(BPIVOT(2)) = CHAINP%BEADENERGY(BPIVOT(2))
             TOTACCEPT2 = TOTACCEPT2+1
          ELSEIF (MVTYPE.EQ.1) THEN
             TOTACCEPT1 = TOTACCEPT1 + 1
          ELSE IF (MVTYPE.EQ.3) THEN
             TOTACCEPT1 = TOTACCEPT1 + 1
             !IF (BPIVOT(1).LT.chainp%NPT) THEN
             !   PREVCHAINP%BEADENERGY(BPIVOT(1)+1) = CHAINP%BEADENERGY(BPIVOT(1)+1)
             !ENDIF
          ENDIF

          PREVENERGY = PREVENERGY + DELE                
          PREVCHAINP%POS = CHAINP%POS
          PREVCHAINP%UVEC = CHAINP%UVEC
          !PREVCHAINP%BEADENERGY(BPIVOT(1)) = CHAINP%BEADENERGY(BPIVOT(1))      
          PREVCHAINP%BEADENERGY = CHAINP%BEADENERGY
          PREVCHAINP%FORCEENERGY = CHAINP%FORCEENERGY

       ELSE ! restore old coordinates and energies
          CHAINP%POS = PREVCHAINP%POS
          CHAINP%UVEC = PREVCHAINP%UVEC

          !CHAINP%BEADENERGY(BPIVOT(1)) = PREVCHAINP%BEADENERGY(BPIVOT(1))
          CHAINP%FORCEENERGY = PREVCHAINP%FORCEENERGY
          !IF (MVTYPE.EQ.2) THEN
          !   CHAINP%BEADENERGY(BPIVOT(2)) = PREVCHAINP%BEADENERGY(BPIVOT(2))
          !ELSEIF (MVTYPE.EQ.3.AND.BPIVOT(1).LT.chainp%NPT) THEN
          !   CHAINP%BEADENERGY(BPIVOT(1)+1) = PREVCHAINP%BEADENERGY(BPIVOT(1)+1)
          !          ENDIF
          CHAINP%BEADENERGY = PREVCHAINP%BEADENERGY
       ENDIF

       DR = CHAINP%POS(:,chainp%NPT) - CHAINP%POS(:,1)
       R2CUR = DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3)

       IF (STEP.GT.INITSTEPS.AND.MOD(STEP,STATSTEPS).EQ.0) THEN
          ! update statistics
          NCT = NCT + 1          
          R2AVG = R2AVG + R2CUR
          CALL GETCHAINRG(CHAINP,RG)
          RGAVG = RGAVG + RG
          TANCORAVG = TANCORAVG + DOT_PRODUCT(CHAINP%UVEC(:,chainp%NPT),CHAINP%UVEC(:,1))
       ENDIF

       IF (NMV1.GT.0) THEN
          FACC1 =  DBLE(TOTACCEPT1)/NMV1; 
       ENDIF
       IF (NMV2.GT.0) THEN
          FACC2 = DBLE(TOTACCEPT2)/NMV2;
       ENDIF

       IF (MOD(STEP,MCPRINTFREQ).EQ.0) THEN
          PRINT '(A,I15,4G20.10,1x,L1,1x,2I5)', 'STEP, FACC, R2AVG, R2CUR:', STEP, &
               & FACC1, FACC2, R2AVG/NCT,  PREVENERGY, accept, BPIVOT
       ENDIF
       IF (STEP.GT.INITSTEPS.AND.MOD(STEP,MCOUTPUTFREQ).EQ.0) THEN
          IF (CHAINP%SHEARABLE) THEN
             UU = DOT_PRODUCT(CHAINP%UVEC(:,chainp%NPT-1),CHAINP%UVEC(:,1))
          ELSE
             D1 = CHAINP%POS(:,2)-CHAINP%POS(:,1)
             D2 = CHAINP%POS(:,CHAINP%NPT) - CHAINP%POS(:,CHAINP%NPT-1)
             UU = DOT_PRODUCT(D1,D2)/SQRT(DOT_PRODUCT(D1,D1)*DOT_PRODUCT(D2,D2))
             RHO1 = D1(3)/SQRT(DOT_PRODUCT(D1,D1))
             RHO2= D2(3)/SQRT(DOT_PRODUCT(D2,D2))
             IF (CHAINP%NPT.GE.5) THEN
                D3 = CHAINP%POS(:,3)-CHAINP%POS(:,2)
                D4 = CHAINP%POS(:,4)-CHAINP%POS(:,3)
                RHO3 = D3(3)/SQRT(DOT_PRODUCT(D3,D3))
                RHO4 = D4(3)/SQRT(DOT_PRODUCT(D4,D4))
             ENDIF
          ENDIF
          RU = DOT_PRODUCT(CHAINP%UVEC(:,1),CHAINP%POS(:,CHAINP%NPT)-CHAINP%POS(:,1))
          CALL GETCHAINRG(CHAINP,RG)
          IF (ENERGYFROMBD) THEN
             CALL GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,ENERGY2,.FALSE.)
          ELSE
             CALL GETENERGY(CHAINP,ENERGY2)
          ENDIF
          WRITE(55,*) STEP, FACC1, FACC2, R2AVG/NCT, DR, UU, RG, RU, &
               & ENERGY2, RHO1,RHO2,RHO3,RHO4
          FLUSH(55)
       ENDIF

       IF (DUMPSNAPSHOTS.AND.STEP.GT.INITSTEPS.AND.MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
          SNAPCT = SNAPCT + 1
          IF (OUTPUTBEADWEIGHT) THEN
             ! output weight of each mobile bead by integrating over U and R vectors
             BEADWEIGHT = 0D0
             DO B = 1,CHAINP%NPT
                IF (.NOT.CHAINP%ISFIXED(B)) THEN
                   CALL GETBEADQUINT(CHAINP,B,INTUWEIGHTNPT,INTUWEIGHTNPT,QINT)
                   BEADWEIGHT(1,B) = QINT
                   CALL GETBEADQRINT(CHAINP,B,INTRWEIGHTNPT,QINT)
                   BEADWEIGHT(2,B) = QINT
                ENDIF
             ENDDO

             CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,STEP,SNAPCT.GT.1.OR.STARTSTEP.GT.0,2,BEADWEIGHTPT)
          ELSE
             CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,STEP,SNAPCT.GT.1.OR.STARTSTEP.GT.0)
          ENDIF
       ENDIF

       IF (ADJUSTEVERY.GT.0.AND.MOD(STEP,ADJUSTEVERY).EQ.0) THEN
          ! adjust step range if needed


          ADJUSTED =.FALSE.
          IF (FACC1.GT.FACCTARGET+FACCTOL) THEN
             IF (.NOT.(ARANGE1.GE.2*PI.AND..NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE)) THEN
                ARANGE1 = ARANGE1*ADJUSTSCL; RRANGE1 = RRANGE1*ADJUSTSCL
                ADJUSTED = .TRUE.             
             ENDIF
          ELSE IF (FACC1.LT.FACCTARGET-FACCTOL) THEN
             ARANGE1 = ARANGE1/ADJUSTSCL; RRANGE1 = RRANGE1/ADJUSTSCL
             ADJUSTED = .TRUE.        
          ENDIF
          ARANGE1 = MIN(ARANGE1,2*PI)

          IF (FACC2.GT.FACCTARGET+FACCTOL) THEN
             IF (.NOT.(ARANGE2.GE.2*PI.AND..NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE)) THEN
                ARANGE2 = ARANGE2*ADJUSTSCL; RRANGE2 = RRANGE2*ADJUSTSCL
                ADJUSTED = .TRUE.   
             ENDIF
          ELSE IF (FACC2.LT.FACCTARGET-FACCTOL) THEN
             ARANGE2 = ARANGE2/ADJUSTSCL; RRANGE2 = RRANGE2/ADJUSTSCL
             ADJUSTED = .TRUE.        
          ENDIF
          ARANGE2 = MIN(ARANGE2,2*PI)

          IF (ADJUSTED) PRINT*, 'ADJUSTED RANGES:', ARANGE1, RRANGE1, ARANGE2, RRANGE2
       ENDIF

    ENDDO

    verbose = .true.
    IF (ENERGYFROMBD) THEN
       CALL GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,prevenergy,.FALSE.)
    ELSE
       CALL GETENERGY(CHAINP,PREVENERGY)
    ENDIF
    PRINT*, 'FINAL ENERGY: ', PREVENERGY

    CLOSE(55)
    CALL CLEANUPCHAIN(PREVCHAINP)
  END SUBROUTINE RUNMONTECARLO1CHAIN

   SUBROUTINE LOCALMOVE2bead(CHAINP,CHAINP2,B1,B2,ARANGE,RRANGE,DOROT1,DOSHIFT1,DOROT2,DOSHIFT2)
     ! do a local move of bead B1 and the bead connected to it b2 as well
     ! DOROT and DOSHIFT set whether to rotate and/or shift each bead
    IMPLICIT NONE
    
    TYPE(CHAIN), POINTER :: CHAINP,chainp2
    DOUBLE PRECISION, INTENT(IN) :: ARANGE, RRANGE
    INTEGER, INTENT(IN) :: B1,b2
    LOGICAL, INTENT(IN) :: DOROT1, DOSHIFT1,DOROT2,DOSHIFT2
    DOUBLE PRECISION :: AX(3), ANG, DELR(3), ROTMAT(3,3), TMP(3)
    INTEGER :: I
    
    IF (DOROT1.OR.DOROT2) THEN
       ! get a random axis
       DO I = 1,3
          AX(I) = GRND()
       ENDDO
       AX = AX/SQRT(DOT_PRODUCT(AX,AX))
       
       ! get a random angle
       ANG = GRND()*2*ARANGE-ARANGE
    ENDIF

    ! get a random shift
    IF (DOSHIFT1.OR.DOSHIFT2) THEN
       IF (CHAINP%STRETCHABLE) THEN
          DO I = 1,3
             DELR(I) = GRND()*2*RRANGE - RRANGE
          ENDDO
!       ELSEIF (CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
!          PRINT*, 'STRETCHABLE AND NOT SHEARABLE NOT YET SET UP'
!          STOP 2
       ELSEIF (CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
          PRINT*, 'SHEARABLE AND NOT STRETCHABLE NOT YET SET UP'
          STOP 2
       ELSE
          DELR = 0D0
       ENDIF
    ENDIF

    ! rotate selected bead
    IF (DOROT1) THEN
       CALL ROTANGAX(ANG,AX,CHAINP%UVEC(:,B1),TMP,.TRUE.,ROTMAT)
       CHAINP%UVEC(:,B1) = TMP
    ENDIF
    IF (DOSHIFT1) THEN
       ! shift selected bead
       CHAINP%POS(:,B1) = CHAINP%POS(:,B1) + DELR
    ENDIF

    ! rotate and shift other bead
    IF (DOROT2) THEN
       CALL ROTANGAX(ANG,AX,CHAINP2%UVEC(:,B2),TMP,(.NOT.DOROT1),ROTMAT)
       CHAINP2%UVEC(:,B2) = TMP
    ENDIF
    IF (DOSHIFT2) CHAINP2%POS(:,B2) = CHAINP2%POS(:,B2) + DELR
  END SUBROUTINE LOCALMOVE2BEAD

  SUBROUTINE LOCALMOVE(CHAINP,ARANGE,RRANGE,B,DOROT,DOSHIFT)
    ! do a local move (translation + rotation) of bead B
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: ARANGE, RRANGE
    INTEGER, INTENT(IN) :: B
    LOGICAL, INTENT(IN) :: DOROT, DOSHIFT
    DOUBLE PRECISION :: AX(3), ANG, DELR(3), ROTMAT(3,3), TMP(3)
    INTEGER :: I

    IF (DOROT) THEN
       ! get a random axis
       DO I = 1,3
          AX(I) = GRND()
       ENDDO
       AX = AX/SQRT(DOT_PRODUCT(AX,AX))

       ! get a random angle
       ANG = GRND()*2*ARANGE-ARANGE
    ENDIF

    IF (DOSHIFT) THEN
       ! get a random shift
       IF (CHAINP%STRETCHABLE) THEN
          DO I = 1,3
             DELR(I) = GRND()*2*RRANGE - RRANGE
          ENDDO
!       ELSEIF (CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
!          PRINT*, 'STRETCHABLE AND NOT SHEARABLE NOT YET SET UP'
!          STOP 2
       ELSEIF (CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
          PRINT*, 'SHEARABLE AND NOT STRETCHABLE NOT YET SET UP'
          STOP 2
       ELSE
          DELR = 0D0
       ENDIF
    ENDIF

    IF (DOROT) THEN
       ! rotate selected bead
       CALL ROTANGAX(ANG,AX,CHAINP%UVEC(:,B),TMP,.TRUE.,ROTMAT)
       CHAINP%UVEC(:,B) = TMP
    ENDIF
    IF (DOSHIFT) THEN
       ! shift selected bead
       CHAINP%POS(:,B) = CHAINP%POS(:,B) + DELR
    ENDIF
  END SUBROUTINE LOCALMOVE

  SUBROUTINE CRANKMOVE2(CHAINP,ARANGE,RRANGE,BPIVOT,DELE)
    ! carry out a Monte Carlo move involving rotating (and shifting) chain between 
    ! 2 beads around axis between them
    ! ARANGE is max rotation angle
    ! RRANGE is range of shift

    IMPLICIT NONE

    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: ARANGE, RRANGE
    INTEGER, INTENT(OUT) :: BPIVOT(2)
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: I, B, TMPI
    DOUBLE PRECISION :: ANG, AX(3), DELR(3), R0(3), R(3), PREVE,TMP(3)
    DOUBLE PRECISION :: ROTMAT(3,3)
    LOGICAL :: CALCROTMAT, CLASH
    INTEGER :: STARTB, FIXB1, FIXB2, NFIXSEG, FIXSEG

    !---------
    ! get 2 random beads around which to pivot    
    IF (CHAINP%NFIXBEAD.GT.0) THEN
       ! how many possible segments between fixed beads are there?   
       NFIXSEG = CHAINP%NFIXBEAD-1
       IF (CHAINP%FIXBEAD(CHAINP%NFIXBEAD).LT.CHAINP%NPT) THEN
          NFIXSEG = NFIXSEG + 1
       ENDIF
       IF (CHAINP%FIXBEAD(1).GT.1) THEN
          NFIXSEG = NFIXSEG + 1
       ENDIF

       ! find which segment between fixed beads to work in
       FIXSEG = FLOOR(GRND()*NFIXSEG)+1       
       IF (CHAINP%FIXBEAD(1).GT.1) THEN
          IF (FIXSEG.EQ.1) THEN
             FIXB1 = 1
          ELSE
             FIXB1 = CHAINP%FIXBEAD(FIXSEG-1)
          ENDIF
          IF (FIXSEG.GT.CHAINP%NFIXBEAD) THEN
             FIXB2 = CHAINP%NPT
          ELSE
             FIXB2 = CHAINP%FIXBEAD(FIXSEG)
          ENDIF
       ELSE
          FIXB1 = CHAINP%FIXBEAD(FIXSEG)
          IF (FIXSEG.GE.CHAINP%NFIXBEAD) THEN
             FIXB2 = CHAINP%NPT
          ELSE
             FIXB2 = CHAINP%FIXBEAD(FIXSEG+1)
          ENDIF
       ENDIF

       ! pick the actual pivot beads within the fixed segment, including edges
       BPIVOT = 0
       DO WHILE (BPIVOT(1).EQ.BPIVOT(2))
          DO I = 1,2
             BPIVOT(I) = FLOOR(GRND()*(FIXB2-FIXB1+1))+FIXB1
          ENDDO
       ENDDO
    ELSE
       BPIVOT = 0
       DO WHILE (BPIVOT(1).EQ.BPIVOT(2))
          DO I = 1,2
             BPIVOT(I) = FLOOR(GRND()*CHAINP%NPT)+1
          ENDDO
       ENDDO
    ENDIF
    ! order the pivot beads
    IF (BPIVOT(2).LT.BPIVOT(1)) THEN
       TMPI = BPIVOT(1)
       BPIVOT(1) = BPIVOT(2)
       BPIVOT(2) = TMPI
    ENDIF

    ! get axis between them
    AX = CHAINP%POS(:,BPIVOT(2))-CHAINP%POS(:,BPIVOT(1))
    AX = AX/SQRT(AX(1)*AX(1)+AX(2)*AX(2)+AX(3)*AX(3))

    ! get a random angle
    ANG = GRND()*2*ARANGE-ARANGE

    ! PRINT*, 'TESTX2:', ANG, AX
    ! DO I = 1,chainp%npt
    !    PRINT '(A,I3,8G15.5)', 'TESTX4:', I, CHAINP%POS(:,I), CHAINP%QUATS(I)%W,  & 
    !         & CHAINP%QUATS(I)%X,  CHAINP%QUATS(I)%Y,  CHAINP%QUATS(I)%Z, chainp%beadenergy(i)
    ! ENDDO
    ! get a random shift
    IF (CHAINP%STRETCHABLE) THEN
       DO I = 1,3
          DELR(I) = GRND()*2*RRANGE - RRANGE
       ENDDO
       !    ELSEIF (CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
       !       PRINT*, 'STRETCHABLE AND NOT SHEARABLE NOT YET SET UP'
       !       STOP 2
    ELSEIF (CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
       PRINT*, 'SHEARABLE AND NOT STRETCHABLE NOT YET SET UP'
       STOP 2
    ELSE
       DELR = 0D0
    ENDIF

    ! rotate and shift everything between the pivot beads, including first but not last bead
    ! do not include first bead if it is fixed
    R0 = CHAINP%POS(:,BPIVOT(1))
    IF (CHAINP%ISFIXED(BPIVOT(1))) THEN
       STARTB = BPIVOT(1)+1
    ELSE
       STARTB = BPIVOT(1)
    ENDIF
    DO B = STARTB,BPIVOT(2)-1
       CALCROTMAT = B.EQ.STARTB
       CALL ROTANGAX(ANG,AX,CHAINP%UVEC(:,B),TMP,CALCROTMAT,ROTMAT)
       CHAINP%UVEC(:,B) = TMP

       CALL ROTANGAX(ANG,AX,CHAINP%POS(:,B)-R0,R,.FALSE.,ROTMAT)
       CHAINP%POS(:,B) = R0+R+DELR
    ENDDO

    ! update pivot bead energies
    PREVE = CHAINP%BEADENERGY(STARTB)+ CHAINP%BEADENERGY(BPIVOT(2))+CHAINP%STERICENERGY

    CALL GETBEADENERGY(CHAINP,STARTB,CHAINP%BEADENERGY(STARTB))
    CALL GETBEADENERGY(CHAINP,BPIVOT(2),CHAINP%BEADENERGY(BPIVOT(2)))

    !IF (CHAINP%STERICS) CALL GETSTERICENERGY(CHAINP,CHAINP%STERICENERGY)

    DELE = CHAINP%BEADENERGY(STARTB) + CHAINP%BEADENERGY(BPIVOT(2))+CHAINP%STERICENERGY - PREVE


    IF (STARTB.EQ.1) THEN
       ! end to end extension also changes
       PREVE =  CHAINP%FORCEENERGY     
       CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY)
       DELE = DELE+CHAINP%FORCEENERGY-PREVE
    ENDIF

    IF (CHAINP%STERICS) THEN       
       CALL CHECKCLASHMV2(CHAINP,BPIVOT(1),BPIVOT(2),CLASH)
       IF (CLASH) THEN
          DELE = HUGE(1D0)
       ENDIF
    ENDIF
  END SUBROUTINE CRANKMOVE2

  SUBROUTINE CHECKCLASHMV2(CHAINP,B1,B2,CLASH)
    ! check for a steric clash after a crank move of type 2
    ! where beads B1 through B2 (but not B2 itself) are moved
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B1, B2
    LOGICAL, INTENT(OUT) :: CLASH
    !DOUBLE PRECISION, INTENT(OUT) :: STERENERGY
    INTEGER :: I, J
    DOUBLE PRECISION :: DIFF(3), ND2

   ! STERENERGY = 0D0
    DO I = B1,B2-1
       DO J = 1,B1-1 ! prior stationary beads
          IF (ABS(I-J).LE.CHAINP%STERSKIP) CYCLE
          DIFF = CHAINP%POS(:,I)-CHAINP%POS(:,J)
          ND2 = DOT_PRODUCT(DIFF,DIFF)
          IF (ND2.LT.CHAINP%STERRAD2) THEN
!             STERENERGY = STERENERGY + CHAINP%STERMOD*(CHAINP%STERRAD2-ND2)**2
             CLASH = .TRUE.
             RETURN
          ENDIF
       ENDDO

       DO J = B2,CHAINP%NPT ! subsequent stationary beads
          IF (ABS(I-J).LE.CHAINP%STERSKIP) CYCLE
          DIFF = CHAINP%POS(:,I)-CHAINP%POS(:,J)
          ND2 = DOT_PRODUCT(DIFF,DIFF)
          IF (ND2.LT.CHAINP%STERRAD2) THEN
 !            STERENERGY = STERENERGY + CHAINP%STERMOD*(CHAINP%STERRAD2-ND2)**2
             CLASH = .TRUE.
             RETURN
          ENDIF
       ENDDO
    ENDDO
    
  END SUBROUTINE CHECKCLASHMV2

  SUBROUTINE CRANKMOVE1(CHAINP,ARANGE,RRANGE,BPIVOT,DELE)
    ! carry out a Monte Carlo move involving rotating chain from a random bead B
    ! around an arbitrary axis
    ! ARANGE gives maximal rotation angle
    ! RRANGE gives the range for the position shift in each dimension

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: ARANGE,RRANGE
    INTEGER, INTENT(OUT) :: BPIVOT
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: I, B
    DOUBLE PRECISION :: ANG, AX(3), DELR(3), R0(3), R(3), PREVE, TMP(3)
    DOUBLE PRECISION :: ROTMAT(3,3)
    LOGICAL :: CALCROTMAT, CLASH
    INTEGER :: NPIV1, NPIV2, BEADUPDATE, STARTB, ENDB

    ! get a random axis
    DO I = 1,3
       AX(I) = GRND()
    ENDDO
    AX = AX/SQRT(DOT_PRODUCT(AX,AX))

    ! get a random angle
    ANG = GRND()*2*ARANGE-ARANGE

    ! get a random shift of the chain fragment
    IF (CHAINP%STRETCHABLE) THEN
       DO I = 1,3
          DELR(I) = GRND()*2*RRANGE-RRANGE
       ENDDO
       !    ELSEIF (CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
       !       PRINT*, 'STRETCHABLE AND NOT SHEARABLE NOT YET SET UP'
       !       STOP 2
    ELSEIF (CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
       PRINT*, 'SHEARABLE AND NOT STRETCHABLE NOT YET SET UP'
       STOP 2
    ELSE
       DELR = 0D0
    ENDIF


    ! get a pivot bead at random    
    IF (CHAINP%NFIXBEAD.EQ.0) THEN 
       BPIVOT = FLOOR(GRND()*(CHAINP%NPT-1))+2
    ELSE
       ! number of potential pivot beads at start of chain
       IF (CHAINP%FIXBEAD(1).EQ.1) THEN
          NPIV1 = 0
       ELSE
          NPIV1 = CHAINP%FIXBEAD(1)
       ENDIF
       ! number of potential pivot beads at end of chain
       IF (CHAINP%FIXBEAD(CHAINP%NFIXBEAD).EQ.CHAINP%NPT) THEN
          NPIV2 = 0
       ELSE
          NPIV2 = CHAINP%NPT - CHAINP%FIXBEAD(CHAINP%NFIXBEAD)+1
       ENDIF
      

       B = FLOOR(GRND()*(NPIV1+NPIV2))+1
       IF (B.LE.NPIV1) THEN
          BPIVOT = B
       ELSE
          BPIVOT = B - NPIV1+CHAINP%FIXBEAD(CHAINP%NFIXBEAD)-1
       ENDIF
    ENDIF

    IF (BPIVOT.LT.1.OR.BPIVOT.GT.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN CRANKMOVE1: BAD BPIVOT', BPIVOT, NPIV1, NPIV2
    ENDIF

    IF (CHAINP%NFIXBEAD.GT.0.AND.BPIVOT.LE.NPIV1) THEN
       STARTB = 1
       IF (CHAINP%ISFIXED(BPIVOT)) THEN
          ENDB = BPIVOT - 1
          BEADUPDATE = BPIVOT
       ELSE
          ENDB = BPIVOT
          BEADUPDATE = BPIVOT + 1
       ENDIF       
       IF (BEADUPDATE.GT.CHAINP%NPT) THEN
          PRINT*, 'ERROR IN CRANKMV1:', BPIVOT
          STOP 1
       ENDIF
    ELSE
       IF (CHAINP%ISFIXED(BPIVOT)) THEN
          STARTB = BPIVOT+1
          BEADUPDATE = BPIVOT+1
       ELSE
          STARTB = BPIVOT
          BEADUPDATE = BPIVOT
       ENDIF
       ENDB = CHAINP%NPT  
    ENDIF

    IF (STARTB.LT.0.OR.ENDB.GT.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN CRANKMOVE1: BAD STARTB OR ENDB', STARTB, ENDB
    ENDIF

    ! rotate and shift everything from bead BPIVOT to edge (including BPIVOT orientation)
    R0 = CHAINP%POS(:,BPIVOT)
    DO B = STARTB,ENDB
       CALCROTMAT = B.EQ.STARTB
       CALL ROTANGAX(ANG,AX,CHAINP%UVEC(:,B),TMP,CALCROTMAT,ROTMAT)
       CHAINP%UVEC(:,B) = TMP

       CALL ROTANGAX(ANG,AX,CHAINP%POS(:,B)-R0,R,.FALSE.,ROTMAT)
       CHAINP%POS(:,B) = R0+R+DELR
    ENDDO


    ! update pivot bead energy (for segment before pivot bead)    
    PREVE = CHAINP%BEADENERGY(BEADUPDATE) + CHAINP%FORCEENERGY+CHAINP%STERICENERGY
    CALL GETBEADENERGY(CHAINP,BEADUPDATE,CHAINP%BEADENERGY(BEADUPDATE) )
    CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY)

    ! overall change in energy
    DELE = CHAINP%BEADENERGY(BEADUPDATE) + CHAINP%FORCEENERGY +CHAINP%STERICENERGY- PREVE



    IF (CHAINP%STERICS) THEN
       IF (CHAINP%NFIXBEAD.GT.0) THEN
          PRINT*, 'ERROR: STERICS NOT SET UP WITH FIXBEADS'
          STOP 1
       ENDIF
       CALL CHECKCLASHMV1(CHAINP,BPIVOT,CLASH)
       IF (CLASH) THEN
          ! PRINT*, 'CLASH!'
          DELE = HUGE(1D0)
       ENDIF
    ENDIF


  END SUBROUTINE CRANKMOVE1

    SUBROUTINE CHECKCLASHMV1(CHAINP,B, CLASH)
    ! check for a steric clash after a crank move of type 1
    ! where beads B through end (inclusive) are moved
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B
    LOGICAL, INTENT(OUT) :: CLASH
    INTEGER :: I, J
    DOUBLE PRECISION :: DIFF(3), ND2

    CLASH = .FALSE.
    DO I = B,CHAINP%NPT
       DO J = 1,B-1 ! prior stationary beads
          IF (ABS(I-J).LE.CHAINP%STERSKIP) CYCLE
          DIFF = CHAINP%POS(:,I)-CHAINP%POS(:,J)
          ND2 = DOT_PRODUCT(DIFF,DIFF)
          IF (ND2.LT.CHAINP%STERRAD2) THEN
             CLASH = .TRUE.
             RETURN
          ENDIF
       ENDDO
    ENDDO
    
  END SUBROUTINE CHECKCLASHMV1


  SUBROUTINE GETSTERMV1(CHAINP,B, STERENERGY)
    ! check for a steric clash after a crank move of type 1
    ! where beads B through end (inclusive) are moved
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: B
    DOUBLE PRECISION, INTENT(OUT) :: STERENERGY
    INTEGER :: I, J
    DOUBLE PRECISION :: DIFF(3), ND2

    STERENERGY = 0D0
    DO I = B,CHAINP%NPT
       DO J = 1,B-1 ! prior stationary beads
          IF (ABS(I-J).LE.CHAINP%STERSKIP) CYCLE
          DIFF = CHAINP%POS(:,I)-CHAINP%POS(:,J)
          ND2 = DOT_PRODUCT(DIFF,DIFF)
          IF (ND2.LT.CHAINP%STERRAD2) THEN
             STERENERGY = STERENERGY + CHAINP%STERMOD*(ND2-CHAINP%STERRAD2)**2
          ENDIF
       ENDDO
    ENDDO
    
  END SUBROUTINE GETSTERMV1

END MODULE MONTECARLO
