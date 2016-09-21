MODULE BROWNDYN
  USE GENUTIL
  ! utilities for running BD simulations with the chain
  IMPLICIT NONE
  ! extra printing for this chain
  INTEGER :: VERBOSECHAIN = 18160
  LOGICAL :: VERBOSE = .FALSE.

CONTAINS
  SUBROUTINE RUNBROWNDYN(CHAINLIST,NCHAIN,NSTEP,DELT,KT,OUTFILE,RUNGEKUTTA,DOBROWN)
    ! run brownian dynamics for NSTEP steps for a bunch of chains in parallel
    ! using time step DELT and temperature KT
    ! Periodically print out energy and end-to-end vector for all chains
    ! RUNGEKUTTA = 1 for euler method or 4 for 4-th order runge-kutta
    ! if DOBROWN=false, do not include brownian forces
    USE CHAINUTIL, ONLY : CHAIN, OUTPUTSNAPSHOT, RESETLINKERS
    USE KEYS,ONLY: BDPRINTEVERY, BDPRINTLOG, STRESSFILE, SNAPSHOTFILE, &
         & SNAPSHOTEVERY, DUMPSNAPSHOTS, BRCRELAX, LOOPRAD, TRACKLOOPING, LOOPFILE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NCHAIN, NSTEP
    DOUBLE PRECISION, INTENT(IN) :: DELT, KT
    CHARACTER(LEN=*) :: OUTFILE
    INTEGER, INTENT(IN) :: RUNGEKUTTA
    LOGICAL :: DOBROWN
    TYPE(CHAIN), INTENT(IN),TARGET :: CHAINLIST(NCHAIN)
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: STEP, C, B
    DOUBLE PRECISION :: DR(3), ENERGY, COM(3), MEANSTRESSCORR, STDSTRESSCORR
    DOUBLE PRECISION :: STRESSXY(6), STRESSXY0(NCHAIN,6),STRESSCORR(NCHAIN)
    INTEGER :: NEXTPRINT
    LOGICAL :: START(nchain), FIRSTSNAP, reset, DONELOOP(NCHAIN)
    DOUBLE PRECISION, ALLOCATABLE :: TMPARRAY(:,:),TMPARRAY2(:,:) 
    DOUBLE PRECISION :: LOOPRAD2

    allocate(TMPARRAY(3,CHAINLIST(1)%NPT), TMPARRAY2(3,CHAINLIST(1)%NPT))

    PRINT*, 'Running BD simulation with ', NCHAIN, 'chains, for ', NSTEP, 'steps'
    START = .TRUE.
    FIRSTSNAP = .TRUE.

    IF (BDPRINTLOG) THEN
       NEXTPRINT = 1
    ENDIF

    OPEN(88,FILE=STRESSFILE,STATUS='UNKNOWN')
    write(88,*) NCHAIN, DELT, KT, CHAINLIST(1)%NPT

    OPEN(55,FILE=OUTFILE,STATUS='UNKNOWN')
    write(55,*) NCHAIN, DELT, KT  

    OPEN(77,FILE=SNAPSHOTFILE,STATUS='UNKNOWN')
    WRITE(77,*) 'X', NCHAIN,DELT,CHAINLIST(1)%NPT
    CLOSE(77)
    CHAINP=>CHAINLIST(1)
    CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,0,.true.)

    DONELOOP = .FALSE.
    LOOPRAD2 = LOOPRAD*LOOPRAD
    IF (TRACKLOOPING) THEN
       OPEN(44,FILE=LOOPFILE,STATUS='UNKNOWN')
       WRITE(44,*) NCHAIN, DELT, CHAINLIST(1)%NPT, LOOPRAD
    ENDIF

    IF (TRACKLOOPING) THEN
       DO C = 1,NCHAIN
          CHAINP=>CHAINLIST(C)    
          DR = CHAINP%POS(:,CHAINP%NPT)-CHAINP%POS(:,1)
          IF (DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3).LT.LOOPRAD2) THEN
             print*, 'Chain ', C, ' has successfully looped at time ', 0
             WRITE(44,*) C,0,0D0, DR
             FLUSH(44)
             DONELOOP(C)=.TRUE.
          ENDIF
       ENDDO
    ENDIF

    DO STEP=1,NSTEP
       IF (ALL(DONELOOP)) THEN
          PRINT*, 'ALL CHAINS HAVE LOOPED', STEP 
          EXIT
       ENDIF
       DO C = 1,NCHAIN          
          ! stop propagating this chain if it has already looped
          IF (DONELOOP(C)) CYCLE

          CHAINP=>CHAINLIST(C)        

          IF (.NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
             !reset all linkers to ground-state
             !this avoids numerical errors that push the lengths away from
             !the desired values
             CALL RESETLINKERS(CHAINP,0.005D0,RESET)

             !IF (RESET) PRINT*, 'resetting linker lengths:', STEP, C
          ENDIF

          ! print*, 'testx1:', sqrt(sum((chainp%pos(:,2)-chainp%pos(:,1))**2))
          IF ((BDPRINTLOG.AND.STEP.EQ.NEXTPRINT).OR.&
               & (.NOT.BDPRINTLOG.AND.MOD(STEP,NINT(BDPRINTEVERY)).EQ.0)) THEN
             IF (RUNGEKUTTA.EQ.4) THEN
                !                VERBOSE = C.EQ.VERBOSECHAIN
                CALL LANGEVINSTEPRK4(CHAINP,DELT,ENERGY,DOBROWN,STRESSXY)
                !                IF (VERBOSE) PRINT*, 'TESTX3:', ENERGY
             ELSE IF (RUNGEKUTTA.EQ.1) THEN
                CALL LANGEVINSTEP(CHAINP,DELT,KT,ENERGY,STRESSXY)
             ELSE
                PRINT*, 'RUNGEKUTTA must be 1 or 4', RUNGEKUTTA
                stop 1
             END IF

             ! IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
             !    CALL GETCHAINFORCEINT(CHAINP,TMPARRAY,TMPARRAY2,ENERGY,.FALSE.)         
             ! ELSE
             !    CALL GETBEADRODFORCE(CHAINP,TMPARRAY,BRCRELAX,TMPARRAY2,ENERGY)
             ! ENDIF

             ! IF (.NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
             !    ! reset all linkers to groundstat
             !    ! this avoids numerical errors that push the lengths away from
             !    ! the desired values
             !    CALL RESETLINKERS(CHAINP)
             ! ENDIF

             DR = CHAINP%POS(:,CHAINP%NPT)-CHAINP%POS(:,1)
             COM(1) = SUM(CHAINP%POS(1,:))
             COM(2) = SUM(CHAINP%POS(2,:))
             COM(3) = SUM(CHAINP%POS(3,:))
             COM = COM/CHAINP%NPT         
             IF (START(C)) THEN
                STRESSXY0(C,:) = STRESSXY
                START(C) = .FALSE.
             ENDIF

             STRESSCORR(c) = SUM(STRESSXY*STRESSXY0(C,:))/6
             ! if (stresscorr(c).gt.1000) then
             !    print*, c,stresscorr(c)
             !    do b = 1,chainp%npt
             !       print*, chainp%pos(:,b)
             !    enddo
             !    print*, 'testx2'
             !    do b = 1,chainp%npt
             !       print*, chainp%uvec(:,b)
             !    enddo
             !    stop 1
             ! endif
             !STRESSCORR(c) = STRESSXY(1)*STRESSXY0(C,1)



             IF (C.EQ.1) PRINT*, 'STEP, ENERGY:', STEP, ENERGY, DR
             IF (.NOT.CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
                CHAINP%UVEC(:,1) = CHAINP%POS(:,2)-CHAINP%POS(:,1)              
                IF (CHAINP%NPT.GE.3) THEN
                   COM = CHAINP%POS(:,3)-CHAINP%POS(:,2)
                ENDIF
             ENDIF
             WRITE(55,*) STEP, C, ENERGY, DR, COM, CHAINP%UVEC(:,1)!,STRESSXY    
             FLUSH(55)
          ELSE
             IF (RUNGEKUTTA.EQ.4) THEN
                CALL LANGEVINSTEPRK4(CHAINP,DELT,ENERGY,DOBROWN)
             ELSE IF (RUNGEKUTTA.EQ.1) THEN
                CALL LANGEVINSTEP(CHAINP,DELT,KT,ENERGY)
             ELSE
                PRINT*, 'RUNGEKUTTA must be 1 or 4', RUNGEKUTTA
                stop 1
             END IF
          ENDIF

          IF (TRACKLOOPING) THEN
             DR = CHAINP%POS(:,CHAINP%NPT)-CHAINP%POS(:,1)
             IF (DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3).LT.LOOPRAD2) THEN
                print*, 'Chain ', C, ' has successfully looped at time ', STEP, STEP*DELT
                WRITE(44,*) C,STEP, STEP*DELT, DR
                FLUSH(44)
                DONELOOP(C)=.TRUE.
             ENDIF
          ENDIF

          IF (DUMPSNAPSHOTS.AND.MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
             CALL OUTPUTSNAPSHOT(CHAINP,SNAPSHOTFILE,STEP,.true.)
             FIRSTSNAP=.FALSE.
             !DO B = 1,CHAINP%NPT
             !   WRITE(77,*) STEP, C, B, CHAINP%POS(:,B)
             !ENDDO
          ENDIF

       ENDDO
       IF ((BDPRINTLOG.AND.STEP.EQ.NEXTPRINT).OR.&
            & (.NOT.BDPRINTLOG.AND.MOD(STEP,NINT(BDPRINTEVERY)).EQ.0)) THEN
          MEANSTRESSCORR = SUM(STRESSCORR)/NCHAIN
          STDSTRESSCORR = SQRT(SUM((STRESSCORR-MEANSTRESSCORR)**2)/NCHAIN)     
          WRITE(88,*) STEP,MEANSTRESSCORR, STDSTRESSCORR
          FLUSH(88)
       END IF
       IF (BDPRINTLOG.AND.STEP.EQ.NEXTPRINT) THEN
          NEXTPRINT = CEILING(NEXTPRINT * EXP(BDPRINTEVERY*LOG(10D0)))
          IF (NEXTPRINT.GT.NSTEP) NEXTPRINT = NSTEP
       ENDIF
    ENDDO
    CLOSE(55)
    CLOSE(88)
    IF (TRACKLOOPING) CLOSE(44)
    !    CLOSE(77)

    DEALLOCATE(TMPARRAY,TMPARRAY2)
  END SUBROUTINE RUNBROWNDYN

  SUBROUTINE LANGEVINSTEPRK4(CHAINP,DELT,ENERGY,DOBROWN,STRESSXY)
    ! propagate forward using a fourth-order Runge-Kutta method

    USE KEYS, ONLY : FIXBEAD1,FIXBEADMID,BRCRELAX
    USE MT19937, ONLY : RNORM, grnd
    USE CHAINUTIL, ONLY : CHAIN, OBSTACLE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(OBSTACLE), POINTER :: OBP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    LOGICAL, INTENT(IN) :: DOBROWN
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: STRESSXY(6)
    DOUBLE PRECISION :: POSFORCE(3,CHAINP%NPT), UFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION :: TENSBRFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION :: BROWNPOS(3,CHAINP%NPT), BROWNU(3,CHAINP%NPT)
    INTEGER :: C,M, B, I, J,CT
    DOUBLE PRECISION :: BROWNF, F, S2DTR(CHAINP%NPT), S2DTU(CHAINP%NPT), S2DTO, RHO, PHI,ST,CP,SP
    DOUBLE PRECISION :: XVEC(3), YVEC(3), POS0(3,CHAINP%NPT),UVEC0(3,CHAINP%NPT)
    DOUBLE PRECISION :: BROWNF1, BROWNF2, TORQUE(3), COM(3)
    DOUBLE PRECISION, DIMENSION(3,CHAINP%NPT) :: K1POS,K1UVEC,K2POS,K2UVEC,K3POS,K3UVEC,K4POS,K4UVEC,POSFORCE0
    DOUBLE PRECISION :: CRELAX

    IF (FIXBEAD1.OR.FIXBEADMID) THEN
       PRINT*, 'FIXBEAD1 AND FIXBEADMID not yet set up with runge kutta', FIXBEAD1, FIXBEADMID
       STOP 1
    ENDIF

    ! coefficient for the additional force to keep bond lengths constrained
    IF (.NOT.CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
       CRELAX = MINVAL(CHAINP%FRICTR(1:CHAINP%NPT))*BRCRELAX/DELT 
    ENDIF

    S2DTR = SQRT(2*CHAINP%FRICTR(1:CHAINP%NPT)/DELT)
    S2DTU = SQRT(2*CHAINP%FRICTU(1:CHAINP%NPT)/DELT)

    POS0 = CHAINP%POS; UVEC0 = CHAINP%UVEC
    IF (PRESENT(STRESSXY)) THEN
       STRESSXY = 0D0            
    ENDIF

    IF (DOBROWN) THEN
       ! get the brownian forces
       BROWNU = 0D0; BROWNPOS = 0D0   

       DO B = 1,CHAINP%NPT
          ! translational brownian force
          BROWNPOS(1,B) = RNORM()*S2DTR(B)
          BROWNPOS(2,B) = RNORM()*S2DTR(B)
          BROWNPOS(3,B) = RNORM()*S2DTR(B)

          IF (CHAINP%SHEARABLE) THEN
             ! rotational brownian force
             BROWNF1 = RNORM()*S2DTU(B)
             BROWNF2 = RNORM()*S2DTU(B)
             IF (CHAINP%UVEC(2,B).EQ.0D0.AND.CHAINP%UVEC(3,B).EQ.0D0) THEN
                XVEC = (/0D0,1D0,0D0/)
                YVEC = (/0D0,0D0,1D0/)
             ELSE
                CALL CROSS_PRODUCT(CHAINP%UVEC(:,B),(/1D0,0D0,0D0/),XVEC)
                CALL CROSS_PRODUCT(CHAINP%UVEC(:,B),XVEC,YVEC)
             ENDIF
             BROWNU(:,B) = BROWNF1*XVEC + BROWNF2*YVEC        
          ENDIF
       ENDDO
    ELSE
       BROWNPOS = 0D0; BROWNU = 0D0;
    ENDIF


    ! OPEN(UNIT=44,FILE='brownpos.out')
    ! DO B = 1,CHAINP%NPT
    !    WRITE(44,*) BROWNPOS(:,B)
    ! ENDDO
    ! close(44)

    ! --------- 1ST RK STEP---------------

    ! Get the deterministic forces
    IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
       CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.TRUE.)  
       ! IF (VERBOSE) THEN
       !    PRINT*, 'TESTX4:' 
       !    DO B = 1,CHAINP%NPT
       !       PRINT*, B, POSFORCE(:,B)
       !    ENDDO
       ! ENDIF
    ELSE
       ! also project out the components of the brownian force 
       ! that would change segment lengths
!       print*, 'testx2:', brownpos(:,5), s2dtr
       CALL GETBEADRODFORCE(CHAINP,BROWNPOS,CRELAX,POSFORCE,ENERGY,.true.)
    ENDIF

    IF (CHAINP%SHEARABLE) THEN
       ! Get just the perpendicular component of the U forces
       DO B = 1,CHAINP%NPT        
          UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)
          CALL CROSS_PRODUCT(UVEC0(:,B),UFORCE(:,B),TORQUE)
          CALL CROSS_PRODUCT(TORQUE,UVEC0(:,B),UFORCE(:,B))
          K1UVEC(:,B) = (BROWNU(:,B) + UFORCE(:,B))/CHAINP%FRICTU(B)
       ENDDO
    ENDIF

    DO B = 1,CHAINP%NPT
       K1POS(:,B) = (BROWNPOS(:,B) + POSFORCE(:,B))/CHAINP%FRICTR(B)
    ENDDO

    POSFORCE0 = POSFORCE;


    ! propagate forward
    CHAINP%POS = POS0 + DELT/2*K1POS    
    
!    print*, 'testx1:', chainp%pos(:,3)
    IF (CHAINP%SHEARABLE) THEN
       CHAINP%UVEC = UVEC0 + DELT/2*K1UVEC
       ! renormalize UVECs
       DO B = 1,CHAINP%NPT
          CALL NORMALIZE(CHAINP%UVEC(:,B))
       ENDDO
    ENDIF

    ! print*, 'testx5:'
    ! do b = 1,chainp%npt
    !    print*, b, chainp%pos(:,b)
    ! enddo
    ! stop 1
    !--------- 2ND RK STEP--------------
    ! Get the deterministic forces
    IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
       CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.TRUE.)  
    ELSE
       CALL GETBEADRODFORCE(CHAINP,BROWNPOS,CRELAX,POSFORCE,ENERGY,.FALSE.)
    ENDIF

    IF (CHAINP%SHEARABLE) THEN
       ! Get just the perpendicular component of the U forces
       DO B = 1,CHAINP%NPT        
          UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)
          CALL CROSS_PRODUCT(UVEC0(:,B),UFORCE(:,B),TORQUE)
          CALL CROSS_PRODUCT(TORQUE,UVEC0(:,B),UFORCE(:,B))
          K2UVEC(:,B) = (BROWNU(:,B) + UFORCE(:,B))/CHAINP%FRICTU(B)
       ENDDO
    ENDIF

    DO B = 1,CHAINP%NPT
       K2POS(:,B) = (BROWNPOS(:,B) + POSFORCE(:,B))/CHAINP%FRICTR(B)
    ENDDO

    ! propagate forward
    CHAINP%POS = POS0 + DELT/2*K2POS
    IF (CHAINP%SHEARABLE) THEN
       CHAINP%UVEC = UVEC0 + DELT/2*K2UVEC
       ! renormalize UVECs
       DO B = 1,CHAINP%NPT
          CALL NORMALIZE(CHAINP%UVEC(:,B))
       ENDDO
    ENDIF

    ! --------- 3RD RK STEP-----------------
    ! Get the deterministic forces
    IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
       CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.TRUE.)  
    ELSE
       CALL GETBEADRODFORCE(CHAINP,BROWNPOS,CRELAX,POSFORCE,ENERGY,.FALSE.)
    ENDIF

    IF (CHAINP%SHEARABLE) THEN
       ! Get just the perpendicular component of the U forces
       DO B = 1,CHAINP%NPT        
          UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)
          CALL CROSS_PRODUCT(UVEC0(:,B),UFORCE(:,B),TORQUE)
          CALL CROSS_PRODUCT(TORQUE,UVEC0(:,B),UFORCE(:,B))
          K3UVEC(:,B) = (BROWNU(:,B) + UFORCE(:,B))/CHAINP%FRICTU(B)
       ENDDO
    ENDIF

    DO B = 1,CHAINP%NPT
       K3POS(:,B) = (BROWNPOS(:,B) + POSFORCE(:,B))/CHAINP%FRICTR(B)
    ENDDO


    ! propagate forward
    CHAINP%POS = POS0 + DELT*K3POS
    IF (CHAINP%SHEARABLE) THEN
       CHAINP%UVEC = UVEC0 + DELT*K3UVEC
       ! renormalize UVECs
       DO B = 1,CHAINP%NPT
          CALL NORMALIZE(CHAINP%UVEC(:,B))
       ENDDO
    ENDIF


    ! ---------- 4TH RK STEP------------------
    ! Get the deterministic forces
    IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
       CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.TRUE.)  
    ELSE
       CALL GETBEADRODFORCE(CHAINP,BROWNPOS,CRELAX,POSFORCE,ENERGY,.FALSE.)
    ENDIF

    IF (CHAINP%SHEARABLE) THEN
       ! Get just the perpendicular component of the U forces
       DO B = 1,CHAINP%NPT        
          UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)
          CALL CROSS_PRODUCT(UVEC0(:,B),UFORCE(:,B),TORQUE)
          CALL CROSS_PRODUCT(TORQUE,UVEC0(:,B),UFORCE(:,B))
           K4UVEC(:,B) = (BROWNU(:,B) + UFORCE(:,B))/CHAINP%FRICTU(B)
       ENDDO       
    ENDIF

    DO B = 1,CHAINP%NPT
       K4POS(:,B) = (BROWNPOS(:,B) + POSFORCE(:,B))/CHAINP%FRICTR(B)
    ENDDO

    ! propagate forward
    CHAINP%POS = POS0 + DELT/6*(K1POS+2*K2POS+2*K3POS+K4POS)
    IF (CHAINP%SHEARABLE) THEN
       CHAINP%UVEC = UVEC0 + DELT/6*(K1UVEC+2*K2UVEC+2*K3UVEC+K4UVEC)
       ! renormalize UVECs
       DO B = 1,CHAINP%NPT
          CALL NORMALIZE(CHAINP%UVEC(:,B))
       ENDDO
    ENDIF

    ! calculate the stress tensor (XY component) on the beads
    IF (PRESENT(STRESSXY)) THEN
       IF (CHAINP%STRETCHABLE.OR.CHAINP%SHEARABLE) THEN
          CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.TRUE.)  
       ELSE
          CALL GETBEADRODFORCE(CHAINP,BROWNPOS,CRELAX,POSFORCE,ENERGY,.FALSE.,TENSBRFORCE)
       ENDIF
       COM(1) = SUM(CHAINP%POS(:,1))/CHAINP%NPT
       COM(2) = SUM(CHAINP%POS(:,2))/CHAINP%NPT
       COM(3) = SUM(CHAINP%POS(:,3))/CHAINP%NPT
       CT = 0
       DO I = 1,3
          DO J = 1,3
             IF (I.EQ.J) CYCLE
             CT = CT + 1
             IF (.NOT.CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN        
               STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*BROWNPOS(J,:) &                
                     & + (chainp%pos(I,:)-COM(I))*(POSFORCE(J,:)-0.5*TENSBRFORCE(J,:)))
               ! STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*BROWNPOS(J,:) &                
               !      & + (chainp%pos(I,:)-COM(I))*(POSFORCE(J,:)-TENSBRFORCE(J,:)))
               !STRESSXY(CT) = STRESSXY(CT)  - SUM(&
               !      &(chainp%pos(I,:)-COM(I))*(posforce(j,:)+BROWNPOS(J,:)))
               ! STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*(BROWNPOS(J,:)+TENSBRFORCE(J,:)) &
               !      & + (chainp%pos(I,:)-COM(I))*(POSFORCE(J,:)-TENSBRFORCE(J,:)))
             ELSE IF (.NOT.CHAINP%SHEARABLE) THEN
                STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*BROWNPOS(J,:) &
                     & + (chainp%pos(I,:)-COM(I))*(POSFORCE(J,:)-UFORCE(J,:)))
             ELSE
                STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*BROWNPOS(J,:) &
                     & + (chainp%pos(I,:)-COM(I))*POSFORCE(J,:))
             ENDIF
            ! STRESSXY(CT) = STRESSXY(CT)  + SUM((chainp%pos(I,:))*(POSFORCE(J,:)))
          ENDDO
       ENDDO
      
       ! DO B = 1,CHAINP%NPT       
       !    !STRESSXY = STRESSXY - (POS0(1,B)-COM(1))*BROWNPOS(2,B) - (POS0(1,B)-COM(1))*POSFORCE(2,B)
       !    !STRESSXY = STRESSXY  - (0.5*(CHAINP%POS(1,B)-POS0(1,B))*BROWNPOS(2,B) + (chainp%POS(1,B)-COM(1))*POSFORCE(2,B))
       !     DO I = 1,3
       !        DO J = 1,3
       !           IF (I.EQ.J) CYCLE
       !           STRESSXY = STRESSXY  - (0.5*(CHAINP%POS(1,B)-POS0(1,B))*BROWNPOS(2,B) + (POS0(1,B)-COM(1))*POSFORCE(2,B))/6
       !        ENDDO
       !     ENDDO
       ! END DO
    ENDIF

  END SUBROUTINE LANGEVINSTEPRK4

  SUBROUTINE LANGEVINSTEP(CHAINP,DELT,KT,ENERGY,STRESSXY)
    ! take a single step via overdamped Langevin dynamics    
   !!! ! use midpoint formula (eg: Grassia & Hinch, 1996)
    ! based on Tao et al, 2005
    ! for orientation vectors, project all forces to perpendicular
    ! and renormalize at each step
    ! returns energy from the midpoint of the step (multiplied by kT)
    USE KEYS, ONLY : FIXBEAD1,FIXBEADMID
    USE MT19937, ONLY : RNORM, grnd
    USE CHAINUTIL, ONLY : CHAIN, OBSTACLE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(OBSTACLE), POINTER :: OBP
    DOUBLE PRECISION, INTENT(IN) :: KT,DELT
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: STRESSXY(6)
    DOUBLE PRECISION :: POSFORCE(3,CHAINP%NPT), UFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION :: BROWNPOS(3,CHAINP%NPT), BROWNU(3,CHAINP%NPT)
    INTEGER :: C,M, B, I ,J ,CT
    DOUBLE PRECISION :: BROWNF, F, S2DTR(CHAINP%NPT), S2DTU(CHAINP%NPT), S2DTO, RHO, PHI,ST,CP,SP
    DOUBLE PRECISION :: XVEC(3), YVEC(3), POS0(3,CHAINP%NPT),UVEC0(3,CHAINP%NPT)
    DOUBLE PRECISION :: BROWNF1, BROWNF2, TORQUE(3), COM(3)

    PRINT*, 'LANGEVINSTEP SUBROUTINE HAS NOT BEEN TESTED IN A WHILE AND MAY NEED DEBUGGING'
    stop 1

    S2DTR = SQRT(2*KT/CHAINP%FRICTR(1:CHAINP%NPT)*DELT)
    S2DTU = SQRT(2*KT/CHAINP%FRICTU(1:CHAINP%NPT)*DELT)

    !PRINT*, 'TESTX1:', S6DTR, S4DTU

    ! Get the deterministic forces
    CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.true.)
    !CALL GETCHAINFORCEGAUSS(CHAINP,POSFORCE,ENERGY)
    ENERGY = ENERGY*KT; POSFORCE = POSFORCE*KT; UFORCE = UFORCE*KT

    POS0 = CHAINP%POS; UVEC0 = CHAINP%UVEC
    BROWNU = 0D0; BROWNPOS = 0D0
    IF (PRESENT(STRESSXY)) THEN
       STRESSXY = 0D0       
    ENDIF
    
    DO B = 1,CHAINP%NPT
       ! translational brownian force

       BROWNPOS(1,B) = RNORM()*S2DTR(B)
       BROWNPOS(2,B) = RNORM()*S2DTR(B)
       BROWNPOS(3,B) = RNORM()*S2DTR(B)

       ! BROWNF = RNORM()*S6DTR ! magnitude
       ! RHO = GRND()*2-1; ST = SQRT(1-RHO*RHO)
       ! PHI = GRND()*2*PI
       ! BROWNPOS(1,B) = BROWNF*ST*COS(PHI)
       ! BROWNPOS(2,B) = BROWNF*ST*SIN(PHI)
       ! BROWNPOS(3,B) = BROWNF*RHO

       IF (CHAINP%SHEARABLE) THEN
          ! rotational brownian force
          BROWNF1 = RNORM()*S2DTU(B)
          BROWNF2 = RNORM()*S2DTU(B)
          IF (CHAINP%UVEC(2,B).EQ.0D0.AND.CHAINP%UVEC(3,B).EQ.0D0) THEN
             XVEC = (/0D0,1D0,0D0/)
             YVEC = (/0D0,0D0,1D0/)
          ELSE
             CALL CROSS_PRODUCT(CHAINP%UVEC(:,B),(/1D0,0D0,0D0/),XVEC)
             CALL CROSS_PRODUCT(CHAINP%UVEC(:,B),XVEC,YVEC)
          ENDIF
          !PHI = GRND()*2*PI; CP = COS(PHI); SP = SIN(PHI)
          !BROWNU(:,B) = BROWNF*(CP*XVEC+SP*YVEC)       
          BROWNU(:,B) = BROWNF1*XVEC + BROWNF2*YVEC

          ! get just the perpendicular component of the U forces
          UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)
          CALL CROSS_PRODUCT(UVEC0(:,B),UFORCE(:,B),TORQUE)
          CALL CROSS_PRODUCT(TORQUE,UVEC0(:,B),UFORCE(:,B))
       ENDIF

       ! move to midpoint positions       
       !CHAINP%POS(:,B) = POS0(:,B) + BROWNPOS(:,B)/2 + DELT/2*POSFORCE(:,B)/CHAINP%FRICTR
       !CHAINP%UVEC(:,B) = UVEC0(:,B) + BROWNU(:,B)/2 + DELT/2*UFORCE(:,B)/CHAINP%FRICTU
       !IF (.not.B.EQ.1) THEN
       ! IF (.not.B.EQ.(CHAINP%NPT+1)/2) THEN
       IF (.NOT.((FIXBEAD1.AND.B.EQ.1).or.(FIXBEADMID.AND.B.EQ.(CHAINP%NPT+1)/2))) THEN
          CHAINP%POS(:,B) = POS0(:,B) + BROWNPOS(:,B) + DELT*POSFORCE(:,B)/CHAINP%FRICTR(B)          
       ENDIF
       IF (CHAINP%SHEARABLE.and.B.LT.CHAINP%NPT) THEN
          CHAINP%UVEC(:,B) = UVEC0(:,B) + BROWNU(:,B) + DELT*UFORCE(:,B)/CHAINP%FRICTU(B)
          ! renormalize orientation vector
          CALL NORMALIZE(CHAINP%UVEC(:,B))
       ENDIF

       !SHEARXY = SHEARXY - 0.5*(CHAINP%POS(1,B)-POS0(1,B))*BROWNPOS(2,B) - CHAINP%POS0(1,B)*POSFORCE(2,B)
  !      IF (PRESENT(STRESSXY)) THEN
! !          STRESSXY = STRESSXY - (POS0(1,B)-COM(1))*BROWNPOS(2,B) - (POS0(1,B)-COM(1))*POSFORCE(2,B)
!           STRESSXY = STRESSXY  - 0.5*(CHAINP%POS(1,B)-POS0(1,B))*BROWNPOS(2,B) - POS0(1,B)*POSFORCE(2,B)
!        ENDIF
       !STRESSCURV = STRESSCURV - POS0(1,B)*BENDFORCE(2,B)
    END DO
    
    IF (PRESENT(STRESSXY)) THEN
       CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY,.true.)  
       COM(1) = SUM(CHAINP%POS(:,1))/CHAINP%NPT
       COM(2) = SUM(CHAINP%POS(:,2))/CHAINP%NPT
       COM(3) = SUM(CHAINP%POS(:,3))/CHAINP%NPT
       CT = 0
       DO I = 1,3
          DO J = 1,3
             IF (I.EQ.J) CYCLE
             CT = CT + 1
              STRESSXY(CT) = STRESSXY(CT)  - SUM(0.5*(CHAINP%POS(I,:)-POS0(I,:))*BROWNPOS(J,:)*CHAINP%FRICTR(:)/DELT &
                   & + (POS0(I,:)-COM(I))*POSFORCE(J,:))
            ! STRESSXY(CT) = STRESSXY(CT)  - SUM((POS0(I,:)-COM(I))*(POSFORCE(J,:)+brownpos(j,:)*chainp%frictr/delt))
          ENDDO
       ENDDO
    ENDIF
    ! Get the deterministic forces at the midpoint
    !CALL GETCHAINFORCEGAUSS(CHAINP,POSFORCE,ENERGY)
    ! print*, 'testx0'
    ! CALL GETCHAINFORCEINT(CHAINP,POSFORCE,UFORCE,ENERGY)
    ! ENERGY = ENERGY*KT; POSFORCE = POSFORCE*KT; UFORCE = UFORCE*KT

    ! DO B = 1,CHAINP%NPT
    !    ! get just the perpendicular component of the U forces
    !    UFORCE(:,B) = UFORCE(:,B)-DOT_PRODUCT(UFORCE(:,B),UVEC0(:,B))*UVEC0(:,B)

    !    ! move to final positions
    !    CHAINP%POS(:,B) = POS0(:,B) + BROWNPOS(:,B) + DELT*POSFORCE(:,B)/CHAINP%FRICTR
    !    CHAINP%UVEC(:,B) = UVEC0(:,B) + BROWNU(:,B) + DELT*UFORCE(:,B)/CHAINP%FRICTU
    !    ! renormalize orientation vector
    !    CALL NORMALIZE(CHAINP%UVEC(:,B))
    ! ENDDO

  END SUBROUTINE LANGEVINSTEP

  SUBROUTINE GETCHAINFORCEINT(CHAINP,RFORCE,UFORCE,ENERGY,GETFORCE)
    USE CHAINUTIL, ONLY : CHAIN
    USE KEYS, ONLY : LOGRTERM, GAUSSIANCHAIN
    ! get intrinsic force on each bead of the chain, as well as intrinsic energy
    ! if GETFORCE=false only the energy is calculated
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: RFORCE(3,CHAINP%NPT),UFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    LOGICAL, INTENT(IN) :: GETFORCE
    INTEGER :: S,I,J
    DOUBLE PRECISION :: EXTRAFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION :: POS1(3), POS2(3), UVEC1(3), UVEC2(3), DU(3), TMP
    DOUBLE PRECISION :: DR(3), DRU, COEFF, NU
    DOUBLE PRECISION :: DCRD(6), DRPERP(3,6), RPERP(3)
    DOUBLE PRECISION :: EBEND, ESTRETCH, ESHEAR, ECOUPLE, ECONST    
    DOUBLE PRECISION :: LP,GAM,EPAR,EPERP,EC, NU1, NU2, CONSTMOD
    DOUBLE PRECISION :: DVEC1(3), DVEC2(3), NDV1, NDV2, NDV12
    DOUBLE PRECISION :: LOGMAT(2),TMPMAT(2),BMINUS
    DOUBLE PRECISION :: LOGMATFORCE(3,CHAINP%NPT,2), tmp2, DTMPNDV(3,3), LOGMATFORCEPREV(3,CHAINP%NPT,2)

    IF (GAUSSIANCHAIN) THEN
       ! gaussian force components only
       UFORCE = 0D0
       CALL GETCHAINFORCEGAUSS(CHAINP,RFORCE,ENERGY)
       RETURN
    ELSEIF (.NOT.CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
       ! force for a bead-rod chain with inextensible segments,
       ! including constraint tensions and the metric pseudopoential
       PRINT*, 'CANNOT USE GETCHAINFORCEINT WITH BEAD-ROD MODEL'
       STOP 1
       !UFORCE = 0D0
       !CALL GETBEADRODFORCE(CHAINP,RFORCE,ENERGY)
       !RETURN
    ENDIF

    EBEND = 0D0; ESTRETCH = 0D0; ESHEAR = 0D0; ECOUPLE = 0D0; ECONST = 0D0
    RFORCE = 0D0; UFORCE = 0D0; EXTRAFORCE = 0D0
    LOGMATFORCE = 0D0    

    CONSTMOD = CHAINP%CONSTMOD
    LOGMAT = (/1D0,1D0/)
    DO S = 1,CHAINP%NPT-1 ! Go through each segment
       LP = CHAINP%LP(S)
       EPAR = CHAINP%EPAR(S)
       EPERP = CHAINP%EPERP(S)
       EC = CHAINP%EC(S)
       GAM = CHAINP%GAM(S)

       POS1 = CHAINP%POS(:,S)
       POS2 = CHAINP%POS(:,S+1)
       UVEC1 = CHAINP%UVEC(:,S)
       UVEC2 = CHAINP%UVEC(:,S+1)

       ! ! stretch energy
       COEFF = EPAR/(2*CHAINP%DELS(S));       
       DR = POS2-POS1
       IF (CHAINP%SHEARABLE) THEN
          DRU = DOT_PRODUCT(DR,UVEC1)
          TMP = DRU - CHAINP%DELS(S)*GAM
          IF (LOGRTERM) THEN
             ESTRETCH = ESTRETCH + COEFF*TMP**2 - 2*LOG(DRu)
          ELSE
             ESTRETCH = ESTRETCH + COEFF*TMP**2
          ENDIF

          IF (GETFORCE) THEN
             RFORCE(:,S+1) = RFORCE(:,S+1)+2*COEFF*TMP*UVEC1
             RFORCE(:,S) = RFORCE(:,S)-2*COEFF*TMP*UVEC1
             UFORCE(:,S) = UFORCE(:,S)+2*COEFF*TMP*DR
          ENDIF
       ELSE
          DRU = SQRT(DOT_PRODUCT(DR,DR))
          TMP = DRU - CHAINP%DELS(S)*GAM
          !print*, 'testx1:', epar/chainp%dels(s), chainp%dels(s)*gam 
          !          PRINT*, 'TESTX1:', DRU
          IF (LOGRTERM) THEN
             ESTRETCH = ESTRETCH + COEFF*TMP**2 + 2*LOG(DRU) 
             IF (VERBOSE) PRINT*, 'TESTX5:', DRU, LOG(DRU)

             ! update additional matrix term
             TMP2 = EXP(-2*DRU*EPAR*GAM)
             LOGMAT(2) = LOGMAT(2)*TMP2   
             IF (GETFORCE) THEN
                LOGMATFORCE(:,:,2) = LOGMATFORCE(:,:,2)*TMP2
                LOGMATFORCE(:,S+1,2) = LOGMATFORCE(:,S+1,2) - LOGMAT(2)*2*EPAR*GAM*DR/DRU
                LOGMATFORCE(:,S,2) = LOGMATFORCE(:,S,2) + LOGMAT(2)*2*EPAR*GAM*DR/DRU             

                !print*, 'testx1:', dot_product(dr,dr), 2*log(dru)
                RFORCE(:,S+1) = RFORCE(:,S+1) + 2*COEFF*TMP/DRU*DR + 2/DRU**2*DR 
                RFORCE(:,S) = RFORCE(:,S) - 2*COEFF*TMP/DRU*DR - 2/DRU**2*DR 
                !EXTRAFORCE(:,S+1) = EXTRAFORCE(:,S+1)  + 2/DRU**2*DR
                !EXTRAFORCE(:,S) = EXTRAFORCE(:,S)  - 2/DRU**2*DR
             ENDIF
          ELSE
             ESTRETCH = ESTRETCH + COEFF*TMP**2
             IF (GETFORCE) THEN
                RFORCE(:,S+1) = RFORCE(:,S+1) + 2*COEFF*TMP/DRU*DR
                RFORCE(:,S) = RFORCE(:,S) - 2*COEFF*TMP/DRU*DR
             ENDIF
          ENDIF

       ENDIF

       ! bend energy
       IF (S.LT.CHAINP%NPT) THEN
          IF (CHAINP%SHEARABLE) THEN
              DU = UVEC2-UVEC1
             ! !print*, 'testx1', s, du
             ! COEFF = LP/(2*CHAINP%DELS(S))
             ! EBEND = EBEND + COEFF*DOT_PRODUCT(DU,DU)

             ! IF (GETFORCE) THEN
             !    UFORCE(:,S+1) = UFORCE(:,S+1) +2*COEFF*DU
             !    UFORCE(:,S) = UFORCE(:,S) - 2*COEFF*DU               
             ! ENDIF
             !print*, 'testx2:', s, uforce(:,s)
             EBEND = EBEND + LP/(CHAINP%DELS(S))*(1-DOT_PRODUCT(UVEC1,UVEC2))
             IF (GETFORCE) THEN
               UFORCE(:,S+1) = UFORCE(:,S+1) -LP/(CHAINP%DELS(S))*UVEC1
               UFORCE(:,S) = UFORCE(:,S) - LP/(CHAINP%DELS(S)) *UVEC2
             ENDIF
          ELSE     
             IF (S.EQ.CHAINP%NPT-1) CYCLE
             DVEC1 = CHAINP%POS(:,S+1)-CHAINP%POS(:,S)
             DVEC2 = CHAINP%POS(:,S+2)-CHAINP%POS(:,S+1)
             NDV1 = SQRT(DOT_PRODUCT(DVEC1,DVEC1))
             NDV2 = SQRT(DOT_PRODUCT(DVEC2,DVEC2))
             NDV12 = NDV1*NDV2

             COEFF = LP/CHAINP%DELS(S)
             TMP = DOT_PRODUCT(DVEC1,DVEC2)
             !print*, 'testx1:', coeff, tmp/ndv12, COEFF*(1-TMP/NDV12)
             EBEND = EBEND + COEFF*(1-TMP/NDV12)

             IF (LOGRTERM) THEN
                BMINUS = EXP(-2*COEFF*TMP/NDV12)
                TMPMAT(1) = LOGMAT(1)+LOGMAT(2)*BMINUS
                TMPMAT(2) = LOGMAT(1)*BMINUS + LOGMAT(2)                

                IF (GETFORCE) THEN
                   ! derivatives of TMP/NDV12 wrt S+2,S,S+1 coords
                   DTMPNDV(:,3) = DVEC1/NDV12 - TMP*DVEC2/(NDV2**3*NDV1)
                   DTMPNDV(:,1) = - DVEC2/NDV12 + TMP*DVEC1/(NDV1**3*NDV2)
                   DTMPNDV(:,2) = - DVEC1/NDV12 + TMP*DVEC2/(NDV2**3*NDV1) &
                        & + DVEC2/NDV12 - TMP*DVEC1/(NDV1**3*NDV2)

                   LOGMATFORCEPREV = LOGMATFORCE
                   LOGMATFORCE(:,:,1) = LOGMATFORCE(:,:,1) + LOGMATFORCEPREV(:,:,2)*BMINUS
                   LOGMATFORCE(:,:,2) = LOGMATFORCE(:,:,2) + LOGMATFORCEPREV(:,:,1)*BMINUS
                   DO I = 0,2
                      LOGMATFORCE(:,S+I,1) = LOGMATFORCE(:,S+I,1) &
                           & - LOGMAT(2)*2*COEFF*BMINUS*DTMPNDV(:,I+1)
                      LOGMATFORCE(:,S+I,2) = LOGMATFORCE(:,S+I,2) &
                           & - LOGMAT(1)*2*COEFF*BMINUS*DTMPNDV(:,I+1)
                   ENDDO
                ENDIF
                LOGMAT = TMPMAT
             ENDIF
             IF (GETFORCE) THEN
                RFORCE(:,S+2) = RFORCE(:,S+2) - COEFF*(DVEC1/NDV12 - TMP/NDV1/NDV2**3*DVEC2)
                RFORCE(:,S+1) = RFORCE(:,S+1) - COEFF*(-DVEC1/NDV12 + TMP/NDV1/NDV2**3*DVEC2 + DVEC2/NDV12 - TMP/NDV2/NDV1**3*DVEC1)
                RFORCE(:,S) = RFORCE(:,S) - COEFF*(-DVEC2/NDV12 + TMP/NDV2/NDV1**3*DVEC1)
             ENDIF
          ENDIF
       ENDIF


       IF (CHAINP%SHEARABLE) THEN
          ! ! shear energy
          COEFF = EPERP/(2*CHAINP%DELS(S));

          RPERP = DR - DRU*UVEC1
          ESHEAR = ESHEAR + COEFF*DOT_PRODUCT(RPERP,RPERP)

          IF (GETFORCE) THEN
             ! deriv of R_perp wrt coordinates dr, U
             DRPERP = 0D0
             DO I = 1,3
                DRPERP(I,I) = 1D0
                DRPERP(I,I+3) = -DRU
                DO J = 1,3
                   DRPERP(I,J) = DRPERP(I,J) - UVEC1(I)*UVEC1(J)
                   DRPERP(I,J+3) = DRPERP(I,J+3) -DR(J)*UVEC1(I)     
                ENDDO
             ENDDO

             DO J = 1,6
                DCRD(J) = DOT_PRODUCT(2*COEFF*RPERP,DRPERP(:,J))          
             ENDDO

             RFORCE(:,S+1) = RFORCE(:,S+1)+ DCRD(1:3)
             RFORCE(:,S) = RFORCE(:,S)-DCRD(1:3)
             UFORCE(:,S) = UFORCE(:,S)+ DCRD(4:6)
          ENDIF

          ! bend-shear coupling
          IF (S.LT.CHAINP%NPT-1) THEN
             COEFF = EC/CHAINP%DELS(S)
             ECOUPLE = ECOUPLE + COEFF*DOT_PRODUCT(DU,RPERP)

             IF (GETFORCE) THEN
                DO J = 1,6
                   DCRD(J) = COEFF*DOT_PRODUCT(DU,DRPERP(:,J))
                ENDDO

                RFORCE(:,S+1) = RFORCE(:,S+1)+DCRD(1:3)
                RFORCE(:,S) = RFORCE(:,S)-DCRD(1:3)
                UFORCE(:,S)=UFORCE(:,S) + DCRD(4:6)-COEFF*RPERP
                UFORCE(:,S+1) = UFORCE(:,S+1) +COEFF*RPERP
             ENDIF
          ENDIF
       ENDIF

    END DO

    !ENERGY = ECONST; FORCE = FCONST;

    IF (CHAINP%SHEARABLE) THEN
       ENERGY = EBEND + ESTRETCH+ESHEAR+ECOUPLE !+ECONST
    ELSE
       ENERGY = EBEND + ESTRETCH
       IF (LOGRTERM) THEN
          TMP = LOGMAT(1)+LOGMAT(2)
          ENERGY = ENERGY - LOG(TMP)       
          IF (GETFORCE) THEN
             !EXTRAFORCE = EXTRAFORCE - (LOGMATFORCE(:,:,1)+LOGMATFORCE(:,:,2))/TMP
             RFORCE = RFORCE - (LOGMATFORCE(:,:,1)+LOGMATFORCE(:,:,2))/TMP
          ENDIF
       ENDIF
       !       ENERGY = ESTRETCH
    ENDIF

    IF (GETFORCE) THEN
       RFORCE = - RFORCE; 
       IF (CHAINP%SHEARABLE) THEN
          UFORCE = -UFORCE       
       ELSE
          ! IF (LOGRTERM) THEN
          !    UFORCE = -EXTRAFORCE
          ! ELSE
             UFORCE = 0D0
          !ENDIF
       ENDIF
    ENDIF
    !ENERGY = LOGMAT(1)+LOGMAT(2)
    !RFORCE = LOGMATFORCE(:,:,1)+LOGMATFORCE(:,:,2)
  END SUBROUTINE GETCHAINFORCEINT

  SUBROUTINE GETCHAINFORCEGAUSS(CHAINP,RFORCE,ENERGY)
    ! get chain energy and forces using just a gaussian model with EPAR as the moduluss
    USE CHAINUTIL, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: RFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    INTEGER :: S,I,J
    DOUBLE PRECISION :: POS1(3), POS2(3), UVEC1(3), UVEC2(3), DU(3), TMP
    DOUBLE PRECISION :: DR(3), DRU, COEFF, NU
    DOUBLE PRECISION :: DCRD(6), DRPERP(3,6), RPERP(3)
    DOUBLE PRECISION :: EBEND, ESTRETCH, ESHEAR, ECOUPLE, ECONST    
    DOUBLE PRECISION :: LP,GAM,EPAR,EPERP,EC, NDR
    
    ENERGY = 0D0
    RFORCE = 0D0; 

    DO S = 1,CHAINP%NPT-1 ! Go through each segment
       EPAR = CHAINP%EPAR(S)

       POS1 = CHAINP%POS(:,S)
       POS2 = CHAINP%POS(:,S+1)

       DR = POS2-POS1
       NDR = DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3)

       ENERGY = ENERGY + EPAR/2/CHAINP%LS(1)*NDR
       RFORCE(:,S+1) = RFORCE(:,S+1)+EPAR/CHAINP%LS(1)*DR
       RFORCE(:,S) = RFORCE(:,S)-EPAR/CHAINP%LS(1)*DR
    ENDDO

    RFORCE = -RFORCE
  END SUBROUTINE GETCHAINFORCEGAUSS

  SUBROUTINE GETBEADRODFORCE(CHAINP,BROWNFORCE,CRELAX,RFORCE,ENERGY,PROJBROWN,TENSBRFORCE)
    ! get chain energy and forces using a bead-rod model with constant
    ! chain segment lengths (for now assumed the same throughout)
    ! includes the constraint tensions and the metric pseudo-potential force
    ! BROWNFORCE contains precalculated brownian forces for use in getting tension forces 
    ! CRELAX is the coefficient for the additional force to keep bond lengths numerically constrained
    ! RFORCE, on exit, contains the potential, tension, and pseudo-potential forces
    ! if PROJBROWN is true, alter BROWNFORCE so that on exit the components
    ! that would alter the segment lengths are projected out
    ! TENSBRFORCE is the component of the tension force that arises from the brownian force (after projection if projection is being done)

    USE CHAINUTIL, ONLY : CHAIN
    USE KEYS, ONLY : USEPSEUDOFORCE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(INout) :: BROWNFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION, INTENT(IN) :: CRELAX
    LOGICAL, INTENT(IN) :: PROJBROWN
    DOUBLE PRECISION, INTENT(OUT) :: RFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: TENSBRFORCE(3,CHAINP%NPT)
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY
    INTEGER :: S,I,J, NPT, NMAT, INFO
    DOUBLE PRECISION :: TDET(0:CHAINP%NPT-1),BDET(CHAINP%NPT)
    double precision, DIMENSION(CHAINP%NPT-1) :: DIAG, LOWDIAG, UPDIAG, UPDIAGSAVE, &
         & GINVOFF,ndr!GINVDIAG
    DOUBLE PRECISION :: DR(3,CHAINP%NPT-1), LAM(CHAINP%NPT-1),LAMBR(CHAINP%NPT-1)
    DOUBLE PRECISION :: DETER, KAPPA, KEFF, DERVP(3), DERVM(3),DERVC(3),NDR2   
    DOUBLE PRECISION :: TOTFORCE(3,CHAINP%NPT), L0, RFORCESAVE(3,CHAINP%NPT)
    DOUBLE PRECISION :: DRTMP1(3), DRTMP2(3), NDRTMP1, NDRTMP2

    ENERGY = 0D0; RFORCE = 0D0
    IF (PRESENT(TENSBRFORCE)) THEN
       TENSBRFORCE = 0D0
    ENDIF

    ! set up the tridiagonal matrix for solving constraint forces
    NPT= CHAINP%NPT
    NMAT = NPT - 1 ! size of matrix = number of segmetns
    DO I = 1,NMAT ! cycle through segments
       ! segment vector
       DR(:,I) = CHAINP%POS(:,I+1)-CHAINP%POS(:,I)
       NDR2 = DOT_PRODUCT(DR(:,I),DR(:,I)) 
       NDR(I) = SQRT(NDR2)
       !PRINT*, 'TESTX3:', I, NDR(I)
       DIAG(I) = 2*NDR2
       !print*, 'testx3:', i, diag(i)
       IF (I.GT.1) THEN
          LOWDIAG(I) = -DOT_PRODUCT(DR(:,I-1),DR(:,I))
          !print*, 'testx4:', i, lowdiag(i)
       ENDIF
    ENDDO
    UPDIAG(1:NMAT-1)=LOWDIAG(2:NMAT)
    LOWDIAG(1) = 0D0; UPDIAG(CHAINP%NPT-1) = 0D0

    IF (USEPSEUDOFORCE) THEN
       ! calculate the pseudopotential forces; algorithm from Pasquali & Morse 2002
       ! first get, top and bottom matrix determinants
       TDET(0) = 1 
       TDET(1) = DIAG(1)
       BDET(NPT) = 1
       BDET(NPT-1) = DIAG(NPT-1)
       DO I = 1,NPT-2
          TDET(I+1) = DIAG(I+1)*TDET(I) - LOWDIAG(I+1)**2*TDET(I-1)
          BDET(NPT-I-1) = DIAG(NPT-I-1)*BDET(NPT-I) - LOWDIAG(NPT-I)**2*BDET(NPT-I+1)      
       ENDDO
       ! overall matrix determinant
       DETER = TDET(NPT-1)

       !PRINT*, 'TESTX2:', dr(:,2)
       !STOP 1
       ! get the diagonal and off-diagonal inverse G matrix elements
       ! inv(G)_i-1,i
       DO I = 1,NMAT
          IF (I.GT.1) THEN
             GINVOFF(I) = -LOWDIAG(I)*TDET(I-2)*BDET(I+1)/DETER
          ENDIF
          !GINVDIAG(I) = TDET(I-1)*BDET(I+1)/DETER
          !print*, 'testx5:', i, ginv(i)
       ENDDO
    ENDIF

    ! and the bend potential and pseudopotential energy and forces
    KAPPA = CHAINP%LP(1)/CHAINP%DELS(1)

    !ENERGY = 0.5*LOG(DETER)       
    !RFORCE(:,2) = RFORCE(:,2)-2*GINVDIAG(1)*DR(:,1)
    !RFORCE(:,1) = RFORCE(:,1)+2*GINVDIAG(1)*DR(:,1)
    DO I = 2,CHAINP%NPT-1 ! cycle over joints
       ! potential and pseudopotential energies
       ! DRTMP1 = CHAINP%POS(:,I)-CHAINP%POS(:,I-1);
       ! DRTMP2 = CHAINP%POS(:,I+1)-CHAINP%POS(:,I);
       ! NDRTMP1 = SQRT(SUM(DRTMP1**2))
       ! NDRTMP2 = SQRT(SUM(DRTMP2**2))
       ! !PRINT*, 'TESTX2RHO: ', I, SUM(DRTMP1*DRTMP2)/NDRTMP1/NDRTMP2
       ! ENERGY = ENERGY + KAPPA*(1-SUM(DRTMP1*DRTMP2)/NDRTMP1/NDRTMP2)
       ENERGY = ENERGY + KAPPA*(1+LOWDIAG(I)/(NDR(I-1)*NDR(I))) 

       ! upper bead derivatives
       DERVP = DR(:,I-1)/NDR(I-1)/NDR(I) + LOWDIAG(I)/NDR(I-1)/NDR(I)**3*DR(:,I)
       ! lower bead derivatives
       DERVM = -DR(:,I)/NDR(I-1)/NDR(I) - LOWDIAG(I)/NDR(I)/NDR(I-1)**3*DR(:,I-1)
       ! center bead derivatives
       DERVC = -(DERVP+DERVM)


       ! effective rigidity
       IF (USEPSEUDOFORCE) THEN
          RFORCE(:,I+1) = RFORCE(:,I+1) +KAPPA*DERVP + GINVOFF(I)*DR(:,I-1)
          RFORCE(:,I) = RFORCE(:,I)+KAPPA*DERVC +GINVOFF(I)*(-DR(:,I-1)+DR(:,I))
          RFORCE(:,I-1) = RFORCE(:,I-1)+KAPPA*DERVM -GINVOFF(I)*DR(:,I)
          ! if (present(tensbrforce)) then
          !    tensbrFORCE(:,I+1) = tensbrFORCE(:,I+1) +KAPPA*DERVP 
          !    tensbrFORCE(:,I) = tensbrFORCE(:,I)+KAPPA*DERVC 
          !    tensbrFORCE(:,I-1) = tensbrFORCE(:,I-1)+KAPPA*DERVM 
          ! endif
       ELSE
          RFORCE(:,I+1) = RFORCE(:,I+1) +KAPPA*DERVP
          RFORCE(:,I) = RFORCE(:,I)+KAPPA*DERVC 
          RFORCE(:,I-1) = RFORCE(:,I-1)+KAPPA*DERVM 
       ENDIF


       !RFORCE(:,I+1) = RFORCE(:,I+1) + GINVOFF(I)*DR(:,I-1)!-2*GINVDIAG(I)*DR(:,I)
       !RFORCE(:,I) = RFORCE(:,I)+GINVOFF(I)*(-DR(:,I-1)+DR(:,I))!+2*GINVDIAG(I)*DR(:,I)
       !RFORCE(:,I-1) = RFORCE(:,I-1)-GINVOFF(I)*DR(:,I)
    ENDDO
  
    RFORCESAVE = RFORCE
    
    ! IF (PRESENT(TENSBRFORCE)) THEN
    !    TENSBRFORCE = RFORCE
    ! ENDIF

    ! find the tension forces arising specifically from the brownian forces    
    IF (PRESENT(TENSBRFORCE).OR.PROJBROWN) THEN
       DO I = 1,NMAT       
          ! keep track of the tensions arising from brownian forces separately
          LAMBR(I) = DOT_PRODUCT(DR(:,I),brownFORCE(:,I+1)-brownFORCE(:,I))
       ENDDO

       ! PRINT*, 'TESTX0:'
       !  DO I = 1,NPT-1
       !     PRINT*, LOWDIAG(I),DIAG(I),UPDIAG(I)
       !  ENDDO
       UPDIAGSAVE = UPDIAG
       CALL DGTSL(NMAT,LOWDIAG,DIAG,UPDIAG,LAMBR,INFO)
       UPDIAG = UPDIAGSAVE
       IF (INFO.NE.0) THEN
          PRINT*, 'ERROR IN SOLVING FOR BROWNIAN TENSION FORCES', info          
          STOP 1
       ENDIF

       IF (PROJBROWN) THEN
          ! project out the tension parts of the brownian force
          BROWNFORCE(:,1) = BROWNFORCE(:,1) + LAMBR(1)*DR(:,1)
          DO I = 2,NPT-1
             BROWNFORCE(:,I) = BROWNFORCE(:,I) + LAMBR(I)*DR(:,I)-LAMBR(I-1)*DR(:,I-1)
          ENDDO
          BROWNFORCE(:,NPT) = BROWNFORCE(:,NPT) - LAMBR(NPT-1)*DR(:,NPT-1)
       ENDIF
    ENDIF
   
    ! find the overall tension forces
    ! other forces include the bend, metric, brownian forces, and some 
    ! additional force to maintain bond lengths numerically
    L0 = CHAINP%DELS(1)*CHAINP%GAM(1)
    !    print*, 'TESTX1:', CRELAX, L0, ndr(1)
    TOTFORCE = RFORCE + BROWNFORCE
    DO I = 1,NMAT
       LAM(I) = DOT_PRODUCT(DR(:,I),TOTFORCE(:,I+1)-TOTFORCE(:,I))+ CRELAX*(NDR(I)-L0)
    ENDDO

    UPDIAGSAVE = UPDIAG
    CALL DGTSL(NMAT,LOWDIAG,DIAG,UPDIAG,LAM,INFO)
    UPDIAG = UPDIAGSAVE
    IF (INFO.NE.0) THEN
       PRINT*, 'ERROR IN SOLVING FOR TENSION FORCES'
       STOP 1
    ENDIF

    RFORCE(:,1) = RFORCE(:,1) + LAM(1)*DR(:,1)
    DO I = 2,NPT-1
       RFORCE(:,I) = RFORCE(:,I) + LAM(I)*DR(:,I)-LAM(I-1)*DR(:,I-1)
    ENDDO
    RFORCE(:,NPT) = RFORCE(:,NPT) - LAM(NPT-1)*DR(:,NPT-1)

    ! PRINT*, 'TESTX2:', DOT_PRODUCT(RFORCE(:,2)+BROWNFORCE(:,2)-RFORCE(:,1)-BROWNFORCE(:,1),&
    !      & DR(:,1))/SQRT(SUM(DR(:,1)**2))    

    ! tension components from just the brownian forces
    IF (PRESENT(TENSBRFORCE)) THEN
       TENSBRFORCE = 0D0
       IF (.NOT.PROJBROWN) THEN
          TENSBRFORCE(:,1) = TENSBRFORCE(:,1) + LAMBR(1)*DR(:,1)
          DO I = 2,NPT-1
             TENSBRFORCE(:,I) = TENSBRFORCE(:,I) + LAMBR(I)*DR(:,I)-LAMBR(I-1)*DR(:,I-1)
          ENDDO
          TENSBRFORCE(:,NPT) = TENSBRFORCE(:,NPT) - LAMBR(NPT-1)*DR(:,NPT-1)
       ENDIF
       !PRINT*, 'TESTX3:', RFORCESAVE(2,1), RFORCE(2,1)-RFORCESAVE(2,1), BROWNFORCE(2,1), TENSBRFORCE(2,1)

       ! --------- all the tension forces ---------
      ! TENSBRFORCE = RFORCE-RFORCESAVE
    ENDIF

    
  END SUBROUTINE GETBEADRODFORCE

END MODULE BROWNDYN
