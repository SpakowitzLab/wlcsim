MODULE REDISC
  ! utilities for rediscretizing the chain on the fly
  USE CHAINUTIL, ONLY : CHAIN
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: PARAMDATA(:,:)
  INTEGER :: NPARAMDATA 
  LOGICAL :: PARAMDATASET = .FALSE.


CONTAINS
   ! for two segments, criterion to add a bead to both:
    ! imagine all the possible position of a center bead on one of them (intersection of two balls of radius LS/2). 
    !Circumscribe this set by a cylinder for ease of computation 
    !(height = L-d, radius = sqrt(L^2-d^2)/2, where d is distance btwn points; 
    ! cylinder axis is vector btwn points
    ! if these potential cylinders for the two segments overlap, then need to 
    ! add center point for each of the two segments
    ! for removing points, check every even bead for possible removal in the same way
    ! in a given run of this function, at most one bead is added to each segment
    ! and at most every other bead is removed
    ! removal happens before addition
  

    SUBROUTINE REDISCADD(CHAINP,DELE)
    ! rediscretize chain, adding beads where potential steric conflicts exist
    ! set of possible positions is circumscribed by cylinder for simplicity of calculation
   
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: S1,S2, NADD, CH
    LOGICAL :: HASCONFLICT(CHAINP%NPT), CYLSET(CHAINP%NPT), INTERSECT
    DOUBLE PRECISION :: LS1, LS2, LS3, LS4, PT1(3), PT2(3), PT3(3), PT4(3)
    DOUBLE PRECISION :: DIFF(3), NDIFF, TMPDELE
    INTEGER :: ADDSEG(CHAINP%NPT)
    DOUBLE PRECISION :: HCYL(CHAINP%NPT), RCYL(CHAINP%NPT), AXCYL(3,CHAINP%NPT),CENTCYL(3,CHAINP%NPT)
    DOUBLE PRECISION :: CUMLEN1, CUMLEN2


    ! for all pairs of segments, check for possible conflict
    
    NADD = 0
    HASCONFLICT = .FALSE.
    CYLSET = .FALSE.

    CUMLEN1 = 0D0
    DO S1 = 1,CHAINP%NPT-1
       PT1 = CHAINP%POS(:,S1); PT2 = CHAINP%POS(:,S1+1)
       LS1 = CHAINP%LS(S1)/2+CHAINP%STERRAD; LS2 = LS1
       
       CUMLEN1 = CUMLEN1 + CHAINP%LS(S1)

       
       CUMLEN2 = 0D0
       DO S2 = 1,CHAINP%NPT
          IF (S2.GT.1) CUMLEN2 = CUMLEN2 + CHAINP%LS(S2-1)
          IF (S2.EQ.S1.OR.S2.EQ.S1+1) CYCLE
          IF (ABS((CUMLEN1-CHAINP%LS(S1)/2)-(CUMLEN2)).LE.2*CHAINP%STERRAD) CYCLE

          PT3 = CHAINP%POS(:,S2); 
          LS3 = CHAINP%STERRAD

          CALL INTERSECT3SPHERE(PT1,PT2,PT3,LS1,LS2,LS3,INTERSECT)

          ! IF (S1.EQ.13.AND.S2.EQ.5) THEN
          !    PRINT*, 'TESTX0:', PT1, PT2
          !    PRINT*, 'TESTX1:', INTERSECT,CHAINP%LS(S1)/2, SQRT(DOT_PRODUCT(PT2-PT1,PT2-PT1))
          ! ENDIF

          IF (INTERSECT) THEN
             !PRINT*, 'TESTX ADDBEAD CONFLICT:', S1,S2

             IF (CHAINP%LS(S1)/2.GT.CHAINP%MINSEGLEN) THEN
                NADD = NADD + 1
                ADDSEG(NADD) = S1
                HASCONFLICT(S1) = .TRUE.
             ENDIF
             
             EXIT
          ENDIF

       ENDDO
       
    ENDDO

    !PRINT*, 'TESTX3:', ADDSEG(1:NADD)

    ! Add beads to all selected segments
    DELE = 0
    DO S1 = 1,NADD     
       !PRINT*, 'TESTX4:', S1, ADDSEG(S1), nadd
       CALL ADDBEAD(CHAINP,ADDSEG(S1),TMPDELE) 
       DO S2 = S1+1,NADD
          IF (ADDSEG(S2).GT.ADDSEG(S1)) ADDSEG(S2) = ADDSEG(S2)+1
       ENDDO
!       ADDSEG(S1+1:NADD) = ADDSEG(S1+1:NADD)+1
       DELE = DELE + TMPDELE
    ENDDO

  END SUBROUTINE REDISCADD
  
  SUBROUTINE REDISCADDOLD(CHAINP,DELE)
    ! rediscretize chain, adding beads where potential steric conflicts exist
    ! set of possible positions is circumscribed by cylinder for simplicity of calculation
   
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: S1,S2, NADD, CH
    LOGICAL :: HASCONFLICT(CHAINP%NPT), CYLSET(CHAINP%NPT), INTERSECT
    DOUBLE PRECISION :: LS1, LS2, LS3, LS4, PT1(3), PT2(3), PT3(3), PT4(3)
    DOUBLE PRECISION :: DIFF(3), NDIFF, TMPDELE
    INTEGER :: ADDSEG(CHAINP%NPT)
    DOUBLE PRECISION :: HCYL(CHAINP%NPT), RCYL(CHAINP%NPT), AXCYL(3,CHAINP%NPT),CENTCYL(3,CHAINP%NPT)
    DOUBLE PRECISION :: CUMLEN1, CUMLEN2


    ! for all pairs of segments, check for possible conflict
    
    NADD = 0
    HASCONFLICT = .FALSE.
    CYLSET = .FALSE.

    CUMLEN1 = 0D0
    DO S1 = 1,CHAINP%NPT-1
       PT1 = CHAINP%POS(:,S1); PT2 = CHAINP%POS(:,S1+1)
       LS1 = CHAINP%LS(S1)/2+CHAINP%STERRAD; LS2 = LS1
       
       CUMLEN1 = CUMLEN1 + CHAINP%LS(S1)

       IF (HASCONFLICT(S1)) CYCLE

       CUMLEN2 = 0D0
       DO S2 = 1,CHAINP%NPT-1
          CUMLEN2 = CUMLEN2 + CHAINP%LS(S2)
          IF (s1.eq.s2) CYCLE
          !IF (S2.EQ.S1) CYCLE
          IF (ABS((CUMLEN1-CHAINP%LS(S1)/2)-(CUMLEN2-CHAINP%LS(S2)/2)).LE.2*CHAINP%STERRAD) CYCLE


          PT3 = CHAINP%POS(:,S2); PT4 = CHAINP%POS(:,S2+1)
          LS3 = CHAINP%LS(S2)/2+CHAINP%STERRAD; LS4 = LS3

          CALL CHECKCONFLICT(PT1,PT2,LS1,LS2,PT3,PT4,LS3,LS4,&
               & CYLSET(S1),HCYL(S1),RCYL(S1),AXCYL(:,S1),CENTCYL(:,S1),&
               & CYLSET(S2), HCYL(S2), RCYL(S2), AXCYL(:,S2), CENTCYL(:,S2), INTERSECT)
          !print*, 'testx2:', s1, s2, intersect

          IF (INTERSECT) THEN
             !PRINT*, 'TESTX ADDBEAD CONFLICT:', S1,S2

             IF (CHAINP%LS(S1)/2.GT.2*CHAINP%STERRAD.AND.CHAINP%LS(S1)/2.GT.CHAINP%MINSEGLEN) THEN
                NADD = NADD + 1
                ADDSEG(NADD) = S1
                HASCONFLICT(S1) = .TRUE.
             ENDIF

             IF (.NOT.HASCONFLICT(S2).AND.CHAINP%LS(S2)/2.GT.2*CHAINP%STERRAD&
                  & .AND.CHAINP%LS(S2)/2.GT.CHAINP%MINSEGLEN) THEN
                NADD = NADD + 1
                ADDSEG(NADD) = S2
                HASCONFLICT(S2) = .TRUE.
             ENDIF
             EXIT
          ENDIF

       ENDDO
       
    ENDDO

    !PRINT*, 'TESTX3:', ADDSEG(1:NADD)

    ! Add beads to all selected segments
    DELE = 0
    DO S1 = 1,NADD     
       !PRINT*, 'TESTX4:', S1, ADDSEG(S1)
       CALL ADDBEAD(CHAINP,ADDSEG(S1),TMPDELE) 
       DO S2 = S1+1,NADD
          IF (ADDSEG(S2).GT.ADDSEG(S1)) ADDSEG(S2) = ADDSEG(S2)+1
       ENDDO
!       ADDSEG(S1+1:NADD) = ADDSEG(S1+1:NADD)+1
       DELE = DELE + TMPDELE
    ENDDO

  END SUBROUTINE REDISCADDOLD

  SUBROUTINE REDISCREMOVE(CHAINP,DELE)
    ! rediscretize chain by removing beads where no potential steric conflicts exist
    ! a conflict is if all possible positions of the bead could sterically clash with current position of some other bead
    ! will not remove neighboring bead
   
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: B, C, NRM, CH
    LOGICAL :: HASCONFLICT(CHAINP%NPT), CYLSET(CHAINP%NPT), INTERSECT, WILLREMOVE(CHAINP%NPT)
    DOUBLE PRECISION :: LS1, LS2, LS3, LS4, PT1(3), PT2(3), PT3(3), PT4(3)
    DOUBLE PRECISION :: DIFF(3), NDIFF, TMPDELE, CUMLEN1, CUMLEN2
    INTEGER :: RMBEADS(CHAINP%NPT)
    DOUBLE PRECISION :: HCYL(CHAINP%NPT), RCYL(CHAINP%NPT), AXCYL(3,CHAINP%NPT),CENTCYL(3,CHAINP%NPT)



    ! for all beads, check if they can be removed
    HASCONFLICT = .FALSE.
    WILLREMOVE = .FALSE.
    CYLSET = .FALSE.
    NRM = 0

    CUMLEN1 = CHAINP%LS(1)
    DO B = 2,CHAINP%NPT-1          
       CUMLEN1 = CUMLEN1 + CHAINP%LS(B)
       ! previous bead is going to be removed so keep this one
       IF (WILLREMOVE(B-1)) CYCLE

       PT1 = CHAINP%POS(:,B-1); PT2 = CHAINP%POS(:,B+1)
       LS1 = CHAINP%LS(B-1)+CHAINP%STERRAD; LS2 = CHAINP%LS(B)+CHAINP%STERRAD
             
       CUMLEN2 = 0
       DO C = 1,CHAINP%NPT
          IF (C.GT.1) CUMLEN2 = CUMLEN2 + CHAINP%LS(C-1)
          IF (abs(b-c).le.1) CYCLE
          IF (ABS(CUMLEN1-CUMLEN2).LT.2*CHAINP%STERRAD) CYCLE

          PT3 = CHAINP%POS(:,C); 
          LS3 = CHAINP%STERRAD; 
          CALL INTERSECT3SPHERE(PT1,PT2,PT3,LS1,LS2,LS3,INTERSECT)
          IF (INTERSECT) THEN
             !PRINT*, 'TESTX CONFLICT:', B,C, LS1-CHAINP%STERRAD, LS2-CHAINP%STERRAD, LS3
             HASCONFLICT(B) = .TRUE.;
             EXIT
          ENDIF

       ENDDO

       IF (.NOT.HASCONFLICT(B).AND. CHAINP%LS(B-1)+CHAINP%LS(B).LE.CHAINP%MAXSEGLEN) THEN
          ! bead has no potential steric conflicts and can be removed
          !PRINT*, 'NO CONFLICT:', B
          NRM = NRM +1
          RMBEADS(NRM) = B
          WILLREMOVE(B) = .TRUE.
       ENDIF
    ENDDO

    ! DO B = 2,CHAINP%NPT
    !    DIFF = CHAINP%POS(:,B)-CHAINP%POS(:,B-1)
    !    PRINT*, 'TESTX0:', B, SQRT(DOT_PRODUCT(DIFF,DIFF)),CHAINP%LS(B)
    ! ENDDO

    ! Remove all the beads marked for removal
    DELE = 0
    DO B = 1,NRM
       CALL REMOVEBEAD(CHAINP,RMBEADS(B),TMPDELE)       
       RMBEADS(B+1:NRM) = RMBEADS(B+1:NRM)-1
       DELE = DELE + TMPDELE
    ENDDO

   
  END SUBROUTINE REDISCREMOVE

  SUBROUTINE INTERSECT3SPHERE(P1,P2,P3,RAD1,RAD2,RAD3,INTERSECT)
    ! check whether 3 spheres centered at points P1,P2,P3
    ! with radii RAD1, RAD2, RAD3 intersect
    ! WARNING: for now this just checks that each pair of balls intersects
    ! this is not sufficient for there to be a common intersection
    USE GENUTIL, ONLY : NORMALIZE
     
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: P1(3), P2(3), P3(3), RAD1,RAD2,RAD3
    LOGICAL, INTENT(OUT) :: INTERSECT
    DOUBLE PRECISION :: DIFF12(3), DIFF13(3), DIFF23(3), EX(3), EY(3)
    DOUBLE PRECISION :: D12, D13, D23, DVAL, IVAL, JVAL, XVAL, YVAL

    ! first make sure each pair intersects
    DIFF12=P2-P1
    D12 = DOT_PRODUCT(DIFF12,DIFF12)
    IF (D12.GT.(RAD1+RAD2)**2) THEN
       INTERSECT = .FALSE.
       RETURN
    ENDIF

    DIFF23 = P3-P2
    D23 = DOT_PRODUCT(DIFF23,DIFF23)
    IF (D23.GT.(RAD2+RAD3)**2) THEN
       INTERSECT = .FALSE.
       RETURN
    ENDIF

    DIFF13 = P3-P1
    D13 = DOT_PRODUCT(DIFF13,DIFF13)
    IF (D13.GT.(RAD1+RAD3)**2) THEN
       INTERSECT = .FALSE.
       RETURN
    ENDIF

    INTERSECT = .TRUE.
    ! calculate possible intersection, using algorithm from wikipedia
    ! DVAL = SQRT(D12)
    ! EX = DIFF12/DVAL
    ! IVAL = DOT_PRODUCT(EX,DIFF13)
    ! EY = DIFF13-IVAL*EX
    ! CALL NORMALIZE(EY)
    ! JVAL = DOT_PRODUCT(EY,DIFF13)

    ! XVAL = (RAD1**2-RAD2**2+DVAL**2)/(2*DVAL)
    ! YVAL = (RAD1**2-RAD3**2+IVAL**2+JVAL**2)/(2*JVAL) - IVAL/JVAL*XVAL
    
    ! INTERSECT = (RAD1**2 - XVAL**2-YVAL**2).GT.0
    
  END SUBROUTINE INTERSECT3SPHERE

  SUBROUTINE REDISCREMOVEOLD(CHAINP,DELE)
    ! rediscretize chain by removing beads where no potential steric conflicts exist
    ! for each even bead:
    ! for each other bead,  consider all the possible positions given the neighbor beads
    ! if any steric overlap is possible, then there is a conflict 
    ! if no conflict remove the even bead.
    ! set of possible positions is circumscribed by cylinder for simplicity of calculation

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER :: B, C, NRM, CH
    LOGICAL :: HASCONFLICT(CHAINP%NPT), CYLSET(CHAINP%NPT), INTERSECT, WILLREMOVE(CHAINP%NPT)
    DOUBLE PRECISION :: LS1, LS2, LS3, LS4, PT1(3), PT2(3), PT3(3), PT4(3)
    DOUBLE PRECISION :: DIFF(3), NDIFF, TMPDELE, CUMLEN1, CUMLEN2
    INTEGER :: RMBEADS(CHAINP%NPT)
    DOUBLE PRECISION :: HCYL(CHAINP%NPT), RCYL(CHAINP%NPT), AXCYL(3,CHAINP%NPT),CENTCYL(3,CHAINP%NPT)



    ! 1) for all beads, check if they can be removed
    HASCONFLICT = .FALSE.
    WILLREMOVE = .FALSE.
    CYLSET = .FALSE.
    NRM = 0

    CUMLEN1 = CHAINP%LS(1)
    DO B = 2,CHAINP%NPT-1   
       
       CUMLEN1 = CUMLEN1 + CHAINP%LS(B)
       IF (HASCONFLICT(B)) CYCLE ! previous conflict already found
       ! previous bead is going to be removed so keep this one
       IF (WILLREMOVE(B-1)) CYCLE

       PT1 = CHAINP%POS(:,B-1); PT2 = CHAINP%POS(:,B+1)
       LS1 = CHAINP%LS(B-1)+CHAINP%STERRAD; LS2 = CHAINP%LS(B)+CHAINP%STERRAD
       
      ! PRINT*, 'TESTX1:', B, CHAINP%LS(B-1), CHAINP%STERRAD, ls1

       ! currently not checking the edge beads; should probably fix this at some point
       ! Also, many pairs of beads get checked twice the way this 
       !is currently set up; should fix this at some point too 
       ! (probably need to save entire matrix of possible conflicts)
       CUMLEN2 = CHAINP%LS(1)
       DO C = 2,CHAINP%NPT-1          
          CUMLEN2 = CUMLEN2 + CHAINP%LS(C)
          IF (ABS(B-C).LE.2) CYCLE
          IF (ABS(CUMLEN1-CUMLEN2).LT.2*CHAINP%STERRAD) CYCLE

          PT3 = CHAINP%POS(:,C-1); PT4 = CHAINP%POS(:,C+1)
          LS3 = CHAINP%LS(C-1)+CHAINP%STERRAD; LS4 = CHAINP%LS(C)+CHAINP%STERRAD
          CALL CHECKCONFLICT(PT1,PT2,LS1,LS2,PT3,PT4,LS3,LS4,&
               & CYLSET(B),HCYL(B),RCYL(B),AXCYL(:,B),CENTCYL(:,B),&
               & CYLSET(C), HCYL(C), RCYL(C), AXCYL(:,C), CENTCYL(:,C), INTERSECT)          
          IF (INTERSECT) THEN
             !PRINT*, 'TESTX CONFLICT:', B,C, LS1, LS2, LS3, LS4
             HASCONFLICT(B) = .TRUE.; HASCONFLICT(C) = .TRUE.
             EXIT
          ENDIF

       ENDDO

       IF (.NOT.HASCONFLICT(B).AND. CHAINP%LS(B-1)+CHAINP%LS(B).LE.CHAINP%MAXSEGLEN) THEN
          ! bead has no potential steric conflicts and can be removed
          !PRINT*, 'NO CONFLICT:', B
          NRM = NRM +1
          RMBEADS(NRM) = B
          WILLREMOVE(B) = .TRUE.
       ENDIF
    ENDDO

    ! Remove all the beads marked for removal
    DELE = 0
    DO B = 1,NRM
       CALL REMOVEBEAD(CHAINP,RMBEADS(B),TMPDELE)       
       RMBEADS(B+1:NRM) = RMBEADS(B+1:NRM)-1
       DELE = DELE + TMPDELE
    ENDDO

   
  END SUBROUTINE REDISCREMOVEOLD

  SUBROUTINE CHECKCONFLICT(PT1,PT2,LS1,LS2,PT3,PT4,LS3,LS4,&
       & CYLSET1,HCYL1,RCYL1,AXCYL1,CENTCYL1,&
       & CYLSET2, HCYL2, RCYL2, AXCYL2, CENTCYL2, CONFLICT)
    ! check whether two segments have a potential conflict
    ! between point A at max distance LS1 from PT1 & LS2 from PT2
    ! and point B at max distance LS3 from PT3 & LS4 from PT4
    ! CYLSET says whether cylinder info has already been set or not
    ! will be set if isn't already
    USE CYLINDERUTIL, ONLY : CYLINDERINTERSECT
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PT1(3), PT2(3), LS1, LS2, PT3(3), PT4(3), LS3, LS4
    DOUBLE PRECISION, INTENT(INOUT) :: HCYL1, RCYL1, HCYL2, RCYL2, AXCYL1(3), AXCYL2(3), CENTCYL1(3), CENTCYL2(3)
    LOGICAL, INTENT(INOUT) :: CYLSET1, CYLSET2
    LOGICAL, INTENT(OUT) :: CONFLICT
    DOUBLE PRECISION :: DIFF(3), NDIFF, RAD13, RAD23, RAD14, RAD24

    CONFLICT = .FALSE.
    RAD13 = (LS1+LS3)**2; RAD23 = (LS2+LS3)**2
    RAD14 = (LS1+LS4)**2; RAD24 = (LS2+LS4)**2
    
    ! do some preliminary checking 
    ! if any pair of the spheres can't intersect then no intersection possible
    DIFF = PT3-PT1
    NDIFF = DOT_PRODUCT(DIFF,DIFF)         
    IF (NDIFF.GT.RAD13) RETURN
    
    DIFF = PT4 - PT1
    NDIFF = DOT_PRODUCT(DIFF,DIFF)          
    IF (NDIFF.GT.RAD14) RETURN

    DIFF = PT3-PT2
    NDIFF = DOT_PRODUCT(DIFF,DIFF)          
    IF (NDIFF.GT.RAD23) RETURN

    DIFF = PT4-PT2
    NDIFF = DOT_PRODUCT(DIFF,DIFF)          
    IF (NDIFF.GT.RAD24) RETURN
          
    ! define the enveloping cylinders          
    IF (.NOT.CYLSET1) THEN
       CALL ENVELOPECYLINDER(PT1,PT2,LS1,LS2,&
            & HCYL1,RCYL1,AXCYL1,CENTCYL1)
       CYLSET1 = .TRUE.
       IF (HCYL1.LT.0) THEN
          PRINT*, 'ERROR IN CHECKCONFLICT: impossible point exists', PT1, PT2, LS1, LS2
          STOP 2
       ENDIF
    ENDIF

    IF (.NOT.CYLSET2) THEN
       CALL ENVELOPECYLINDER(PT3,PT4,LS3,LS4,&
                  & HCYL2,RCYL2,AXCYL2,CENTCYL2)
       CYLSET2 = .TRUE.
       IF (HCYL2.LT.0) THEN
          PRINT*, 'ERROR IN REDISCMV1: impossible point exists', PT3, PT4, LS3, LS4
          STOP 2
       ENDIF
    ENDIF

    ! Check if enveloping cylinders intersect
    CONFLICT = CYLINDERINTERSECT(RCYL1,RCYL2,HCYL1,HCYL2,&
         & CENTCYL1,AXCYL1,CENTCYL2,AXCYL2)
    
    
  END SUBROUTINE CHECKCONFLICT

  SUBROUTINE ENVELOPECYLINDER(PT1, PT2, RAD1, RAD2, HCYL, RCYL, AXCYL, CENTCYL)
    ! for a particular pair of points and radii
    ! define a cylinder enveloping the intersection of two filled spheres
    ! centered at those points
    ! HCYL = Height of cylinder
    ! RCYL = radius of cylinder
    ! AXCYL = normalized axis of cylinder
    ! CENTCYL = center of cylinder
    ! returns negative HCYL if the two spheres do not intersect
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PT1(3), PT2(3), RAD1, RAD2
    DOUBLE PRECISION, INTENT(OUT) :: HCYL, RCYL, AXCYL(3), CENTCYL(3)
    DOUBLE PRECISION :: DIFF(3), D, DEL1, DEL2

    DIFF = PT2-PT1
    D = SQRT(DOT_PRODUCT(DIFF,DIFF))          
    IF (D.GT.RAD1+RAD2) THEN
       HCYL = -1
       RETURN
    ENDIF

    DEL1 = MIN(RAD1,RAD2)
    DEL2 = MAX(RAD1,RAD2)
    IF (DEL2.GT.DEL1+D) THEN ! one sphere inside another
       HCYL = 2*DEL1; RCYL = DEL1
       AXCYL = DIFF/D
       IF (RAD1.LT.RAD2) THEN
          CENTCYL = PT1
       ELSE
          CENTCYL = PT2
       ENDIF
    ELSE
       HCYL = DEL1+DEL2-D
       RCYL = SQRT(DEL1**2 - ((DEL1**2-DEL2**2+D**2)/(2*D))**2)
       AXCYL = DIFF/D
       CENTCYL = PT1+(DEL1-DEL2+D)/2*AXCYL
    ENDIF

  END SUBROUTINE ENVELOPECYLINDER

  SUBROUTINE ADDBEAD(CHAINP,BEAD,DELE)
    ! add an additional bead after the specified BEAD (cannot be last one)
    ! returns consequent change in energy

    USE GENUTIL
    USE CHAINUTIL, ONLY : GETBEADENERGY, OUTPUTSNAPSHOT, GETSTERICENERGY
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: BEAD
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    INTEGER, PARAMETER :: MAXNTRY = 1000000
    DOUBLE PRECISION :: ENEW(2), EPREV, INTERP(6)
    INTEGER :: NPT, IND, TRY, B, FAILFINITEXT, NCLASH
    DOUBLE PRECISION :: XAX(3), YAX(3), UVEC(3), POS(3), UREL(3)
    DOUBLE PRECISION :: ST, RHO, PHI, R1, R2, R3,U, LOGU, ND, ND2
    DOUBLE PRECISION :: DEL, LP, EC, GAM, EPAR, EPH, ECHECK
    DOUBLE PRECISION :: DU1(3), DU2(3), DPERP(3), DPAR, DVEC(3), DIFF(3)
    LOGICAL :: CLASH
    DOUBLE PRECISION :: R1P, R2P, R3P, CVAL, COEFF(4)
    DOUBLE PRECISION :: DIFFPT(3), DIST, RAD1, RAD2, H1, H2, VOL1, VOL2, FRAC1
    DOUBLE PRECISION :: LOGV, XAX2(3), YAX2(3), TMP, TMPSCL, CONST1, CONST2, TMP2, RAD, R1S, R2S, R3S
    LOGICAL :: SUCCESS
    DOUBLE PRECISION :: THETA, AX(3), UVEC1(3), UVEC2(3), ROTMAT(3,3)

    CALL OUTPUTSNAPSHOT(CHAINP,'startaddbead.out',0,.false.)

    IF (CHAINP%NPT.GE.CHAINP%MAXNPT) THEN
       PRINT*, 'ERROR IN ADDBEAD: cannot add new bead because too many beads', CHAINP%NPT, CHAINP%MAXNPT
       STOP 1
    ENDIF

    IF (.NOT.PARAMDATASET) THEN
       PRINT*, 'ERROR IN ADDBEAD: parameter data array has not been set.'
       STOP 1
    ENDIF

    IF (BEAD.LE.0.OR.BEAD.GE.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN ADDBEAD: cannot add after edge bead', BEAD
       STOP 1
    ENDIF

    du1 = chainp%pos(:,bead+1)-chainp%pos(:,bead)
    !print*, 'testx1:', bead, chainp%ls(bead), chainp%gam(bead), sqrt(dot_product(du1,du1))

    EPREV = CHAINP%BEADENERGY(BEAD+1)+CHAINP%STERICENERGY
    

    ! add bead to chain arrays (only those used for MC, not COORDS)
    NPT = CHAINP%NPT
    CHAINP%POS(:,BEAD+2:NPT+1) = CHAINP%POS(:,BEAD+1:NPT)
    CHAINP%UVEC(:,BEAD+2:NPT+1) = CHAINP%UVEC(:,BEAD+1:NPT)    
    CHAINP%BEADENERGY(BEAD+2:NPT+1) = CHAINP%BEADENERGY(BEAD+1:NPT)    
    CHAINP%LS(BEAD+2:NPT) = CHAINP%LS(BEAD+1:NPT-1)
    CHAINP%LP(BEAD+2:NPT) = CHAINP%LP(BEAD+1:NPT-1)
    CHAINP%GAM(BEAD+2:NPT) = CHAINP%GAM(BEAD+1:NPT-1)
    CHAINP%EPAR(BEAD+2:NPT) = CHAINP%EPAR(BEAD+1:NPT-1)
    CHAINP%EPERP(BEAD+2:NPT) = CHAINP%EPERP(BEAD+1:NPT-1)
    CHAINP%EC(BEAD+2:NPT) = CHAINP%EC(BEAD+1:NPT-1)

    CHAINP%NPT = CHAINP%NPT + 1


    ! reset the parameters 
    CHAINP%LS(BEAD) = CHAINP%LS(BEAD)/2
    CHAINP%LS(BEAD+1) = CHAINP%LS(BEAD)

    CALL INTERPARRAY(PARAMDATA,(/NPARAMDATA,6/),1,CHAINP%LS(BEAD),IND,INTERP)
    IF (IND.LE.0.OR.IND.GE.NPARAMDATA) THEN
       PRINT*, 'ERROR IN ADDBEAD: segment length is out of bounds', CHAINP%LS(BEAD), PARAMDATA(1,1),PARAMDATA(NPARAMDATA,1)
       STOP 1
    ENDIF

    CHAINP%LP(BEAD:BEAD+1) = INTERP(2)
    CHAINP%GAM(BEAD:BEAD+1) = INTERP(3)
    CHAINP%EPAR(BEAD:BEAD+1) = INTERP(4)
    CHAINP%EPERP(BEAD:BEAD+1) = INTERP(5)
    CHAINP%EC(BEAD:BEAD+1) = INTERP(6)

! -----------------------------------------------------
! 
!     ! Rejection sampling to place bead with proper distribution
!     DEL = CHAINP%LS(BEAD); GAM = CHAINP%GAM(BEAD); 
!     EPAR = CHAINP%EPAR(BEAD)
!     EPH = CHAINP%EPERP(BEAD)-CHAINP%EC(BEAD)**2/CHAINP%LP(BEAD)
!     LP = CHAINP%LP(BEAD)
!     !    M = SQRT(2*PI*DEL)**3/EPAR/EPERP**2
!     !    LOGM = LOG(M);

!     IF (CHAINP%UVEC(1,BEAD)**2 + CHAINP%UVEC(2,BEAD)**2.EQ.0) THEN
!        XAX = (/1D0,0D0,0D0/)
!        YAX = (/0D0,1D0,0D0/)
!     ELSE
!        CALL CROSS_PRODUCT(CHAINP%UVEC(:,BEAD),(/0D0,0D0,1D0/),XAX)
!        CALL NORMALIZE(XAX)
!        CALL CROSS_PRODUCT(CHAINP%UVEC(:,BEAD),XAX,YAX)
!     ENDIF

!     ! upper bead in coordinate system of lower bead
!     DIFFPT = CHAINP%POS(:,BEAD+2)-CHAINP%POS(:,BEAD)
!     !R3P = DOT_PRODUCT(DIFF,CHAINP%UVEC(:,BEAD))
!     !R1P = DOT_PRODUCT(DIFF,XAX)
!     !R2P = DOT_PRODUCT(DIFF,YAX)
!     !CVAL = (R1P**2/4 + R2P**2/4)/(2*DEL/EPH) + (R3P-DEL*GAM)**2/4/(2*DEL/EPAR)
        
!     RAD1 = DEL; RAD2 = DEL
!     DIST = SQRT(DOT_PRODUCT(DIFFPT,DIFFPT))
!     IF (DIST.GT.RAD1+RAD2) THEN
!        PRINT*, 'ERROR IN BEAD SAMPLING: beads too far apart', dist, rad1+rad2
!        stop 1
!     ENDIF
       

!     ! heights of spherical caps
!     H1 = (RAD2-RAD1+DIST)*(RAD2+RAD1-DIST)/(2*DIST)
!     H2 = (RAD1-RAD2+DIST)*(RAD1+RAD2-DIST)/(2*DIST)
!     VOL1 = H1**2*(3*RAD1-H1) ! 3/pi times the volume of the cap
!     VOL2 = H2**2*(3*RAD2-H2)
!     LOGV = LOG(VOL1+VOL2)
!     FRAC1 = VOL1/(VOL1+VOL2)
!     CONST1 = -RAD1**2*(RAD1-H1)+(RAD1-H1)**3/3D0
!     CONST2 = -RAD2**2*(RAD2-H2)+(RAD2-H2)**3/3D0    

!     CALL CROSS_PRODUCT(YAX,DIFFPT,XAX2); CALL NORMALIZE(XAX2)
!     CALL CROSS_PRODUCT(DIFFPT,XAX2,YAX2); CALL NORMALIZE(YAX2)

!     FAILFINITEXT = 0; NCLASH = 0
!     DO TRY = 1,MAXNTRY
!        !sample from envelope distrib
!        U = GRND(); LOGU = LOG(U)
!        PHI = GRND()*2*PI
!        RHO = 1D0-ABS(SQRT(DEL/LP)*RNORM())       

!        IF (RHO.LT.-1D0) then
!           !print*, 'testx-1:', try, rho
!           CYCLE
!        endif
!        ST = SQRT(1-RHO**2)
!        PHI = GRND()*2*PI

!        ! uniformly sample from intersection of 2 spheres for possible
!        ! center bead position

!        ! sample position along axis and circle radius
!        TMP = GRND(); TMP2 = GRND()
!        IF (TMP.LT.FRAC1) THEN
!           TMPSCL = TMP/FRAC1          
!           COEFF = (/-1D0/3D0,0D0,RAD1**2, CONST1-H1**2*(3*RAD1-H1)/3*TMPSCL/)
!           CALL CUBICROOT(COEFF,(/RAD1-H1,RAD1/),R3S,SUCCESS)
          
!           RAD = SQRT(RAD1**2 - R3S**2)*SQRT(TMP2)
!        ELSE
!           TMPSCL = (TMP-FRAC1)/(1-FRAC1)
!           COEFF = (/-1D0/3D0,0D0,RAD2**2, CONST2-H2**2*(3*RAD2-H2)/3*TMPSCL/)
!           CALL CUBICROOT(COEFF,(/RAD2-H2,RAD2/),R3S,SUCCESS)
!           RAD = SQRT(RAD2**2 - R3S**2)*SQRT(TMP2)

!           R3S = DIST - R3S
!        ENDIF
!        IF (.NOT.SUCCESS) THEN
!           PRINT*, 'ERROR IN BEAD SAMPLING: no cubic root found', TMP, FRAC1, tmpscl, R3S, RAD1, H1, TMPSCL
!           print*, COEFF
!           print*, rad2, h2, const2
!           STOP 1
!        ENDIF

!        THETA = GRND()*2*PI
!        R1S = RAD*COS(THETA); R2S = RAD*SIN(THETA)
       
!        !R1 = SQRT(DEL/EPH)*RNORM()
!        !R2 = SQRT(DEL/EPH)*RNORM()
!        !R3 = SQRT(DEL/EPAR)*RNORM()+GAM*DEL
!        !R1 = SQRT(2*DEL/EPH)*RNORM()+R1P/2
!        !R2 = SQRT(2*DEL/EPH)*RNORM()+R2P/2
!        !R3 = SQRT(2*DEL/EPAR)*RNORM()+R3P/2

!        !IF (BEAD.EQ.4.AND.ABS(DEL-1).LT.1E-12) THEN
!        !   PRINT*, 'testx5:', (R1-R1P/2)**2 + (R2-R2P/2)**2 + (R3-R3P/2)**2, 2*DEL/EPH, 2*DEL/EPAR
!        !ENDIF

!        !POS  = CHAINP%POS(:,BEAD) + R3*CHAINP%UVEC(:,BEAD)+R1*XAX+R2*YAX
!        POS = R3S*DIFFPT/DIST + R1S*XAX2 + R2S*YAX2
!        R1 = DOT_PRODUCT(POS,XAX)
!        R2 = DOT_PRODUCT(POS,YAX)
!        R3 = DOT_PRODUCT(POS,CHAINP%UVEC(:,BEAD))
       
!        POS = CHAINP%POS(:,BEAD)+POS

!        DVEC = CHAINP%POS(:,BEAD+2)-POS      
!        ND2 = DOT_PRODUCT(DVEC,DVEC); ND = SQRT(ND2);

!        ! Check for finite extension 
!        IF (CHAINP%FINITEXT) THEN
!           IF (R1S**2+R2S**2+R3S**2.GT.DEL**2.OR.ND2.GT.DEL**2) THEN
!              FAILFINITEXT = FAILFINITEXT+1
!             ! print*, 'testx0', try, R1S**2+R2S**2+R3S**2,ND2
!              !stop 1
!              CYCLE
!           ENDIF
!        ENDIF

!        ! u vector relative to previous one
!        UrEl = (/ST*COS(PHI),ST*SIN(PHI),RHO/)
!        DU1 = UREL + (/EC/LP*R1,EC/LP*R2,-1D0/)

!        ! absolute u vector
!        UVEC = UREL(1)*XAX + UREL(2)*YAX + UREL(3)*CHAINP%UVEC(:,BEAD)
!        !UVEC = DOT_PRODUCT(UREL,XAX)*XAX + DOT_PRODUCT(UREL,YAX)*YAX + &
!        !     & DOT_PRODUCT(UREL,CHAINP%UVEC(:,BEAD))*CHAINP%UVEC(:,BEAD)
!        DPAR = DOT_PRODUCT(DVEC,UVEC)
!        DPERP = DVEC - DPAR*UVEC
!        DU2 = CHAINP%UVEC(:,BEAD+2) - UVEC + EC/LP*DPERP

!        ECHECK = -LP/2/DEL*DOT_PRODUCT(DU1,DU1) - LP/2/DEL*DOT_PRODUCT(DU2,DU2)&
!             & -  EPAR/2/DEL*(DPAR-GAM*DEL)**2  - EPH/2/DEL*(ND2 - DPAR*DPAR) &
!             & - EPAR/2/DEL*(R3-GAM*DEL)**2 - EPH/2/DEL*(R1**2+R2**2) &
!             & - LOGV

!          ! IF (BEAD.EQ.3) THEN
!         !  print*, 'TESTX1:', TRY, R1,R2,R3,RHO,ECHECK,LOGU
!         !  PRINT*, 'TESTX2:', -LP/2/DEL*DOT_PRODUCT(DU1,DU1), - LP/2/DEL*DOT_PRODUCT(DU2,DU2)
!         ! PRINT*, 'TESTX3:', - EPAR/2/DEL*(DPAR-GAM*DEL)**2, &
!         !      & - EPH/2/DEL*(ND2 - DPAR*DPAR), - EPAR/2/DEL*(R3-GAM*DEL)**2, &
!         !      & - EPH/2/DEL*(R1**2+R2**2), logv
!         !  ENDIF
       
!        !IF (LOGU.GT.ECHECK) CYCLE

!        ! check for steric clashes with existing beads
!        CLASH = .FALSE.
!        ! IF (CHAINP%STERICS) THEN
!        !    DO B = 1,CHAINP%NPT
!        !       IF (ABS(B-BEAD+1).LT.CHAINP%STERSKIP) CYCLE
!        !       DIFF = CHAINP%POS(:,B)-POS
!        !       ND2 = DOT_PRODUCT(DIFF,DIFF)             
!        !       IF (ND2.LT.CHAINP%STERRAD2) THEN
!        !          NCLASH = NCLASH + 1
!        !          CLASH = .TRUE.; EXIT
!        !       ENDIF
!        !    ENDDO          
!        ! ENDIF      
!        IF (.NOT.CLASH) EXIT

!     ENDDO
!     IF (TRY.GE.MAXNTRY) THEN
!       PRINT*, 'ERROR IN SAMPLING BEAD POSITION! Failed to generate successfull sample.', BEAD, DEL
!       PRINT*, 'FRACTION OF TIME FAILED FINITE EXTENSION:', DBLE(FAILFINITEXT)/MAXNTRY
!       PRINT*, 'FRACTION OF TIME FAILED CLASH:', DBLE(NCLASH)/MAXNTRY
!       CALL OUTPUTSNAPSHOT(CHAINP,'fail.snap.out',0,.FALSE.)
!       !STOP 1
!     ENDIF

!     CHAINP%POS(:,BEAD+1) = POS
!     CHAINP%UVEC(:,BEAD+1) = UVEC

! !    PRINT*, 'TESTX0:', CLASH, TRY
!     CALL MINIMONTECARLO(CHAINP,BEAD+1,1000)
! -------------------------------

    ! interpolate bead position
    CHAINP%POS(:,BEAD+1) = (CHAINP%POS(:,BEAD) + CHAINP%POS(:,BEAD+2))/2

    ! interpolate orientation (rotate halfway from one UVEC to the other
    UVEC1 = CHAINP%UVEC(:,BEAD)
    UVEC2 = CHAINP%UVEC(:,BEAD+2)
    CALL CROSS_PRODUCT(UVEC1,UVEC2,AX)
    CALL NORMALIZE(AX)
    THETA = ACOS(DOT_PRODUCT(UVEC1,UVEC2))
    CALL ROTANGAX(THETA/2,AX,UVEC1,CHAINP%UVEC(:,BEAD+1),.TRUE.,ROTMAT)

!    CALL MINIMONTECARLO(CHAINP,BEAD+1,1000)

    ! recalculate energy for affected beads
    !EPREV = CHAINP%BEADENERGY(BEAD+1)
    CALL GETBEADENERGY(CHAINP,BEAD+1,ENEW(1))
    CALL GETBEADENERGY(CHAINP,BEAD+2,ENEW(2))
    CALL GETSTERICENERGY(CHAINP,CHAINP%STERICENERGY)

    DELE = ENEW(1)+ENEW(2)+CHAINP%STERICENERGY-EPREV-CHAINP%MU

  END SUBROUTINE ADDBEAD

  SUBROUTINE REMOVEBEAD(CHAINP,BEAD,DELE)
    ! coarsen the discretization by removing one bead from the chain    
    ! cannot be an edge bead
    ! returns consequent change in energy
    USE GENUTIL, ONLY : INTERPARRAY
    USE CHAINUTIL, ONLY : GETBEADENERGY
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: BEAD
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    DOUBLE PRECISION :: ENEW, EPREV, INTERP(6),TMP,DIFF(3)
    INTEGER :: NPT, IND
    

    IF (.NOT.PARAMDATASET) THEN
       PRINT*, 'ERROR IN REMOVE BEAD: parameter data array has not been set.'
       STOP 1
    ENDIF

    IF (BEAD.LE.1.OR.BEAD.GE.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN REMOVEBEAD: cannot remove edge bead', BEAD
       STOP 1
    ENDIF

    CALL GETBEADENERGY(CHAINP,BEAD,TMP)
    IF (TMP.GT.1D10) THEN
       DIFF = CHAINP%POS(:,BEAD)-CHAINP%POS(:,BEAD-1)
       PRINT*, 'TESTX2:', BEAD,SQRT(DOT_PRODUCT(DIFF,DIFF)), CHAINP%LS(BEAD)
       STOP 1
    ENDIF

    ! remove bead from chain arrays (only those used for MC, not COORDS)
    NPT = CHAINP%NPT
    CHAINP%POS(:,BEAD:NPT-1) = CHAINP%POS(:,BEAD+1:NPT)
    CHAINP%UVEC(:,BEAD:NPT-1) = CHAINP%UVEC(:,BEAD+1:NPT)
    
    !PRINT*, 'TESTX1:', TMP,CHAINP%BEADENERGY(BEAD), CHAINP%BEADENERGY(BEAD+1)
    EPREV = CHAINP%BEADENERGY(BEAD)+CHAINP%BEADENERGY(BEAD+1)
    CHAINP%BEADENERGY(BEAD:NPT-1) = CHAINP%BEADENERGY(BEAD+1:NPT)
    
    CHAINP%LS(BEAD-1) = CHAINP%LS(BEAD-1)+CHAINP%LS(BEAD)

    CHAINP%LS(BEAD:NPT-2) = CHAINP%LS(BEAD+1:NPT-1)
    CHAINP%LP(BEAD:NPT-2) = CHAINP%LP(BEAD+1:NPT-1)
    CHAINP%GAM(BEAD:NPT-2) = CHAINP%GAM(BEAD+1:NPT-1)
    CHAINP%EPAR(BEAD:NPT-2) = CHAINP%EPAR(BEAD+1:NPT-1)
    CHAINP%EPERP(BEAD:NPT-2) = CHAINP%EPERP(BEAD+1:NPT-1)
    CHAINP%EC(BEAD:NPT-2) = CHAINP%EC(BEAD+1:NPT-1)
    
    CHAINP%LS(NPT-1) = 0; CHAINP%LP(NPT-1) = 0;
    CHAINP%GAM(NPT-1) = 0; CHAINP%EPAR(NPT-1) = 0;
    CHAINP%EPERP(NPT-1) = 0; CHAINP%EC(NPT-1) = 0;
    CHAINP%BEADENERGY(NPT) = 0

    ! decrease number of beads
    CHAINP%NPT = NPT-1

    ! get new parameters
    CALL INTERPARRAY(PARAMDATA,(/NPARAMDATA,6/),1,CHAINP%LS(BEAD-1),IND,INTERP)
    IF (IND.LE.0.OR.IND.GE.NPARAMDATA) THEN
       PRINT*, 'ERROR IN REMOVEBEAD: segment length is out of bounds', CHAINP%LS(BEAD-1), PARAMDATA(1,1),PARAMDATA(NPARAMDATA,1)
       STOP 1
    ENDIF

    CHAINP%LP(BEAD-1) = INTERP(2)
    CHAINP%GAM(BEAD-1) = INTERP(3)
    CHAINP%EPAR(BEAD-1) = INTERP(4)
    CHAINP%EPERP(BEAD-1) = INTERP(5)
    CHAINP%EC(BEAD-1) = INTERP(6)    

    ! recalculate energy for affected bead
    CALL GETBEADENERGY(CHAINP,BEAD,ENEW)
    DELE = ENEW-EPREV+CHAINP%MU
    !PRINT*, 'TESTX2:', ENEW,EPREV,CHAINP%MU
    ! IF (DELE.LT.-1E10) THEN
    !    PRINT*, DELE, BEAD, EPREV, ENEW
    !    STOP 2
    ! ENDIF
  END SUBROUTINE REMOVEBEAD

  SUBROUTINE READPARAMDATA(INFILE)
    ! read in a file containing parameter data
    ! columns in order: del, lp, gam, epar, eperp, ec, plen, err
    ! first 6 columns are saved into global array PARAMDATA

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: INFILE
    DOUBLE PRECISION :: DATAROW(8), TMP1, TMP2
    INTEGER :: C, NROW, ICHECK

    ! find number of lines in the file
    OPEN(UNIT=88,FILE=INFILE,STATUS='OLD')
    NROW = 0
    DO
       READ(88,*,IOSTAT=ICHECK) DATAROW
       IF (ICHECK.LT.0) EXIT
       NROW = NROW + 1
    ENDDO
    NPARAMDATA = NROW

    ALLOCATE(PARAMDATA(NROW,6))
    PARAMDATASET = .TRUE.

    REWIND(88)

    DO C = 1,NROW
       READ(88,*,IOSTAT=ICHECK) PARAMDATA(C,:),TMP1,TMP2
    ENDDO
  END SUBROUTINE READPARAMDATA

  SUBROUTINE CLEANUPDATA
    DEALLOCATE(PARAMDATA)
  END SUBROUTINE CLEANUPDATA

  SUBROUTINE CUBICROOT(COEFF,RANGE,ROOT, SUCCESS)
    ! for a cubic equation Ax^3 + Bx^2 + Cx + D
    ! return the real root within range RANGE if one exists
    ! otherwise SUCCESS=.false.
    USE GENUTIL, ONLY : PI
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: COEFF(4),RANGE(2)
    DOUBLE PRECISION, INTENT(OUT) :: ROOT
    LOGICAL, INTENT(OUT) :: SUCCESS
    DOUBLE PRECISION :: A,B,C
    DOUBLE PRECISION :: Q, R, SQ, THETA, AVAL, BVAL

    A = COEFF(2)/COEFF(1)
    B = COEFF(3)/COEFF(1)
    C = COEFF(4)/COEFF(1)

    ! algorithm from numerical recipes
    Q = (A**2 - 3*B)/9;
    R = (2*A**3 - 9*A*B + 27*c)/54

    IF (R**2.LT.Q**3) THEN
!       print*, 'testx1'
       ! three real roots
       SQ = SQRT(Q)
       THETA = ACOS(R/SQ**3)
       SUCCESS = .TRUE.
       ROOT = -2*SQ*COS(THETA/3)-A/3
       IF (ROOT.GE.RANGE(1).AND.ROOT.LE.RANGE(2)) RETURN
       ROOT = -2*SQ*COS((THETA+2*PI)/3) - A/3
       IF (ROOT.GE.RANGE(1).AND.ROOT.LE.RANGE(2)) RETURN
       ROOT = -2*SQ*COS((THETA-2*PI)/3) - A/3
       IF (ROOT.GE.RANGE(1).AND.ROOT.LE.RANGE(2)) RETURN
       SUCCESS = .FALSE.
    ELSE
      ! print*, 'testx2'
       AVAL = -SIGN(1D0,R)*(ABS(R) + SQRT(R**2-Q**3))**(1D0/3D0)
       IF (AVAL.EQ.0) THEN
          BVAL=0
       ELSE
          BVAL = Q/AVAL
       ENDIF
       ROOT = (AVAL+BVAL) - A/3
       IF (ROOT.GE.RANGE(1).AND.ROOT.LE.RANGE(2)) THEN
          SUCCESS = .TRUE.
       ELSE
          SUCCESS = .FALSE.
       ENDIF
    ENDIF
  END SUBROUTINE CUBICROOT

  SUBROUTINE MINIMONTECARLO(CHAINP,BEAD,NSTEPS)
    ! run a short monte carlo to equilibrate a particular new bead
    USE MT19937, ONLY : GRND
    USE CHAINUTIL, ONLY : GETBEADENERGY, GETFORCEENERGY
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: BEAD, NSTEPS
!    DOUBLE PRECISION, INTENT(FINAL) :: DELEFINAL
    DOUBLE PRECISION :: STARTENERGY, PREVENERGIES(2), PREVPOS(3), PREVUVEC(3), DELE
    INTEGER :: STEP
    LOGICAL :: ACCEPT, STARTCLASH
    DOUBLE PRECISION :: ARANGE, RRANGE, TMP

    ARANGE = 0.1D0
    RRANGE = CHAINP%LS(BEAD)/10

    IF (BEAD.EQ.CHAINP%NPT) THEN
       PRINT*, 'ERROR IN MINIMONTECARLO: not set up for edge bead'
       STOP 1
    ENDIF
    
    STARTCLASH = .FALSE.
    CALL BEADMOVE3(CHAINP,0D0,0D0,BEAD,DELE)
    STARTENERGY = CHAINP%BEADENERGY(BEAD)+CHAINP%BEADENERGY(BEAD+1)+CHAINP%FORCEENERGY+CHAINP%STERICENERGY

    !CALL GETBEADENERGY(CHAINP,BEAD,CHAINP%BEADENERGY(BEAD))
    !CALL GETBEADENERGY(CHAINP,BEAD+1,CHAINP%BEADENERGY(BEAD+1))
    !CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY)
    !STARTENERGY = CHAINP%BEADENERGY(BEAD)+CHAINP%BEADENERGY(BEAD+1)+CHAINP%FORCEENERGY
    PREVENERGIES = CHAINP%BEADENERGY(BEAD:BEAD+1)
    PREVPOS = CHAINP%POS(:,BEAD)
    PREVUVEC = CHAINP%UVEC(:,BEAD)

    DO STEP = 1,NSTEPS
       CALL BEADMOVE3(CHAINP,ARANGE,RRANGE,BEAD,DELE)

       IF (DELE.LT.0) THEN
          ACCEPT = .TRUE.
       ELSE
          TMP = GRND()
          ACCEPT = (TMP.LT.EXP(-DELE))
       ENDIF

       ! IF (STARTCLASH) THEN
       !    ACCEPT = .TRUE.
       !    IF (DELE.LT.HUGE(1D0)) STARTCLASH = .FALSE.
       ! ENDIF

       IF (ACCEPT) THEN
          PREVENERGIES = CHAINP%BEADENERGY(BEAD:BEAD+1)
          PREVPOS = CHAINP%POS(:,BEAD)
          PREVUVEC = CHAINP%UVEC(:,BEAD)
       ELSE
          CHAINP%BEADENERGY(BEAD:BEAD+1) = PREVENERGIES
          CHAINP%POS(:,BEAD) = PREVPOS
          CHAINP%UVEC(:,BEAD) = PREVUVEC
       ENDIF

      ! IF (MOD(STEP,100).EQ.0) THEN
      !    print*, 'step, energy:', STEP, CHAINP%BEADENERGY(BEAD)+CHAINP%BEADENERGY(BEAD+1), dele, accept
      ! ENDIF
    ENDDO

  END SUBROUTINE MINIMONTECARLO

  SUBROUTINE BEADMOVE3(CHAINP,ARANGE,RRANGE,B,DELE)
    ! move and rotate an individual bead
    USE MT19937, ONLY : GRND
    USE CHAINUTIL, ONLY : GETSTERICENERGY, GETBEADENERGY, GETFORCEENERGY, CHAIN
    USE QUAtUTIL, ONLY : ROTANGAX
    IMPLICIT NONE

    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: ARANGE, RRANGE
    INTEGER, INTENT(in) :: B
    DOUBLE PRECISION, INTENT(OUT) :: DELE
    DOUBLE PRECISION :: AX(3), ANG, DELR(3), PREVE, TMP(3),ROTMAT(3,3)
    INTEGER :: I, ALLBEADS(CHAINP%NPT), B2
    LOGICAL :: CLASH
    DOUBLE PRECISION :: DIFF(3), ND2

    ! get a random bead
    !B = FLOOR(GRND()*CHAINP%NPT)+1

    ! get a random axis
    DO I = 1,3
       AX(I) = GRND()
    ENDDO
    AX = AX/SQRT(DOT_PRODUCT(AX,AX))

    ! get a random angle
    ANG = GRND()*2*ARANGE-ARANGE
    
    ! get a random shift
    IF (CHAINP%STRETCHABLE.AND.CHAINP%SHEARABLE) THEN
       DO I = 1,3
          DELR(I) = GRND()*2*RRANGE - RRANGE
       ENDDO
    ELSEIF (CHAINP%STRETCHABLE.AND..NOT.CHAINP%SHEARABLE) THEN
       PRINT*, 'STRETCHABLE AND NOT SHEARABLE NOT YET SET UP'
       STOP 2
    ELSEIF (CHAINP%SHEARABLE.AND..NOT.CHAINP%STRETCHABLE) THEN
       PRINT*, 'SHEARABLE AND NOT STRETCHABLE NOT YET SET UP'
       STOP 2
    ELSE
       DELR = 0D0
    ENDIF

    ! rotate selected bead
    CALL ROTANGAX(ANG,AX,CHAINP%UVEC(:,B),TMP,.TRUE.,ROTMAT)
    CHAINP%UVEC(:,B) = TMP
    ! shift selected bead
    CHAINP%POS(:,B) = CHAINP%POS(:,B) + DELR

    ! update bordering energies
    IF (B.EQ.1) THEN
       PREVE = CHAINP%BEADENERGY(B+1)+CHAINP%FORCEENERGY+CHAINP%STERICENERGY
       CALL GETBEADENERGY(CHAINP,B+1,CHAINP%BEADENERGY(B+1))
       CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY)
       DELE = CHAINP%BEADENERGY(B+1)+CHAINP%FORCEENERGY-PREVE
    ELSEIF (B.EQ.CHAINP%NPT) THEN
       PREVE = CHAINP%BEADENERGY(B)+CHAINP%FORCEENERGY+CHAINP%STERICENERGY
       CALL GETBEADENERGY(CHAINP,B,CHAINP%BEADENERGY(B))
       CALL GETFORCEENERGY(CHAINP,CHAINP%FORCEENERGY)
       DELE = CHAINP%BEADENERGY(B)+CHAINP%FORCEENERGY-PREVE
    ELSE
       PREVE = CHAINP%BEADENERGY(B+1) + CHAINP%BEADENERGY(B)+CHAINP%STERICENERGY
       CALL GETBEADENERGY(CHAINP,B,CHAINP%BEADENERGY(B))
       CALL GETBEADENERGY(CHAINP,B+1,CHAINP%BEADENERGY(B+1))
       DELE = CHAINP%BEADENERGY(B+1) + CHAINP%BEADENERGY(B) - PREVE
    ENDIF
    
    IF (CHAINP%STERICS) CALL GETSTERICENERGY(CHAINP,CHAINP%STERICENERGY)
    DELE = DELE + CHAINP%STERICENERGY

! THEN
!        DO B2 = 1,CHAINP%NPT
!           IF (ABS(B2-B).Le.CHAINP%STERSKIP) CYCLE
!           DIFF = CHAINP%POS(:,B)-CHAINP%POS(:,B2)
!           ND2 = DOT_PRODUCT(DIFF,DIFF)             
!           IF (ND2.LT.CHAINP%STERRAD2) THEN
!              DELE = HUGE(1D0)
!              RETURN
!           ENDIF
!        ENDDO
!     ENDIF

    ! IF (CHAINP%STERICS) THEN
    !    ALLBEADS = (/(I, I=1,CHAINP%NPT)/)
    !    CALL CHECKSTERICCLASH(CHAINP,1,(/B/),CHAINP%NPT,ALLBEADS,CLASH)
    !    IF (CLASH) THEN
    !       DELE = HUGE(1D0)
    !       RETURN
    !    ENDIF
    ! ENDIF

    ! IF (DELE.GT.HUGE(1D0)/100) THEN
    !    PRINT*, 'BAD DELE IN BEADMOVE3:', DELE, CHAINP%BEADENERGY(B), CHAINP%BEADENERGY(B+1)
    !    STOP 1
    ! ENDIF
  END SUBROUTINE BEADMOVE3
END MODULE REDISC
