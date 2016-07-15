MODULE CYLINDERUTIL
  ! utilities for dealing with cylindrical sterics
  USE GENUTIL
  IMPLICIT NONE  

CONTAINS
  LOGICAL FUNCTION CYLINDERINTERSECT(RA,RB,LA,LB,CA,NA,CB,NB)
    ! check if two cylinders intersect
    ! RA, RB are the radii
    ! LA,LB are the cylinder lengths
    ! CA,NA are the center and normalized axis of the 1st cylinder
    ! CB,NB are for 2nd cylinder
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RA, RB, LA, LB
    DOUBLE PRECISION, INTENT(IN) :: CA(3), NA(3), CB(3), NB(3)
    DOUBLE PRECISION :: DIST, CENTA(3),CENTB(3), CENTB2(3), TMIN
    INTEGER :: I, J
    LOGICAL :: INSEG1, INSEG2


    CYLINDERINTERSECT = .FALSE.

    ! get distance between line segments
    CALL LINESEGDIST(NA,CA,NB,CB,LA,LB,DIST,INSEG1,INSEG2)

    IF (DIST.GT.(RA+RB)**2) then
       RETURN ! lower bound is above 0       
    ELSE IF (INSEG1.AND.INSEG2) THEN
       CYLINDERINTERSECT = .TRUE. ! shells intersect
       RETURN
    ENDIF


    DO I = -1,1,2
       ! check if axis B crosses any disc on A
       CENTA = CA + I*NA*LA/2
       CYLINDERINTERSECT = LINEDISCINTERSECT(RA,CENTA,NA,NB,CB,LB)       
       IF (CYLINDERINTERSECT) RETURN

       ! check if axis A crosses any disc on B
       CENTB = CB + I*NB*LB/2
       CYLINDERINTERSECT = LINEDISCINTERSECT(RB,CENTB,NB,NA,CA,LA)       

       IF (CYLINDERINTERSECT) RETURN       

       ! check if any pair of discs intersects
       DO J = -1,1,2
          CENTB2 = CB + J*NB*LB/2         

          CYLINDERINTERSECT =  DISCDISCINTERSECT(RA,CENTA,NA,RB,CENTB2,NB)         
          IF (CYLINDERINTERSECT) RETURN
       ENDDO

       ! check for circle-shell intersections
       ! get distance between circle A and axis B      
       CALL CIRCLELINEDIST(RA,CENTA,NA,NB,CB,DIST,TMIN)
       IF (DIST.LT.RB*RB.AND.ABS(TMIN).LT.LB/2) THEN
          CYLINDERINTERSECT = .TRUE.; RETURN
       ENDIF

       ! get distance between circle B and axis A
       CALL CIRCLELINEDIST(RB,CENTB,NB,NA,CA,DIST,TMIN)
       IF (DIST.LT.RA*RA.AND.ABS(TMIN).LT.LA/2) THEN
          CYLINDERINTERSECT = .TRUE.; RETURN
       ENDIF
    ENDDO

  END FUNCTION CYLINDERINTERSECT

  LOGICAL FUNCTION LINEDISCINTERSECT(RD,CD,ND,M,B,L)
    ! check if a line intersects a disc
    ! RD, CD, ND are radius, center and normal of disc
    ! line is M*T+B; length of line is L
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: CD(3), ND(3), M(3), B(3),L,RD
    DOUBLE PRECISION :: T, BC(3), PT(3)

    ! get point where line passes through disc plane
    BC = B-CD
    T = -DOT_PRODUCT(BC,ND)/DOT_PRODUCT(M,ND)

    LINEDISCINTERSECT = .FALSE.
    IF (ABS(T).LE.L/2) THEN
       ! check distance of intersection from disc center
       PT = M*T + BC
       IF (DOT_PRODUCT(PT,PT).LE.RD*RD) THEN
          LINEDISCINTERSECT = .TRUE.
       ENDIF
    ENDIF

  END FUNCTION LINEDISCINTERSECT

  LOGICAL FUNCTION DISCDISCINTERSECT(RA,CA,NA,RB,CB,NB)
    ! check if two discs intersect
    ! assume na, nb are normalized    
    IMPLICIT NONE
    DOUBLE PRECISION :: RA, RB, CA(3), NA(3), CB(3), NB(3)
    DOUBLE PRECISION, PARAMETER :: TINY = 1D-10
    DOUBLE PRECISION :: CBA(3), V(3), U(3), ST, CT, DIFF(3)
    double precision :: DIST, T,PT1(3)

    DISCDISCINTERSECT = .FALSE.
    CBA = CB-CA   

    IF (DOT_PRODUCT(CBA,CBA).GT.(RA+RB)**2) RETURN

    IF (ABS(ABS(DOT_PRODUCT(NA,NB))-1).LT.TINY) THEN
       ! discs are parallel
       DISCDISCINTERSECT = ABS(DOT_PRODUCT(CBA,NA)).LT.TINY           
    ELSE       
       ! consider case where circle A passes through disc B
       CALL CROSS_PRODUCT(NA,NB,U); U = U/NORM(U)
       CALL CROSS_PRODUCT(NA,U,V); 

       ST = DOT_PRODUCT(CBA,NB)/RA/DOT_PRODUCT(V,NB)        
       IF (ABS(ST).GT.1) RETURN ! circle never passes through plane of disc
       CT = SQRT(1-ST**2)

       IF (CT.EQ.0D0) THEN
          ! circle A only hits disc B at one point
          PT1 = RA*CT*U + RA*ST*V -CBA
          DISCDISCINTERSECT = DOT_PRODUCT(PT1,PT1).LE.RB*RB
          RETURN
       ENDIF

       PT1 = -RA*CT*U + RA*ST*V +CA
       DIFF = 2*RA*CT*U      

       ! check if each relevant endpoint on circle A is within the disc
       ! and otherwise, whether some point on the line between them is within the disc
       IF (DOT_PRODUCT(PT1-CB,PT1-CB).LE.RB*RB) THEN
          DISCDISCINTERSECT = .TRUE.; RETURN          
       ELSEIF (DOT_PRODUCT(PT1+DIFF-CB,PT1+DIFF-CB).LE.RB*RB) THEN          
          DISCDISCINTERSECT = .TRUE.; RETURN
       ELSE
          calL ptlinedist(cb,diff,pt1,DIST,T)          
          IF (DIST.LE.RB*RB.AND.T.LT.1.AND.T.GT.0) THEN
             DISCDISCINTERSECT = .TRUE.; RETURN
          ENDIF
       ENDIF

    ENDIF
  END FUNCTION DISCDISCINTERSECT

  SUBROUTINE PTLINEDIST(PT,M,B,DIST,T)
    ! give the squared distance between a point and a line
    ! and the position on the line L = M*t + B where the distance vector hits
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PT(3), M(3),B(3)
    DOUBLE PRECISION, INTENT(OUT) :: DIST, T
    DOUBLE PRECISION :: DIFF(3)

    T = DOT_PRODUCT(M,PT-B)/DOT_PRODUCT(M,M)
    DIFF = M*T+B-PT
    DIST = DOT_PRODUCT(DIFF,DIFF)

  END SUBROUTINE PTLINEDIST

    SUBROUTINE LINESEGDIST(M,B,U,V,LENA,LENB,DIST,INSEG1,INSEG2)
    ! calculate the minimal squared distance between 2 line segments
    ! line segment 1 has center B, slope M, length LENA
    ! line segment 2 has center V, slope U, length LENB
    ! INSEG is true if the points at minimal separation fall within (not at ends) of each segment
    ! see cylinder-cylinder notes from 10/19/2009

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M(3), B(3),U(3), V(3),LENA,LENB
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    LOGICAL, INTENT(OUT) :: INSEG1,INSEG2
    DOUBLE PRECISION :: BV(3), UU, UM, MM, MBV, TMP, T, S, DIFF(3)
    DOUBLE PRECISION :: Z, CMP, SMM, SUU

    INSEG1 = .FALSE.; INSEG2 = .FALSE.
    BV = B - V
    UU = DOT_PRODUCT(U,U)
    UM = DOT_PRODUCT(U,M)
    MM = DOT_PRODUCT(M,M)
    MBV = DOT_PRODUCT(M,BV)

    IF (UM.EQ.0D0) THEN
       ! lines are perpendicular
       T = MBV/MM
       S = DOT_PRODUCT(U,BV)/UU
    ELSE
       TMP = (UM-UU*MM/UM)       
       IF (ABS(TMP).LT.EPSILON(1D0)) THEN
          ! lines are parallel and do not intersect       
          ! T = -MBV/MM
          ! INSEG2 = .TRUE.
          ! INSEG1 = ABS(T).LT.LENA/2
          ! DIST = DOT_PRODUCT(M*T+BV,M*T+BV)   
          
          SMM = SQRT(MM); SUU = SQRT(UU)
          Z = -MBV/SMM
          CMP = SMM*LENA/2 + SUU*LENB/2

          IF (ABS(Z).LT.CMP) THEN
             ! distance btwn segments is distance btwn lines
             DIST = DOT_PRODUCT(BV,BV)-Z*Z
             INSEG1 = .TRUE.; INSEG2 = .TRUE.
          ELSE ! segments are offset
             T = SIGN(LENA/2,Z)
             S = -SIGN(LENB/2,Z)
             IF (UM.LT.0) THEN ! segments are oppositely oriented
                ! flip which endpoint to use on 2nd segment
                S = -S
             ENDIF
             DIFF = M*T-U*S+BV
             DIST = DOT_PRODUCT(DIFF,DIFF)
          ENDIF
          RETURN
       ENDIF

       T= (UU/UM*MBV-DOT_PRODUCT(BV,U))/TMP
       S = (MM*T + MBV)/UM
    ENDIF

    INSEG1 = ABS(T).LT.LENA/2
    IF (.NOT.INSEG1) THEN   
       T = SIGN(LENA/2,T-LENA/2)     
       CALL PTLINEDIST(M*T + B,U,V,DIST,S)
    ENDIF
    
    INSEG2 = ABS(S).LT.LENB/2
    IF (.NOT.INSEG2) THEN
       S = SIGN(LENB/2,S-LENB/2)
       CALL PTLINEDIST(U*S + V,M,B,DIST,T)
    ENDIF
    
    INSEG1 = ABS(T).LT.LENA/2
    IF (.NOT.INSEG1) THEN
       T = SIGN(LENA/2,T-LENA/2) 
    ENDIF   

    DIFF = M*T + B - U*S - V
    DIST = DOT_PRODUCT(DIFF,DIFF)    
  END SUBROUTINE LINESEGDIST

  SUBROUTINE CIRCLELINEDIST(RA,CA,NA,MV,BV,DIST,TMIN)
    ! minimial squared distance between a circle of radius RA, centered at CA
    ! with normal (normalized) given by NA
    ! and the line L(t) = MV*t + BV
    ! also return the point at which minimal distance happens
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: CA(3), NA(3), MV(3), BV(3), RA
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    DOUBLE PRECISION :: A0,A1,A2,A3,A4,A5,A6,ACOEFF(7)
    DOUBLE PRECISION :: DV(3), EV(3), FV(3)
    DOUBLE PRECISION :: AX, BX, CX, TMIN, FA, FD, DFC, DFA,T
    INTEGER :: TC
    DOUBLE PRECISION :: DBRENT
    EXTERNAL DBRENT
    DOUBLE PRECISION :: C1, B0, B1, B2, B3,B4,D0,D1,D2
    DOUBLE PRECISION :: DSCR, SDSCR    

    DV = BV-CA
    EV = MV - DOT_PRODUCT(NA,MV)*NA
    FV = DV - DOT_PRODUCT(NA,DV)*NA

    A6 = DOT_PRODUCT(MV,MV)
    A5 = 2*DOT_PRODUCT(DV,MV)
    A4 = DOT_PRODUCT(DV,DV) + RA**2
    A3 = -2*RA
    A2 = DOT_PRODUCT(EV,EV)
    A1 = 2*DOT_PRODUCT(EV,FV)
    A0 = DOT_PRODUCT(FV,FV)

    ACOEFF = (/A0,A1,A2,A3,A4,A5,A6/)

    ! bracket the interval
    BX = DOT_PRODUCT(CA-BV,MV)/DOT_PRODUCT(MV,MV)
    AX = BX - RA; CX = BX + RA

    ! -------------
    ! get first minimum of function
    ! -------------
    DIST = DBRENT(AX,BX,CX,CLFUNC,10,ACOEFF,1d-7,TMIN)
    call cLfunc(TMIN,acoeff,FD,dfA)    

    ! get coefficients for checking for 2nd minimum
    C1 = FD-A4
    B4 = A6**2
    B3 = 2*A6*A5
    B2=A5**2 - 2*A6*C1 - A3**2*A2
    B1 = -2*A5*C1-A3**2*A1
    B0 = C1**2 - A3**2*A0    
    
    ! deflate to 2nd degree polynomial by dividing by the double root
    D2 = B4
    D1 = B3 + 2*B4*TMIN
    D0 = B2+3*B4*TMIN**2+2*B3*TMIN

    ! find roots of 2nd degree polynomial (if it exists); use this to bracket
    DSCR = D1**2 - 4*D0*D2   
    IF (DSCR.GT.0D0) THEN
       SDSCR = SQRT(DSCR)
       IF (D2.GT.0) THEN
          AX = (-D1 - SDSCR)/(2*D2)
          CX = (-D1 + SDSCR)/(2*D2)
       ELSE
          AX = (-D1 + SDSCR)/(2*D2)
          CX = (-D1 - SDSCR)/(2*D2)
       ENDIF
       BX = (AX+CX)/2       

       ! get 2nd minimum
       DIST = DBRENT(AX,BX,CX,CLFUNC,10,ACOEFF,1d-7,TMIN)
       
    endif

  END SUBROUTINE CIRCLELINEDIST

  SUBROUTINE CLFUNC(T,ACOEFF,F,DF)
    ! get the squared distance btwn a circle and a line using precalculated coefficients
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: T,ACOEFF(:)
    DOUBLE PRECISION, INTENT(OUT) :: F, DF
    DOUBLE PRECISION :: TMP

    TMP = ACOEFF(3)*T**2 + ACOEFF(2)*T+ACOEFF(1)

    IF (TMP.LT.-100*EPSILON(1D0)) THEN
       PRINT*, 'ERROR IN CLFUNC: bad value inside square root', TMP
       STOP
    ELSE
       TMP = SQRT(MAX(0D0,TMP))
    ENDIF

    F = ACOEFF(7)*T**2 + ACOEFF(6)*T + ACOEFF(5) + ACOEFF(4)*TMP
    DF = 2*ACOEFF(7)*T + ACOEFF(6) + ACOEFF(4)/2*(2*ACOEFF(3)*T + ACOEFF(2))/TMP

  END SUBROUTINE CLFUNC
END MODULE CYLINDERUTIL
