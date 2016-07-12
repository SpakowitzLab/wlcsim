MODULE QUATUTIL
  ! utilities for dealing with quaternions
  ! and other representations of rotation
  ! including euler angles, rotation matrices, and alpha+gamma and z-axis-vector representations

  IMPLICIT NONE
  LOGICAL :: TESTQUAT = .FALSE.
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  ! when z1^2+z2^2 goes below tiny when working with alpha+gamma & zvec coordinates
  ! use the v->0 approximation
  DOUBLE PRECISION, PARAMETER :: NZTINY=1D-14

  ! definition of a quaternion class
  TYPE QUATERNION
     DOUBLE PRECISION :: W, X, Y, Z
  END TYPE QUATERNION

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE QPRODUCT
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE QDIVIDE
  END INTERFACE

CONTAINS
  SUBROUTINE QUAT2SCREW(QUAT,TRANS,HELCRD)
    ! convert from a rotation+translation to screw coordinates
    ! (equivalently to overall helix coordinates)
    ! returns height, angle, radius, orientation of
    ! canonical system relative to helix system (3 euler angles)
    TYPE(QUATERNION), INTENT(IN) :: QUAT
    DOUBLE PRECISION, INTENT(IN) :: TRANS(3)
    DOUBLE PRECISION, INTENT(OUT) :: HELCRD(6)
    DOUBLE PRECISION :: THETA, AX(3), TP(3), ST2, CT2, AXT3, A, B, G, NP, NTP

    ! angle-axis representation from quaternion
    THETA = 2*ACOS(QUAT%W)

    IF (THETA.EQ.0D0) THEN ! translation only
       HELCRD(1) = SQRT(DOT_PRODUCT(TRANS,TRANS))
       HELCRD(2) = 0D0
       HELCRD(3) = 0D0
       HELCRD(4:6) = (/0D0,0D0,0D0/)
       RETURN
    ENDIF

    ! angle of rotation around screw axis
    HELCRD(2) = THETA

    AX = (/QUAT%X,QUAT%Y,QUAT%Z/);
    ST2 = SIN(THETA/2)
    AX = AX/ST2; ! normalize the axis

    ! shift along screw axis
    HELCRD(1) = DOT_PRODUCT(TRANS,AX)

    ! translation in plane perpendicular to axis
    TP = TRANS - HELCRD(1)*AX
    NTP = SQRT(DOT_PRODUCT(TP,TP))

    ! radius of screw
    NP = NTP/(2*ST2)
    HELCRD(3) = NP

    ! cross-product of AX x TP
    AXT3 = AX(1)*TP(2) - AX(2)*TP(1)

    ! euler angles of screw axis coord system relative to canonical
    A = atan2(AX(1),-AX(2))
    B = ACOS(AX(3))
    CT2 = QUAT%W/ST2
    G = ATAN2(-TP(3)-CT2*AXT3,-AXT3+CT2*TP(3))

    ! euler angles of canonical system relative to screw axis
    HELCRD(4:6) = (/PI-G,B,PI-A/)
  END SUBROUTINE QUAT2SCREW

  SUBROUTINE COORDS2QUAT(AG,ZVEC,Q)
    ! convert from an alpha+gamma angle and a vector along the z axis
    ! to a normalized quaternion
    ! note: this doesn't work if beta = pi
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: AG, ZVEC(3)
    TYPE(QUATERNION) :: Q
    DOUBLE PRECISION :: NZ, ZAX(3), ALPHA,BETA,GAMMA

    NZ = SQRT(DOT_PRODUCT(ZVEC,ZVEC))
    ZAX = ZVEC/NZ

    NZ = ZAX(1)**2+ZAX(2)**2

    IF (NZ.EQ.0) THEN
       CALL EULER2QUAT((/AG,0D0,0D0/),Q)
    ELSE
       ALPHA = ATAN2(ZAX(1),-ZAX(2)); GAMMA = AG-ALPHA
       BETA = ACOS(ZAX(3))
       CALL EULER2QUAT((/ALPHA,BETA,GAMMA/),Q)
    ENDIF

  END SUBROUTINE COORDS2QUAT

  SUBROUTINE COORDS2ROTMAT(AG,ZVEC,MAT,DMATAG,DMATZ)
    ! convert from coordinates that include the alpha+gamma euler angle
    ! and a non-normalized vector along the Z axis in canonical reference system
    ! to a rotation matrix
    ! if DMATAG and DMATZ are provided, also get derivatives of all the matrix
    ! components wrt the coordinates

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: AG, ZVEC(3)
    DOUBLE PRECISION, INTENT(OUT) :: MAT(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DMATAG(3,3), DMATZ(3,3,3)
    DOUBLE PRECISION :: NZTOT2, NZTOT, NZ, Z1, Z2, Z3
    DOUBLE PRECISION :: Z11,Z22,Z12,Z13,Z23,CAG,SAG
    DOUBLE PRECISION :: DXDZ(3,3),DYDZ(3,3)

    IF ((PRESENT(DMATAG).AND..NOT.PRESENT(DMATZ)).OR.&
         & (.NOT.PRESENT(DMATAG).AND.PRESENT(DMATZ))) THEN
       PRINT*, 'ERROR IN COORDS2ROTMAT: neither or both of DMATAG and DMATZ must be provided'
       stop 1
    ENDIF

    MAT = 0D0

    ! normalize the zaxis
    NZTOT2 = DOT_PRODUCT(ZVEC,ZVEC); NZTOT = SQRT(NZTOT2)
    MAT(:,3) = ZVEC/NZTOT
    Z1 = MAT(1,3); Z2 = MAT(2,3);
    Z11 = Z1*Z1; Z22 = Z2*Z2;
    NZ = Z11+Z22;
    Z3 = MAT(3,3)
    Z12 = Z1*Z2; Z23 = Z2*Z3; Z13 = Z1*Z3

    ! IF (Z3.LT.-1D0+SQRT(NZTINY)) THEN
    !    PRINT*, 'ERROR IN COORDS2ROTMAT: have hit gimbal lock with Z axis pointing downward. Not set up to deal with this yet.'
    !    STOP 1
    ! ENDIF

    CAG = COS(AG); SAG = SIN(AG)

    !IF (TESTQUAT) PRINT*, 'TESTX1', NZ, NZTINY
    IF (NZ.LT.NZTINY) THEN
       IF (ZVEC(3).LT.0) THEN
          PRINT*, 'PROBLEM IN COORDS2ROTMAT: some z-vector falls along the &
               & negative z axis. This causes gimbal lock. &
               & If working with a single nucleosome you may be able to &
               & avoid this by rotating the entire structure slightly, or &
               & by using RANDSTART to start the linker beads in different positions.'
          stop 1
       endif
       MAT(1,1) = CAG - (Z12*SAG+Z11*CAG)/2
       !if (testquat) print*, 'testx2:', mat(1,1)
       MAT(2,1) = SAG - (Z22*SAG+Z12*CAG)/2

       MAT(1,2) = -SAG+(Z11*SAG-Z12*CAG)/2
       MAT(2,2) = CAG-(Z22*CAG-Z12*SAG)/2
    ELSE
       MAT(1,1) = (Z22*CAG-Z12*SAG+Z3*(Z12*SAG+Z11*CAG))/NZ
       MAT(2,1) = (-Z12*CAG + Z11*SAG+Z3*(Z22*SAG+Z12*CAG))/NZ

       MAT(1,2) = -(Z22*SAG+Z12*CAG+Z3*(-Z12*CAG+Z11*SAG))/NZ
       !if (testquat) print*, 'testx2:', mat(1,2)
       MAT(2,2) = (Z12*SAG+Z11*CAG+Z3*(Z22*CAG-Z12*SAG))/NZ
    ENDIF
    MAT(3,1) = -Z2*SAG-Z1*CAG
    MAT(3,2) = -CAG*Z2+SAG*Z1

    IF (PRESENT(DMATZ)) THEN
       DMATZ = 0D0

       ! derivative of normalized z axis wrt ZVEC coordinates (transposed)
       DMATZ(:,3,1) = -ZVEC(1)*MAT(:,3)/NZTOT2 + (/1D0/NZTOT,0D0,0D0/)
       DMATZ(:,3,2) = -ZVEC(2)*MAT(:,3)/NZTOT2 + (/0D0,1D0/NZTOT,0D0/)
       DMATZ(:,3,3) = -ZVEC(3)*MAT(:,3)/NZTOT2 + (/0D0,0D0,1D0/NZTOT/)

       IF (NZ.LT.NZTINY) THEN

          DXDZ(1,1) = -Z2*SAG/2-Z1*CAG
          DXDZ(1,2) = -Z1*SAG/2
          DXDZ(2,1) = -Z2*CAG/2
          DXDZ(2,2) = -Z2*SAG-Z1*CAG/2
          DXDZ(:,3) = 0D0

          DYDZ(1,1) = Z1*SAG - Z2*CAG/2
          DYDZ(1,2) = -Z1*CAG/2
          DYDZ(2,1) = Z2*SAG/2
          DYDZ(2,2) = -Z2*CAG+Z1*SAG/2
          DYDZ(:,3) = 0D0
       ELSE
          ! derivative wrt normalized z axis
          DXDZ(1,1) = (-Z2*SAG+Z23*SAG+2*Z13*CAG-2*Z1*MAT(1,1))/NZ
          DXDZ(1,2) = (2*Z2*CAG-Z1*SAG+Z13*SAG-2*Z2*MAT(1,1))/NZ
          DXDZ(1,3) = (Z12*SAG+Z11*CAG)/NZ

          DXDZ(2,1) = (-Z2*CAG+2*Z1*SAG+Z23*CAG - 2*MAT(2,1)*Z1)/NZ
          DXDZ(2,2) = (-Z1*CAG+2*Z23*SAG + Z13*CAG-2*MAT(2,1)*Z2)/NZ
          DXDZ(2,3) = (Z22*SAG+Z12*CAG)/NZ

          DYDZ(1,1) = (-Z2*CAG+Z23*CAG-2*Z13*SAG-2*Z1*MAT(1,2))/NZ
          DYDZ(1,2) = (-2*Z2*SAG-Z1*CAG+Z13*CAG-2*Z2*MAT(1,2))/NZ
          DYDZ(1,3) = (Z12*CAG-Z11*SAG)/NZ

          DYDZ(2,1) = (Z2*SAG+2*Z1*CAG-Z23*SAG-2*Z1*MAT(2,2))/NZ
          DYDZ(2,2) = (Z1*SAG+2*Z23*CAG-Z13*SAG-2*Z2*MAT(2,2))/NZ
          DYDZ(2,3) = (Z22*CAG-Z12*SAG)/NZ
       ENDIF
       DXDZ(3,:) = (/-CAG,-SAG,0D0/)
       DYDZ(3,:) = (/SAG,-CAG,0D0/)

       ! derivatives of new axes

       CALL DGEMM('N','N',3,3,3,1D0,DXDZ,3,DMATZ(:,3,:),3,0D0,DMATZ(:,1,:),3)
       CALL DGEMM('N','N',3,3,3,1D0,DYDZ,3,DMATZ(:,3,:),3,0D0,DMATZ(:,2,:),3)


       !CALL DGEMV('N',3,3,1D0,DZ,3,DMATZTMP(2,1,:),1,0D0,DMATZ(2,1,:),1)
       !DMATZ(2,1,:) = DXDZ
    ENDIF

    IF (PRESENT(DMATAG)) THEN
       DMATAG = 0D0

       IF (NZ.LT.NZTINY) THEN
          DMATAG(1,1) = -SAG-(Z12*CAG-Z11*SAG)/2
          DMATAG(2,1) = CAG-(Z22*CAG-Z12*SAG)/2

          DMATAG(1,2) = -CAG+(Z11*CAG+Z12*SAG)/2
          DMATAG(2,2) = -SAG+(Z22*SAG+Z12*CAG)/2
       ELSE
          DMATAG(1,1) = (-Z22*SAG-Z1*Z2*CAG+Z1*Z23*CAG-Z3*Z11*SAG)/NZ
          DMATAG(2,1) = (Z12*SAG+Z11*CAG + Z3*Z22*CAG-Z1*Z23*SAG)/NZ

          DMATAG(1,2) = (-Z22*CAG+Z12*SAG-Z1*Z23*SAG-Z3*Z11*CAG)/NZ
          DMATAG(2,2) = (Z12*CAG-Z11*SAG-Z22*Z3*SAG-Z1*Z23*CAG)/NZ
       ENDIF
       DMATAG(3,1) = -Z2*CAG+Z1*SAG ! dX3/dA
       DMATAG(3,2) = SAG*Z2+CAG*Z1! dY3/dA
    ENDIF
  END SUBROUTINE COORDS2ROTMAT

  ! ----------------- STUFF INVOLVING TREATING QUATERNIONS AS 4-VECTORS -----
  FUNCTION QUAT2QV(Q)
    ! convert a quaternion object to a 4-vector
    IMPLICIT NONE
    DOUBLE PRECISION :: QUAT2QV(4)
    TYPE(QUATERNION), INTENT(IN) :: Q

    QUAT2QV = (/Q%W,Q%X,Q%Y,Q%Z/)
  END FUNCTION QUAT2QV

  SUBROUTINE ROTQV(THETA,AX,QV,DT)
    ! get the quaternion corresponding to rotation around axis AX by angle THETA
    ! as a 4-vector (in QV); optionally, also get the derivative wrt theta
    ! WARNING: AX assumed to be normalized; does not check for this!

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: THETA, AX(3)
    DOUBLE PRECISION, INTENT(OUT) :: QV(4)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DT(4)
    DOUBLE PRECISION :: CT, ST

    CT = COS(THETA/2); ST = SIN(THETA/2);
    QV(1) = CT
    QV(2:4) = ST*AX

    IF (PRESENT(DT)) THEN
       DT(1) = -ST/2
       DT(2:4) = CT/2*AX
    ENDIF

  END SUBROUTINE ROTQV

  SUBROUTINE QVPTMULT(Q,PT,QP,DQ,DPT)
    ! multiply a quaternion by a point in 3-space
    ! returns the result in QP
    ! optionally, returns derivatives wrt the quaternion components in dQ
    ! or wrt point components in DPT
    ! WARNING: no normalization happens here
    ! WARNING: this is  a pretty inefficient way to do things; should fix at some point without resorting to rotation matrices
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Q(4),PT(3)
    DOUBLE PRECISION, INTENT(OUT) :: QP(3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DQ(3,4),DPT(3,3)
    DOUBLE PRECISION :: QN, QINV(4),QTMP(4),QTMP2(4)
    DOUBLE PRECISION :: MAT(3,3),DMAT(3,3,4)
    INTEGER :: I,J
    TYPE(QUATERNION) :: QUAT

    QN = DOT_PRODUCT(Q,Q)
    QUAT%W = Q(1); QUAT%X = Q(2); QUAT%Y = Q(3); QUAT%Z = Q(4)
    IF (PRESENT(DQ)) THEN
       CALL QUAT2ROTMAT(QUAT,MAT,dMAT)
    ELSE
       CALL QUAT2ROTMAT(QUAT,MAT)
    ENDIF
    DO I = 1,3
       QP(I) = DOT_PRODUCT(MAT(I,:),PT)
       IF (PRESENT(DQ)) THEN
          DO J = 1,4
             DQ(I,J) = DOT_PRODUCT(dMAT(I,:,J),PT)
          ENDDO
       ENDIF
    ENDDO

    QP = QP/QN
    IF (PRESENT(DQ)) THEN
       DO I = 1,3
          DO J = 1,4
             DQ(I,J) = (DQ(I,J) - 2*Q(J)*QP(I))/QN
          ENDDO
       ENDDO
    ENDIF

    IF (PRESENT(DPT)) THEN
       DPT =MAT/QN
    END IF
  END SUBROUTINE QVPTMULT

  ! ----------------- STUFF INVOLVING QUATERNION OBJECTS -------------------

  TYPE(QUATERNION) FUNCTION RHQINTERP(Q1,Q2,F)
    ! interpolate from one quaternion to another always in a right-handed sense
    ! (relative to the Q1 z-axis)
    ! F should be between 0 and 1
    TYPE(QUATERNION) :: Q1, Q2
    DOUBLE PRECISION :: F
    DOUBLE PRECISION :: ANG, AX(3), DIR
    TYPE(QUATERNION) :: QREL

    ! relative quaternion for Q2 relative to Q1
    QREL = INVQUAT(Q1)*Q2

    ! angle of rotation
    ANG = ACOS(QREL%W)*2

    ! axis of rotation
    AX = (/QREL%X,QREL%Y,QREL%Z/); AX = AX/SQRT(DOT_PRODUCT(AX,AX))

    DIR = DOT_PRODUCT(AX,QUAT2PT(Q1*PTQUAT((/0D0,0D0,1D0/))/Q1))

    print*, 'testx2:', dir, ang, ax

    IF (DIR.LT.0) THEN
       ! rotation axis points away from quaternion axis, so flip it
       ANG = -ANG
       AX = -AX
    ENDIF

    QREL = ROTQUAT(ANG*F,AX)

    RHQINTERP = Q1*QREL
  END FUNCTION RHQINTERP

  TYPE(QUATERNION) FUNCTION QSLERP(Q1,Q2,F)
    ! Sphreical linear interpolation between two quaternions
    ! F should be between 0 and 1
    TYPE(QUATERNION), INTENT(IN) :: Q1, Q2
    DOUBLE PRECISION, INTENT(IN) :: F
    DOUBLE PRECISION :: QV1(4), QV2(4), ANG, QVANS(4)

    QV1 = (/Q1%W,Q1%X,Q1%Y,Q1%Z/)
    QV2 = (/Q2%W, Q2%X, Q2%Y, Q2%Z/)

    ! angle between them
    ANG = ACOS(DOT_PRODUCT(QV1,QV2))

    QVANS = SIN((1-F)*ANG)/SIN(ANG)*QV1 + SIN(F*ANG)/SIN(ANG)*QV2

    QSLERP%W = QVANS(1); QSLERP%X = QVANS(2); QSLERP%Y = QVANS(3); QSLERP%Z = QVANS(4);
  END FUNCTION QSLERP

  SUBROUTINE QUAT2ROTMAT(Q,MAT,DMAT)
    ! convert a quaternion object to a rotation matrix and optionally return derivatives
    ! NOTE: no normalization
    IMPLICIT NONE
    TYPE(QUATERNION) :: Q
    DOUBLE PRECISION, INTENT(OUT) :: MAT(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DMAT(3,3,4)
    DOUBLE PRECISION :: A,B,C,D,AA,BB,CC,DD,AB,AC,AD,BC,BD,CD

    A = Q%W; B = Q%X; C = Q%Y; D = Q%Z
    AA = A*A; BB = B*B; CC = C*C; DD = D*D
    AB = 2*A*B; AC = 2*A*C; AD = 2*A*D
    BC = 2*B*C; BD = 2*B*D
    CD = 2*C*D

    MAT(1,:) = (/AA+BB-CC-DD,BC-AD,AC+BD/)
    MAT(2,:) = (/AD+BC,AA-BB+CC-DD,CD-AB/)
    MAT(3,:) = (/BD-AC,AB+CD,AA-BB-CC+DD/)

    IF (PRESENT(DMAT)) THEN
       dMAT(1,:,1) = (/A,-D,C/)
       dMAT(1,:,2) = (/B,C,D/)
       dMAT(1,:,3) = (/-C,B,A/)
       dMAT(1,:,4) = (/-D,-A,B/)
       dMAT(2,:,1) = (/D,A,-B/)
       dMAT(2,:,2) = (/C,-B,-A/)
       dMAT(2,:,3) = (/B,C,D/)
       dMAT(2,:,4) = (/A,-D,C/)
       dMAT(3,:,1) = (/-C,B,A/)
       dMAT(3,:,2) = (/D,A,-B/)
       dMAT(3,:,3) = (/-A,D,-C/)
       dMAT(3,:,4) = (/B,C,D/)
       dMAT = dMAT*2
    ENDIF
  END SUBROUTINE QUAT2ROTMAT

  TYPE(QUATERNION) FUNCTION ROTMAT2QUAT(R)
    ! convert from a rotation matrix to a quaternion object
    ! following Deibel, 2006 (but with the rotation matrix transposed)
    ! R(:,3) has the z axis after rotation is applied, etc.
    ! assumes R is orthonormal
    IMPLICIT NONE
    DOUBLE PRECISION :: R(3,3)
    DOUBLE PRECISION :: R11, R22,R33, TMP

    R11 = R(1,1); R22 = R(2,2); R33 = R(3,3)

    IF (R22 .GE.-R33.AND.R11.GE.-R22.AND.R11.GE.-R33) THEN
       TMP = SQRT(1+R11+R22+R33)
       ROTMAT2QUAT%W = TMP/2
       ROTMAT2QUAT%X = (R(3,2)-R(2,3))/(TMP*2)
       ROTMAT2QUAT%Y = (R(1,3)-R(3,1))/(TMP*2)
       ROTMAT2QUAT%Z = (R(2,1)-R(1,2))/(TMP*2)
       !PRINT*, 'TESTX1Q'
    ELSEIF (R22.LE.-R33.AND.R11.GT.R22.AND.R11.GT.R33) THEN
       TMP = SQRT(1+R11-R22-R33)
       ROTMAT2QUAT%W = (R(3,2)-R(2,3))/(TMP*2)
       ROTMAT2QUAT%X = TMP/2
       ROTMAT2QUAT%Y = (R(2,1)+R(1,2))/(TMP*2)
       ROTMAT2QUAT%Z = (R(1,3)+R(3,1))/(TMP*2)
       !PRINT*, 'TESTX2Q'
    ELSEIF (R22.GT.R33.AND.R11.LT.R22.AND.R11.LE.-R33) THEN
       TMP = SQRT(1-R11+R22-R33)
       ROTMAT2QUAT%W = (R(1,3)-R(3,1))/(TMP*2)
       ROTMAT2QUAT%X = (R(2,1)+R(1,2))/(TMP*2)
       ROTMAT2QUAT%Y = TMP/2
       ROTMAT2QUAT%Z = (R(3,2)+R(2,3))/(TMP*2)
       !PRINT*, 'TESTX3Q'
    ELSEIF (R22.LT.R33.AND.R11.LE.-R22.AND.R11.LT.R33) THEN
       TMP = SQRT(1D0-R11-R22+R33)
       ROTMAT2QUAT%W = (R(2,1)-R(1,2))/(TMP*2)
       ROTMAT2QUAT%X = (R(1,3)+R(3,1))/(TMP*2)
       ROTMAT2QUAT%Y = (R(3,2)+R(2,3))/(TMP*2)
       ROTMAT2QUAT%Z = TMP/2
       !PRINT*, 'TESTX4Q'
    ELSE
       PRINT*, 'ERROR IN ROTMAT2QUAT: bad rotation matrix'
       PRINT*, R(1,:)
       PRINT*, R(2,:)
       PRINT*, R(3,:)
       STOP 1
       ROTMAT2QUAT%W = 0; ROTMAT2QUAT%X = 0; ROTMAT2QUAT%Y = 0; ROTMAT2QUAT%Z = 0
    ENDIF
  END FUNCTION ROTMAT2QUAT

  SUBROUTINE EULER2QUAT(EUL,Q,DERV,GETDERV)
    ! convert from Euler angles (z-x-z convention) to a quaternion
    ! if GETDERV is true, also get the derivatives of the quaternion components
    ! with respect to the euler angles (4 rows by 3 columns)
    ! WARNING: since the quaternion representation has more parameters, switching from quaternion
    ! to euler and back again will not always give the exact same quaternion, though it
    ! will give an equivalent one (eg: may invert all components)

    DOUBLE PRECISION, INTENT(IN) :: EUL(3)
    TYPE(QUATERNION), INTENT(OUT) :: Q
    DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: DERV(4,3)
    LOGICAL, INTENT(IN),OPTIONAL :: GETDERV
    DOUBLE PRECISION :: CA,SA,CB,SB,CG,SG,ALPHA,GAMMA,BETA
    LOGICAL :: FLIPBETA

    BETA = ANGLE2PI(EUL(2))
    FLIPBETA = BETA.GT.PI
    ALPHA = EUL(1); GAMMA = EUL(3)
    IF (BETA.GT.PI) THEN
       BETA = 2*PI-BETA
       GAMMA = GAMMA+PI
       ALPHA = ALPHA+PI
    ENDIF

    ALPHA = ANGLE2PI(ALPHA)
    GAMMA = ANGLE2PI(GAMMA)
    !IF (QUATTEST) THEN
    !   PRINT*, 'TESTXE:', EUL
    !   PRINT*, 'A,B,G:',ALPHA,BETA,GAMMA
    !ENDIF

    CA = COS(ALPHA/2); SA = SIN(ALPHA/2)
    CB = COS(BETA/2); SB = SIN(BETA/2)
    CG = COS(GAMMA/2); SG = SIN(GAMMA/2)

    Q%W = CG*CB*CA - SG*CB*SA
    Q%X = SG*SB*SA + CG*SB*CA
    Q%Y = CG*SB*SA - SG*SB*CA
    Q%Z = CG*CB*SA + SG*CB*CA

    IF (PRESENT(DERV)) THEN
       IF (GETDERV) THEN
          DERV(1,:) = 0.5D0*(/-CG*CB*SA - SG*CB*CA, -CG*SB*CA + SG*SB*SA, -SG*CB*CA - CG*CB*SA/)
          DERV(2,:) = 0.5D0*(/SG*SB*CA - CG*SB*SA, SG*CB*SA + CG*CB*CA, CG*SB*SA - SG*SB*CA /)
          DERV(3,:) = 0.5D0*(/CG*SB*CA + SG*SB*SA, CG*CB*SA - SG*CB*CA, -SG*SB*SA - CG*SB*CA/)
          DERV(4,:) = 0.5D0*(/CG*CB*CA - SG*CB*SA, -CG*SB*SA - SG*SB*CA, -SG*CB*SA + CG*CB*CA/)
       ENDIF
       IF (FLIPBETA) THEN
          DERV(:,2) = -DERV(:,2)
       ENDIF
    ENDIF

  END SUBROUTINE EULER2QUAT

  SUBROUTINE QUAT2EULER(Q,EUL)
    ! get the Euler angles (z-x-z convention) corresponding to a unit quaternion
    ! NOTE: the quaternion must already be normalized
    ! currently can't handle gimbal lock
    TYPE(QUATERNION), INTENT(IN) :: Q
    DOUBLE PRECISION, INTENT(OUT) :: EUL(3)
    DOUBLE PRECISION :: DUMMY

    EUL(2) = 1 - 2*(Q%X**2+Q%Y**2)
    IF (EUL(2).GT.1) THEN
       EUL(2) = 0D0
    ELSEIF (EUL(2).LT.-1) THEN
       EUL(2) = PI
    ELSE
       EUL(2) = ACOS(EUL(2))
    ENDIF

    ! deal with the gimbal lock issues
    IF (EUL(2).LE.EPSILON(0D0)) THEN
       EUL(3) = ATAN2(2*(Q%X*Q%Y + Q%W*Q%Z), 1 - 2*(Q%Y**2+Q%Z**2))
       IF (EUL(3).LT.0) EUL(3) = 2*PI+EUL(3)
       EUL(1) = 0D0
       RETURN
    ELSEIF (EUL(2).GE.PI-EPSILON(0D0)) THEN
       EUL(3) = ATAN2(-2*(Q%X*Q%Y + Q%W*Q%Z), 1 - 2*(Q%Y**2+Q%Z**2))
       IF (EUL(3).LT.0) EUL(3) = 2*PI+EUL(3)
       EUL(1) = 0D0
       RETURN
    ENDIF

    DUMMY = Q%W*Q%X-Q%Y*Q%Z

    EUL(1) = ATAN2(Q%W*Q%Y+Q%X*Q%Z, DUMMY)
    IF (EUL(1).LT.0) EUL(1) = 2*PI+EUL(1)

    DUMMY = Q%W*Q%X+Q%Y*Q%Z

    EUL(3) = ATAN2(Q%X*Q%Z-Q%W*Q%Y,DUMMY)
    IF (EUL(3).LT.0) EUL(3) = 2*PI+EUL(3)

  END SUBROUTINE QUAT2EULER

  FUNCTION QUAT2PT(Q)
    ! get the point corresponding to a quaternion
    TYPE(QUATERNION) :: Q
    DOUBLE PRECISION :: QUAT2PT(3)

    QUAT2PT = (/Q%X,Q%Y,Q%Z/)
  END FUNCTION QUAT2PT

  TYPE(QUATERNION) FUNCTION ROTQUAT(THETA,AX)
    ! get the quaternion corresponding to rotation around unit axis AX
    ! by an angle theta
    DOUBLE PRECISION :: THETA, AX(3)
    DOUBLE PRECISION :: T2, ST2

    IF (ABS(SUM(AX**2)-1D0).GT.EPSILON(1d0)*10) THEN
       print*, 'ERROR in ROTQUAT: Axis does not have unit norm'
       STOP 1
    ENDIF

    T2 = THETA/2; ST2 = SIN(THETA/2)
    ROTQUAT%W = COS(T2)
    ROTQUAT%X = ST2*AX(1); ROTQUAT%Y = ST2*AX(2); ROTQUAT%Z=ST2*AX(3)
  END FUNCTION ROTQUAT

  TYPE(QUATERNION) FUNCTION PTQUAT(P)
    ! turn a 3d point into a quaternion
    DOUBLE PRECISION :: P(3)

    PTQUAT%W = 0D0; PTQUAT%X = P(1); PTQUAT%Y = P(2); PTQUAT%Z = P(3)
  END FUNCTION PTQUAT

  TYPE(QUATERNION) FUNCTION INVQUAT(Q)
    ! get the inverse of a quaternion
    IMPLICIT NONE
    TYPE(QUATERNION), INTENT(IN) :: Q
    DOUBLE PRECISION :: QN

    QN = Q%W**2+Q%X**2+Q%Y**2+Q%Z**2
    INVQUAT%W = Q%W/QN; INVQUAT%X = -Q%X/QN; INVQUAT%Y=-Q%Y/QN; INVQUAT%Z=-Q%Z/QN

  END FUNCTION INVQUAT

  TYPE(QUATERNION) FUNCTION QDIVIDE(P,Q)
    ! multiply P by inverse of Q(in that order)
    IMPLICIT NONE
    TYPE(QUATERNION), INTENT(IN) :: P,Q
    TYPE(QUATERNION) :: QINV
    DOUBLE PRECISION :: QN

    ! inverse of the 2nd quaternion
    QN = Q%W**2+Q%X**2+Q%Y**2+Q%Z**2
    QINV%W = Q%W/QN; QINV%X = -Q%X/QN; QINV%Y=-Q%Y/QN; QINV%Z=-Q%Z/QN

    QDIVIDE = P*QINV
  END FUNCTION QDIVIDE

  TYPE(QUATERNION) FUNCTION QPRODUCT(P,Q)
    ! quaternion multiplication
    IMPLICIT NONE
    TYPE(QUATERNION), INTENT(IN) :: P, Q

    QPRODUCT%W = P%W*Q%W - P%X*Q%X - P%Y*Q%Y - P%Z*Q%Z
    QPRODUCT%X = P%W*Q%X + P%X*Q%W + P%Y*Q%Z - P%Z*Q%Y
    QPRODUCT%Y = P%W*Q%Y - P%X*Q%Z + P%Y*Q%W + P%Z*Q%X
    QPRODUCT%Z = P%W*Q%Z + P%X*Q%Y - P%Y*Q%X + P%Z*Q%W

  END FUNCTION QPRODUCT

  ! --------- general angle and euler angle stuff -------------
  SUBROUTINE ROTANGAX(ANG,AX,INVEC,OUTVEC,CALCROTMAT,ROTMAT)
    ! rotate a 3D vector by angle ANG around axis AX
    ! if CALCROTMAT is true, recalculate the rotation matrix
    ! otherwise use the provided one
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: ANG, AX(3), INVEC(3)
    DOUBLE PRECISION, INTENT(OUT) :: OUTVEC(3)
    DOUBLE PRECISION, INTENT(INOUT) :: ROTMAT(3,3)
    LOGICAL, INTENT(IN)  :: CALCROTMAT
    DOUBLE PRECISION :: CT,ST,CT1
    INTEGER :: I

    IF (CALCROTMAT) THEN
       CT = COS(ANG); ST = SIN(ANG)
       CT1 = 1-CT
       ROTMAT(1,:) = (/CT + AX(1)**2*CT1, AX(1)*AX(2)*CT1-AX(3)*ST, AX(1)*AX(3)*CT1+AX(2)*ST/)
       ROTMAT(2,:) = (/AX(2)*AX(1)*CT1+AX(3)*ST,CT+AX(2)**2*CT1,AX(2)*AX(3)*CT1-AX(1)*ST/)
       ROTMAT(3,:) = (/AX(3)*AX(1)*CT1-AX(2)*ST,AX(3)*AX(2)*CT1+AX(1)*ST, CT+AX(3)**2*CT1/)
    ENDIF

    DO I = 1,3
       OUTVEC(I) = DOT_PRODUCT(ROTMAT(I,:),INVEC)
    ENDDO
  END SUBROUTINE ROTANGAX

  SUBROUTINE EUL2ROTMAT(EUL,ROTMAT,DMAT)
    ! get the rotation matrix corresponding to various euler angles
    ! and the appropriate derivatives if DMAT is present
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: EUL(3)
    DOUBLE PRECISION, INTENT(OUT) :: ROTMAT(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DMAT(3,3,3)
    DOUBLE PRECISION :: CA,SA,CB,SB,CG,SG

    CA = COS(EUL(1)); SA = SIN(EUL(1))
    CB = COS(EUL(2)); SB = SIN(EUL(2))
    CG = COS(EUL(3)); SG = SIN(EUL(3))

    ROTMAT(1,:) = (/CA*CG-SA*CB*SG, -CA*SG-SA*CB*CG, SB*SA/)
    ROTMAT(2,:) = (/SA*CG+CA*CB*SG,-SA*SG+CA*CB*CG,-SB*CA/)
    ROTMAT(3,:) = (/SB*SG,SB*CG,CB/)

    IF (PRESENT(DMAT)) THEN
       dMAT(1,:,1) = (/-SA*CG-CA*CB*SG,SA*SG-CA*CB*CG,SB*CA/)
       dMAT(2,:,1) = (/CA*CG-SA*CB*SG,-CA*SG-SA*CB*CG,SB*SA/)
       dMAT(3,:,1) = 0D0

       dMAT(1,:,2) = (/SA*SB*SG,SA*SB*CG,CB*SA/)
       dMAT(2,:,2) = (/-CA*SB*SG,-CA*SB*CG,-CB*CA/)
       dMAT(3,:,2) = (/CB*SG,CB*CG,-SB/)

       dMAT(1,:,3) = (/-CA*SG-SA*CB*CG,-CA*CG+SA*CB*SG,0D0/)
       dMAT(2,:,3) = (/-SA*SG+CA*CB*CG,-SA*CG-CA*CB*SG,0D0/)
       dMAT(3,:,3) = (/SB*CG,-SB*SG,0D0/)
    ENDIF
  END SUBROUTINE EUL2ROTMAT

  SUBROUTINE GETANGLE(IJ,KJ,CST,dCTdIJ,dCTdKJ)
    ! get the angle between three points (I-J-K)
    ! and, optionally, the derivative of that angle
    ! actually this returns the COSINE of the angle and its derivative
    ! IJ = I-J; KJ = K-J
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: IJ(3),KJ(3)
    DOUBLE PRECISION, INTENT(OUT) :: CST
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dCTdIJ(3), dCTdKJ(3)

    DOUBLE PRECISION :: DXI, DYI, DZI, DXJ, DYJ, DZJ
    DOUBLE PRECISION :: RI2, RJ2, RI, RJ, RIR, RJR
    DOUBLE PRECISION :: DXIR,DYIR, DZIR,DXJR,DYJR,DZJR

    DXI=IJ(1); DYI = IJ(2); DZI = IJ(3)
    DXJ = KJ(1); DYJ = KJ(2); DZJ = KJ(3)

    RI2=DXI*DXI+DYI*DYI+DZI*DZI
    RJ2=DXJ*DXJ+DYJ*DYJ+DZJ*DZJ
    RI=SQRT(RI2)
    RJ=SQRT(RJ2)
    RIR=1/RI
    RJR=1/RJ
    DXIR=DXI*RIR
    DYIR=DYI*RIR
    DZIR=DZI*RIR
    DXJR=DXJ*RJR
    DYJR=DYJ*RJR
    DZJR=DZJ*RJR
    CST=DXIR*DXJR+DYIR*DYJR+DZIR*DZJR
    IF (PRESENT(DCTDIJ)) THEN
       dCTdIJ(1)=-(DXIR*CST-DXJR)*RIR
       dCTdIJ(2)=-(DYIR*CST-DYJR)*RIR
       dCTdIJ(3)=-(DZIR*CST-DZJR)*RIR
    ENDIF
    IF (PRESENT(DCTDKJ)) THEN
       dCTdKJ(1)=-(DXJR*CST-DXIR)*RJR
       dCTdKJ(2)=-(DYJR*CST-DYIR)*RJR
       dCTdKJ(3)=-(DZJR*CST-DZIR)*RJR
    ENDIF

  END SUBROUTINE GETANGLE

  SUBROUTINE GETDIHEDRAL(IJ,JK,LK, PHI, dPdIJ, dPdJK, dPdLK)
    ! dihedral angle for 4 atoms I, J, K, L bound in order
    ! IJ = I-J; JK = J-K; LK = L-K
    ! find the dihedral torsion angle; return it in PHI
    ! Also return all derivatives: dP/dIJx, dP/dIJy, dP/dIJz in triplet dPdIJ
    ! same with dPdJK, dPdLK
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) ::  IJ(3), JK(3), LK(3)
    DOUBLE PRECISION, INTENT(OUT) :: PHI
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dPdIJ(3), dPdJK(3), dPdLK(3)
    DOUBLE PRECISION :: DPDJ(3),DPDI(3),DPDK(3),DPDL(3)
    DOUBLE PRECISION :: FX, FY, FZ, GX, GY, GZ, HX, HY, HZ
    DOUBLE PRECISION :: AX, AY, AZ, BX, BY, Bz
    DOUBLE PRECISION :: RF, RG, RH, RF2, RG2, RH2, RFR, RGR, RHR
    DOUBLE PRECISION :: CSTTWO, SNTTWO2, CSTTHREE, SNTTHREE2, SNTTWO2R, SNTTHREE2R
    DOUBLE PRECISION :: RA2, RB2, RA2R, RB2R, RABR, CP
    DOUBLE PRECISION :: MYTX, MYTY, MYTZ, MYSCALAR
    DOUBLE PRECISION :: DUMMY, DUMMY2
    DOUBLE PRECISION :: B1(3), B2(3), B3(3), B12(3), B23(3)
    LOGICAL :: NOCOOR = .FALSE.


    FX=IJ(1)
    FY=IJ(2)
    FZ=IJ(3)
    GX=JK(1)
    GY=JK(2)
    GZ=JK(3)
    HX=LK(1)
    HY=LK(2)
    HZ=LK(3)
    ! A=F x G, B=H x G
    AX=FY*GZ-FZ*GY
    AY=FZ*GX-FX*GZ
    AZ=FX*GY-FY*GX
    BX=HY*GZ-HZ*GY
    BY=HZ*GX-HX*GZ
    BZ=HX*GY-HY*GX
    ! RG=|G|, RGR=1/|G|
    RG2=GX*GX+GY*GY+GZ*GZ
    RG=SQRT(RG2)
    RGR=1/RG
    ! dae for use in evaluating B-matrix
    RF2=FX*FX+FY*FY+FZ*FZ
    RF=SQRT(RF2)
    RFR=1/RF
    RH2=HX*HX+HY*HY+HZ*HZ
    RH=SQRT(RH2)
    RHR=1/RH


    CSTTWO=-(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
    SNTTWO2=1-CSTTWO*CSTTWO
    SNTTWO2R=1/SNTTWO2
    CSTTHREE=(HX*GX+HY*GY+HZ*GZ)*RHR*RGR
    SNTTHREE2=1-CSTTHREE*CSTTHREE
    SNTTHREE2R=1/SNTTHREE2

    RA2=AX*AX+AY*AY+AZ*AZ
    RB2=BX*BX+BY*BY+BZ*BZ
    RA2R=1/RA2
    RB2R=1/RB2
    RABR=SQRT(RA2R*RB2R)

    PHI =  ATAN2(-RG*(FX*BX+FY*BY+FZ*BZ),AX*BX+AY*BY+AZ*BZ)


    IF (PRESENT(DPDIJ).OR.PRESENT(DPDJK).OR.PRESENT(DPDLK)) THEN
       DUMMY=RFR*RFR*RGR*SNTTWO2R
       dPdI = (/-AX*DUMMY, -AY*DUMMY, -AZ*DUMMY/)
       DUMMY=RFR*RFR*RGR*RGR*SNTTWO2R*(RG-RF*CSTTWO)
       DUMMY2=RHR*RGR*RGR*SNTTHREE2R*CSTTHREE
       dPdJ = (/AX*DUMMY-BX*DUMMY2, AY*DUMMY-BY*DUMMY2, AZ*DUMMY-BZ*DUMMY2/)
       DUMMY=RHR*RHR*RGR*SNTTHREE2R
       dPdL = (/BX*DUMMY,BY*DUMMY,BZ*DUMMY/)
    ENDIF
    IF (PRESENT(DPDIJ)) DPDIJ = DPDI
    IF (PRESENT(DPDLK)) DPDLK = DPDL
    IF (PRESENT(DPDJK)) DPDJK = DPDI + DPDJ


  END SUBROUTINE GETDIHEDRAL

  DOUBLE PRECISION FUNCTION ANGLE0(ANGLE)
    ! convert an angle to one between +/- pi
    ! (so keeping it as close as possible to zero)
    IMPLICIT NONE
    DOUBLE PRECISION :: ANGLE

    ANGLE0 = ANGLE2PI(ANGLE)
    IF (ANGLE0.GT.PI) THEN
       ANGLE0 = ANGLE0 - 2*PI
    ENDIF
  END FUNCTION ANGLE0

  DOUBLE PRECISION FUNCTION ANGLE2PI(ANGLE)
    ! convert an angle to one between 0 and 2pi by adding or subtracting multiples of 2pi
    IMPLICIT NONE
    DOUBLE PRECISION :: ANGLE
    INTEGER :: N2PI

    N2PI = INT(ANGLE/(2*PI))
    IF (ANGLE.LT.0) THEN
       ANGLE2PI = ANGLE+(-N2PI+1)*2*PI
    ELSE
       ANGLE2PI = ANGLE - N2PI*2*PI
    ENDIF
  END FUNCTION ANGLE2PI

END MODULE QUATUTIL
