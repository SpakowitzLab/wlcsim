module QUATUTIL
  ! utilities for dealing with quaternions
  ! and other representations of rotation
  ! including euler angles, rotation matrices, and alpha + gamma and z-axis-vector representations

  implicit none
  LOGICAL :: TESTQUAT = .FALSE.
  real(dp), PARAMETER :: PI = 3.141592653589793D0
  ! when z1^2 + z2^2 goes below tiny when working with alpha + gamma & zvec coordinates
  ! use the v->0 approximation
  real(dp), PARAMETER :: NZTinY = 1D-14

  ! definition of a quaternion class
  TYPE QUATERNION
     real(dp) :: W, X, Y, Z
  END TYPE QUATERNION

  inTERFACE OPERATOR (*)
     module PROCEDURE QPRODUCT
  END inTERFACE

  inTERFACE OPERATOR (/)
     module PROCEDURE QDIVIDE
  END inTERFACE

contains
  subroutine QUAT2SCREW(QUAT,TRANS,HELCRD)
    ! convert from a rotation + translation to screw coordinates
    ! (equivalently to overall helix coordinates)
    ! returns height, angle, radius, orientation of
    ! canonical system relative to helix system (3 euler angles)
    TYPE(QUATERNION), intent(in) :: QUAT
    real(dp), intent(in) :: TRANS(3)
    real(dp), intent(out) :: HELCRD(6)
    real(dp) :: THETA, AX(3), TP(3), ST2, CT2, AXT3, A, B, G, NP, NTP

    ! angle-axis representation from quaternion
    THETA = 2*ACOS(QUAT%W)

    if (THETA == 0D0) then ! translation only
       HELCRD(1) = SQRT(doT_PRODUCT(TRANS,TRANS))
       HELCRD(2) = 0D0
       HELCRD(3) = 0D0
       HELCRD(4:6) = (/0D0,0D0,0D0/)
       RETURN
    ENDif

    ! angle of rotation around screw axis
    HELCRD(2) = THETA

    AX = (/QUAT%X,QUAT%Y,QUAT%Z/);
    ST2 = Sin(THETA/2)
    AX = AX/ST2; ! normalize the axis

    ! shift along screw axis
    HELCRD(1) = doT_PRODUCT(TRANS,AX)

    ! translation in plane perpendicular to axis
    TP = TRANS - HELCRD(1)*AX
    NTP = SQRT(doT_PRODUCT(TP,TP))

    ! radius of screw
    NP = NTP/(2*ST2)
    HELCRD(3) = NP

    ! cross-product of AX x TP
    AXT3 = AX(1)*TP(2) - AX(2)*TP(1)

    ! euler angles of screw axis coord system relative to canonical
    A = atan2(AX(1),-AX(2))
    B = ACOS(AX(3))
    CT2 = QUAT%W/ST2
    G = ATAN2(-TP(3)-CT2*AXT3,-AXT3 + CT2*TP(3))

    ! euler angles of canonical system relative to screw axis
    HELCRD(4:6) = (/PI-G,B,PI-A/)
  END subroutine QUAT2SCREW

  subroutine COORDS2QUAT(AG,ZVEC,Q)
    ! convert from an alpha + gamma angle and a vector along the z axis
    ! to a normalized quaternion
    ! note: this doesn't work if beta = pi
    implicit none
    real(dp), intent(in) :: AG, ZVEC(3)
    TYPE(QUATERNION) :: Q
    real(dp) :: NZ, ZAX(3), ALPHA,BETA,GAMMA

    NZ = SQRT(doT_PRODUCT(ZVEC,ZVEC))
    ZAX = ZVEC/NZ

    NZ = ZAX(1)**2 + ZAX(2)**2

    if (NZ == 0) then
       CALL EULER2QUAT((/AG,0D0,0D0/),Q)
    else
       ALPHA = ATAN2(ZAX(1),-ZAX(2)); GAMMA = AG-ALPHA
       BETA = ACOS(ZAX(3))
       CALL EULER2QUAT((/ALPHA,BETA,GAMMA/),Q)
    ENDif

  END subroutine COORDS2QUAT

  subroutine COORDS2ROTMAT(AG,ZVEC,MAT,DMATAG,DMATZ)
    ! convert from coordinates that include the alpha + gamma euler angle
    ! and a non-normalized vector along the Z axis in canonical reference system
    ! to a rotation matrix
    ! if DMATAG and DMATZ are provided, also get derivatives of all the matrix
    ! components wrt the coordinates

    implicit none
    real(dp), intent(in) :: AG, ZVEC(3)
    real(dp), intent(out) :: MAT(3,3)
    real(dp), intent(out), OPTIONAL :: DMATAG(3,3), DMATZ(3,3,3)
    real(dp) :: NZTOT2, NZTOT, NZ, Z1, Z2, Z3
    real(dp) :: Z11,Z22,Z12,Z13,Z23,CAG,SAG
    real(dp) :: DXDZ(3,3),DYDZ(3,3)

    if ((PRESENT(DMATAG).AND..NOT.PRESENT(DMATZ)).OR.&
         & (.NOT.PRESENT(DMATAG).AND.PRESENT(DMATZ))) then
       PRinT*, 'ERROR in COORDS2ROTMAT: neither or both of DMATAG and DMATZ must be provided'
       stop 1
    ENDif

    MAT = 0D0

    ! normalize the zaxis
    NZTOT2 = doT_PRODUCT(ZVEC,ZVEC); NZTOT = SQRT(NZTOT2)
    MAT(:,3) = ZVEC/NZTOT
    Z1 = MAT(1,3); Z2 = MAT(2,3);
    Z11 = Z1*Z1; Z22 = Z2*Z2;
    NZ = Z11 + Z22;
    Z3 = MAT(3,3)
    Z12 = Z1*Z2; Z23 = Z2*Z3; Z13 = Z1*Z3

    ! if (Z3 < -1D0 + SQRT(NZTinY)) then
    !    PRinT*, 'ERROR in COORDS2ROTMAT: have hit gimbal lock with Z axis pointing downward. Not set up to deal with this yet.'
    !    STOP 1
    ! ENDif

    CAG = COS(AG); SAG = Sin(AG)

    !if (TESTQUAT) PRinT*, 'TESTX1', NZ, NZTinY
    if (NZ < NZTinY) then
       if (ZVEC(3) < 0) then
          PRinT*, 'PROBLEM in COORDS2ROTMAT: some z-vector falls along the &
               & negative z axis. This causes gimbal lock. &
               & If working with a single nucleosome you may be able to &
               & avoid this by rotating the entire structure slightly, or &
               & by using RANDSTART to start the linker beads in different positions.'
          stop 1
       endif
       MAT(1,1) = CAG - (Z12*SAG + Z11*CAG)/2
       !if (testquat) print*, 'testx2:', mat(1,1)
       MAT(2,1) = SAG - (Z22*SAG + Z12*CAG)/2

       MAT(1,2) = -SAG + (Z11*SAG-Z12*CAG)/2
       MAT(2,2) = CAG-(Z22*CAG-Z12*SAG)/2
    else
       MAT(1,1) = (Z22*CAG-Z12*SAG + Z3*(Z12*SAG + Z11*CAG))/NZ
       MAT(2,1) = (-Z12*CAG + Z11*SAG + Z3*(Z22*SAG + Z12*CAG))/NZ

       MAT(1,2) = -(Z22*SAG + Z12*CAG + Z3*(-Z12*CAG + Z11*SAG))/NZ
       !if (testquat) print*, 'testx2:', mat(1,2)
       MAT(2,2) = (Z12*SAG + Z11*CAG + Z3*(Z22*CAG-Z12*SAG))/NZ
    ENDif
    MAT(3,1) = -Z2*SAG-Z1*CAG
    MAT(3,2) = -CAG*Z2 + SAG*Z1

    if (PRESENT(DMATZ)) then
       DMATZ = 0D0

       ! derivative of normalized z axis wrt ZVEC coordinates (transposed)
       DMATZ(:,3,1) = -ZVEC(1)*MAT(:,3)/NZTOT2 + (/1D0/NZTOT,0D0,0D0/)
       DMATZ(:,3,2) = -ZVEC(2)*MAT(:,3)/NZTOT2 + (/0D0,1D0/NZTOT,0D0/)
       DMATZ(:,3,3) = -ZVEC(3)*MAT(:,3)/NZTOT2 + (/0D0,0D0,1D0/NZTOT/)

       if (NZ < NZTinY) then

          DXDZ(1,1) = -Z2*SAG/2-Z1*CAG
          DXDZ(1,2) = -Z1*SAG/2
          DXDZ(2,1) = -Z2*CAG/2
          DXDZ(2,2) = -Z2*SAG-Z1*CAG/2
          DXDZ(:,3) = 0D0

          DYDZ(1,1) = Z1*SAG - Z2*CAG/2
          DYDZ(1,2) = -Z1*CAG/2
          DYDZ(2,1) = Z2*SAG/2
          DYDZ(2,2) = -Z2*CAG + Z1*SAG/2
          DYDZ(:,3) = 0D0
       else
          ! derivative wrt normalized z axis
          DXDZ(1,1) = (-Z2*SAG + Z23*SAG + 2*Z13*CAG-2*Z1*MAT(1,1))/NZ
          DXDZ(1,2) = (2*Z2*CAG-Z1*SAG + Z13*SAG-2*Z2*MAT(1,1))/NZ
          DXDZ(1,3) = (Z12*SAG + Z11*CAG)/NZ

          DXDZ(2,1) = (-Z2*CAG + 2*Z1*SAG + Z23*CAG - 2*MAT(2,1)*Z1)/NZ
          DXDZ(2,2) = (-Z1*CAG + 2*Z23*SAG + Z13*CAG-2*MAT(2,1)*Z2)/NZ
          DXDZ(2,3) = (Z22*SAG + Z12*CAG)/NZ

          DYDZ(1,1) = (-Z2*CAG + Z23*CAG-2*Z13*SAG-2*Z1*MAT(1,2))/NZ
          DYDZ(1,2) = (-2*Z2*SAG-Z1*CAG + Z13*CAG-2*Z2*MAT(1,2))/NZ
          DYDZ(1,3) = (Z12*CAG-Z11*SAG)/NZ

          DYDZ(2,1) = (Z2*SAG + 2*Z1*CAG-Z23*SAG-2*Z1*MAT(2,2))/NZ
          DYDZ(2,2) = (Z1*SAG + 2*Z23*CAG-Z13*SAG-2*Z2*MAT(2,2))/NZ
          DYDZ(2,3) = (Z22*CAG-Z12*SAG)/NZ
       ENDif
       DXDZ(3,:) = (/-CAG,-SAG,0D0/)
       DYDZ(3,:) = (/SAG,-CAG,0D0/)

       ! derivatives of new axes

       CALL DGEMM('N','N',3,3,3,1D0,DXDZ,3,DMATZ(:,3,:),3,0D0,DMATZ(:,1,:),3)
       CALL DGEMM('N','N',3,3,3,1D0,DYDZ,3,DMATZ(:,3,:),3,0D0,DMATZ(:,2,:),3)


       !CALL DGEMV('N',3,3,1D0,DZ,3,DMATZTMP(2,1,:),1,0D0,DMATZ(2,1,:),1)
       !DMATZ(2,1,:) = DXDZ
    ENDif

    if (PRESENT(DMATAG)) then
       DMATAG = 0D0

       if (NZ < NZTinY) then
          DMATAG(1,1) = -SAG-(Z12*CAG-Z11*SAG)/2
          DMATAG(2,1) = CAG-(Z22*CAG-Z12*SAG)/2

          DMATAG(1,2) = -CAG + (Z11*CAG + Z12*SAG)/2
          DMATAG(2,2) = -SAG + (Z22*SAG + Z12*CAG)/2
       else
          DMATAG(1,1) = (-Z22*SAG-Z1*Z2*CAG + Z1*Z23*CAG-Z3*Z11*SAG)/NZ
          DMATAG(2,1) = (Z12*SAG + Z11*CAG + Z3*Z22*CAG-Z1*Z23*SAG)/NZ

          DMATAG(1,2) = (-Z22*CAG + Z12*SAG-Z1*Z23*SAG-Z3*Z11*CAG)/NZ
          DMATAG(2,2) = (Z12*CAG-Z11*SAG-Z22*Z3*SAG-Z1*Z23*CAG)/NZ
       ENDif
       DMATAG(3,1) = -Z2*CAG + Z1*SAG ! dX3/dA
       DMATAG(3,2) = SAG*Z2 + CAG*Z1! dY3/dA
    ENDif
  END subroutine COORDS2ROTMAT

  ! ----------------- STUFF inVOLVinG TREATinG QUATERNIONS AS 4-VECTORS -----
  FUNCTION QUAT2QV(Q)
    ! convert a quaternion object to a 4-vector
    implicit none
    real(dp) :: QUAT2QV(4)
    TYPE(QUATERNION), intent(in) :: Q

    QUAT2QV = (/Q%W,Q%X,Q%Y,Q%Z/)
  END FUNCTION QUAT2QV

  subroutine ROTQV(THETA,AX,QV,DT)
    ! get the quaternion corresponding to rotation around axis AX by angle THETA
    ! as a 4-vector (in QV); optionally, also get the derivative wrt theta
    ! WARNinG: AX assumed to be normalized; does not check for this!

    implicit none
    real(dp), intent(in) :: THETA, AX(3)
    real(dp), intent(out) :: QV(4)
    real(dp), intent(out), OPTIONAL :: DT(4)
    real(dp) :: CT, ST

    CT = COS(THETA/2); ST = Sin(THETA/2);
    QV(1) = CT
    QV(2:4) = ST*AX

    if (PRESENT(DT)) then
       DT(1) = -ST/2
       DT(2:4) = CT/2*AX
    ENDif

  END subroutine ROTQV

  subroutine QVPTMULT(Q,PT,QP,DQ,DPT)
    ! multiply a quaternion by a point in 3-space
    ! returns the result in QP
    ! optionally, returns derivatives wrt the quaternion components in dQ
    ! or wrt point components in DPT
    ! WARNinG: no normalization happens here
    ! WARNinG: this is  a pretty inefficient way to do things; should fix at some point without resorting to rotation matrices
    implicit none
    real(dp), intent(in) :: Q(4),PT(3)
    real(dp), intent(out) :: QP(3)
    real(dp), intent(out), OPTIONAL :: DQ(3,4),DPT(3,3)
    real(dp) :: QN, QinV(4),QTMP(4),QTMP2(4)
    real(dp) :: MAT(3,3),DMAT(3,3,4)
    integer :: I,J
    TYPE(QUATERNION) :: QUAT

    QN = doT_PRODUCT(Q,Q)
    QUAT%W = Q(1); QUAT%X = Q(2); QUAT%Y = Q(3); QUAT%Z = Q(4)
    if (PRESENT(DQ)) then
       CALL QUAT2ROTMAT(QUAT,MAT,dMAT)
    else
       CALL QUAT2ROTMAT(QUAT,MAT)
    ENDif
    do I = 1,3
       QP(I) = doT_PRODUCT(MAT(I,:),PT)
       if (PRESENT(DQ)) then
          do J = 1,4
             DQ(I,J) = doT_PRODUCT(dMAT(I,:,J),PT)
          ENDdo
       ENDif
    ENDdo

    QP = QP/QN
    if (PRESENT(DQ)) then
       do I = 1,3
          do J = 1,4
             DQ(I,J) = (DQ(I,J) - 2*Q(J)*QP(I))/QN
          ENDdo
       ENDdo
    ENDif

    if (PRESENT(DPT)) then
       DPT =MAT/QN
    END if
  END subroutine QVPTMULT

  ! ----------------- STUFF inVOLVinG QUATERNION OBJECTS -------------------

  TYPE(QUATERNION) FUNCTION RHQinTERP(Q1,Q2,F)
    ! interpolate from one quaternion to another always in a right-handed sense
    ! (relative to the Q1 z-axis)
    ! F should be between 0 and 1
    TYPE(QUATERNION) :: Q1, Q2
    real(dp) :: F
    real(dp) :: ANG, AX(3), DIR
    TYPE(QUATERNION) :: QREL

    ! relative quaternion for Q2 relative to Q1
    QREL = inVQUAT(Q1)*Q2

    ! angle of rotation
    ANG = ACOS(QREL%W)*2

    ! axis of rotation
    AX = (/QREL%X,QREL%Y,QREL%Z/); AX = AX/SQRT(doT_PRODUCT(AX,AX))

    DIR = doT_PRODUCT(AX,QUAT2PT(Q1*PTQUAT((/0D0,0D0,1D0/))/Q1))

    print*, 'testx2:', dir, ang, ax

    if (DIR < 0) then
       ! rotation axis points away from quaternion axis, so flip it
       ANG = -ANG
       AX = -AX
    ENDif

    QREL = ROTQUAT(ANG*F,AX)

    RHQinTERP = Q1*QREL
  END FUNCTION RHQinTERP

  TYPE(QUATERNION) FUNCTION QSLERP(Q1,Q2,F)
    ! Sphreical linear interpolation between two quaternions
    ! F should be between 0 and 1
    TYPE(QUATERNION), intent(in) :: Q1, Q2
    real(dp), intent(in) :: F
    real(dp) :: QV1(4), QV2(4), ANG, QVANS(4)

    QV1 = (/Q1%W,Q1%X,Q1%Y,Q1%Z/)
    QV2 = (/Q2%W, Q2%X, Q2%Y, Q2%Z/)

    ! angle between them
    ANG = ACOS(doT_PRODUCT(QV1,QV2))

    QVANS = Sin((1-F)*ANG)/Sin(ANG)*QV1 + Sin(F*ANG)/Sin(ANG)*QV2

    QSLERP%W = QVANS(1); QSLERP%X = QVANS(2); QSLERP%Y = QVANS(3); QSLERP%Z = QVANS(4);
  END FUNCTION QSLERP

  subroutine QUAT2ROTMAT(Q,MAT,DMAT)
    ! convert a quaternion object to a rotation matrix and optionally return derivatives
    ! NOTE: no normalization
    implicit none
    TYPE(QUATERNION) :: Q
    real(dp), intent(out) :: MAT(3,3)
    real(dp), intent(out), OPTIONAL :: DMAT(3,3,4)
    real(dp) :: A,B,C,D,AA,BB,CC,DD,AB,AC,AD,BC,BD,CD

    A = Q%W; B = Q%X; C = Q%Y; D = Q%Z
    AA = A*A; BB = B*B; CC = C*C; DD = D*D
    AB = 2*A*B; AC = 2*A*C; AD = 2*A*D
    BC = 2*B*C; BD = 2*B*D
    CD = 2*C*D

    MAT(1,:) = (/AA + BB-CC-DD,BC-AD,AC + BD/)
    MAT(2,:) = (/AD + BC,AA-BB + CC-DD,CD-AB/)
    MAT(3,:) = (/BD-AC,AB + CD,AA-BB-CC + DD/)

    if (PRESENT(DMAT)) then
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
    ENDif
  END subroutine QUAT2ROTMAT

  TYPE(QUATERNION) FUNCTION ROTMAT2QUAT(R)
    ! convert from a rotation matrix to a quaternion object
    ! following Deibel, 2006 (but with the rotation matrix transposed)
    ! R(3,:) has the z axis after rotation is applied, etc.
    ! assumes R is orthonormal
    implicit none
    real(dp) :: R(3,3)
    real(dp) :: R11, R22,R33, TMP

    R11 = R(1,1); R22 = R(2,2); R33 = R(3,3)

    if (R22 >= -R33.AND.R11 >= -R22.AND.R11 >= -R33) then
       TMP = SQRT(1 + R11 + R22 + R33)
       ROTMAT2QUAT%W = TMP/2
       ROTMAT2QUAT%X = (R(2,3)-R(3,2))/(TMP*2)
       ROTMAT2QUAT%Y = (R(3,1)-R(1,3))/(TMP*2)
       ROTMAT2QUAT%Z = (R(1,2)-R(2,1))/(TMP*2)
       !PRinT*, 'TESTX1Q'
    elseif (R22 <= -R33.AND.R11 > R22.AND.R11 > R33) then
       TMP = SQRT(1 + R11-R22-R33)
       ROTMAT2QUAT%W = (R(2,3)-R(3,2))/(TMP*2)
       ROTMAT2QUAT%X = TMP/2
       ROTMAT2QUAT%Y = (R(1,2) + R(2,1))/(TMP*2)
       ROTMAT2QUAT%Z = (R(3,1) + R(1,3))/(TMP*2)
       !PRinT*, 'TESTX2Q'
    elseif (R22 > R33.AND.R11 < R22.AND.R11 <= -R33) then
       TMP = SQRT(1-R11 + R22-R33)
       ROTMAT2QUAT%W = (R(3,1)-R(1,3))/(TMP*2)
       ROTMAT2QUAT%X = (R(1,2) + R(2,1))/(TMP*2)
       ROTMAT2QUAT%Y = TMP/2
       ROTMAT2QUAT%Z = (R(2,3) + R(3,2))/(TMP*2)
       !PRinT*, 'TESTX3Q'
    elseif (R22 < R33.AND.R11 <= -R22.AND.R11 < R33) then
       TMP = SQRT(1D0-R11-R22 + R33)
       ROTMAT2QUAT%W = (R(1,2)-R(2,1))/(TMP*2)
       ROTMAT2QUAT%X = (R(3,1) + R(1,3))/(TMP*2)
       ROTMAT2QUAT%Y = (R(2,3) + R(3,2))/(TMP*2)
       ROTMAT2QUAT%Z = TMP/2
       !PRinT*, 'TESTX4Q'
    else
       PRinT*, 'ERROR in ROTMAT2QUAT: bad rotation matrix'
       PRinT*, R(:,1)
       PRinT*, R(:,2)
       PRinT*, R(:,3)
       STOP 1
       ROTMAT2QUAT%W = 0; ROTMAT2QUAT%X = 0; ROTMAT2QUAT%Y = 0; ROTMAT2QUAT%Z = 0
    ENDif
  END FUNCTION ROTMAT2QUAT

  subroutine EULER2QUAT(EUL,Q,DERV,GETDERV)
    ! convert from Euler angles (z-x-z convention) to a quaternion
    ! if GETDERV is true, also get the derivatives of the quaternion components
    ! with respect to the euler angles (4 rows by 3 columns)
    ! WARNinG: since the quaternion representation has more parameters, switching from quaternion
    ! to euler and back again will not always give the exact same quaternion, though it
    ! will give an equivalent one (eg: may invert all components)

    real(dp), intent(in) :: EUL(3)
    TYPE(QUATERNION), intent(out) :: Q
    real(dp), intent(out),OPTIONAL :: DERV(4,3)
    LOGICAL, intent(in),OPTIONAL :: GETDERV
    real(dp) :: CA,SA,CB,SB,CG,SG,ALPHA,GAMMA,BETA
    LOGICAL :: FLIPBETA

    BETA = ANGLE2PI(EUL(2))
    FLIPBETA = BETA > PI
    ALPHA = EUL(1); GAMMA = EUL(3)
    if (BETA > PI) then
       BETA = 2*PI-BETA
       GAMMA = GAMMA + PI
       ALPHA = ALPHA + PI
    ENDif

    ALPHA = ANGLE2PI(ALPHA)
    GAMMA = ANGLE2PI(GAMMA)
    !if (QUATTEST) then
    !   PRinT*, 'TESTXE:', EUL
    !   PRinT*, 'A,B,G:',ALPHA,BETA,GAMMA
    !ENDif

    CA = COS(ALPHA/2); SA = Sin(ALPHA/2)
    CB = COS(BETA/2); SB = Sin(BETA/2)
    CG = COS(GAMMA/2); SG = Sin(GAMMA/2)

    Q%W = CG*CB*CA - SG*CB*SA
    Q%X = SG*SB*SA + CG*SB*CA
    Q%Y = CG*SB*SA - SG*SB*CA
    Q%Z = CG*CB*SA + SG*CB*CA

    if (PRESENT(DERV)) then
       if (GETDERV) then
          DERV(1,:) = 0.5D0*(/-CG*CB*SA - SG*CB*CA, -CG*SB*CA + SG*SB*SA, -SG*CB*CA - CG*CB*SA/)
          DERV(2,:) = 0.5D0*(/SG*SB*CA - CG*SB*SA, SG*CB*SA + CG*CB*CA, CG*SB*SA - SG*SB*CA /)
          DERV(3,:) = 0.5D0*(/CG*SB*CA + SG*SB*SA, CG*CB*SA - SG*CB*CA, -SG*SB*SA - CG*SB*CA/)
          DERV(4,:) = 0.5D0*(/CG*CB*CA - SG*CB*SA, -CG*SB*SA - SG*SB*CA, -SG*CB*SA + CG*CB*CA/)
       ENDif
       if (FLIPBETA) then
          DERV(:,2) = -DERV(:,2)
       ENDif
    ENDif

  END subroutine EULER2QUAT

  subroutine QUAT2EULER(Q,EUL)
    ! get the Euler angles (z-x-z convention) corresponding to a unit quaternion
    ! NOTE: the quaternion must already be normalized
    ! currently can't handle gimbal lock
    TYPE(QUATERNION), intent(in) :: Q
    real(dp), intent(out) :: EUL(3)
    real(dp) :: DUMMY

    EUL(2) = 1 - 2*(Q%X**2 + Q%Y**2)
    if (EUL(2) > 1) then
       EUL(2) = 0D0
    elseif (EUL(2) < -1) then
       EUL(2) = PI
    else
       EUL(2) = ACOS(EUL(2))
    ENDif

    ! deal with the gimbal lock issues
    if (EUL(2) <= EPSILON(0D0)) then
       EUL(3) = ATAN2(2*(Q%X*Q%Y + Q%W*Q%Z), 1 - 2*(Q%Y**2 + Q%Z**2))
       if (EUL(3) < 0) EUL(3) = 2*PI + EUL(3)
       EUL(1) = 0D0
       RETURN
    elseif (EUL(2) >= PI-EPSILON(0D0)) then
       EUL(3) = ATAN2(-2*(Q%X*Q%Y + Q%W*Q%Z), 1 - 2*(Q%Y**2 + Q%Z**2))
       if (EUL(3) < 0) EUL(3) = 2*PI + EUL(3)
       EUL(1) = 0D0
       RETURN
    ENDif

    DUMMY = Q%W*Q%X-Q%Y*Q%Z

    EUL(1) = ATAN2(Q%W*Q%Y + Q%X*Q%Z, DUMMY)
    if (EUL(1) < 0) EUL(1) = 2*PI + EUL(1)

    DUMMY = Q%W*Q%X + Q%Y*Q%Z

    EUL(3) = ATAN2(Q%X*Q%Z-Q%W*Q%Y,DUMMY)
    if (EUL(3) < 0) EUL(3) = 2*PI + EUL(3)

  END subroutine QUAT2EULER

  FUNCTION QUAT2PT(Q)
    ! get the point corresponding to a quaternion
    TYPE(QUATERNION) :: Q
    real(dp) :: QUAT2PT(3)

    QUAT2PT = (/Q%X,Q%Y,Q%Z/)
  END FUNCTION QUAT2PT

  TYPE(QUATERNION) FUNCTION ROTQUAT(THETA,AX)
    ! get the quaternion corresponding to rotation around unit axis AX
    ! by an angle theta
    real(dp) :: THETA, AX(3)
    real(dp) :: T2, ST2

    if (ABS(SUM(AX**2)-1D0) > EPSILON(1d0)*10) then
       print*, 'ERROR in ROTQUAT: Axis does not have unit norm'
       STOP 1
    ENDif

    T2 = THETA/2; ST2 = Sin(THETA/2)
    ROTQUAT%W = COS(T2)
    ROTQUAT%X = ST2*AX(1); ROTQUAT%Y = ST2*AX(2); ROTQUAT%Z = ST2*AX(3)
  END FUNCTION ROTQUAT

  TYPE(QUATERNION) FUNCTION PTQUAT(P)
    ! turn a 3d point into a quaternion
    real(dp) :: P(3)

    PTQUAT%W = 0D0; PTQUAT%X = P(1); PTQUAT%Y = P(2); PTQUAT%Z = P(3)
  END FUNCTION PTQUAT

  TYPE(QUATERNION) FUNCTION inVQUAT(Q)
    ! get the inverse of a quaternion
    implicit none
    TYPE(QUATERNION), intent(in) :: Q
    real(dp) :: QN

    QN = Q%W**2 + Q%X**2 + Q%Y**2 + Q%Z**2
    inVQUAT%W = Q%W/QN; inVQUAT%X = -Q%X/QN; inVQUAT%Y = -Q%Y/QN; inVQUAT%Z = -Q%Z/QN

  END FUNCTION inVQUAT

  TYPE(QUATERNION) FUNCTION QDIVIDE(P,Q)
    ! multiply P by inverse of Q(in that order)
    implicit none
    TYPE(QUATERNION), intent(in) :: P,Q
    TYPE(QUATERNION) :: QinV
    real(dp) :: QN

    ! inverse of the 2nd quaternion
    QN = Q%W**2 + Q%X**2 + Q%Y**2 + Q%Z**2
    QinV%W = Q%W/QN; QinV%X = -Q%X/QN; QinV%Y = -Q%Y/QN; QinV%Z = -Q%Z/QN

    QDIVIDE = P*QinV
  END FUNCTION QDIVIDE

  TYPE(QUATERNION) FUNCTION QPRODUCT(P,Q)
    ! quaternion multiplication
    implicit none
    TYPE(QUATERNION), intent(in) :: P, Q

    QPRODUCT%W = P%W*Q%W - P%X*Q%X - P%Y*Q%Y - P%Z*Q%Z
    QPRODUCT%X = P%W*Q%X + P%X*Q%W + P%Y*Q%Z - P%Z*Q%Y
    QPRODUCT%Y = P%W*Q%Y - P%X*Q%Z + P%Y*Q%W + P%Z*Q%X
    QPRODUCT%Z = P%W*Q%Z + P%X*Q%Y - P%Y*Q%X + P%Z*Q%W

  END FUNCTION QPRODUCT

  ! --------- general angle and euler angle stuff -------------
  subroutine ROTANGAX(ANG,AX,inVEC,outVEC,CALCROTMAT,ROTMAT)
    ! rotate a 3D vector by angle ANG around axis AX
    ! if CALCROTMAT is true, recalculate the rotation matrix
    ! otherwise use the provided one
    implicit none
    real(dp), intent(in) :: ANG, AX(3), inVEC(3)
    real(dp), intent(out) :: outVEC(3)
    real(dp), intent(inout) :: ROTMAT(3,3)
    LOGICAL, intent(in)  :: CALCROTMAT
    real(dp) :: CT,ST,CT1
    integer :: I

    if (CALCROTMAT) then
       CT = COS(ANG); ST = Sin(ANG)
       CT1 = 1-CT
       ROTMAT(1,:) = (/CT + AX(1)**2*CT1, AX(1)*AX(2)*CT1-AX(3)*ST, AX(1)*AX(3)*CT1 + AX(2)*ST/)
       ROTMAT(2,:) = (/AX(2)*AX(1)*CT1 + AX(3)*ST,CT + AX(2)**2*CT1,AX(2)*AX(3)*CT1-AX(1)*ST/)
       ROTMAT(3,:) = (/AX(3)*AX(1)*CT1-AX(2)*ST,AX(3)*AX(2)*CT1 + AX(1)*ST, CT + AX(3)**2*CT1/)
    ENDif

    do I = 1,3
       outVEC(I) = doT_PRODUCT(ROTMAT(I,:),inVEC)
    ENDdo
  END subroutine ROTANGAX

  subroutine EUL2ROTMAT(EUL,ROTMAT,DMAT)
    ! get the rotation matrix corresponding to various euler angles
    ! and the appropriate derivatives if DMAT is present
    implicit none
    real(dp), intent(in) :: EUL(3)
    real(dp), intent(out) :: ROTMAT(3,3)
    real(dp), intent(out), OPTIONAL :: DMAT(3,3,3)
    real(dp) :: CA,SA,CB,SB,CG,SG

    CA = COS(EUL(1)); SA = Sin(EUL(1))
    CB = COS(EUL(2)); SB = Sin(EUL(2))
    CG = COS(EUL(3)); SG = Sin(EUL(3))

    ROTMAT(1,:) = (/CA*CG-SA*CB*SG, -CA*SG-SA*CB*CG, SB*SA/)
    ROTMAT(2,:) = (/SA*CG + CA*CB*SG,-SA*SG + CA*CB*CG,-SB*CA/)
    ROTMAT(3,:) = (/SB*SG,SB*CG,CB/)

    if (PRESENT(DMAT)) then
       dMAT(1,:,1) = (/-SA*CG-CA*CB*SG,SA*SG-CA*CB*CG,SB*CA/)
       dMAT(2,:,1) = (/CA*CG-SA*CB*SG,-CA*SG-SA*CB*CG,SB*SA/)
       dMAT(3,:,1) = 0D0

       dMAT(1,:,2) = (/SA*SB*SG,SA*SB*CG,CB*SA/)
       dMAT(2,:,2) = (/-CA*SB*SG,-CA*SB*CG,-CB*CA/)
       dMAT(3,:,2) = (/CB*SG,CB*CG,-SB/)

       dMAT(1,:,3) = (/-CA*SG-SA*CB*CG,-CA*CG + SA*CB*SG,0D0/)
       dMAT(2,:,3) = (/-SA*SG + CA*CB*CG,-SA*CG-CA*CB*SG,0D0/)
       dMAT(3,:,3) = (/SB*CG,-SB*SG,0D0/)
    ENDif
  END subroutine EUL2ROTMAT

  subroutine GETANGLE(IJ,KJ,CST,dCTdIJ,dCTdKJ)
    ! get the angle between three points (I-J-K)
    ! and, optionally, the derivative of that angle
    ! actually this returns the COSinE of the angle and its derivative
    ! IJ = I-J; KJ = K-J
    implicit none
    real(dp), intent(in) :: IJ(3),KJ(3)
    real(dp), intent(out) :: CST
    real(dp), intent(out), OPTIONAL :: dCTdIJ(3), dCTdKJ(3)

    real(dp) :: DXI, DYI, DZI, DXJ, DYJ, DZJ
    real(dp) :: RI2, RJ2, RI, RJ, RIR, RJR
    real(dp) :: DXIR,DYIR, DZIR,DXJR,DYJR,DZJR

    DXI = IJ(1); DYI = IJ(2); DZI = IJ(3)
    DXJ = KJ(1); DYJ = KJ(2); DZJ = KJ(3)

    RI2 = DXI*DXI + DYI*DYI + DZI*DZI
    RJ2 = DXJ*DXJ + DYJ*DYJ + DZJ*DZJ
    RI = SQRT(RI2)
    RJ = SQRT(RJ2)
    RIR = 1/RI
    RJR = 1/RJ
    DXIR = DXI*RIR
    DYIR = DYI*RIR
    DZIR = DZI*RIR
    DXJR = DXJ*RJR
    DYJR = DYJ*RJR
    DZJR = DZJ*RJR
    CST = DXIR*DXJR + DYIR*DYJR + DZIR*DZJR
    if (PRESENT(DCTDIJ)) then
       dCTdIJ(1) = -(DXIR*CST-DXJR)*RIR
       dCTdIJ(2) = -(DYIR*CST-DYJR)*RIR
       dCTdIJ(3) = -(DZIR*CST-DZJR)*RIR
    ENDif
    if (PRESENT(DCTDKJ)) then
       dCTdKJ(1) = -(DXJR*CST-DXIR)*RJR
       dCTdKJ(2) = -(DYJR*CST-DYIR)*RJR
       dCTdKJ(3) = -(DZJR*CST-DZIR)*RJR
    ENDif

  END subroutine GETANGLE

  subroutine GETDIHEDRAL(IJ,JK,LK, PHI, dPdIJ, dPdJK, dPdLK)
    ! dihedral angle for 4 atoms I, J, K, L bound in order
    ! IJ = I-J; JK = J-K; LK = L-K
    ! find the dihedral torsion angle; return it in PHI
    ! Also return all derivatives: dP/dIJx, dP/dIJy, dP/dIJz in triplet dPdIJ
    ! same with dPdJK, dPdLK
    implicit none
    real(dp), intent(in) ::  IJ(3), JK(3), LK(3)
    real(dp), intent(out) :: PHI
    real(dp), intent(out), OPTIONAL :: dPdIJ(3), dPdJK(3), dPdLK(3)
    real(dp) :: DPDJ(3),DPDI(3),DPDK(3),DPDL(3)
    real(dp) :: FX, FY, FZ, GX, GY, GZ, HX, HY, HZ
    real(dp) :: AX, AY, AZ, BX, BY, Bz
    real(dp) :: RF, RG, RH, RF2, RG2, RH2, RFR, RGR, RHR
    real(dp) :: CSTTWO, SNTTWO2, CSTTHREE, SNTTHREE2, SNTTWO2R, SNTTHREE2R
    real(dp) :: RA2, RB2, RA2R, RB2R, RABR, CP
    real(dp) :: MYTX, MYTY, MYTZ, MYSCALAR
    real(dp) :: DUMMY, DUMMY2
    real(dp) :: B1(3), B2(3), B3(3), B12(3), B23(3)
    LOGICAL :: NOCOOR = .FALSE.


    FX = IJ(1)
    FY = IJ(2)
    FZ = IJ(3)
    GX = JK(1)
    GY = JK(2)
    GZ = JK(3)
    HX = LK(1)
    HY = LK(2)
    HZ = LK(3)
    ! A = F x G, B = H x G
    AX = FY*GZ-FZ*GY
    AY = FZ*GX-FX*GZ
    AZ = FX*GY-FY*GX
    BX = HY*GZ-HZ*GY
    BY = HZ*GX-HX*GZ
    BZ = HX*GY-HY*GX
    ! RG = |G|, RGR = 1/|G|
    RG2 = GX*GX + GY*GY + GZ*GZ
    RG = SQRT(RG2)
    RGR = 1/RG
    ! dae for use in evaluating B-matrix
    RF2 = FX*FX + FY*FY + FZ*FZ
    RF = SQRT(RF2)
    RFR = 1/RF
    RH2 = HX*HX + HY*HY + HZ*HZ
    RH = SQRT(RH2)
    RHR = 1/RH


    CSTTWO = -(FX*GX + FY*GY + FZ*GZ)*RFR*RGR
    SNTTWO2 = 1-CSTTWO*CSTTWO
    SNTTWO2R = 1/SNTTWO2
    CSTTHREE = (HX*GX + HY*GY + HZ*GZ)*RHR*RGR
    SNTTHREE2 = 1-CSTTHREE*CSTTHREE
    SNTTHREE2R = 1/SNTTHREE2

    RA2 = AX*AX + AY*AY + AZ*AZ
    RB2 = BX*BX + BY*BY + BZ*BZ
    RA2R = 1/RA2
    RB2R = 1/RB2
    RABR = SQRT(RA2R*RB2R)

    PHI =  ATAN2(-RG*(FX*BX + FY*BY + FZ*BZ),AX*BX + AY*BY + AZ*BZ)


    if (PRESENT(DPDIJ).OR.PRESENT(DPDJK).OR.PRESENT(DPDLK)) then
       DUMMY = RFR*RFR*RGR*SNTTWO2R
       dPdI = (/-AX*DUMMY, -AY*DUMMY, -AZ*DUMMY/)
       DUMMY = RFR*RFR*RGR*RGR*SNTTWO2R*(RG-RF*CSTTWO)
       DUMMY2 = RHR*RGR*RGR*SNTTHREE2R*CSTTHREE
       dPdJ = (/AX*DUMMY-BX*DUMMY2, AY*DUMMY-BY*DUMMY2, AZ*DUMMY-BZ*DUMMY2/)
       DUMMY = RHR*RHR*RGR*SNTTHREE2R
       dPdL = (/BX*DUMMY,BY*DUMMY,BZ*DUMMY/)
    ENDif
    if (PRESENT(DPDIJ)) DPDIJ = DPDI
    if (PRESENT(DPDLK)) DPDLK = DPDL
    if (PRESENT(DPDJK)) DPDJK = DPDI + DPDJ


  END subroutine GETDIHEDRAL

  real(dp) FUNCTION ANGLE0(ANGLE)
    ! convert an angle to one between +/- pi
    ! (so keeping it as close as possible to zero)
    implicit none
    real(dp) :: ANGLE

    ANGLE0 = ANGLE2PI(ANGLE)
    if (ANGLE0 > PI) then
       ANGLE0 = ANGLE0 - 2*PI
    ENDif
  END FUNCTION ANGLE0

  real(dp) FUNCTION ANGLE2PI(ANGLE)
    ! convert an angle to one between 0 and 2pi by adding or subtracting multiples of 2pi
    implicit none
    real(dp) :: ANGLE
    integer :: N2PI

    N2PI = inT(ANGLE/(2*PI))
    if (ANGLE < 0) then
       ANGLE2PI = ANGLE + (-N2PI + 1)*2*PI
    else
       ANGLE2PI = ANGLE - N2PI*2*PI
    ENDif
  END FUNCTION ANGLE2PI

END module QUATUTIL
