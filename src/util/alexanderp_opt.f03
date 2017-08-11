
! This subroutine calculates the Alexander polynomial A(t) of a closed
! polygonal space curve evaluated at t = -1

subroutine ALEXANDERP(R,N,DELTA,Cross,CrossSize,NCross)
  use params, only : dp
  implicit none

  !inPUT VARIABLES
  integer N                     ! Number of points in space curve
  real(dp) R(3,N)       !Space curve
  integer CrossSize            !Size of cross matrix. Overallocate to avoid resizing
  !outPUT VARIABLES
  integer  DELTA
  !inTERMEDIATE VARIABLES
  real(dp), ALLOCATABLE ::  A(:,:) ! Alexander polynomial matrix evaluated at t = -1
  real(dp) NV(3)        ! normal vector for projection
  real(dp) RP(3,N)      ! projection of R onto plane defined by NV
  real(dp) RdoTN(N)        ! Dot product of R and NV
  integer Ncross                ! number of crossings in projection
  integer Ndegen                ! number of crossings for a given segment
  real(dp) Cross(CrossSize,6)        ! matrix of cross indices and coordinates
  integer I,J,K,IP1,JP1,KP1           ! iteration indices
  real(dp) smax,tmax    ! length of segments in the projection
  real(dp) ui(3),uj(3)  ! tangent vectors of segments in the projection
  real(dp) udot         ! dot product of tangents in the projection
  real(dp) t
  real(dp) sint,tint    ! intersection coordinates for projected segments
  real(dp) srmax,trmax  ! maximum true length of segment (non-projected)
  real(dp) srint,trint  ! intersection coordinates in unprojected coordinates
  real(dp) uri(3),urj(3) ! tangent vectors of unprojected segments
  real(dp) DRI(3),DRJ(3) !Displacement vectors of unprojected segments
  real(dp) thetai,thetaj ! angle between real segment and projected segment
  integer, ALLOCATABLE :: over_ind(:) !Vector of over_pass indices (index of over-pass of Nth crossing)
  real(dp) delta_double ! real(dp) form of delta. To be converted to integer
  integer index,inD

  !Set the normal vector for the plane of projection. Currently set to parallel to z axis
  NV = 0.
  NV(3) = 1.

  !Calculate the projection of R onto the projection plane
  do I = 1,N
     RdoTN(I) = R(1,I)*NV(1) + R(2,I)*NV(2) + R(3,I)*NV(3)
  ENDdo


  !Set-up matrix containing information on crossings.
  Cross = 0.
  Cross(:,5) = 1
  !Calculate the projection of the curve into the plane with normal NV

  do I = 1,N
     RP(:,I) = R(:,I)-RdoTN(I)*NV
  ENDdo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Determine all crossings in the projection of the curve
  !Current algorithm uses a brute-force search for all intersections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Ncross = 0

  do I = 1,N

     do J = 1,I-2
        !Calculate succeeding segments IP1 and JP1
        if (I == N) then
           IP1 = 1
        else
           IP1 = I + 1
        ENDif

        if (J == N) then
           JP1 = 1
        else
           JP1 = J + 1
        ENDif

        !Skip calculation for adjacent segments
        if (I == J.OR.J == IP1.OR.I == JP1) then
           GOTO 10
        ENDif

        !Calculate lengths of segments in the projection and the tangents
        smax = SQRT(SUM((RP(:,IP1)-RP(:,I))**2))
        tmax = SQRT(SUM((RP(:,JP1)-RP(:,J))**2))
        ui = (RP(:,IP1)-RP(:,I))/smax
        uj = ((RP(:,JP1)-RP(:,J)))/tmax
        udot = ui(1)*uj(1) + ui(2)*uj(2) + ui(3)*uj(3)

        !If segments are parallel, continue to next segment
        if (udot == 1..OR.udot == -1.) then
           GOTO 10
        ENDif

        !Compute the point of intersection
        tint = (RP(2,j)-RP(2,i)-(ui(2)/ui(1))*(RP(1,j)-RP(1,i)))/((ui(2)*uj(1)/ui(1))-uj(2))
        sint = (RP(1,j)-RP(1,i) + uj(1)*tint)/ui(1)

        !If the intersection point is within length of segments, count as an intersection
        if (sint >= 0.AND.sint < smax.AND.tint >= 0.AND.tint < tmax) then
           !Determine if this is an undercrossing (RI under RJ) or overcrossing

           !Compute lengths and tangents  of true segments (not projected)
           srmax = SQRT(SUM((R(:,IP1)-R(:,I))**2))
           trmax = SQRT(SUM((R(:,JP1)-R(:,J))**2))
           DRI = R(:,IP1)-R(:,I)
           DRJ = R(:,JP1)-R(:,J)
           uri = DRI/srmax
           urj = DRJ/trmax
           !Calculate the angle between the real segment and the projection
           thetai = ATAN(DRI(3)/smax)
           thetaj = ATAN(DRJ(3)/tmax)
           !Calculate point of intersection in the projection in terms of true length
           srint = sint/cos(thetai)
           trint = tint/cos(thetaj)

           !Determine whether this is an undercrossing or an overcrossing.
           !Save the indices appropriately (the index of the undercrossing segment
           !must come first

           if (R(3,i) + uri(3)*srint<R(3,j) + urj(3)*trint) then
              Ncross = Ncross + 1
              Cross(Ncross,1) = i;
              Cross(Ncross,2) = j;
              Cross(Ncross,3) = sint;
              Cross(Ncross,4) = tint;
              Cross(Ncross,6) = -uj(2)*ui(1) + ui(2)*uj(1); !cross(ui,uj), + for RH, - for LH
           else
              Ncross = Ncross + 1
              Cross(Ncross,1) = j
              Cross(Ncross,2) = i
              Cross(Ncross,3) = tint
              Cross(Ncross,4) = sint
              Cross(Ncross,6) = -(-uj(2)*ui(1) + ui(2)*uj(1)); !cross(ui,uj), + for RH, - for LH

           ENDif
        ENDif

10      continue

     ENDdo

  ENDdo

  !Sort the undercrossings according to order of occurrence. The undercrossings that come first
  !as you traverse the chain along its contour length come first. Hence, sort based on the first column of Cross
  CALL bubble_sort(Cross(1:Ncross,:),Ncross,6,1)

  !Sort under-crossings of same segment with respect to order of occurrence (w.r.t sint)

  index = 1

  do while (index < Ncross)
     !Determine the number of degenerate crossings for the current segment index
     Ndegen = 1
     do while (nint(Cross(index + Ndegen,1)) == nint(Cross(index,1)))
        Ndegen = Ndegen + 1
     ENDdo
     if (Ndegen > 1) then
        CALL bubble_sort(Cross(index:index + Ndegen-1,:),Ndegen,6,3)
     ENDif
     Cross(index:index + Ndegen-1,5) = Ndegen
     index = index + Ndegen

  ENDdo

  !Construct vector of over-pass indices according to indexing described by Vologodskii
  !The element in the Nth row is the index of the segment that overpasses the Nth crossing

  allocate(over_ind(Ncross))

  do I = 1,Ncross
     !get original polymer index of overpass
     J = Cross(I,2)
     !Special cases: j compes before first crossing or
     !after the last crossing
     if (J < nint(Cross(1,1))) then
        over_ind(I) = 1
     elseif (J > nint(Cross(Ncross,1))) then
        over_ind(I) = Ncross + 1
     ENDif

     !Sum over all crossings to determine which undercrossing this lies between
     do K = 1,Ncross
        !If J lies between cross K and cross K + 1, then it is segment K + 1
        if (J > nint(Cross(K,1)).AND.J < nint(Cross(K + 1,1))) then
           over_ind(I) = K + 1
           GOTO 20
           !If J = K, then segment j contains undercrossings
           ! then need to determine where overpass lies relative to undercrossing
        elseif (J == nint(CROSS(K,1))) then
           t = Cross(I,4)
           Ndegen = Cross(K,5)
           !Special case: t is before the first under-pass or after the last under-pass
           !of segment j
           if (t <= Cross(K,3)) then
              over_ind(I) = K
              GOTO 20
           elseif (t >= Cross(K + Ndegen-1,3)) then
              OVER_inD(I) = K + Ndegen
              GOTO 20
              !Otherwise, determine which under-crossings t lies between
           else
              inD = 1
              !loop over all degenerate undercrossings
              do while (inD < Ndegen)
                 !if t lies between the s of undercrossing k + ind-1 and the next,
                 !then this over_pass has a new index of k + ind in the re-indexing
                 !scheme
                 if (t > Cross(K + inD-1,3).AND.t < CROSS(k + inD,3)) then
                    over_ind(I) = K + inD
                    GOTO 20
                 ENDif
                 inD = inD + 1
              ENDdo
           ENDif
        ENDif

     ENDdo
20   continue

  ENDdo


  !Calculate the Alexander matrix evaluated at t = -1
  ! Note that the Alexander matrix is correct only up to
  ! a factor of +-1. Since the alexander polynomial evaluated
  ! at t = -1 is always positive, take the absolute value of the
  ! determinant


  allocate(A(Ncross,Ncross))
  A = 0.

  do K = 1,Ncross
     KP1 = K + 1
     I = over_ind(K)
     !Periodic index conditions
     if (K == Ncross) then
        KP1 = 1
     ENDif

     if (I >= Ncross + 1) then
        I = 1
     ENDif

     !calculate values of the matrix
     if (I == K.OR.I == KP1) then
        A(K,K) = -1.
        A(K,KP1) = 1.
     else
        A(K,K) = 1.
        A(K,KP1) = 1.
        A(K,I) = -2.
     ENDif

  ENDdo

  !Calculate the determinant of the matrix with one row and one column removed


  !If A has one crossing or less, it is the trivial knot

  if (Ncross <= 1) then
     delta_double = 1.
  else
     CALL abs_determinant(A(1:Ncross-1,1:Ncross-1),NCross-1,delta_double)
  ENDif

  delta = nint(delta_double)

  !De-allocate arrays

  DEallocate(A)
  DEallocate(over_ind)

  RETURN

END subroutine ALEXANDERP
