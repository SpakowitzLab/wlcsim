!This subroutine calculates the alexander polynomial of a chain represented by a set of
!discrete beads after a crankshaft move has been performed. The routine updates the values
!in the Cross matrix that change during the move. This subroutine currently only works if there
!is one polymer

subroutine alexanderp_crank(wlc_p,R,Delta,Cross,CrossSize,NCross,IT1,IT2,DIB)
  use params, only : dp,wlcsim_params
  implicit none
  type(wlcsim_params), intent(in) :: wlc_p
  !implicit none
  !inPUT VARIABLES
 ! integer N                     ! Number of points in space curve
  integer CrossSize             !Size of the cross matrix (larger than the total number of crossings to prevent reallocation)
  real(dp) R(3,wlc_p%NB)       !Space curve
  real(dp) Cross(CrossSize,6)        !Matrix of cross indices and coordinates
  real(dp) CrossNew(CrossSize,6)
  integer Ncross                !Total number of crossings
  integer NCrossNew
  integer  IT1,IT2               !Indices of beads that bound the segment rotated during crankshaft move
  integer DIB                   !Number of segments in the segment rotated
  integer DIO                   !Number of segments outside the segment rotated
  !outPUT VARIABLES
  integer  DELTA

  !inTERMEDIATE VARIABLES
  integer IT1N,IT2N             !New indicies for segments that bound the crankshaft move. Necessary to capture part that changes when IT1 = IT2
  real(dp), ALLOCATABLE ::  A(:,:) ! Alexander polynomial matrix evaluated at t = -1
  real(dp) NV(3)        ! normal vector for projection
  real(dp) RP(3,wlc_p%NB)      ! projection of R onto plane defined by NV
  real(dp) RdoTN(wlc_p%NB)        ! Dot product of R and NV
  integer Ndegen                ! number of crossings for a given segment
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
  integer II,IO,IIP1,IOP1
  LOGICAL Copy

  !Performance testing for development
  real(dp) TIME1
  real(dp) TIME2
  real(dp) DT_PRUNE
  real(dp) DT_inTERSECT


  !Set the normal vector for the plane of projection. Currently set to parallel to z axis
  NV = 0.
  NV(3) = 1.

  !Calculate the projection of R onto the projection plane
  do I = 1,wlc_p%NB
     RdoTN(I) = R(1,I)*NV(1) + R(2,I)*NV(2) + R(3,I)*NV(3)
  ENDdo

  !Calculate the projection of the curve into the plane with normal NV

  do I = 1,wlc_p%NB
     RP(:,I) = R(:,I)-RdoTN(I)*NV
  ENDdo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Update Cross matrix by removing all instances in which segments that are in the rotated
  !segment appear. This is done by copying to a new matrix only the crossings that don't
  !involve the portion of the chain that was moved
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !If IT1 = IT2, then the two segments that sandwich IT1 are changed, so the indices need to be updated
  !to calculate the part of the Cross matrix that changes
  IT1N = IT1
  IT2N = IT2

  if (IT2N == IT1N) then
     if (DIB == 0) then
        if (IT1N /= 1.AND.IT1N /= wlc_p%Nb) then
           IT2N = IT1N + 1
           IT1N = IT1N-1
        elseif (IT1N == 1) then
           IT2N = IT1N + 1
           IT1N = wlc_p%NB
        elseif (IT1N == wlc_p%Nb) then
           IT2N = 1
           IT1N = wlc_p%Nb-1
        ENDif
     ENDif

     if (DIB == wlc_p%Nb) then
        if (IT1N /= 1.AND.IT1N /= wlc_p%Nb) then
           IT2N = IT1N-1
           IT1N = IT1N + 1
        elseif (IT1N == 1) then
           IT1N = IT1N + 1
           IT2N = wlc_p%Nb
        elseif (IT1N == wlc_p%Nb) then
           IT2N =wlc_p%Nb-1
           IT1N = 1
        ENDif
     ENDif

  ENDif

  DIO = wlc_p%Nb-DIB
  II = IT1N
  inD = 1
  NCrossNew = 0
  CrossNew = 0.
  !Copy the old cross matrix to the new cross matrix. Only copy rows in which intersections do not involve
  !one of the segments moved

  do J = 1,NCross
     Copy = .TRUE.
     II = IT1N
     do I = 1,DIB
        if (II ==wlc_p%Nb + 1) then
           II = 1
        ENDif
        if (II == nint(Cross(J,1)).OR.II == nint(Cross(J,2))) then
           Copy = .FALSE.
           EXIT
        ENDif
        II = II + 1
     ENDdo

     if (Copy) then
        CrossNew(inD,:) = Cross(J,:)
        inD = inD + 1
        NCrossNew = NCrossNew + 1
     ENDif
  ENDdo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Determine all crossings involving the segment rotated. Current method using a brute force search for all
  !intersections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CrossNew(:,5) = 1.
  II = IT1N
  NCross = NCrossNew
  Cross = CrossNew
  !Loop over all segments within the portion of the chain that was moved

  do I = 1,DIB

     if(II == wlc_p%Nb + 1) then
        II = 1
        IIP1 = 2
     elseif (II == wlc_p%Nb) then
        IIP1 = 1
     else
        IIP1 = II + 1
     ENDif

     !Loop over all segments outside of the portion of the chain that was moved
     IO = IT2N
     do J = 1,DIO
        if(IO == wlc_p%Nb + 1) then
           IO = 1
           IOP1 = 2
        elseif (IO == wlc_p%Nb) then
           IOP1 = 1
        else
           IOP1 = IO + 1
        ENDif
        !Skip calculation for adjacent segments
        if (II == IOP1.OR.IO == IIP1.OR.II == IO) then
           GOTO 10
        ENDif

        !Calculate lengths of segments in the projection and the tangents
        smax = SQRT(SUM((RP(:,IIP1)-RP(:,II))**2))
        tmax = SQRT(SUM((RP(:,IOP1)-RP(:,IO))**2))
        ui = (RP(:,IIP1)-RP(:,II))/smax
        uj = ((RP(:,IOP1)-RP(:,IO)))/tmax
        udot = ui(1)*uj(1) + ui(2)*uj(2) + ui(3)*uj(3)

        !If segments are parallel, continue to next segment
        if (udot == 1..OR.udot == -1.) then
           GOTO 10
        ENDif

        !Compute the point of intersection
        tint = (RP(2,IO)-RP(2,II)-(ui(2)/ui(1))*(RP(1,IO)-RP(1,II)))/((ui(2)*uj(1)/ui(1))-uj(2))
        sint = (RP(1,IO)-RP(1,II) + uj(1)*tint)/ui(1)

        !If the intersection point is within length of segments, count as an intersection
        if (sint >= 0.AND.sint < smax.AND.tint >= 0.AND.tint < tmax) then
           !Determine if this is an undercrossing (RI under RJ) or overcrossing

           !Compute lengths and tangents  of true segments (not projected)
           srmax = SQRT(SUM((R(:,IIP1)-R(:,II))**2))
           trmax = SQRT(SUM((R(:,IOP1)-R(:,IO))**2))
           DRI = R(:,IIP1)-R(:,II)
           DRJ = R(:,IOP1)-R(:,IO)
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

           if (R(3,II) + uri(3)*srint<r(3,IO) + urj(3)*trint) then
              Ncross = Ncross + 1
              Cross(Ncross,1) = II;
              Cross(Ncross,2) = IO;
              Cross(Ncross,3) = sint;
              Cross(Ncross,4) = tint;
              Cross(Ncross,6) = -uj(2)*ui(1) + ui(2)*uj(1); !cross(ui,uj), + for RH, - for LH
           else
              Ncross = Ncross + 1
              Cross(Ncross,1) = IO
              Cross(Ncross,2) = II
              Cross(Ncross,3) = tint
              Cross(Ncross,4) = sint
              Cross(Ncross,6) = -(-uj(2)*ui(1) + ui(2)*uj(1)); !cross(ui,uj), + for RH, - for LH

           ENDif
        ENDif

10      continue
        IO = IO + 1
     ENDdo

     !Loop over all segments inside the portion that was moved, up to the segment that lies before the segment
     !immediately adjacent to the II segment in the outermost loop

     IO = IT1N

     do J = 1,I-2
        if(IO == wlc_p%Nb + 1) then
           IO = 1
           IOP1 = 2
        elseif (IO == wlc_p%Nb) then
           IOP1 = 1
        else
           IOP1 = IO + 1
        ENDif
        !Skip calculation for adjacent segments
        if (II == IOP1.OR.IO == IIP1.OR.II == IO) then
           GOTO 20
        ENDif

        !Calculate lengths of segments in the projection and the tangents
        smax = SQRT(SUM((RP(:,IIP1)-RP(:,II))**2))
        tmax = SQRT(SUM((RP(:,IOP1)-RP(:,IO))**2))
        ui = (RP(:,IIP1)-RP(:,II))/smax
        uj = ((RP(:,IOP1)-RP(:,IO)))/tmax
        udot = ui(1)*uj(1) + ui(2)*uj(2) + ui(3)*uj(3)

        !If segments are parallel, continue to next segment
        if (udot == 1..OR.udot == -1.) then
           GOTO 20
        ENDif

        !Compute the point of intersection
        tint = (RP(2,IO)-RP(2,II)-(ui(2)/ui(1))*(RP(1,IO)-RP(1,II)))/((ui(2)*uj(1)/ui(1))-uj(2))
        sint = (RP(1,IO)-RP(1,II) + uj(1)*tint)/ui(1)

        !If the intersection point is within length of segments, count as an intersection
        if (sint >= 0.AND.sint < smax.AND.tint >= 0.AND.tint < tmax) then
           !Determine if this is an undercrossing (RI under RJ) or overcrossing

           !Compute lengths and tangents  of true segments (not projected)
           srmax = SQRT(SUM((R(:,IIP1)-R(:,II))**2))
           trmax = SQRT(SUM((R(:,IOP1)-R(:,IO))**2))
           DRI = R(:,IIP1)-R(:,II)
           DRJ = R(:,IOP1)-R(:,IO)
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

           if (R(3,II) + uri(3)*srint<r(3,IO) + urj(3)*trint) then
              Ncross = Ncross + 1
              Cross(Ncross,1) = II;
              Cross(Ncross,2) = IO;
              Cross(Ncross,3) = sint;
              Cross(Ncross,4) = tint;
              Cross(Ncross,6) = -uj(2)*ui(1) + ui(2)*uj(1); !cross(ui,uj), + for RH, - for LH
           else
              Ncross = Ncross + 1
              Cross(Ncross,1) = IO
              Cross(Ncross,2) = II
              Cross(Ncross,3) = tint
              Cross(Ncross,4) = sint
              Cross(Ncross,6) = -(-uj(2)*ui(1) + ui(2)*uj(1)); !cross(ui,uj), + for RH, - for LH

           ENDif
        ENDif

20      continue
        IO = IO + 1
     ENDdo

     II = II + 1
  ENDdo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Continue with the alexander polynomial calculation as usual. The remainder is no different
  !from the subroutine that calclates the alexander polynomial from scratch (i.e. w/o a Cross matrix
  !from the previous move.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
           GOTO 30
           !If J = K, then segment j contains undercrossings
           ! then need to determine where overpass lies relative to undercrossing
        elseif (J == nint(CROSS(K,1))) then
           t = Cross(I,4)
           Ndegen = Cross(K,5)
           !Special case: t is before the first under-pass or after the last under-pass
           !of segment j
           if (t <= Cross(K,3)) then
              over_ind(I) = K
              GOTO 30
           elseif (t >= Cross(K + Ndegen-1,3)) then
              OVER_inD(I) = K + Ndegen
              GOTO 30
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
                    GOTO 30
                 ENDif
                 inD = inD + 1
              ENDdo
           ENDif
        ENDif

     ENDdo
30   continue

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

  !Deallocate Arrays

  DEallocate(A)
  DEallocate(over_ind)

  RETURN

END subroutine alexanderp_crank
