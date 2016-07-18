      SUBROUTINE INSERTION_SORT(N,A)
      IMPLICIT NONE
      INTEGER N,I,J
      DOUBLE PRECISION A(N),X
      DO 30 I=2,N
      X=A(I)
      J=I
10    J=J-1
      IF(J.EQ.0 .OR. A(J).LE.X) GO TO 20
      A(J+1)=A(J)
      GO TO 10
20    A(J+1)=X
30    CONTINUE
      END

      SUBROUTINE CHECK_COLLISIONS(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      INTEGER NT, COL_TYPE
      DOUBLE PRECISION FPT_DIST, TIME
      DOUBLE PRECISION R(NT,3), HAS_COLLIDED(NT, NT)
      if (COL_TYPE.EQ.0) then
         return
      else if (COL_TYPE.EQ.1) then
         call CHECK_COLLISIONS_BRUTE(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      else if (COL_TYPE.EQ.2) then
         call CHECK_COLLISIONS_KD(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      else if (COL_TYPE.EQ.3) then
         call CHECK_COLLISIONS_BB(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      end if
      END

      SUBROUTINE CHECK_COLLISIONS_BRUTE(R, NT, HAS_COLLIDED, FPT_DIST, &
          TIME, COL_TYPE)
      INTEGER NT, COL_TYPE
      DOUBLE PRECISION FPT_DIST, TIME
      DOUBLE PRECISION R(NT,3), HAS_COLLIDED(NT, NT)
!     Check if the particles have collided
          DO 140 K1 = 1, NT
              DO 150 K2 = 1, NT
                  IF (HAS_COLLIDED(K1,K2).LT.0.0d0 .AND. K1.NE.K2 &
                        .AND. abs(R(K1,1) - R(K2,1)) < FPT_DIST &
                        .AND. abs(R(K1,2) - R(K2,2)) < FPT_DIST &
                        .AND. abs(R(K1,3) - R(K2,3)) < FPT_DIST) THEN
                     HAS_COLLIDED(K1,K2) = TIME
                  END IF
150           CONTINUE
140       CONTINUE

      END

      SUBROUTINE CHECK_COLLISIONS_KD(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create, kdtree2_r_nearest_around_point

      INTEGER NT, COL_TYPE, NFOUND, NALLOC, K1, K2, I
      DOUBLE PRECISION FPT_DIST, TIME
      DOUBLE PRECISION R(NT,3), HAS_COLLIDED(NT, NT)
      type(kdtree2), pointer :: col_tree
      type(kdtree2_result), allocatable :: kd_results(:)

      col_tree => kdtree2_create(R, rearrange = .true., sort = .false.)
      do K1 = 1,NT
         call kdtree2_r_nearest_around_point(col_tree, idxin = K1, &
               correltime = 1, r2 = FPT_DIST, nfound = NFOUND, nalloc = NALLOC, &
               results = kd_results)
         do I = 1,NFOUND
             K2 = kd_results(I)%idx
             if (HAS_COLLIDED(K1,K2) .LT. 0) then
                 HAS_COLLIDED(K1,K2) = TIME
             endif
         enddo
      enddo
      END

      SUBROUTINE CHECK_COLLISIONS_BB(R, NT, HAS_COLLIDED, FPT_DIST, TIME, COL_TYPE)
      INTEGER NT, COL_TYPE
      DOUBLE PRECISION FPT_DIST, TIME
      DOUBLE PRECISION R(NT,3), HAS_COLLIDED(NT, NT)
         WRITE(*,*) "Not yet implemented!"
         STOP 1
      END
