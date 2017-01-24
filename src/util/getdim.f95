!This subroutine calculates the center of mass and radius of gyrations of
!a set of polymers.

!Brad Krajina
!Written 4-27-14

SUBROUTINE getdim(N,NP,NT,R,RCOM,RCOMSQ,RGYRSQ,R2,R4,R6,DR)

  INTEGER, INTENT(IN) :: N !Number of beads per chain
  INTEGER, INTENT(IN) :: NP !Number of chains
  INTEGER, INTENT(IN) :: NT !Total number of beads
  DOUBLE PRECISION, INTENT(IN) :: R(NT,3) !Bead positions
  DOUBLE PRECISION, INTENT(OUT) ::  RCOM(NP,3) ! Center of mass matrix
  DOUBLE PRECISION, INTENT(OUT) ::  RCOMSQ(NP) ! Center of mass squared
  DOUBLE PRECISION, INTENT(OUT) ::  RGYRSQ(NP) ! Radius of gyration squared
  DOUBLE PRECISION X(N,3) ! Bead positions of particular polymer
  DOUBLE PRECISION  XSQ(NP) !Sum of sqaures of position vectors
  DOUBLE PRECISION XGYR(N,3) ! Bead positions relative to center of mass
  DOUBLE PRECISION R2(NP) !Squared end-to-end displacement
  DOUBLE PRECISION R4(NP) !4th order moment
  DOUBLE PRECISION R6(NP) !6th order moment
  DOUBLE PRECISION DR(NP) !Magnitude of end-to-end distance

  INTEGER IP !Polymer Index
  INTEGER IB1 ! First bead on polymer (total index)
  INTEGER IBN ! Final bead on polymer (total index)

  RCOM=0.0d0
  RCOMSQ=0.0d0
  RGYRSQ=0.0d0

  !Calculate center of mass and radius of gyration squared

  DO IP=1,NP
     IB1=N*(IP-1)+1
     IBN=N*IP
     X(:,1)=R(IB1:IBN,1)
     X(:,2)=R(IB1:IBN,2)
     X(:,3)=R(IB1:IBN,3)
     RCOM(IP,:)=SUM(X,DIM = 1)/REAL(N)
     XGYR(:,1)=X(:,1)-RCOM(IP,1)
     XGYR(:,2)=X(:,2)-RCOM(IP,2)
     XGYR(:,3)=X(:,3)-RCOM(IP,3)
     RGYRSQ(IP)=SUM(XGYR**2)/N
     R2(IP)=SUM((R(IB1,:)-R(IBN,:))**2)
     R4(IP)=SUM((R(IB1,:)-R(IBN,:))**4)
     R6(IP)=SUM((R(IB1,:)-R(IBN,:))**4)
     DR(IP)=SQRT(SUM((R(IB1,:)-R(IBN,:))**2))

  ENDDO


  !Calculate square of center of mass for each polymer

  RCOMSQ=SUM(RCOM**2,DIM = 2)



  RETURN

END SUBROUTINE








