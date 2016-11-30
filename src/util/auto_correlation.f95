!This subroutine calculates the auto-correlation in
!the entries of a vector spaced by a index difference, delta
!a vector of auto-correlations is generated, beginning with the auto-correlation
!at delta_min and ending with the autocorrelation at delta_max

SUBROUTINE auto_correlation(v,N,delta,auto)

  !Input variables
  DOUBLE PRECISION v(N) !vector
  INTEGER N             !length of vector
  INTEGER, INTENT(IN) :: delta         !spacing between elements
  !Subroutine variables
  DOUBLE PRECISION v_avg  !average value of v
  DOUBLE PRECISION var    !variance in v
  DOUBLE PRECISION vv(N-delta)  !vector of products of fluctions

  !Output variables
  DOUBLE PRECISION auto   !autocorrelation

 ! ALLOCATE(vv(N-delta))

  v_avg=SUM(v)/N
  var=SUM((v-v_avg)**2)/N

  DO I=1,N-delta
     vv(I)=(v(I)-v_avg)*(v(I+delta)-v_avg)

  ENDDO

  auto=SUM(vv)/((N-delta)*var)

RETURN
ENDSUBROUTINE auto_correlation


