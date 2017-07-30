!This subroutine calculates the auto-correlation in
!the entries of a vector spaced by a index difference, delta
!a vector of auto-correlations is generated, beginning with the auto-correlation
!at delta_min and ending with the autocorrelation at delta_max

subroutine auto_correlation(v,N,delta,auto)

  use params, only : dp

  !Input variables
  real(dp) v(N) !vector
  integer N             !length of vector
  integer, intent(in) :: delta         !spacing between elements
  !Subroutine variables
  real(dp) v_avg  !average value of v
  real(dp) var    !variance in v
  real(dp) vv(N-delta)  !vector of products of fluctions

  !Output variables
  real(dp) auto   !autocorrelation

 ! allocate(vv(N-delta))

  v_avg = SUM(v)/N
  var = SUM((v-v_avg)**2)/N

  do I = 1,N-delta
     vv(I) = (v(I)-v_avg)*(v(I + delta)-v_avg)

  ENDdo

  auto = SUM(vv)/((N-delta)*var)

RETURN
ENDsubroutine auto_correlation


