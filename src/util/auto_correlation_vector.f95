!This subroutine generates a vector of auto-correlations of elements
!in a vector v. Each element in the autocorrelation vector corresponds
!to a particular delta (in the index) between successive elements of the vector
!for which the autocorrelation is computed. The auto_correlation is computed,
!beginning with delta=delta_min and ending with delta=delta_max

SUBROUTINE auto_correlation_vector(v,N,delta_min,delta_max,auto_vector)

  !Input variables
  DOUBLE PRECISION v(N)
  INTEGER N
  INTEGER delta_min
  INTEGER delta_max


  !Intermediate variables

  INTEGER Ndelta
  INTEGER  delta
  DOUBLE PRECISION auto

  !output variables

  DOUBLE PRECISION auto_vector(delta_max-delta_min+1)

  Ndelta=delta_max-delta_min+1

  DO I=0,Ndelta-1
     delta=delta_min+I
     CALL auto_correlation(v,N,delta,auto)
     auto_vector(I+1)=auto
  ENDDO

RETURN

ENDSUBROUTINE
