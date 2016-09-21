!        Program: Example Inx.f90
!
!        Brian Geelen.   Created: november 2001.
!
!        Example use of Module Bessf_In_x                          
!
	PROGRAM example_Inx

	USE Bessf_In_x
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: i15 = SELECTED_REAL_KIND(15,300)
	REAL (KIND = i15)  :: x, Bess_In

	INTEGER, PARAMETER  :: i9 = SELECTED_INT_KIND(9)
    INTEGER (KIND = i9) :: n, loop_1  
	
    n = -1
    x =  1.0D0
		

        DO loop_1 = 1, 21

            n = n + 1
             
	        CALL Bessf_In(x, n, Bess_In) 
    
	        		
		    WRITE(*,100)x, n, Bess_In
	    100 FORMAT(2x, F5.1, 4x, I3, 4x, ES16.9)		
	
	    END DO
	
	STOP
	END PROGRAM example_Inx


	MODULE Bessf_In_x
!
!   Module Program: Bessf_In.f90
!   FN: Bessf_In    Version 1; Created: november 2001.
!
!   Purpose: evaluation of the modifed Bessel function of the first kind,
!   In(x), where n, an integer value, over the interval 0 to 100, and
!   x, is a real positive argument, over the interval 1 to 100. Calculation 
!   procedure is based on convergent power series and asymptotic series.
!
!   Accuracy: returns precision of at least 9 significant decimal places.
!	
!   Reference: Geelen, B., "Accurate solution for the modified Bessel 
!   function of the first kind", Advances in Engineering Software, No 23, 
!   1995, pp. 105-109.
!
!   Uses apart subroutines for I0(x), I1(x) and In(x) where n = or > 2.
!
!   Brian Geelen; PB10416, 6000GK Weert, The Netherlands. Email: bwave@iae.nl
!
!   Example output result listing:
!   x        n               In(x)
!   1.0      0     1.266065878E+00
!   1.0      1     5.651591040E-01
!   1.0      2     1.357476698E-01
!   1.0      3     2.216842492E-02
!   1.0     20     3.966835986E-25
!   1.0    100     8.473674008E-189
!
!***
!***
!***
!   
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Bessf_In
    CONTAINS
	
    SUBROUTINE Bessf_In(x, n, Bess_In) 
	
    INTEGER, PARAMETER  :: i15 = SELECTED_REAL_KIND(15,300)
    INTEGER, PARAMETER  :: i9  = SELECTED_INT_KIND(9)
    
    REAL (KIND = i15), INTENT(IN)   :: x
    REAL (KIND = i15), INTENT(OUT)  :: Bess_In
    INTEGER (KIND = i9), INTENT(IN) :: n

        IF (n < 1) THEN

            CALL Bessf_I0(x, Bess_In)
            
        ELSEIF( n < 2) THEN

            CALL Bessf_I1(x, Bess_In)
            
        ELSE

            CALL Bessf_In_gt2(x, n, Bess_In)
        
        ENDIF

    END SUBROUTINE Bessf_In
!***
!******************************************************************************
!***
	SUBROUTINE Bessf_I0(x, Bess_I0)
!
!   Subroutine Program: Bessf_I0.f90
!   FN: Bessf_I0_ver1                  Created: november 2001.
!
!   Version 1; Brian Geelen.
!
!   Purpose: evaluation of the modifed Bessel function, I0(x),
!   where x, is a positive real positive argument, over the interval 1 to 100.
!
!   Accuracy: returns precision of at least 9 significant decimal places.
!	
!   Convergent power series solution for x < 14, asymptotic series for x > 14.
!
!   
!	precision and variable declarations

    IMPLICIT NONE
    INTEGER, PARAMETER  :: i15 = SELECTED_REAL_KIND(15,300)
    REAL (KIND = i15)   :: ak, px, an, series_delta, series_rslt,    &
                           odd, seq_num, ser_num, fact, b1, b2, b3

	INTEGER, PARAMETER  :: i9 = SELECTED_INT_KIND(9)
	INTEGER (KIND = i9) :: loop_1, loop_2

    REAL (KIND = i15), PARAMETER ::	e = 2.718281828459045D0,         &
	                   pi = 3.141592653589793D0, error = 0.5D-12
	
    REAL (KIND = i15), INTENT(IN)  :: x
    REAL (KIND = i15), INTENT(OUT) :: Bess_I0
       
!   select power series or asympotic series solution

	   expression_select: IF (x < 14.0) THEN              ! power series for x < 14
		
          ak = 0.0D0; an = 1.0D0; series_rslt = 0.0D0     ! set variables
 
            series_solution: DO loop_1 = 1, 50      

                ak             = ak + 2.0D0
                px             = x **(ak)
                an             = an * (ak * ak)
                series_delta   = (px / an)
                series_rslt    = series_rslt + (px / an)
                Bess_I0        = 1.0D0 + series_rslt

				IF (series_delta < error) EXIT        ! exit if error term satisfied       
				
			END DO series_solution		              
	
	   ELSE		!  asymptotic series solution for x > 14

          odd   = 1.0D0; seq_num = 1.0D0; ser_num  = 1.0D0 ! set variables
          fact  = 2.0D0; an      = 2.0D0

          b1            = (e **x) / DSQRT(2.0D0 * pi * x)
          b2            = 8.0D0 * x
          series_rslt   = 1.0D0 + (1.0D0 / b2) + (9.0D0 / (2.0D0 * (b2**2.0D0) ) )

            asymptotic_series: DO loop_2 = 1,20

                odd             = odd + 2.0D0
                seq_num         = seq_num * (odd** 2.0D0) 
                ser_num         = seq_num * ( (odd + 2.0D0) **2.0D0)

                an              = an + 1.0D0
                fact            = fact * an
                b3              = (8.0D0 * x) **an

                series_delta    = (ser_num / (fact * b3) )
                series_rslt     = series_rslt + (ser_num / (fact * b3) )
                Bess_I0         = b1 * series_rslt
			
            END DO asymptotic_series		 
						
	   END IF expression_select			     

	END SUBROUTINE Bessf_I0
!***
!******************************************************************************
!***
    SUBROUTINE Bessf_I1(x, Bess_I1)
!
!   Subroutine Program: Bessf_I1
!   FN: Bessf_I1_ver1.f90              Created: november 2001.
!
!   Version 1; Brian Geelen.
!
!   Purpose: evaluation of the modified Bessel function, I1(x),
!   where x, is a real positive argument, over the interval 1 to 100
!
!   Accuracy: returns precision of at least 9 significant decimal places.
!  
!   Convergent power series solution for x < 14, asymptotic series for x > 14.
!
!   
!   precision and variable declarations

    IMPLICIT NONE
    INTEGER, PARAMETER  :: i15 = SELECTED_REAL_KIND(15,300)
    REAL (KIND = i15)   :: ak, an, fact1, sx, a2, fact2, series_delta, &
                           odd, seq_num, ser_num, c1, c2, c3, series
    
    INTEGER, PARAMETER  :: i9 = SELECTED_INT_KIND(9)
    INTEGER (KIND = i9) :: loop_1, loop_2

   	REAL (KIND = i15), PARAMETER ::	e = 2.718281828459045D0,           &
	                   pi = 3.141592653589793D0, error = 0.5D-12
	
    REAL (KIND = i15), INTENT(IN)  :: x
    REAL (KIND = i15), INTENT(OUT) :: Bess_I1

!   select power series or asympotic series solution

	   expression_select: IF (x < 14.0) THEN      ! power series solution for x < 14
		          		              
          ak = 1.0D0; an = 0.0D0; fact1 = 1.0D0   ! set variables
		      
          Bess_I1 = (x / 2.0D0)
            
            series_solution: DO loop_1 = 1,31
				
                ak    = ak + 2.0D0
                sx    = x **(ak)

                a2    = 2.0D0**(ak)
				
                an    = an + 1.0D0
                fact1 = fact1 * an
                fact2 = fact1 * (an + 1.0D0)
								
                Bess_I1       = ( sx / (a2 * fact1 * fact2 ) ) +  Bess_I1 
                series_delta  = ( sx / (a2 * fact1 * fact2 ) )
             
                IF (series_delta < error) EXIT    ! exit if error term satisfied 
			   
            END DO series_solution
		         
       ELSE      !  asymptotic series solution for x > 14

          odd    = 1.0D0; seq_num = 1.0D0; ser_num = 1.0D0 ! set variables
          an     = 2.0D0; fact1 = 2.0D0

          c1     = (e **x) / DSQRT(2.0D0 * pi * x)
		  c2     = 8.0D0 * x
		  series = 1.0D0 - (3.0D0 / c2) - (15.0D0 / (2.0D0 * (c2**2.0D0) ))
          
            asymptotic_series: DO loop_2 = 1,29

                odd      = odd + 2.0D0
                seq_num  = (odd * odd) * seq_num
                ser_num  = seq_num * (odd + 2.0D0) * (odd + 4.0D0)
				
                an       = an + 1.0D0
                fact1    = fact1 * an
				        
                c3       = (8.0D0 * x) **an
                series   = series - ( ser_num / ( fact1 * c3 ))
				
                Bess_I1  = series * c1

			END DO asymptotic_series		  
            						
	   END IF expression_select			     

	END SUBROUTINE Bessf_I1
!***
!******************************************************************************
!***
	Subroutine Bessf_In_gt2(x, n, Bess_In) 
!
!   Subroutine Program: Bessf_In.f90
!   FN: Bessf_In_GT2_ver1              Created: november 2001.
!
!   Version 1; Brian Geelen.
!
!   Purpose: evaluation of the modifed Bessel function, In(x),
!   where n, a integer value, is = or > 2, and x, a positive real
!   argument, over the interval 1 to 100.
!
!   Accuracy: returns precision of at least 9 significant decimal places.
!	
!   Employs a downward recurrence approach with nornmalisation to I0(x).
!	External call made to subroutine Bessf_I0
!   to solve the modifed Bessel function I0(x).
!
!
!	precision and variable declarations

	IMPLICIT NONE
	INTEGER, PARAMETER  :: i15 = SELECTED_REAL_KIND(15,300)
	REAL (KIND = i15)   :: a1, a2, a3, a4, Bess_I0

	INTEGER, PARAMETER  :: i9 = SELECTED_INT_KIND(9)
	INTEGER (KIND = i9) :: j1, m, k, iexp

	REAL (KIND = i15), INTENT(IN)  :: x
	REAL (KIND = i15), INTENT(OUT) :: Bess_In

	INTEGER (KIND = i9), INTENT(IN) :: n
	
	iexp = maxexponent(x)/2           ! initialisation parameter
	
	      
    a1 = 1.0D0; a2 = 0.0D0; a3 = 2.0D0 / x
    Bess_In = 0.0D0

    j1 = 200

    m = 2 * ( (n + INT (SQRT (real (j1 * n)) )))
    
		 recursion: DO k = m, 1, -1   ! perform downward recursion
			
			 a4 = a2 + k * a1 * a3 
			 a2 = a1
			 a1 = a4
				
				 IF (exponent(a1) > iexp) THEN
						
						 Bess_In = scale(Bess_In, -iexp) 
						 a1      = scale(a1, -iexp)
						 a2      = scale(a2, -iexp)

				 END IF
						 IF (k == n) Bess_In = a2 
	
		 END DO recursion

!   solve Bessel function I0(x) with subroutine call
!   and perform normalisation to I0(x)

    CALL Bessf_I0(x, Bess_I0)       

    Bess_In = Bess_In * Bess_I0 / a1 
      
    END Subroutine Bessf_In_gt2
!***
!******************************************************************************
!***
END MODULE Bessf_In_x
