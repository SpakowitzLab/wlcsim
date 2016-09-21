SUBROUTINE Bessf_I0(x, Bess_I0)
! --------------------------------------------
! EFK: copied from bessfin.f90
! downloaded from: http://www.netlib.org/a/sf/bessfin.f90
! --------------------------------------------
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
