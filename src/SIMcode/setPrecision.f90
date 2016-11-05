! This is for universally setting the precison of constants
Module setPrecision 
    IMPLICIT NONE
    INTEGER, Parameter :: dp = SELECTED_REAL_KIND(15,307)
end module
