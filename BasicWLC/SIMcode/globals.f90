module globals
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public :: pi, dp

    integer, parameter :: dp = REAL64
    real(dp), parameter :: pi = 3.1415926535897931_dp

endmodule
