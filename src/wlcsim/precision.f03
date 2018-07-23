module precision
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: IEEE_ARITHMETIC
    implicit none
    public
    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over selected_real_kind(15,307)
                                      ! only available as of fortran 2008
    real(dp), parameter :: eps = 0.00000001_dp
    real(dp), parameter :: epsApprox = 0.001_dp  ! for compairison but allows for more ronding error
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
end module precision
