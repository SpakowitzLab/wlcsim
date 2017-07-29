module precision
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: IEEE_ARITHMETIC
    implicit none
    public
    ! precision of simulations
    integer, parameter :: dp = real64 ! preferred over selected_real_kind(15,307)
                                      ! only available as of fortran 2008
end module precision
