module sim_params
    implicit none

! code-specific parameters
    integer, parameter :: dp = selected_real_kind(15, 307)

! simulation parameters
    real(dp), parameter :: lp = 1.0_dp
    real(dp), parameter :: l = 100.0_dp
    real(dp), parameter :: lbox = 30.0_dp
    real(dp), parameter :: rend = 0.0_dp
    real(dp), parameter :: vhc = 0.0_dp
    integer, parameter :: n = 101
    integer, parameter :: np = 1
    real(dp), parameter :: tf = 1.0_dp
    integer, parameter :: indmax = 100
    real(dp), parameter :: dt = 0.01_dp
    integer, parameter :: frmfile = 0
    integer, parameter :: brown = 1
    integer, parameter :: inton = 0
    integer, parameter :: logtime = 0
    integer, parameter :: ninit = 1000000
    integer, parameter :: nstep = 0
    real(dp), parameter :: fpt_dist = 2.0_dp
    integer, parameter :: col_type = 1

! universal constants
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197169_dp
    real(dp), parameter :: e = 2.7182818284590452353602874713526624977572471_dp
end module
