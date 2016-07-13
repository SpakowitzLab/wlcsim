! file: test_sort.f
program test_sort
implicit none
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)
real(dp) :: r_min(3) = (/ 0.2_dp, 0.1_dp, 0.3_dp /)
integer :: ind(3) = (/ 1, 2, 3 /)
write(*,*) r_min
write(*,*) ind
call sortp_1r8(3, ind, r_min)
write(*,*) r_min
write(*,*) ind

end program
