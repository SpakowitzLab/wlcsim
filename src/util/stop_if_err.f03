subroutine stop_if_err(err, msg)
#ifdef f2003
    use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
#else
#define stderr 0
#endif
    implicit none
    integer err
    character(len=*), intent(in) :: msg

    if (err /= 0) then
        write(stderr,*) msg
        stop
    endif

end subroutine stop_if_err
