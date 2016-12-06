subroutine stop_if_err(err, msg)
    implicit none
    integer err
    character(len=*), intent(in) :: msg

    if (err /= 0) then
        print *, msg
        stop
    endif

end subroutine stop_if_err
