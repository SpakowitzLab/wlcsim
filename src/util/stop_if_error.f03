subroutine stop_if_error(err, msg)
    use params, only: MAXFILENAMELEN
    implicit none
    integer err
    character(MAXFILENAMELEN) msg

    if (err /= 0) then
        print *, msg
        stop
    endif

end subroutine stop_if_error
