subroutine stop_if_err(err, msg)
    use params, only: MAXFILENAMELEN
    implicit none
    integer err
    character(MAXFILENAMELEN) msg

    if (err /= 0) then
        print *, msg
        stop
    endif

end subroutine stop_if_err
