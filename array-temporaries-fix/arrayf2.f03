module arrayf
contains
    function f(R)
        real, intent(in) :: R(3)
        real f(3)
        f = R + 1
    end function
end module
