
module mymod
contains
   function fun()
        implicit none
        real, dimension(4) :: fun
        fun(1)=0.000
        fun(2)=0.005
        fun(3)=0.005
        fun(4)=0.000
    end function fun
    subroutine sub(input)
        implicit none
        real, intent(inout) :: input(4)
        input(1)=input(1)+0.000
        input(2)=input(2)+0.005
        input(3)=input(3)+0.005
        input(4)=input(4)+0.000
    end subroutine
end module
program main

    use mymod, only : fun, sub

    integer i,nloops
    real timeval, tfinal
    real, dimension(4) :: output
    real, dimension(4) :: temp
    output = 1.0

    nloops =1000000
    ! vvvvvvv Doesn't cause warning vvvvv
    call cpu_time(timeval)
    do i=1,nloops
        call sub(output)
    enddo
    call cpu_time(tfinal)
    print*, "subroutine took",tfinal-timeval

    ! vvvvvvv Causes warning vvvvvvvv
    call cpu_time(timeval)
    do i=1,nloops
        output = fun()
        output = output + temp
    enddo
    call cpu_time(tfinal)
    print*, "function took",tfinal-timeval

    ! vvvvvvv Doesn't cause warning vvvvvv
    output = fun()

    ! vvvvvvv Doesn't cause warning vvvvv
    call cpu_time(timeval)
    do i=1,nloops
        call sub(output)
    enddo
    call cpu_time(tfinal)
    print*, "subroutine took",tfinal-timeval

    PRINT *, output

end program main

