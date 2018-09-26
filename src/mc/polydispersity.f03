#include "../defines.inc"
module polydispersity

implicit none
integer, allocatable, dimension(:):: chainID
integer, allocatable, dimension(:):: firstBead

contains

function first_bead_of_chain(IP) result(IB)
    implicit none
    integer, intent(in) :: IP
    integer IB
    if (WLC_P__POLY_DISP_TYPE == "None") then
        IB = (IP-1)*WLC_P__NB+1
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        IB = firstBead(IP)
    endif
end function

function last_bead_of_chain(IP) result(IB)
    implicit none
    integer, intent(in) :: IP
    integer IB
    if (WLC_P__POLY_DISP_TYPE == "None") then
        IB = IP*WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        IB = firstBead(IP+1)-1
    endif
end function

subroutine setup_polydispersity()
    implicit none
    integer length, IB, IP

    if (WLC_P__POLY_DISP_TYPE == "None") then
        continue
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        allocate(firstBead(WLC_P__NP+1))
        open (unit = 1, file = "input/polyLengths", status = 'OLD')
        firstBead(1) = 1
        IB = 1
        do IP = 1,WLC_P__NP
            read(1,"(I7)") length
            firstBead(IP+1) = firstBead(IP)+length
        enddo
        close(1)

        if (WLC_P__NT .ne. firstBead(WLC_P__NP+1) - 1) then
            print*, "Error, wrong number of total beads"
            print*, WLC_P__NT, ".ne.", firstBead(WLC_P__NP+1)
            stop
        endif

        allocate(chainID(WLC_P__NT))
        do IP = 1,WLC_P__NP
            do IB = firstBead(IP),firstBead(IP+1) - 1
                chainID(IB)=IP
            enddo
        enddo
    endif

end subroutine

function length_of_chain(IP) result(length)
    implicit none
    integer, intent(in) :: IP
    integer :: length


    if (WLC_P__POLY_DISP_TYPE == "None") then
        length = WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        length = firstBead(IP+1)-firstBead(IP)
    endif

end function

function length_of_chain_containing(I) result(length)
    implicit none
    integer, intent(in) :: I
    integer IP
    integer length

    if (WLC_P__POLY_DISP_TYPE == "None") then
        length = WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        IP = chainID(I)
        length = firstBead(IP+1)-firstBead(IP)
    endif
end function

function chain_ID(I) result(ID)
    implicit none
    integer, intent(in) :: I
    integer ID
    if (WLC_P__POLY_DISP_TYPE == "None") then
        ID = (I-1)/WLC_P__NB + 1
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        ID = chainID(I)
    endif
end function

function leftmost_from(I) result(left_end)
    implicit none
    integer, intent(in) :: I
    integer left_end
    if (WLC_P__POLY_DISP_TYPE == "None") then
        left_end = ((I-1)/WLC_P__NB)*WLC_P__NB + 1
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        left_end = firstBead(chainID(I))
    endif
end function

function rightmost_from(I) result(right_end)
    implicit none
    integer, intent(in) :: I
    integer right_end
    if (WLC_P__POLY_DISP_TYPE == "None") then
        right_end = ((I-1)/WLC_P__NB + 1)*WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        right_end = firstBead(chainID(I) + 1) - 1
    endif
end function

function get_IB(I) result(IB)
    implicit none
    integer, intent(in) :: I
    integer IB
    if (WLC_P__POLY_DISP_TYPE == "None") then
        IB = I - ((I-1)/WLC_P__NB)*WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        IB = I - firstBead(chainID(I)) + 1
    endif
end function

function get_I(IB,IP) result(I)
    implicit none
    integer, intent(in) :: IB
    integer, intent(in) :: IP
    integer I
    if (WLC_P__POLY_DISP_TYPE == "None") then
        I = IB + (IP-1)*WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        I = IB + firstBead(IP) - 1
    endif
end function

function get_IP(I) result(IP)
    implicit none
    integer, intent(in) :: I
    integer IP
    if (WLC_P__POLY_DISP_TYPE == "None") then
        IP = (I-1)/WLC_P__NB + 1
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        IP = chainID(I)
    endif
end function

function is_left_end(I) result(output)
    implicit none
    integer, intent(in) :: I
    logical output

    if (WLC_P__POLY_DISP_TYPE == "None") then
        output = MOD(I,WLC_P__NB) == 1
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        if (I == 1) then
            output = .True.
            return
        endif
        if (chainID(I-1) .ne. chainID(I)) then
            output = .True.
        else
            output = .False.
        endif
    endif
end function

function is_right_end(I) result(output)
    implicit none
    integer, intent(in) :: I
    logical output

    if (WLC_P__POLY_DISP_TYPE == "None") then
        output = MOD(I,WLC_P__NB) == 0
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        if (I == WLC_P__NT) then
            output = .True.
            return
        endif
        if (chainID(I+1) .ne. chainID(I)) then
            output = .True.
        else
            output = .False.
        endif
    endif
end function

! Test whether beads A and B are from the same chain.
! If ether A or B is not a bead return .False.
function are_on_same_chain(A,B) result(output)
    implicit none
    integer, intent(in) :: A
    integer, intent(in) :: B
    logical output
    if (A>WLC_P__NT .or. B>WLC_P__NT) then
        output = .False.
        return
    elseif (A<1 .or. B<1) then
        output = .False.
        return
    endif

    if (WLC_P__POLY_DISP_TYPE == "None") then
        output = ( (A-1)/WLC_P__NB == (B-1)/WLC_P__NB )
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        output = ( chainID(A) == chainID(B) )
    endif
end function

function n_mono_per_poly(IP) result(NMPP)
    implicit none
    integer, intent(in) :: IP
    integer NMPP

    if (WLC_P__POLY_DISP_TYPE == "None") then
        NMPP = WLC_P__NMPP ! a.k.a. WLC_P__NB/WLC_P__NBPM
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        NMPP = length_of_chain(IP)/WLC_P__NBPM
        if (MOD(length_of_chain(IP),WLC_P__NBPM) .ne. 0) then
            print*, "Warming: monomers don't fit in polymer"
        endif
    endif
end function

function max_chain_length() result(high)
    implicit none
    integer high
    integer IP

    if (WLC_P__POLY_DISP_TYPE == "None") then
        high = WLC_P__NB
    elseif (WLC_P__POLY_DISP_TYPE == "FromFile") then
        high = firstBead(2)-firstBead(1)
        do IP=2,WLC_P__NP
            high = max(high,firstBead(IP+1)-firstBead(IP))
        enddo
    endif
end function

end module
