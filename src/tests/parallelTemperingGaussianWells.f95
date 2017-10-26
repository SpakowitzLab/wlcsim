Subroutine simpleSim(rand_stat)

use mersenne_twister
use params
implicit none
type(random_stat), intent(inout) :: rand_stat ! state of random number generator
Type(wlcsim_params) wlc_p
TYPE(MCData) wlc_d
real(dp) r(8)
real(dp) rp(8)
integer term
real urnd(1) ! single random number
real(dp) energy, prob, test
character*16 fileName
character*32 fullName
logical isfile
integer istep
real(dp) chiOld, xchiOld, EchiOld

fileName = 'input/params'
call wlcsim_params_setParams(mc,fileName)
call wlcsim_params_allocate(wlc_p,wlc_d)
call PT_override(wlc_p,wlc_d)

do term = 1,8
    r(term) = 0.0_dp
enddo

ISTEP = -1
wlc_p%inD = 1
do while (wlc_p%inD.le.wlc_p%inDMAX)
do ISTEP = 1,wlc_p%NSTEP
    if (abs(wlc_p%ECHI-wlc_p%CHI*(r(1)**2)).gt.0.001) then
        print*, "-----------------------"
        print*, "Energy mismatch"
        print*, "wlc_p%ECHI = ",wlc_p%ECHI
        print*, "wlc_p%CHI*r(1)**2",wlc_p%CHI*r(1)**2
        print*, "inD", wlc_p%inD, "  ISTEP", ISTEP
        print*, "r1", r(1), " rp(1)", rp(1)
        print*, "chiOld", chiOld, " CHI", wlc_p%CHI
        print*, "xchiOld", xchiOld, " wlc_p%x_chi",wlc_p%x_chi
        print*, "EChiOld",EchiOld," Echi", wlc_p%Echi
        stop 1
    endif
    if (wlc_p%ECHI.lt.-0.00000001_dp) then
        print*, "--------------"
        print*, "negitive EChi"
        print*, "wlc_p%ECHI", wlc_p%ECHI
        print*, "inD", wlc_p%inD, "  ISTEP", ISTEP
        print*, "r(1)", r(1), "rp(1)", rp(1)
        stop 1
    endif
    wlc_p%ECHI = wlc_p%CHI*(r(1)**2)
    ! Make Move
    do term = 1,8
        if (mod(istep,8) + 1.eq.term) then
            call random_number(urnd,rand_stat)
            rp(term) = r(term) + urnd(1)-0.5_dp
        else
            rp(term) = r(term)
        endif
    enddo

    ! Calculate change in energy

    wlc_p%dx_chi=   (rp(1)**2-r(1)**2)
    wlc_p%dx_mu=    (rp(2)**2-r(2)**2)
    wlc_p%dx_Field= (rp(3)**2-r(3)**2)
    wlc_p%dx_couple = -(rp(4)**2-r(4)**2)
    wlc_p%dx_kap=   (rp(5)**2-r(5)**2)

    wlc_p%DEBind = wlc_p%mu*(rp(2)**2-r(2)**2) + (rp(2)-r(2))

    wlc_p%dx_chi = wlc_p%dx_chi*wlc_p%CHI_ON
    wlc_p%dx_couple = wlc_p%dx_couple*wlc_p%Couple_ON
    wlc_p%dx_Kap = wlc_p%dx_Kap*wlc_p%KAP_ON

    wlc_p%DEChi = wlc_p%Chi*        wlc_p%dx_chi
    wlc_p%DECouple = wlc_p%HP1_Bind*wlc_p%dx_couple
    wlc_p%DEKap = wlc_p%Kap*        wlc_p%dx_Kap
    wlc_p%DEField = wlc_p%h_A*      wlc_p%dx_Field

    ! wlc_p%deelas(1) = wlc_p%Para(1)*(rp(6)**2-r(6)**2)
    ! wlc_p%deelas(2) = wlc_p%Para(2)*(rp(7)**2-r(7)**2)
    ! wlc_p%deelas(3) = wlc_p%Para(3)*(rp(8)**2-r(8)**2)

    ! accept or reject
    ENERGY = wlc_p%DEELAS(1) + wlc_p%DEELAS(2) + wlc_p%DEELAS(3) &
                         +wlc_p%DEKap + wlc_p%DECouple + wlc_p%DEChi + wlc_p%DEBind + wlc_p%DEField
    PROB = exp(-ENERGY)
    call random_number(urnd,rand_stat)
    TEST = urnd(1)
    if (TEST <= PROB) then
        r = rp
        wlc_p%EBind = wlc_p%EBind + wlc_p%DEBind
        wlc_p%x_mu = wlc_p%x_mu + wlc_p%dx_mu
        wlc_p%EELAS(1) = wlc_p%EELAS(1) + wlc_p%DEELAS(1)
        wlc_p%EELAS(2) = wlc_p%EELAS(2) + wlc_p%DEELAS(2)
        wlc_p%EELAS(3) = wlc_p%EELAS(3) + wlc_p%DEELAS(3)
        wlc_p%ECouple = wlc_p%ECouple + wlc_p%DECouple
        wlc_p%EKap = wlc_p%EKap + wlc_p%DEKap
        wlc_p%EChi = wlc_p%EChi + wlc_p%DEChi
        wlc_p%EField = wlc_p%EField + wlc_p%DEField
        wlc_p%x_Couple = wlc_p%x_couple + wlc_p%dx_couple
        wlc_p%x_kap = wlc_p%x_Kap + wlc_p%dx_kap
        wlc_p%x_chi = wlc_p%x_chi + wlc_p%dx_chi
        wlc_p%x_field = wlc_p%x_field + wlc_p%dx_field

    endif
    if (abs(wlc_p%EChi-wlc_p%x_chi*wlc_p%Chi).gt.0.0000001_dp) then
        print*, "~~~~~~~~~~~~~~~"
        print*, "Error. wlc_p%Echi", wlc_p%Echi," wlc_p%x_chi*wlc_p%Chi",wlc_p%x_chi*wlc_p%Chi
        stop 1
    endif
    chiOld = wlc_p%Chi
    xchiOld = wlc_p%x_chi
    EchiOld = wlc_p%EChi
    if ((mod(ISTEP,4)).eq.0) then
        call replicaExchange(mc)
    ENDif


    ! output
    if ((mod(istep + wlc_p%NSTEP*(wlc_p%inD-1),1000).eq.0).and.(wlc_p%inD.gt.wlc_p%indEndRepAdapt)) then
        ! Record position
        fullName = 'data/r' // (wlc_p%repSufix)
        fullName = trim(fullName)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (UNIT = 1, FILE = fullName, STATUS = 'OLD', POSITION = "append")
        else
            open (UNIT = 1, FILE = fullName, STATUS = 'NEW')

            fullName = 'data/cof' // (wlc_p%repSufix)
            inquire(file = fullName, exist = isfile)
            if (.not.isfile) then
                open (UNIT = 2, FILE = fullName, STATUS = 'NEW')
                write(2,*), wlc_p%chi, wlc_p%mu, wlc_p%h_A, wlc_p%HP1_Bind,wlc_p%kap,&
                            ! wlc_p%Para(1),wlc_p%Para(2), wlc_p%Para(3),
                            wlc_p%id
                close(2)
            endif
        endif
        write(1,*) r,wlc_p%ind, wlc_p%NSTEP*(wlc_p%inD-1) + ISTEP
        Close(1)
    endif

enddo
print*, wlc_p%inD
wlc_p%inD = wlc_p%inD + 1
enddo
end subroutine
