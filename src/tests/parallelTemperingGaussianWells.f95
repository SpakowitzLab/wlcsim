Subroutine simpleSim(rand_stat)

use mersenne_twister
use params
implicit none
type(random_stat), intent(inout) :: rand_stat ! state of random number generator
Type(wlcsim_params) mc
TYPE(MCData) md
double precision r(8)
double precision rp(8)
integer term
real urnd(1) ! single random number
double precision energy, prob, test
character*16 fileName
character*32 fullName
logical isfile
integer istep
double precision chiOld, xchiOld, EchiOld

fileName='input/params'
call wlcsim_params_setParams(mc,fileName)
call wlcsim_params_allocate(mc,md)
call PT_override(mc,md)

do term=1,8
    r(term)=0.0_dp
enddo

ISTEP=-1
mc%IND=1
do while (mc%IND.le.mc%INDMAX)
do ISTEP=1,mc%NSTEP
    if (abs(mc%ECHI-mc%CHI*(r(1)**2)).gt.0.001) then
        print*, "-----------------------"
        print*, "Energy mismatch"
        print*, "mc%ECHI=",mc%ECHI
        print*, "mc%CHI*r(1)**2",mc%CHI*r(1)**2
        print*, "IND", mc%IND, "  ISTEP", ISTEP
        print*, "r1", r(1), " rp(1)", rp(1)
        print*, "chiOld", chiOld, " CHI", mc%CHI
        print*, "xchiOld", xchiOld, " mc%x_chi",mc%x_chi
        print*, "EChiOld",EchiOld," Echi", mc%Echi
        stop 1
    endif
    if (mc%ECHI.lt.-0.00000001_dp) then
        print*, "--------------"
        print*, "negitive EChi"
        print*, "mc%ECHI", mc%ECHI
        print*, "IND", mc%IND, "  ISTEP", ISTEP
        print*, "r(1)", r(1), "rp(1)", rp(1)
        stop 1
    endif
    mc%ECHI=mc%CHI*(r(1)**2)
    ! Make Move
    do term=1,8
        if (mod(istep,8)+1.eq.term) then
            call random_number(urnd,rand_stat)
            rp(term)=r(term)+urnd(1)-0.5_dp
        else
            rp(term)=r(term)
        endif
    enddo

    ! Calculate change in energy

    mc%dx_chi=   (rp(1)**2-r(1)**2)
    mc%dx_mu=    (rp(2)**2-r(2)**2)
    mc%dx_Field= (rp(3)**2-r(3)**2)
    mc%dx_couple=-(rp(4)**2-r(4)**2)
    mc%dx_kap=   (rp(5)**2-r(5)**2)

    mc%DEBind=mc%mu*(rp(2)**2-r(2)**2) + (rp(2)-r(2))

    mc%dx_chi=mc%dx_chi*mc%CHI_ON
    mc%dx_couple=mc%dx_couple*mc%Couple_ON
    mc%dx_Kap=mc%dx_Kap*mc%KAP_ON

    mc%DEChi=mc%Chi*        mc%dx_chi
    mc%DECouple=mc%HP1_Bind*mc%dx_couple
    mc%DEKap=mc%Kap*        mc%dx_Kap
    mc%DEField=mc%h_A*      mc%dx_Field

    mc%deelas(1)=mc%Para(1)*(rp(6)**2-r(6)**2)
    mc%deelas(2)=mc%Para(2)*(rp(7)**2-r(7)**2)
    mc%deelas(3)=mc%Para(3)*(rp(8)**2-r(8)**2)

    ! accept or reject
    ENERGY=mc%DEELAS(1)+mc%DEELAS(2)+mc%DEELAS(3) &
                         +mc%DEKap+mc%DECouple+mc%DEChi+mc%DEBind+mc%DEField
    PROB=exp(-ENERGY)
    call random_number(urnd,rand_stat)
    TEST=urnd(1)
    if (TEST.LE.PROB) then
        r=rp
        mc%EBind=mc%EBind+mc%DEBind
        mc%x_mu=mc%x_mu+mc%dx_mu
        mc%EELAS(1)=mc%EELAS(1)+mc%DEELAS(1)
        mc%EELAS(2)=mc%EELAS(2)+mc%DEELAS(2)
        mc%EELAS(3)=mc%EELAS(3)+mc%DEELAS(3)
        mc%ECouple=mc%ECouple+mc%DECouple
        mc%EKap=mc%EKap+mc%DEKap
        mc%EChi=mc%EChi+mc%DEChi
        mc%EField=mc%EField+mc%DEField
        mc%x_Couple=mc%x_couple+mc%dx_couple
        mc%x_kap=mc%x_Kap+mc%dx_kap
        mc%x_chi=mc%x_chi+mc%dx_chi
        mc%x_field=mc%x_field+mc%dx_field

    endif
    if (abs(mc%EChi-mc%x_chi*mc%Chi).gt.0.0000001_dp) then
        print*, "~~~~~~~~~~~~~~~"
        print*, "Error. mc%Echi", mc%Echi," mc%x_chi*mc%Chi",mc%x_chi*mc%Chi
        stop 1
    endif
    chiOld=mc%Chi
    xchiOld=mc%x_chi
    EchiOld=mc%EChi
    IF ((mod(ISTEP,4)).eq.0) THEN
        call replicaExchange(mc)
    ENDIF


    ! output
    if ((mod(istep+mc%NSTEP*(mc%IND-1),1000).eq.0).and.(mc%IND.gt.mc%indEndRepAdapt)) then
        ! Record position
        fullName='data/r' // (mc%repSufix)
        fullName=trim(fullName)
        inquire(file = fullName, exist=isfile)
        if (isfile) then
            OPEN (UNIT = 1, FILE = fullName, STATUS = 'OLD', POSITION="append")
        else
            OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')

            fullName='data/cof' // (mc%repSufix)
            inquire(file = fullName, exist=isfile)
            if (.not.isfile) then
                OPEN (UNIT = 2, FILE = fullName, STATUS = 'NEW')
                write(2,*), mc%chi, mc%mu, mc%h_A, mc%HP1_Bind,mc%kap,&
                            mc%Para(1),mc%Para(2), mc%Para(3), mc%id
                close(2)
            endif
        endif
        write(1,*) r,mc%ind, mc%NSTEP*(mc%IND-1)+ISTEP
        Close(1)
    endif

enddo
print*, mc%IND
mc%IND=mc%IND+1
enddo
end subroutine
