#include "../defines.inc"
module nucleosome
!   --------------------------------------------------------------
!
!    This module is designed to handle the various geometrical
!    considerations of nucleosomes.
!
!   --------------------------------------------------------------

    use precision, only: dp, pi
    implicit none
    real(dp), parameter :: basePairsPerTurn = 10.5_dp
    real(dp), dimension(10,WLC_P__MAX_BACEPAIRS_PER_BEAD) :: multiParams
    real(dp), dimension(3,3,147) :: nucleosomeROT
    real(dp), dimension(3,147) :: nucleosomeTran

    private :: nucleosomeROT, nucleosomeTran
contains

subroutine nucleosomeProp(Uin,Vin,Rin,linkBP,wrapBP,Uout,Vout,Rout)
!Calculates the final position and orientation
!Rout,Uout,Vout
!based on the the incomming position and orientation
!Rin, Uin, Vin
!for a nucleosome with wrapBP bace paris of DNA wrapped around it.
!The intrisic rotation of a linkBP long linker that follows is also included.
    use vector_utils, only: cross
    implicit none
    real(dp), intent(in), dimension(3) :: Uin ! Entry direction along DNA
    real(dp), intent(in), dimension(3) :: Vin ! Entry tangent vector
    real(dp), intent(in), dimension(3) :: Rin ! Position of entry
    integer, intent(in) :: linkBP
    integer, intent(in) :: wrapBP
    real(dp), intent(out), dimension(3) :: Uout ! Exit position
    real(dp), intent(out), dimension(3) :: Vout ! Exit tangent vector
    real(dp), intent(out), dimension(3) :: Rout ! Exit position

    real(dp), dimension(3,3) :: linkRot
    real(dp), dimension(3,3) :: mtrx
    real(dp), parameter :: angle = 2*pi/basePairsPerTurn ! intrinsic rotation/bp

    mtrx(:,1) = Vin
    mtrx(:,2) = cross(Uin,Vin)
    mtrx(:,3) = Uin

    Rout = Rin + MATMUL(mtrx,nucleosomeTran(:,wrapBP))

    linkRot(1,:) = [cos(angle*linkBP),-sin(angle*linkBP), 0.0_dp]
    linkRot(2,:) = [sin(angle*linkBP), cos(angle*linkBP), 0.0_dp]
    linkRot(3,:) = [           0.0_dp,            0.0_dp, 1.0_dp]

    mtrx = MATMUL(MATMUL(mtrx,nucleosomeROT(:,:,wrapBP)),linkRot)

    Uout = mtrx(:,3)/norm2(mtrx(:,3))
    Vout = mtrx(:,1)/norm2(mtrx(:,1))
    return

end subroutine nucleosomeProp

function nucleosome_energy(RP1,R,UP1,U,VP1,V,linkBP,wrapBP)
!Caltulate the bending energy of a nucleosome
!at position R and orientation U, V
!with wrapBP of DNA wrapped around it followed by a liner linkBP
!bace pairs long that ends at RP1, UP1, and VP1.
!Assumes a SSWLCWT linker.
    use MC_wlc, only: E_SSWLCWT
    real(dp), intent(in), dimension(3) :: R ! R of bead i
    real(dp), intent(in), dimension(3) :: RP1 ! R of bead i+1
    real(dp), intent(in), dimension(3) :: U ! U of bead i
    real(dp), intent(in), dimension(3) :: UP1 ! U of bead i+1
    real(dp), intent(in), dimension(3) :: V ! V of bead i
    real(dp), intent(in), dimension(3) :: VP1 ! V of bead i+1
    integer, intent(in) :: linkBP
    integer, intent(in) :: wrapBP
    real(dp) nucleosome_energy(4)

    real(dp) Rtemp(3)
    real(dp) Utemp(3)
    real(dp) Vtemp(3)


    call nucleosomeProp(U,V,R,linkBP,wrapBP,Utemp,Vtemp,Rtemp)

    nucleosome_energy =  E_SSWLCWT(RP1,Rtemp,UP1,Utemp,VP1,Vtemp, &
                   multiParams(1,linkBP),&   ! EB
                   multiParams(2,linkBP),&   ! EPAR
                   multiParams(3,linkBP),&   ! EPERP
                   multiParams(5,linkBP),&   ! ETA
                   multiParams(4,linkBP),&   ! GAM
                   multiParams(9,linkBP))   ! etwist
end function nucleosome_energy

subroutine get_params(i,EB,EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype)
!Return parameters for a linker that is i bace pairs long
    implicit none
    integer, intent(in) :: i ! Number of bace pairs in linder
    real(dp), intent(out) :: EB ! Bending modulus
    real(dp), intent(out) :: EPAR ! Stretch modulus
    real(dp), intent(out) :: EPERP ! Shear modulus
    real(dp), intent(out) :: GAM ! ground state segment length
    real(dp), intent(out) :: ETA ! bend-shear coupling
    real(dp), intent(out) :: XIR,XIU,sigma
    real(dp), intent(out) :: etwist ! Twist coefficient
    real(dp), intent(out) :: simtype ! WLC, SSWLC, or Gaussian Chain

        EB     = multiParams(1,i)
        EPAR   = multiParams(2,i)
        EPERP  = multiParams(3,i)
        GAM    = multiParams(4,i)
        ETA    = multiParams(5,i)
        XIR    = multiParams(6,i)
        XIU    = multiParams(7,i)
        sigma  = multiParams(8,i)
        etwist = multiParams(9,i)
        simtype = int(multiParams(10,i)+0.1_dp) ! the +0.1 is to prevent round down due to lack of precision
end subroutine get_params

! ----------------------------------------------------------------------
!
!          Look up elastic constants for linkers,
!          Look up translation and rotation for nucleosomes.
!
! ----------------------------------------------------------------------
subroutine setup_nucleosome_constants()
    use precision, only: nan
    use MC_wlc, only: calc_elastic_constants
    implicit none
    integer i,j, simtype
    real(dp) DEL
    real(dp) EB, EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist, dt

    multiParams=nan
    do i = 2,WLC_P__MAX_BACEPAIRS_PER_BEAD
        DEL = WLC_P__LENGTH_PER_BP*i/WLC_P__LP
        call calc_elastic_constants(DEL,WLC_P__LP,WLC_P__LT,EB,EPAR,GAM,XIR,EPERP,ETA,XIU,DT,&
                                SIGMA,ETWIST,simtype)
        if (simtype .ne. 2) then
            print*, "Error!  nucleosomes only setup for SSWLCWT"
            print*, "Simtype ",simtype, " encountered"
            stop
        endif
        multiParams(1,i) = EB
        multiParams(2,i) = EPAR
        multiParams(3,i) = EPERP
        multiParams(4,i) = GAM
        multiParams(5,i) = ETA
        multiParams(6,i) = XIR
        multiParams(7,i) = XIU
        multiParams(8,i) = sigma
        multiParams(9,i) = etwist
        multiParams(10,i) = real(simtype,dp)

    enddo

    open (UNIT = 5, FILE = "input/nucleosomeR", STATUS = "OLD")
    do i = 1,147
        do j = 1,3
            read(5,*) nucleosomeROT(j,:,148-i)
        enddo
    enddo
    close(5)

    if (WLC_P__INCLUDE_NUC_TRANS) then
        open (UNIT = 5, FILE = "input/nucleosomeT", STATUS = "OLD")
        do i = 1,147
            read(5,*) nucleosomeTran(:,148-i)
        enddo
        close(5)
    else
        nucleosomeTran = 0.0_dp
    endif
end subroutine setup_nucleosome_constants

subroutine loadNucleosomePositions(wlc_nucleosomeWrap,wlc_basepairs)
    use precision, only: nan
    implicit none
    integer, intent(out) :: wlc_nucleosomeWrap(WLC_P__NT)
    integer, intent(out) :: wlc_basepairs(WLC_P__NT)

    ! In the future you can set up code here to choose nucleosome spacing
    wlc_nucleosomeWrap = 147
    wlc_basepairs = WLC_P__LL


end subroutine

end module nucleosome
