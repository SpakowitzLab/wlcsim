#include "../defines.inc"
!   --------------------------------------------------------------
!
!    This module is designed to handle the various geometrical
!    considerations of nucleosomes.
!
!   --------------------------------------------------------------

module nucleosome
    use precision, only: dp, pi
    implicit none
    real(dp), parameter :: basePairsPerTurn = 10.5_dp
    real(dp), dimension(10,WLC_P__MAX_BACEPAIRS_PER_BEAD) :: multiParams
    real(dp), dimension(3,3,147) :: nucleosomeROT
    real(dp), dimension(3,147) :: nucleosomeTran

    private :: nucleosomeROT, nucleosomeTran
contains

! -------------------------------------------------------------
!
!  nucleosomeProp calculates the final position and orientation
!  Rout,Uout,Vout
!  based on the the incomming position and orientation
!  Rin, Uin, Vin
!  for a nucleosome with wrapBP bace paris of DNA wrapped around it.
!  The intrisic rotation of a liner that follows.
!
! -------------------------------------------------------------
subroutine nucleosomeProp(Uin,Vin,Rin,linkBP,wrapBP,Uout,Vout,Rout)
    use vector_utils, only: cross
    implicit none
    real(dp), intent(in), dimension(3) :: Uin
    real(dp), intent(in), dimension(3) :: Vin
    real(dp), intent(in), dimension(3) :: Rin
    integer, intent(in) :: linkBP
    integer, intent(in) :: wrapBP
    real(dp), intent(out), dimension(3) :: Uout
    real(dp), intent(out), dimension(3) :: Vout
    real(dp), intent(out), dimension(3) :: Rout

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


!  ------------------------------------------------------
!
!  nucleosome_energy caltulate the bending energy of a nucleosome
!  at position R and orientation U, V
!  with wrapBP of DNA wrapped around it followed by a liner linkBP
!  bace pairs long that ends at RP1, UP1, and VP1.
!  It does this assuming a SSWLCWT linker.
!
!  ------------------------------------------------------
function nucleosome_energy(RP1,R,UP1,U,VP1,V,linkBP,wrapBP)
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


! ---------------------------------------------------------------------
!  Return parameters for a linker i bace pairs long
! ---------------------------------------------------------------------
subroutine get_params(i,EB,EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype)
    implicit none
    integer, intent(in) :: i
    real(dp), intent(out) :: EB, EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype

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
    ! sterics testing !
    use LineLineIntersection, only: LineLineIntersectionTestOverlapA1B1, LineLineIntersectionTestOverlapA1B2, &
      LineLineIntersectionTestOverlapA2B1, LineLineIntersectionTestOverlapA2B2, LineLineIntersectionTestSameLine, &
      LineLineIntersectionTestParallelOverlapA1B1, LineLineIntersectionTestParallelOverlapA1B2, &
      LineLineIntersectionTestParallelOverlapA2B1, LineLineIntersectionTestParallelOverlapA2B2, &
      LineLineIntersectionTestParallelA1B1, LineLineIntersectionTestParallelA1B2, LineLineIntersectionTestParallelA2B1, &
      LineLineIntersectionTestParallelA2B2, LineLineIntersectionTestIntersectA1, LineLineIntersectionTestIntersectA2, &
      LineLineIntersectionTestIntersectB1, LineLineIntersectionTestIntersectB2, LineLineIntersectionTestIntersectMiddle, &
      LineLineIntersectionTestIntersectProjectionCollideZ, LineLineIntersectionTestIntersectProjectionNoCollideZ, &
      LineLineIntersectionTestIntersectProjectionCollideY, LineLineIntersectionTestIntersectProjectionNoCollideY, &
      LineLineIntersectionTestIntersectProjectionCollideX, LineLineIntersectionTestIntersectProjectionNoCollideX, &
      LineLineIntersectionTestIntersectProjection
    use SphereLineIntersection, only: SphereLineIntersectionTestLineInside, SphereLineIntersectionTestLineInsideEdgeA1, &
      SphereLineIntersectionTestLineInsideEdgeA2, SphereLineIntersectionTestLineTangent, &
      SphereLineIntersectionTestLineOutsideA1, SphereLineIntersectionTestLineOutsideA2, &
      SphereLineIntersectionTestLineOutsideBoth, SphereLineIntersectionTestLineCloseA1, &
      SphereLineIntersectionTestLineCloseA2
    use SphereSphereIntersection, only: SphereSphereIntersectionTestAinB, SphereSphereIntersectionTestBinA, &
      SphereSphereIntersectionTestTangent, SphereSphereIntersectionTestOverlap, SphereSphereIntersectionTestNoOverlap
    implicit none
    integer, intent(out) :: wlc_nucleosomeWrap(WLC_P__NT)
    integer, intent(out) :: wlc_basepairs(WLC_P__NT)
    real(dp), parameter :: L_in_bp = WLC_P__L/WLC_P__LENGTH_PER_BP
    integer, parameter :: nNucs = nint((L_in_bp-WLC_P__LL)/(147+WLC_P__LL)) ! overhang LL 
    real(dp) discretization, num_link_beads
    real(dp) discretization_overhang, num_link_beads_overhang
    integer iter, off_discretization, off_discretization_overhang
    integer, parameter :: outFileUnit = 99
    LOGICAL isfile
    character(100) :: filename = 'data/discretization'

    ! In the future you can set up code here to choose nucleosome spacing
    print*, nNucs, WLC_P__NB, WLC_P__LL
    if (WLC_P__INCLUDE_NUC_TRANS) then
        if (WLC_P__INCLUDE_DISCRETIZE_LINKER) then 
            if (WLC_P__NT /= WLC_P__NB) then
                print*, "oops havent set up for multipolymer sims yet"
                stop
            endif
            ! figure out main discretization scheme
            discretization = WLC_P__LL/((WLC_P__NB-2-nNucs)/(nNucs+1)+1)
            call discretizationScheme(discretization, discretization, num_link_beads, &
                    off_discretization)
            ! figure out overhang discretization scheme
            discretization_overhang = WLC_P__LL/ ((WLC_P__NB - ((num_link_beads-1)*(nNucs-1) + nNucs)) / 2)
            call discretizationScheme(discretization_overhang, discretization_overhang, num_link_beads_overhang, &
                off_discretization_overhang)
            ! print for sanity check
            print*, discretization, num_link_beads, off_discretization
            print*, discretization_overhang, num_link_beads_overhang, off_discretization_overhang
            ! set first overhang
            iter = 1
            wlc_basepairs(iter) = off_discretization_overhang
            wlc_nucleosomeWrap(iter:iter+num_link_beads_overhang-1) = 1
            iter = iter + 1
            wlc_basepairs(iter:iter+num_link_beads_overhang-2) = discretization_overhang  
            ! set middle beads
            iter = iter + num_link_beads_overhang - 1
            do while (iter <= WLC_P__NT-num_link_beads_overhang)
                wlc_nucleosomeWrap(iter) = 147
                wlc_basepairs(iter) = off_discretization
                iter = iter + 1
                if (iter + num_link_beads - 2 <= WLC_P__NT - num_link_beads) then
                    wlc_basepairs(iter:iter+num_link_beads-2) = discretization  
                    wlc_nucleosomeWrap(iter:iter+num_link_beads-2) = 1 
                    iter = iter + num_link_beads - 1
                endif
            enddo
            ! set last overhang
            wlc_basepairs(iter-1) = off_discretization_overhang
            wlc_nucleosomeWrap(iter:iter+num_link_beads_overhang-1) = 1
            wlc_basepairs(iter:iter+num_link_beads_overhang-1) = discretization_overhang
            ! set last wlc_basepairs to 0 as reminder that this is not an actual extension
            wlc_basepairs(WLC_P__NT) = 0

            ! testing sterics here !
            if(WLC_P__CYLINDRICAL_CHAIN_EXCLUSION) then
                call LineLineIntersectionTestOverlapA1B1()
                call LineLineIntersectionTestOverlapA1B2()
                call LineLineIntersectionTestOverlapA2B1()
                call LineLineIntersectionTestOverlapA2B2()
                call LineLineIntersectionTestSameLine()
                call LineLineIntersectionTestParallelOverlapA1B1()
                call LineLineIntersectionTestParallelOverlapA1B2()
                call LineLineIntersectionTestParallelOverlapA2B1()
                call LineLineIntersectionTestParallelOverlapA2B2()
                call LineLineIntersectionTestParallelA1B1()
                call LineLineIntersectionTestParallelA1B2()
                call LineLineIntersectionTestParallelA2B1()
                call LineLineIntersectionTestParallelA2B2()
                call LineLineIntersectionTestIntersectA1()
                call LineLineIntersectionTestIntersectA2()
                call LineLineIntersectionTestIntersectB1()
                call LineLineIntersectionTestIntersectB2()
                call LineLineIntersectionTestIntersectMiddle()
                call LineLineIntersectionTestIntersectProjectionCollideZ()
                call LineLineIntersectionTestIntersectProjectionNoCollideZ()
                call LineLineIntersectionTestIntersectProjectionCollideY()
                call LineLineIntersectionTestIntersectProjectionNoCollideY()
                call LineLineIntersectionTestIntersectProjectionCollideX()
                call LineLineIntersectionTestIntersectProjectionNoCollideX()
                call LineLineIntersectionTestIntersectProjection()
                print*, "SUCCESS: successful completion of all 25 line-line collision unit tests"
                call SphereLineIntersectionTestLineInside()
                call SphereLineIntersectionTestLineInsideEdgeA1()
                call SphereLineIntersectionTestLineInsideEdgeA2()
                call SphereLineIntersectionTestLineTangent()
                call SphereLineIntersectionTestLineOutsideA1()
                call SphereLineIntersectionTestLineOutsideA2()
                call SphereLineIntersectionTestLineOutsideBoth()
                call SphereLineIntersectionTestLineCloseA1()
                call SphereLineIntersectionTestLineCloseA2()
                print*, "SUCCESS: successful completion of all 9 sphere-line collision unit tests"
                call SphereSphereIntersectionTestAinB()
                call SphereSphereIntersectionTestBinA()
                call SphereSphereIntersectionTestTangent()
                call SphereSphereIntersectionTestOverlap()
                call SphereSphereIntersectionTestNoOverlap()
                print*, "SUCCESS: successful completion of all 5 sphere-sphere collision unit tests"
            endif
        endif
    else
        wlc_nucleosomeWrap = 147
        wlc_basepairs = WLC_P__LL
    endif

    ! save discretization state
    inquire(file = filename, exist = isfile)
    if (isfile) then
        open (unit = outFileUnit, file = filename, status ='OLD', POSITION = "append")
    else
        open (unit = outFileUnit, file = filename, status = 'NEW')
    endif
        write(outFileUnit,*) wlc_nucleosomeWrap
        write(outFileUnit,*) wlc_basepairs
    close(outFileUnit)

end subroutine

subroutine discretizationScheme(discretizationIN, discretization, num_link_beads, off_discretization)
    implicit none
    real(dp), intent(in) :: discretizationIN
    real(dp), intent(out) :: discretization
    real(dp), intent(out) :: num_link_beads
    integer, intent(out) :: off_discretization

    real(dp), parameter:: threshold = 0.0001
    real(dp) off_link_beads
    integer iter

    discretization = discretizationIN
    ! figure out discretization scheme
        num_link_beads = WLC_P__LL/discretization
        if (modulo(num_link_beads,1.0) >= threshold) then 
            print*, "choose a bead value to discretize linker DNA at an integer number"
            stop
        endif
        off_discretization = WLC_P__LL-floor(discretization)*num_link_beads
        off_link_beads = num_link_beads
        do while ( modulo(discretization,1.0) >= threshold .OR. (off_discretization < 2 .AND. off_link_beads > 0))
            discretization = floor(discretization)
            off_link_beads = off_link_beads - 1
            off_discretization = WLC_P__LL-discretization*off_link_beads
        enddo
        if (off_discretization > 2*discretization) then
            print*, off_discretization, discretization
            print*, "offset discretization is twice the size of the actual discretization &
                    &(not sure if this is actually an issue but jic)"
            stop
        else if (off_discretization == 0) then
            off_discretization = discretization
        endif
        if (off_link_beads < 1) then
            print*, "lower discretization length or change linker/fragment length pairing"
            stop
        endif
end subroutine discretizationScheme

end module nucleosome
