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
    real(dp), intent(in) :: linkBP
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
    real(dp), intent(in) :: linkBP
    integer, intent(in) :: wrapBP
    real(dp) nucleosome_energy(4)
    real(dp) EB, EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype


    real(dp) Rtemp(3)
    real(dp) Utemp(3)
    real(dp) Vtemp(3)


    call nucleosomeProp(U,V,R,linkBP,wrapBP,Utemp,Vtemp,Rtemp)

    !  need to interpolate between basepairs if fractional
    call get_params(linkBP,EB,EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype)

    nucleosome_energy =  E_SSWLCWT(RP1,Rtemp,UP1,Utemp,VP1,Vtemp, &
                   EB,&!multiParams(1,linkBP),&   ! EB
                   EPAR,&!multiParams(2,linkBP),&   ! EPAR
                   EPERP,&!multiParams(3,linkBP),&   ! EPERP
                   ETA,&!multiParams(5,linkBP),&   ! ETA
                   GAM,&!multiParams(4,linkBP),&   ! GAM
                   ETWIST)!multiParams(9,linkBP))   ! etwist
end function nucleosome_energy

!  ------------------------------------------------------
!
! internucleosomal energy (i.e. faces attracted via harmonic
! oscillator), need spring constant and preferred distance
!
!  ------------------------------------------------------
function internucleosome_energy(RI,RJ,UI,UJ,VI,VJ)
    use MC_wlc, only: E_SSWLCWT
    use vector_utils, only: cross
    use energies, only: energyOf, internucleosome_
    real(dp), intent(in), dimension(3) :: RI ! R of nuc i
    real(dp), intent(in), dimension(3) :: RJ! R of nuc j
    real(dp), intent(in), dimension(3) :: UI ! U of nuc i
    real(dp), intent(in), dimension(3) :: UJ ! U of nuc j
    real(dp), intent(in), dimension(3) :: VI ! V of nuc i
    real(dp), intent(in), dimension(3) :: VJ ! V of nuc j
    real(dp), parameter :: tau_faceface = 0.83
    real(dp), parameter :: e_faceface = 4.6
    real(dp), parameter :: tau_faceside = 0.35
    real(dp), parameter :: e_faceside = 1.5
    real(dp), parameter :: tau_sideside = 1.2
    real(dp), parameter :: e_sideside = 2.0
    real(dp), dimension(3), parameter :: center = [4.8455, -2.4445, 0.6694]
    real(dp), dimension(3,3) :: mtrxI, mtrxJ
    real(dp), dimension(3) :: polyI, faceI, faceItop, faceIbot
    real(dp), dimension(3) :: polyJ, faceJ, faceJtop, faceJbot
    real(dp), dimension(4,3) :: faceDistList
    real(dp), dimension(4) :: distList
    real(dp), dimension(3) :: distS, distC, dist
    real(dp) cospsi, costhetaS, costhetaC
    integer i, indList(1)
    real(dp) internucleosome_energy

    ! initialiaze
    internucleosome_energy = 0

    ! construct matrices
    mtrxI(:,1) = VI
    mtrxI(:,2) = cross(UI,VI)
    mtrxI(:,3) = UI
    mtrxJ(:,1) = VJ
    mtrxJ(:,2) = cross(UJ,VJ)
    mtrxJ(:,3) = UJ

    ! center of nucs
    polyI = RI + MATMUL(mtrxI, center)
    polyJ = RJ +  MATMUL(mtrxJ, center)

    ! construct face I normal vector (pointing up through face)
    faceItop = RI + MATMUL(mtrxI, center+[0.0_dp,WLC_P__NUCLEOSOME_HEIGHT/2,0.0_dp])
    faceIbot = RI + MATMUL(mtrxI, center+[0.0_dp,-WLC_P__NUCLEOSOME_HEIGHT/2,0.0_dp])
    faceI = faceItop-faceIbot
    ! construct face J normal vector (pointing up through face)
    faceJtop = RJ + MATMUL(mtrxJ, center+[0.0_dp,WLC_P__NUCLEOSOME_HEIGHT/2,0.0_dp])
    faceJbot = RJ + MATMUL(mtrxJ, center+[0.0_dp,-WLC_P__NUCLEOSOME_HEIGHT/2,0.0_dp])
    faceJ = faceJtop-faceJbot
    cospsi = dot_product(faceI/norm2(faceI),faceJ/norm2(faceJ))

    ! list of combinatorial face attractions
    faceDistList(1,:) = faceItop - faceJbot 
    faceDistList(2,:) = faceItop - faceJtop 
    faceDistList(3,:) = faceIbot - faceJbot 
    faceDistList(4,:) = faceIbot - faceJtop 

    ! find the closest faces to define the face oritentation vector
    do i = 1,4
        distList(i) = norm2(faceDistList(i,:))
    enddo
    indList = minloc(distList)
    distS = faceDistList(indList(1),:)
    costhetaS = dot_product(distS/norm2(distS),faceJ/norm2(faceJ))!-cospsi*faceJ/(abs(cospsi)*norm2(faceJ)))
    ! face-face (histone-histone attraction)
    if (norm2(distS) <= tau_faceface) then 
        internucleosome_energy = internucleosome_energy &
                - e_faceface*(cospsi**2)*(costhetaS**2)/tau_faceface
    else
        internucleosome_energy = internucleosome_energy &
                - e_faceface*(cospsi**2)*(costhetaS**2)/norm2(distS)
    endif
    distC = polyI-polyJ
    costhetaC = dot_product(distC/norm2(distC),faceJ/norm2(faceJ))!-cospsi*faceJ/(abs(cospsi)*norm2(faceJ)))
    ! face-side (histone-DNA attraction)
    dist = norm2(distC)-WLC_P__NUCLEOSOME_HEIGHT/2-WLC_P__NUCLEOSOME_RADIUS
    if (norm2(dist) <= tau_faceside) then 
        internucleosome_energy = internucleosome_energy &
                - e_faceside*(1-cospsi*2)*(costhetaC**2)/tau_faceside
    else
        internucleosome_energy = internucleosome_energy &
                - e_faceside*(1-cospsi**2)*(costhetaC**2)/norm2(dist)
    endif
    ! side-side (DNA-DNA attraction)
    dist = norm2(distC)-2*WLC_P__NUCLEOSOME_RADIUS
    if (norm2(dist) <= tau_sideside) then 
        internucleosome_energy = internucleosome_energy &
                - e_sideside*(1-costhetaC**2)/tau_sideside
    else
        internucleosome_energy = internucleosome_energy &
                - e_sideside*(1-costhetaC**2)/norm2(dist)
    endif

end function internucleosome_energy

! ---------------------------------------------------------------------
!  Return parameters for a linker i bace pairs long
! ---------------------------------------------------------------------
subroutine get_params(i,EB,EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype)
    implicit none
    real(dp), intent(in) :: i
    real(dp), intent(out) :: EB, EPAR,EPERP,GAM,ETA,XIR,XIU,sigma,etwist,simtype
    integer indUp, indDown
    real(dp) ratio, offratio

        ! interpolate between points
        indDown = floor(i)
        indUp = ceiling(i)
        ratio = i/indUp
        offratio = 1-ratio

        EB     = ratio*multiParams(1,indUp) + offratio*multiParams(1,indDown)
        EPAR   = ratio*multiParams(2,indUp) + offratio*multiParams(2,indDown)
        EPERP  = ratio*multiParams(3,indUp) + offratio*multiParams(3,indDown)
        GAM    = ratio*multiParams(4,indUp) + offratio*multiParams(4,indDown)
        ETA    = ratio*multiParams(5,indUp) + offratio*multiParams(5,indDown)
        XIR    = ratio*multiParams(6,indUp) + offratio*multiParams(6,indDown)
        XIU    = ratio*multiParams(7,indUp) + offratio*multiParams(7,indDown)
        sigma  = ratio*multiParams(8,indUp) + offratio*multiParams(8,indDown)
        etwist = ratio*multiParams(9,indUp) + offratio*multiParams(9,indDown)
        simtype = int(ratio*multiParams(10,indUp) + offratio*multiParams(10,indDown)+0.1_dp) ! the +0.1 is to prevent round down due to lack of precision
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
    use GJKAlgorithm, only: GJK, sameShapeTest, noIntersectX, intersectX, tangentX, runtimeTest5, runtimeTest6, &
                            noIntersectY, intersectY, tangentY, runtimeTest1, runtimeTest3, &
                            noIntersectZ, intersectZ, tangentZ, runtimeTest2, runtimeTest4
    use polydispersity, only: first_bead_of_chain, last_bead_of_chain
    implicit none
    integer, intent(out) :: wlc_nucleosomeWrap(WLC_P__NT)
    real(dp), intent(out) :: wlc_basepairs(WLC_P__NT)
    real(dp), parameter :: L_in_bp = WLC_P__L/WLC_P__LENGTH_PER_BP
    integer, parameter :: nNucs = nint((L_in_bp-WLC_P__LL)/(147+WLC_P__LL)) ! assuming all octasomes
    real(dp) discretization, off_discretization, num_link_beads
    integer iter, i

    ! choose nucleosome spacing
    print*, WLC_P__L0
    print*, nNucs, WLC_P__NB, WLC_P__LL
    if (WLC_P__INCLUDE_NUC_TRANS) then
        if (WLC_P__INCLUDE_DISCRETIZE_LINKER) then 
            ! figure out main discretization scheme
            discretization = WLC_P__LL/((WLC_P__NB-2-nNucs)/(nNucs+1)+1)
            call discretizationScheme(discretization, num_link_beads, off_discretization)
            ! print for sanity check
            print*, discretization, num_link_beads, off_discretization
            do i = 1, WLC_P__NP
                ! set first linker
                iter = first_bead_of_chain(i)
                wlc_basepairs(iter) = off_discretization
                wlc_nucleosomeWrap(iter:iter+num_link_beads-1) = 1
                iter = iter + 1
                wlc_basepairs(iter:iter+num_link_beads-2) = discretization
                ! set middle beads
                iter = iter + num_link_beads - 1
                do while (iter <= last_bead_of_chain(i)-num_link_beads)
                    wlc_nucleosomeWrap(iter) = 147
                    wlc_basepairs(iter) = off_discretization
                    iter = iter + 1
                    if (iter + num_link_beads - 2 <= last_bead_of_chain(i) - num_link_beads) then
                        wlc_basepairs(iter:iter+num_link_beads-2) = discretization  
                        wlc_nucleosomeWrap(iter:iter+num_link_beads-2) = 1 
                        iter = iter + num_link_beads - 1
                    endif
                enddo
                ! set last last linker
                wlc_basepairs(iter-1) = off_discretization
                wlc_nucleosomeWrap(iter:iter+num_link_beads-1) = 1
                wlc_basepairs(iter:iter+num_link_beads-1) = discretization
                ! set last wlc_basepairs to 0 as reminder that this is not an actual extension
                wlc_basepairs(last_bead_of_chain(i)) = 0
            enddo
            ! testing sterics here !
            if(WLC_P__GJK_STERICS) then
            do iter = 1, 10000 ! check to make sure GJK is not stochastic
                call sameShapeTest()
                call noIntersectX()
                call intersectX()
                call tangentX()
                call noIntersectY()
                call intersectY()
                call tangentY()
                call noIntersectZ()
                call intersectZ()
                call tangentZ()
                call runtimeTest1()
                call runtimeTest2()
                call runtimeTest3()
                call runtimeTest4()
                call runtimeTest5()
                call runtimeTest6()
            enddo
                print*, "SUCCESS: successful completion of all GJK collision unit tests"
            endif
        endif
    else
        wlc_nucleosomeWrap = 147
        wlc_basepairs = WLC_P__LL
    endif
    print*, wlc_basepairs
    print*, wlc_nucleosomeWrap

end subroutine

subroutine discretizationScheme(discretization, num_link_beads, off_discretization)
    implicit none
    real(dp), intent(inout) :: discretization
    real(dp), intent(out) :: num_link_beads
    real(dp), intent(out) :: off_discretization

    real(dp), parameter:: threshold = 0.0001
    real(dp) off_link_beads
    integer iter

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
