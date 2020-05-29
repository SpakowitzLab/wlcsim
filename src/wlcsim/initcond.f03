#include "../defines.inc"
!! ---------------------------------------------------------------*

!
!     This subroutine sets the initial condition for a chain within
!     the capsid
!
!     Andrew Spakowitz
!     Written 4-16-04
!
!     Updated by Quinn in 2016
!
subroutine initcond(R,U,NT,NP,FRMFILE,rand_stat, wlc_p)
! values from wlcsim_data
use params, only: wlc_V, wlc_ExplicitBindingPair, wlc_basepairs, wlc_nucleosomeWrap, &
                  wlc_network_start_index

!use mt19937, only : grnd, init_genrand, rnorm, mt, mti
use mersenne_twister
use params, only: dp, pi, wlcsim_params,  nan
use vector_utils, only: randomUnitVec, random_perp, cross
use nucleosome, only: nucleosomeProp, multiParams
use polydispersity, only: length_of_chain
implicit none

type(wlcsim_params), intent(in) :: wlc_p

integer NP,NT           ! Number of beads
real(dp) R(3,NT)  ! Bead positions
real(dp) U(3,NT)  ! Tangent vectors
real(dp) LBOX(3)  ! Box edge length
integer I,J,IB            ! Index Holders
LOGICAL FRMFILE           ! Is conformation in file?
!real(dp) RMin
real(dp) R0(3)

!     Varibles for type 2

real(dp) Uold(3) ! save previous direction
real(dp) Vold(3) ! save previous direction
real(dp) Rold(3) ! save previous position
real(dp) theta   ! random angle
real(dp) z       ! random z position
real(dp) rr      ! random radial position
LOGICAL search           ! for exiting while loop
real(dp) test(3) ! test position for inside confinment
real(dp) Rc      ! radius of confinement
integer ii !for testing

!   temporary variables
real(dp) mag    ! magnitude of U for reload, or of U when smoothing

!      Random number generator initiation
type(random_stat) rand_stat
real(dp) urand(3)
logical in_confinement

real(dp) center(3)
real(dp) RloopVec(3)
real(dp) perpVec(3)
real(dp) trash(3)
integer length
integer otherEnd
real(dp) nloops

LBOX(1)=WLC_P__LBOX_X
LBOX(2)=WLC_P__LBOX_Y
LBOX(3)=WLC_P__LBOX_Z

!     Setup the choice parameters

if (FRMFILE)then
   open (UNIT = 5, FILE = 'input/r0', STATUS = 'OLD')
   Do I = 1,NT
      READ(5,*) R(1,I),R(2,I),R(3,I)
   enddo
   CLOSE(5)

   open (UNIT = 5, FILE = 'input/u0', STATUS = 'OLD')
   Do I = 1,NT
      READ(5,*) U(1,I),U(2,I),U(3,I)
       mag = sqrt(U(1,I)**2 + U(2,I)**2 + U(3,I)**2)
       U(1,I) = U(1,I)/mag
       U(2,I) = U(2,I)/mag
       U(3,I) = U(3,I)/mag
   enddo
   CLOSE(5)

   if (WLC_P__LOCAL_TWIST) then
       open (UNIT = 5, FILE = 'input/v0', STATUS = 'OLD')
       Do I = 1,NT
          READ(5,*) wlc_V(1,I),wlc_V(2,I),wlc_V(3,I)
           wlc_V(:,I) = wlc_V(:,I)/norm2(wlc_V(:,I))
       enddo
   endif
   CLOSE(5)
   return
endif
!     Fix the initial condition
if (WLC_P__INITCONDTYPE.eq.'lineInYFromOrigin') then
! straight line in y direction starting at origin, equilibrium bond lengths
    iB = 1
    do i = 1,nP
        do j = 1,length_of_chain(I)
            R(1,iB) = 0.0_dp
            R(2,iB) = wlc_p%GAM*j
            R(3,iB) = 0.0_dp
            U(1,iB) = 0.0_dp
            U(2,iB) = 1.0_dp
            U(3,iB) = 0.0_dp
            if (WLC_P__LOCAL_TWIST) then
                wlc_V(:,IB) = [1.0_dp, 0.0_dp, 0.0_dp]
            endif
            iB = iB + 1
        enddo
    enddo

else if (WLC_P__INITCONDTYPE.eq.'lineInY') then
! staight line in y direction with random starting position

    IB = 1
    do I = 1,NP
       call random_number(urand,rand_stat)
       R0(1) = urand(1)*LBOX(1)
       R0(2) = urand(2)*LBOX(2)
       R0(3) = urand(3)*LBOX(3)

       length = length_of_chain(I)
       do J = 1,length
          R(1,IB) = R0(1)
          R(2,IB) = R0(2) + wlc_p%GAM*(J - length/2.0_dp - 0.5_dp) ! center on box
          R(3,IB) = R0(3)
          U(1,IB) = 0.0_dp
          U(2,IB) = 1.0_dp
          U(3,IB) = 0.0_dp
          if (WLC_P__LOCAL_TWIST) then
              wlc_V(:,IB) = [1.0_dp, 0.0_dp, 0.0_dp]
          endif
          IB = IB + 1

       enddo
    enddo

else if (WLC_P__INITCONDTYPE.eq.'randomLineSlitInZBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! slit boundary in z direction

    IB = 1
    do  I = 1,NP
        call random_number(urand,rand_stat)
        Rold(1) = urand(1)*LBOX(1)
        Rold(2) = urand(2)*LBOX(2)
        Rold(3) = urand(3)*LBOX(3)
        call randomUnitVec(Uold,rand_stat)
        if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)

        length = length_of_chain(I)
        do J = 1,length
           search = .TRUE.
           ii = 0
           do while(search)
                ii = ii + 1
                if (ii.gt.100) then
                    print*,'stuck in loop'
                    print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                    print*,'test = ',test(1),test(2),test(3)
                    exit
                endif
                test(1) = Rold(1) + Uold(1)*wlc_p%GAM
                test(2) = Rold(2) + Uold(2)*wlc_p%GAM
                test(3) = Rold(3) + Uold(3)*wlc_p%GAM
                search = .not. in_confinement(test, 1, 1, 1)
                if (search) then
                     call randomUnitVec(Uold,rand_stat)
                     if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
                endif
           enddo
           R(1,IB) = test(1)
           R(2,IB) = test(2)
           R(3,IB) = test(3)
           Rold(1) = test(1)
           Rold(2) = test(2)
           Rold(3) = test(3)
           U(1,IB) = Uold(1)
           U(2,IB) = Uold(2)
           U(3,IB) = Uold(3)
           if (WLC_P__LOCAL_TWIST) wlc_V(:,IB) = Vold

           IB = IB + 1
        enddo
    enddo
else if (WLC_P__INITCONDTYPE.eq.'randomLineCubeBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! square boundary

    IB = 1
    do  I = 1,NP
       call random_number(urand,rand_stat)
       if (WLC_P__BOUNDARY_TYPE == "ExtendBinsPast") then
           Rold(1) = urand(1)*(LBOX(1)-2*WLC_P__DBIN)+WLC_P__DBIN
           Rold(2) = urand(2)*(LBOX(2)-2*WLC_P__DBIN)+WLC_P__DBIN
           Rold(3) = urand(3)*(LBOX(3)-2*WLC_P__DBIN)+WLC_P__DBIN
       else
           Rold(1) = urand(1)*LBOX(1)
           Rold(2) = urand(2)*LBOX(2)
           Rold(3) = urand(3)*LBOX(3)
       endif
       call randomUnitVec(Uold,rand_stat)
       if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)

       length = length_of_chain(I)
       do J = 1,length
          search = .TRUE.
          ii = 0
          do while(search)
               ii = ii + 1
               if (ii.gt.100) then
                   print*,'stuck in loop'
                   print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                   print*,'test = ',test(1),test(2),test(3)
                   exit
               endif
               test(1) = Rold(1) + Uold(1)*wlc_p%GAM
               test(2) = Rold(2) + Uold(2)*wlc_p%GAM
               test(3) = Rold(3) + Uold(3)*wlc_p%GAM
               search = .not. in_confinement(test, 1, 1, 1)
               if (search) then
                    call randomUnitVec(Uold,rand_stat)
                    if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
               endif
          enddo
          R(1,IB) = test(1)
          R(2,IB) = test(2)
          R(3,IB) = test(3)
          Rold(1) = test(1)
          Rold(2) = test(2)
          Rold(3) = test(3)
          U(1,IB) = Uold(1)
          U(2,IB) = Uold(2)
          U(3,IB) = Uold(3)
          if (WLC_P__LOCAL_TWIST) wlc_V(:,IB) = Vold

          IB = IB + 1
       enddo
    enddo

else if (WLC_P__INITCONDTYPE.eq.'randomLineOutsideOfSphere') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! internalshpere boundary
    Rc = WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp ! use LBOX as radius
    IB = 1
    do  I = 1,NP
       search = .TRUE.
       do while(search)
           call random_number(urand,rand_stat)
           Rold(1) = urand(1)*LBOX(1)
           Rold(2) = urand(2)*LBOX(2)
           Rold(3) = urand(3)*LBOX(3)
           search = .not. in_confinement(Rold, 1, 1, 1)
       enddo
       call randomUnitVec(Uold,rand_stat)
       if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
       length = length_of_chain(I)
       do J = 1,length
           search = .TRUE.
           ii=0
           do while(search)
               test(1) = Rold(1) + Uold(1)*wlc_p%GAM
               test(2) = Rold(2) + Uold(2)*wlc_p%GAM
               test(3) = Rold(3) + Uold(3)*wlc_p%GAM
               search = .not. in_confinement(test, 1, 1, 1)
               if (search) then
                    ii = ii + 1
                    if (ii.gt.100) then
                        print*,'stuck in loop'
                        print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                        print*,'test = ',test(1),test(2),test(3)
                        print*,'diameter',WLC_P__CONFINEMENT_SPHERE_DIAMETER
                        print*,'center',WLC_P__LBOX_X/2.0_dp,&
                                        WLC_P__LBOX_Y/2.0_dp,&
                                        WLC_P__LBOX_Z/2.0_dp
                        stop
                    endif
                    call randomUnitVec(Uold,rand_stat)
                    if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
               endif
           enddo
           R(1,IB) = test(1)
           R(2,IB) = test(2)
           R(3,IB) = test(3)
           Rold(1) = test(1)
           Rold(2) = test(2)
           Rold(3) = test(3)
           U(1,IB) = Uold(1)
           U(2,IB) = Uold(2)
           U(3,IB) = Uold(3)
           if (WLC_P__LOCAL_TWIST) wlc_V(:,IB) = Vold
           IB = IB + 1
       enddo ! loop to N
    enddo ! loop to np
else if (WLC_P__INITCONDTYPE.eq.'randomLineSphereBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! shpere boundary
    ! radius of LBox/2 centered at LBox/2
    Rc = WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0_dp ! use LBOX as radius
    IB = 1
    do  I = 1,NP
       call random_number(urand,rand_stat)
       theta = urand(1)*2.0_dp*PI
       z = urand(2)*2.0_dp-1.0_dp
       rr = Rc*urand(3)  ! should have an r**2 from jacobian
       Rold(1) = sqrt(1.0_dp-z*z)*cos(theta)*rr + WLC_P__LBOX_X/2.0_dp
       Rold(2) = sqrt(1.0_dp-z*z)*sin(theta)*rr + WLC_P__LBOX_Y/2.0_dp
       Rold(3) = z*rr + WLC_P__LBOX_Z/2.0_dp
       call randomUnitVec(Uold,rand_stat)
       if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
       length = length_of_chain(I)
       do J = 1,length
           search = .TRUE.
           ii=0
           do while(search)
               test(1) = Rold(1) + Uold(1)*wlc_p%GAM
               test(2) = Rold(2) + Uold(2)*wlc_p%GAM
               test(3) = Rold(3) + Uold(3)*wlc_p%GAM
               search = .not. in_confinement(test, 1, 1, 1)
               if (search) then
                    ii = ii + 1
                    if (ii.gt.100) then
                        print*,'stuck in loop'
                        print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                        print*,'test = ',test(1),test(2),test(3)
                        print*,'diameter',WLC_P__CONFINEMENT_SPHERE_DIAMETER
                        print*,'center',WLC_P__LBOX_X/2.0_dp,&
                                        WLC_P__LBOX_Y/2.0_dp,&
                                        WLC_P__LBOX_Z/2.0_dp
                        stop
                    endif
                    call randomUnitVec(Uold,rand_stat)
                    if (WLC_P__LOCAL_TWIST) call random_perp(Uold,Vold,trash,rand_stat)
               endif
           enddo
           R(1,IB) = test(1)
           R(2,IB) = test(2)
           R(3,IB) = test(3)
           Rold(1) = test(1)
           Rold(2) = test(2)
           Rold(3) = test(3)
           U(1,IB) = Uold(1)
           U(2,IB) = Uold(2)
           U(3,IB) = Uold(3)
           if (WLC_P__LOCAL_TWIST) wlc_V(:,IB) = Vold
           IB = IB + 1
       enddo ! loop to N
    enddo ! loop to np
else if (WLC_P__INITCONDTYPE.eq.'randomlyDistributeBeadsInSphere') then
    ! randomly distribute beads in shereical confinement
    do IB = 1,NT
        search = .true.
        do while(search)
             call random_number(urand,rand_stat)
             test(1) = urand(1)*LBox(1)
             test(2) = urand(2)*LBox(2)
             test(3) = urand(3)*LBox(3)
             search = .not. in_confinement(test, 1, 1, 1)
        enddo
        R(1,IB) = test(1)
        R(2,IB) = test(2)
        R(3,IB) = test(3)
        U(1,IB) = 0.00_dp
        U(2,IB) = 1.00_dp
        U(3,IB) = 0.00_dp
        if (WLC_P__LOCAL_TWIST) wlc_V(:,IB) = [1.0_dp, 0.0_dp, 0.0_dp]
    enddo
else if (WLC_P__INITCONDTYPE == 'ring') then
    IB = 1
    do  I = 1,NP
        call random_number(urand,rand_stat)
        R0(1) = urand(1)*LBOX(1)
        R0(2) = urand(2)*LBOX(1)
        R0(3) = urand(3)*LBOX(1)
        length = length_of_chain(I)
        do  J = 1,length
             R(1,IB) = R0(1) + ((wlc_p%GAM*length)/(2*PI))*Cos(J*2.0_dp*PI/length)
             R(2,IB) = R0(2) + ((wlc_p%GAM*length)/(2*PI))*Sin(J*2.0_dp*PI/length)
             R(3,IB) = 0.0_dp
             U(1,IB) = -Sin(J*2.0_dp*PI/length)
             U(2,IB) = Cos(J*2.0_dp*PI/length)
             U(3,IB) = 0.0_dp;
             if (WLC_P__LOCAL_TWIST) then
                wlc_V(:,IB) = [0.0_dp, 0.0_dp, 1.0_dp]
             endif
             IB = IB + 1
        ENDdo
    ENDdo
elseif (WLC_P__INITCONDTYPE == 'multiRing') then
    IB = 1
    center(1) = WLC_P__LBOX_X/2.0_dp
    center(2) = WLC_P__LBOX_Y/2.0_dp
    center(3) = WLC_P__LBOX_Z/2.0_dp
    do while (IB .le. WLC_P__NT)
        if (IB == WLC_P__NT) then
            R(:,IB) = center
            call randomUnitVec(RloopVec,rand_stat)
            U(:,IB) = RloopVec
            if (WLC_P__LOCAL_TWIST) call random_perp(U(:,IB),wlc_V(:,IB),trash,rand_stat)
            exit
        endif
        otherEnd = IB+1
        do
            if (otherEnd == WLC_P__NT) exit
            if (WLC_P__NETWORK) then
                if (wlc_network_start_index(otherEnd) /= &
                    wlc_network_start_index(otherEnd+1)) exit
            else
                if (wlc_ExplicitBindingPair(otherEnd) /= -1) exit
            endif
            otherEnd=otherEnd+1
        enddo

        call randomUnitVec(RloopVec,rand_stat)
        !if (otherEnd == IB) then
        !    R(:,I) = center
        !    U(:,I) = RloopVec
        !    if (WLC_P__LOCAL_TWIST) call random_perp(U(:,I),wlc_V(:,I),trash,rand_stat)
        !endif
        call random_perp(RloopVec,perpVec,trash,rand_stat)

        length = ((otherEnd-IB)*wlc_p%GAM/(2.0_dp*PI))
        nloops = real(ceiling(length/(WLC_P__CONFINEMENT_SPHERE_DIAMETER*0.25_dp)),dp)
        length = length/nloops

        do I = IB,otherEnd
            R(:,I) = center + length*( &
                     + cos(2.0_dp*PI*nloops*(I-IB)/real(otherEnd-IB,dp))*RloopVec &
                     + sin(2.0_dp*PI*nloops*(I-IB)/real(otherEnd-IB,dp))*perpVec &
                     - RloopVec)
            U(:,I) = cos(2.0_dp*PI*nloops*(I-IB)/real(otherEnd-IB,dp))*perpVec &
                    -sin(2.0_dp*PI*nloops*(I-IB)/real(otherEnd-IB,dp))*RloopVec
            if (WLC_P__LOCAL_TWIST) wlc_V(:,I) = trash
        enddo
        IB=otherEnd+1
    enddo
elseif (WLC_P__INITCONDTYPE == 'nucleosome') then
    if (WLC_P__NEIGHBOR_BINS) then 
        R(1,1) = WLC_P__LBOX_X/2
        R(2,1) = WLC_P__LBOX_Y/2
        R(3,1) = WLC_P__LBOX_Z/2
    else
        R(1,1) = 0.0_dp
        R(2,1) = 0.0_dp
        R(3,1) = 0.0_dp
    endif
    U(1,1) = 0.0_dp
    U(2,1) = 0.0_dp
    U(3,1) = 1.0_dp
    wlc_V(1,1) = 1.0_dp
    wlc_V(2,1) = 0.0_dp
    wlc_V(3,1) = 0.0_dp

    ! NP EDIT : initialization is at odds with how energy is calculated 
    ! (i.e i-1 rotation to i vs i rotation to i+1, intuitvely the same physics but not the same results)
    do IB=1,WLC_P__NT-1
        ! Rotation (and translation) due to nucleosome
        call nucleosomeProp(U(:,IB), wlc_V(:,IB), R(:,IB), &
                            wlc_basepairs(IB),wlc_nucleosomeWrap(IB), &
                            U(:,IB+1), wlc_V(:,IB+1), R(:,IB+1))
        ! Translation due to zero-enery linker
        R(:,IB+1) = R(:,IB+1) + U(:,IB+1)*WLC_P__LENGTH_PER_BP*wlc_basepairs(IB)
    enddo
else if (WLC_P__INITCONDTYPE == 'WormlikeChain') then
    call effective_wormlike_chain_init(R, U, NT, wlc_p, rand_stat)
else if (WLC_P__INITCONDTYPE == 'randomWalkWithBoundary') then
    call gaus_init(R, U, NT, wlc_p, rand_stat)
else
    print*, "Unknown version of chain initialization WLC_P__INITCONDTYPE....."
    stop 1
endif

RETURN
END

subroutine wlc_init(R, U, NB, EPS, l0, rand_stat)
        ! takes R(3,NB) with R(:,1) preset and makes a WLC given EPS
    use mersenne_twister, only : random_number, random_stat
    use params, only : dp, pi
    use vector_utils, only: cross

    implicit none

    integer, intent(in) :: NB
    real(dp), intent(in) :: EPS, l0
    real(dp), intent(inout) :: R(3,NB), U(3,NB)
    type(random_stat), intent(inout) :: rand_stat

    integer J
    real(dp) N1(3), N2(3), z, theta
    real(dp) urand(3)

    do J = 2,NB

         call random_number(urand,rand_stat)
         theta = urand(1)*2.0_dp*pi
         z = (1.0_dp/EPS)*log(2.0_dp*sinh(EPS)*urand(2)+exp(-EPS))

         N1 = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
         N1 = N1 - dot_product(N1, U(:,J-1))*U(:,J-1)
         N1 = N1/norm2(N1)

         N2 = cross(N1, U(:,J-1))
         N2 = N2/norm2(N2)

         U(:,J) = sqrt(1-z*z)*(cos(theta)*N1 + sin(theta)*N2 + z*U(:,J-1))
         U(:,J) = U(:,J)/norm2(U(:,J))

         if (WLC_P__LOCAL_TWIST) then
             print*, "wlc chain initialization is not implimented for local twist"
             stop
         endif

         R(:,J) = R(:,J-1) + l0*U(:,J)
     enddo
end subroutine wlc_init

subroutine effective_wormlike_chain_init(R, U, NT, wlc_p, rand_stat)
    use mersenne_twister
    use params, only : wlcsim_params, dp, max_wlc_l0, maxWlcDelta
    use vector_utils, only: randomUnitVec
    use polydispersity, only: length_of_chain
    implicit none
    integer, intent(in) :: nt
    type(wlcsim_params), intent(in) :: wlc_p
    real(dp), intent(out) :: R(3,nt), U(3,nt)
    type(random_stat), intent(inout) :: rand_stat

    integer IB, NgB, i, j
    real(dp) l0, eps
    real(dp) urand(3)
    real(dp), dimension(:,:), allocatable :: tmpR, tmpU

    if (maxWlcDelta < wlc_p%DEL) then
        print *, "You can use gaussian chain-based initialization since your beads are so far apart, stopping..."
        stop 1
    end if
    if (WLC_P__LOCAL_TWIST) then
        print*, "wlc chain initialization is not implimented for local twist"
        stop
    endif
    IB = 1
    if (wlc_p%DEL > max_wlc_l0) then
        NgB = ceiling(wlc_p%DEL/max_wlc_l0) + 1
    else
        NgB = 1 + 1
    endif
    allocate(tmpR(3,NgB))
    allocate(tmpU(3,NgB))
    l0 = wlc_p%DEL/(NgB - 1)
    EPS = WLC_P__LP/l0 ! bending rigidity for wormlike chain
    do I = 1,WLC_P__NP
        ! uniformly first bead inside box
        call random_number(urand,rand_stat)
        R(1,IB) = urand(1)*WLC_P__LBOX_X
        R(2,IB) = urand(2)*WLC_P__LBOX_Y
        R(3,IB) = urand(3)*WLC_P__LBOX_Z
        ! uniformly from unit sphere first tan vec
        call randomUnitVec(U(:,IB),rand_stat)
        IB = IB + 1
        do J = 2,length_of_chain(I)
            tmpR(:,1) = R(:,IB-1)
            tmpU(:,1) = U(:,IB-1)
            call wlc_init(tmpR, tmpU, NgB, EPS, l0, rand_stat)
            R(:,IB) = tmpR(:,NgB)
            U(:,IB) = tmpU(:,NgB)
            IB = IB + 1
        enddo
    enddo
    deallocate(tmpR)
    deallocate(tmpU)
end subroutine effective_wormlike_chain_init

subroutine gaus_init(R, U, NT, wlc_p, rand_stat)
    ! values from wlcsim_data
    use params, only: wlc_V
    use params, only: wlcsim_params, dp
    use mersenne_twister
    use vector_utils, only: random_perp
    use polydispersity, only: length_of_chain

    implicit none
    integer, intent(in) :: NT
    type(wlcsim_params), intent(in) :: wlc_p
    real(dp), intent(out) :: R(3,NT), U(3,NT)
    type(random_stat), intent(inout) :: rand_stat
    real(dp) urand(3)
    real(dp) :: init_e2e(3)
    integer i, j, ib, length
    real(dp) trash(3)


    ! init_e2e makes easy to add fixed distances between specific beads in future
    call random_number(urand,rand_stat)
    if (WLC_P__RING) then
        init_e2e = 0.0_dp
        ib = 1
        do i = 1,WLC_P__NP
            length = length_of_chain(i)
            call make_rw_fix_end2end(R(:,IB:IB+length-1), length, init_e2e, wlc_p, rand_stat)
            IB = IB + length
        enddo
    else
        ib = 1
        do i = 1, WLC_P__NP
            length = length_of_chain(i)
            call make_rw_with_boundary(R(:,IB:IB+length-1), length, wlc_p, rand_stat)
            IB = IB + length
        enddo

    end if
    IB = 1
    do I = 1,WLC_P__NP
        U(1,IB) = R(1,IB + 1) - R(1,IB)
        U(2,IB) = R(2,IB + 1) - R(2,IB)
        U(3,IB) = R(3,IB + 1) - R(3,IB)
        U(:,IB) = U(:,IB)/norm2(U(:,IB))
        if (WLC_P__LOCAL_TWIST) call random_perp(U(:,IB),wlc_V(:,IB),trash,rand_stat)
        IB = IB + 1
        do J = 2,length_of_chain(I)-1
            U(1,IB) = R(1,IB + 1) - R(1,IB - 1)
            U(2,IB) = R(2,IB + 1) - R(2,IB - 1)
            U(3,IB) = R(3,IB + 1) - R(3,IB - 1)
            U(:,IB) = U(:,IB)/norm2(U(:,IB))
            if (WLC_P__LOCAL_TWIST) call random_perp(U(:,IB),wlc_V(:,IB),trash,rand_stat)
            IB = IB + 1
        enddo
        U(1,IB) = R(1,IB) - R(1,IB - 1)
        U(2,IB) = R(2,IB) - R(2,IB - 1)
        U(3,IB) = R(3,IB) - R(3,IB - 1)
        U(:,IB) = U(:,IB)/norm2(U(:,IB))
        if (WLC_P__LOCAL_TWIST) call random_perp(U(:,IB),wlc_V(:,IB),trash,rand_stat)
        IB = IB + 1
     enddo
end subroutine gaus_init

subroutine make_rw_fix_end2end(R, NB, e2e, wlc_p, rand_stat)
    ! for a GC, we should have
    ! SIGMA = sqrt(2.0_dp*WLC_P__LP*WLC_P__L/3.0_dp/WLC_P__NB)

    use params, only : dp, wlcsim_params
    use mersenne_twister

    implicit none

    type(wlcsim_params), intent(in) :: wlc_p
    type(random_stat) rand_stat  ! state of random number chain
    real(dp), intent(in) :: e2e(3)
    integer, intent(in) :: nb
    real(dp), intent(out) :: R(3,nb)
    integer :: j
    real(dp) :: actual_e2e(3)

    call make_rw_with_boundary(R, NB, wlc_p, rand_stat)
    actual_e2e = R(:,NB) - R(:,1)
    do J = 2,NB
        R(:,J) = R(:,J) - actual_e2e*(J-1)/(NB-1)
        R(:,J) = R(:,J) + e2e*(J-1)/(NB-1)
    enddo
end subroutine

subroutine make_rw_with_boundary(R, NB, wlc_p, rand_stat)
    ! for a GC, we should have
    ! SIGMA = sqrt(2.0_dp*WLC_P__LP*WLC_P__L/3.0_dp/WLC_P__NB)

    use params, only : dp, wlcsim_params
    use mersenne_twister

    implicit none

    type(wlcsim_params), intent(in) :: wlc_p
    type(random_stat) rand_stat  ! state of random number chain
    integer, intent(in) :: nb
    real(dp), intent(out) :: R(3,nb)
    integer :: ib, j
    real(dp) nrand(3)
    logical in_confinement, is_inside_boundary


    ! initialize as if it were a gaussian chain, first bead at zero
    !TODO add function to src/util/confinement to make this unif in confinement
    IB = 1
    R(1,IB) = 0.0_dp
    R(2,IB) = 0.0_dp
    R(3,IB) = 0.0_dp
    IB = IB + 1
    do J = 2,NB
        is_inside_boundary = .False.
        do while (.not. is_inside_boundary)
            call random_gauss(nrand, rand_stat)
            R(1,IB) = R(1,IB-1) + wlc_p%SIGMA*nrand(1)
            R(2,IB) = R(2,IB-1) + wlc_p%SIGMA*nrand(2)
            R(3,IB) = R(3,IB-1) + wlc_p%SIGMA*nrand(3)
            is_inside_boundary = in_confinement(R, NB, IB, IB)
        enddo
        IB = IB + 1
    enddo
end subroutine
